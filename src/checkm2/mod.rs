//! # CheckM2 1.1.0, ported
//!
//! Bacterial/archaeal completeness and contamination as CheckM2 computes
//! them:
//!
//! 1. **Genes** — Prodigal-style calls with translation tables 4 and 11
//!    (`single` mode from 100 kb, `meta` below, `meta` again if `single`
//!    finds nothing); table 4 wins when its coding density is > 0.7 and more
//!    than 0.05 above table 11's. Proteins keep Prodigal's terminal `*`.
//! 2. **Features** (21,241, CheckM2's order) — amino-acid counts, total
//!    protein length (with stops) and CDS count; per-KO protein counts from a
//!    DIAMOND-compatible search (`--max-target-seqs 1 --query-cover 80
//!    --subject-cover 80 --id 30 --evalue 1e-5`) against CheckM2's
//!    `uniref100.KO.1.dmnd`; pathway and category presence fractions; module
//!    sums of KO counts over the module definition length.
//! 3. **Models** — the general LightGBM completeness and contamination
//!    models; the specific (neural network) completeness model on
//!    MinMax-scaled features; CheckM2's choice between the two from the mean
//!    completeness, the protein-length ratio and the cosine similarity to its
//!    5,300 reference genomes; the divergence note.
//! 4. **Output** — `quality_report.tsv` with CheckM2's columns.

pub mod nn;

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::common::genes::{call_genes_with, CallOptions};
use crate::common::pyfmt::{np_round, py_float};
use crate::common::{Contig, Domain, GeneSet, QcRecord, Quality};
use crate::gbm::Booster;

const FEATURE_ORDERING: &str = include_str!("../../assets/checkm2/feature_ordering.json");
const MODULE_DEFS: &str = include_str!("../../assets/checkm2/module_definitions.json");
const PATH_CAT_MAP: &str = include_str!("../../assets/checkm2/kegg_path_category_mapping.json");

/// CheckM2's feature layout.
pub struct FeatureSpace {
    pub metadata_order: Vec<String>,
    pub ko_genes: Vec<String>,
    pub ko_index: HashMap<String, usize>,
    pub pathways: Vec<String>,
    /// KO indices per pathway (one per mapping row).
    pathway_kos: Vec<Vec<usize>>,
    pub categories: Vec<String>,
    category_kos: Vec<Vec<usize>>,
    pub modules: Vec<String>,
    /// KO indices per module (definition entries that are KO features,
    /// duplicates kept) and the full definition length.
    module_kos: Vec<(Vec<usize>, usize)>,
}

impl FeatureSpace {
    pub fn len(&self) -> usize {
        self.metadata_order.len() + self.ko_genes.len() + self.pathways.len() + self.modules.len() + self.categories.len()
    }
    pub fn load() -> Result<Self> {
        let fo: serde_json::Value = serde_json::from_str(FEATURE_ORDERING).context("parsing feature_ordering.json")?;
        let md: serde_json::Value = serde_json::from_str(MODULE_DEFS).context("parsing module_definitions.json")?;
        let pcm: serde_json::Value = serde_json::from_str(PATH_CAT_MAP).context("parsing kegg_path_category_mapping.json")?;
        let list = |v: &serde_json::Value| -> Vec<String> {
            v.as_array().map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect()).unwrap_or_default()
        };
        let metadata_order = list(&fo["Metadata"]);
        let ko_genes = list(&fo["KO_Genes"]);
        let pathways = list(&fo["KO_Pathways"]);
        let categories = list(&fo["KO_Categories"]);
        let modules = list(&fo["KO_Modules"]);
        let ko_index: HashMap<String, usize> = ko_genes.iter().enumerate().map(|(i, k)| (k.clone(), i)).collect();

        // mapping rows in JSON (= DataFrame) order
        let mut rows: Vec<(usize, &str, &str, &str)> = Vec::new();
        if let Some(obj) = pcm["Kegg_ID"].as_object() {
            for (idx, kid) in obj {
                let n: usize = idx.parse().unwrap_or(usize::MAX);
                rows.push((
                    n,
                    kid.as_str().unwrap_or(""),
                    pcm["KO_Pathways"].get(idx).and_then(|v| v.as_str()).unwrap_or(""),
                    pcm["KO_Categories"].get(idx).and_then(|v| v.as_str()).unwrap_or(""),
                ));
            }
        }
        rows.sort_by_key(|r| r.0);
        let mut by_path: HashMap<&str, Vec<usize>> = HashMap::new();
        let mut by_cat: HashMap<&str, Vec<usize>> = HashMap::new();
        for (_, kid, p, c) in &rows {
            if let Some(&i) = ko_index.get(*kid) {
                by_path.entry(p).or_default().push(i);
                by_cat.entry(c).or_default().push(i);
            }
        }
        let pathway_kos = pathways.iter().map(|p| by_path.get(p.as_str()).cloned().unwrap_or_default()).collect();
        let category_kos = categories.iter().map(|c| by_cat.get(c.as_str()).cloned().unwrap_or_default()).collect();
        let module_kos = modules
            .iter()
            .map(|m| {
                let def: Vec<String> = md.get(m).map(list).unwrap_or_default();
                let n = def.len();
                (def.iter().filter_map(|k| ko_index.get(k).copied()).collect(), n)
            })
            .collect();
        Ok(FeatureSpace {
            metadata_order,
            ko_genes,
            ko_index,
            pathways,
            pathway_kos,
            categories,
            category_kos,
            modules,
            module_kos,
        })
    }

    /// The 22 metadata features from proteins (with Prodigal's `*`): counts
    /// of the 20 amino acids, total length and CDS count.
    pub fn metadata(&self, proteins: &[Vec<u8>]) -> Vec<f64> {
        let mut aa = [0u64; 256];
        let mut total = 0u64;
        for p in proteins {
            total += p.len() as u64;
            for &b in p {
                aa[b as usize] += 1;
            }
        }
        self.metadata_order
            .iter()
            .map(|col| match col.as_str() {
                "AALength" => total as f64,
                "CDS" => proteins.len() as f64,
                c if c.len() == 1 => aa[c.as_bytes()[0] as usize] as f64,
                _ => 0.0,
            })
            .collect()
    }

    /// KO counts (per `ko_genes` index) from KO labels; labels that are not
    /// CheckM2 features are ignored, as CheckM2 does.
    pub fn ko_counts<'a, I: IntoIterator<Item = (&'a str, u32)>>(&self, kos: I) -> Vec<u32> {
        let mut v = vec![0u32; self.ko_genes.len()];
        for (k, n) in kos {
            if let Some(&i) = self.ko_index.get(k) {
                v[i] += n;
            }
        }
        v
    }

    /// The full feature vector.
    pub fn build_vector(&self, metadata: &[f64], ko: &[u32]) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.len());
        v.extend_from_slice(metadata);
        v.extend(ko.iter().map(|&c| c as f64));
        let presence = |kos: &[usize]| -> f64 {
            if kos.is_empty() {
                return 0.0;
            }
            kos.iter().filter(|&&i| ko[i] > 0).count() as f64 / kos.len() as f64
        };
        for kos in &self.pathway_kos {
            v.push(presence(kos));
        }
        for (kos, n) in &self.module_kos {
            let s: u64 = kos.iter().map(|&i| ko[i] as u64).sum();
            v.push(if *n == 0 { 0.0 } else { s as f64 / *n as f64 });
        }
        for kos in &self.category_kos {
            v.push(presence(kos));
        }
        v
    }
}

#[cfg_attr(not(feature = "embed-models"), allow(dead_code))]
fn inflate_text(z: &[u8]) -> Result<String> {
    use std::io::Read;
    let mut s = String::new();
    flate2::read::ZlibDecoder::new(z).read_to_string(&mut s).context("decompressing embedded model")?;
    Ok(s)
}

/// CheckM2's LightGBM models.
pub struct Models {
    pub comp: Booster,
    pub cont: Booster,
}

impl Models {
    pub fn load(comp_path: &str, cont_path: &str) -> Result<Self> {
        Ok(Models {
            comp: Booster::load(comp_path).context("loading completeness model")?,
            cont: Booster::load(cont_path).context("loading contamination model")?,
        })
    }

    #[cfg(feature = "embed-models")]
    pub fn load_embedded() -> Result<Self> {
        const COMP: &[u8] = include_bytes!("../../assets/checkm2/models/general_model_COMP.gbm.z");
        const CONT: &[u8] = include_bytes!("../../assets/checkm2/models/model_CONT.gbm.z");
        Ok(Models {
            comp: Booster::parse(&inflate_text(COMP)?).context("parsing embedded completeness model")?,
            cont: Booster::parse(&inflate_text(CONT)?).context("parsing embedded contamination model")?,
        })
    }

    pub fn load_embedded_opt() -> Option<Result<Self>> {
        #[cfg(feature = "embed-models")]
        {
            Some(Self::load_embedded())
        }
        #[cfg(not(feature = "embed-models"))]
        {
            None
        }
    }

    /// General completeness (clipped to [0, 100]) and contamination (≥ 0).
    pub fn predict(&self, features: &[f64]) -> (f64, f64) {
        (self.comp.predict(features).clamp(0.0, 100.0), self.cont.predict(features).max(0.0))
    }
}

pub const GENERAL: &str = "Gradient Boost (General Model)";
pub const SPECIFIC: &str = "Neural Network (Specific Model)";

/// A CheckM2 prediction.
#[derive(Clone, Debug)]
pub struct Prediction {
    pub general: f64,
    /// The network's completeness (float32 like Keras), if available.
    pub specific: Option<f32>,
    pub contamination: f64,
    /// The chosen completeness.
    pub completeness: f64,
    pub model: &'static str,
    pub cosine: Option<f64>,
    pub notes: String,
}

/// `cosine_decider` (modelPostprocessing.py).
pub fn choose_model(general: f64, specific: f64, cosine: f64, aa_length: f64) -> (f64, &'static str) {
    let novelty = general / (cosine * cosine);
    let mean = (general + specific) / 2.0;
    let aa_ratio = aa_length / mean;
    let pick = |limit: f64| if novelty < limit { (specific, SPECIFIC) } else { (general, GENERAL) };
    if mean < 55.0 && aa_ratio < 1500.0 {
        (general, GENERAL)
    } else if mean > 90.0 {
        pick(160.0)
    } else if mean > 80.0 || mean > 70.0 {
        pick(165.0)
    } else if mean > 60.0 {
        pick(170.0)
    } else if mean > 50.0 || mean > 40.0 {
        pick(175.0)
    } else {
        (specific, SPECIFIC)
    }
}

/// `__flag_divergent_predictions`.
pub fn divergence_note(general: f64, specific: f64) -> String {
    let diff = (general - specific).abs();
    if specific < 50.0 || diff < 25.0 {
        "None".into()
    } else {
        format!("Low confidence prediction - substantial ({}%) disagreement between completeness prediction models", diff as i64)
    }
}

/// Predict completeness/contamination for one feature vector.
pub fn predict(features: &[f64], models: &Models, nn: Option<&nn::SpecificModel>) -> Prediction {
    let (general, contamination) = models.predict(features);
    match nn {
        Some(m) => {
            let scaled = m.scale_features(features);
            let specific = m.predict(&scaled);
            let cosine = m.max_cosine(&scaled);
            let (completeness, model) = choose_model(general, specific as f64, cosine, features[20]);
            Prediction {
                general,
                specific: Some(specific),
                contamination,
                completeness,
                model,
                cosine: Some(cosine),
                notes: divergence_note(general, specific as f64),
            }
        }
        None => Prediction {
            general,
            specific: None,
            contamination,
            completeness: general,
            model: GENERAL,
            cosine: None,
            notes: "Specific model not available in this build; general model used".into(),
        },
    }
}

/// Gene calls chosen the way CheckM2 chooses them.
#[derive(Clone, Debug)]
pub struct GenomeCalls {
    pub genes: GeneSet,
    pub ttable: u32,
    pub coding_density: f64,
    /// `single` or `meta`.
    pub procedure: &'static str,
}

/// CheckM2's coding density: per contig, the union of `[start, end]` gene
/// intervals in a mask as long as the last coding base (so that base itself
/// is not counted), summed and divided by the genome size.
pub fn coding_density(contigs: &[Contig], genes: &GeneSet) -> f64 {
    let total: usize = contigs.iter().map(|c| c.len()).sum();
    if total == 0 {
        return 0.0;
    }
    let per = genes.per_contig(contigs.len());
    let mut coding = 0usize;
    for idx in per {
        if idx.is_empty() {
            continue;
        }
        let last = idx.iter().map(|&i| genes.genes[i].end).max().unwrap_or(0);
        let mut mask = vec![false; last];
        for &i in &idx {
            let g = &genes.genes[i];
            for m in mask.iter_mut().take(g.end + 1).skip(g.begin) {
                *m = true;
            }
        }
        coding += mask.iter().filter(|&&b| b).count();
    }
    coding as f64 / total as f64
}

/// Call genes with CheckM2's procedure and translation-table choice.
pub fn call_genes(contigs: &[Contig], threads: usize, ttable: Option<u32>) -> Result<GenomeCalls> {
    use dynamite::model::Mode;
    let total: usize = contigs.iter().map(|c| c.len()).sum();
    if total < 100_000 {
        // Prodigal's meta mode ignores -g: both tables give the same calls
        let genes = call_genes_with(
            contigs,
            Domain::Bacteria,
            threads,
            CallOptions { mode: Some(Mode::Meta), closed: false, code: None },
        )?;
        let cd = coding_density(contigs, &genes);
        return Ok(GenomeCalls { genes, ttable: ttable.unwrap_or(11), coding_density: cd, procedure: "meta" });
    }
    let tables: Vec<u32> = match ttable {
        Some(t) => vec![t],
        None => vec![4, 11],
    };
    let mut calls: Vec<(u32, GeneSet, f64, &'static str)> = Vec::new();
    for t in tables {
        let mut genes = call_genes_with(
            contigs,
            Domain::Bacteria,
            threads,
            CallOptions { mode: Some(Mode::Bacteria), closed: false, code: Some(t) },
        )?;
        let mut procedure = "single";
        if genes.genes.is_empty() {
            genes = call_genes_with(
                contigs,
                Domain::Bacteria,
                threads,
                CallOptions { mode: Some(Mode::Meta), closed: false, code: Some(t) },
            )?;
            procedure = "meta";
        }
        let cd = coding_density(contigs, &genes);
        calls.push((t, genes, cd, procedure));
    }
    let best = match ttable {
        Some(t) => t,
        None => {
            let cd = |t: u32| calls.iter().find(|c| c.0 == t).map(|c| c.2).unwrap_or(0.0);
            if cd(4) - cd(11) > 0.05 && cd(4) > 0.7 {
                4
            } else {
                11
            }
        }
    };
    let (t, genes, cd, procedure) = calls.into_iter().find(|c| c.0 == best).unwrap();
    Ok(GenomeCalls { genes, ttable: t, coding_density: cd, procedure })
}

/// CheckM2's per-genome statistics columns.
#[derive(Clone, Debug)]
pub struct GenomeStats {
    pub ttable: u32,
    pub coding_density: f64,
    pub contig_n50: f64,
    pub average_gene_length: f64,
    pub genome_size: usize,
    pub gc: f64,
    pub total_coding_sequences: usize,
    pub total_contigs: usize,
    pub max_contig_length: usize,
}

/// CheckM2's N50: the length-weighted median contig length (mean of the
/// two middle values for an even total).
pub fn checkm2_n50(lengths: &[usize]) -> f64 {
    let total: usize = lengths.iter().sum();
    if total == 0 {
        return 0.0;
    }
    let mut ls = lengths.to_vec();
    ls.sort_unstable();
    let at = |k: usize| -> usize {
        let mut acc = 0usize;
        for &l in &ls {
            acc += l;
            if k < acc {
                return l;
            }
        }
        *ls.last().unwrap()
    };
    if total % 2 == 0 {
        (at(total / 2 - 1) + at(total / 2)) as f64 / 2.0
    } else {
        at(total / 2) as f64
    }
}

pub fn genome_stats(contigs: &[Contig], calls: &GenomeCalls) -> GenomeStats {
    let lengths: Vec<usize> = contigs.iter().map(|c| c.len()).collect();
    let (mut gc, mut at) = (0u64, 0u64);
    for c in contigs {
        for &b in &c.seq {
            match b {
                b'G' | b'C' => gc += 1,
                b'A' | b'T' => at += 1,
                _ => {}
            }
        }
    }
    let n = calls.genes.genes.len();
    let total_len: u64 = calls.genes.genes.iter().map(|g| (g.aa.len() + g.stop as usize) as u64).sum();
    GenomeStats {
        ttable: calls.ttable,
        coding_density: calls.coding_density,
        contig_n50: checkm2_n50(&lengths),
        average_gene_length: if n == 0 { f64::NAN } else { total_len as f64 / n as f64 },
        genome_size: lengths.iter().sum(),
        gc: gc as f64 / (at + gc + 1) as f64,
        total_coding_sequences: n,
        total_contigs: lengths.len(),
        max_contig_length: lengths.iter().copied().max().unwrap_or(0),
    }
}

/// KO annotation of proteins against CheckM2's DIAMOND database
/// (`.dmnd` read in place, or a FASTA copy): KO = the part of the hit's id
/// after the first `~`; returns protein counts per KO label.
pub fn annotate_kos(proteins: &[(String, Vec<u8>)], db: &Path) -> Result<(HashMap<String, u32>, crate::protsearch::Stats)> {
    let (hits, stats) = crate::protsearch::search_path(proteins, db, &crate::protsearch::Params::checkm2())?;
    let mut counts: HashMap<String, u32> = HashMap::new();
    for h in hits {
        if let Some((_, ko)) = h.subject.split_once('~') {
            *counts.entry(ko.to_string()).or_insert(0) += 1;
        }
    }
    Ok((counts, stats))
}

/// KO counts from a user file: one `K#####` per line (one per protein), or
/// `K##### <count>`; any other columns (gene ids, `ko:` prefixes) are
/// tolerated.
pub fn parse_ko_file(text: &str) -> HashMap<String, u32> {
    let is_ko = |t: &str| t.len() == 6 && t.starts_with('K') && t[1..].bytes().all(|b| b.is_ascii_digit());
    let mut m: HashMap<String, u32> = HashMap::new();
    for line in text.lines() {
        let toks: Vec<&str> = line.split_whitespace().map(|t| t.trim_start_matches("ko:")).collect();
        let Some(pos) = toks.iter().position(|t| is_ko(t)) else { continue };
        let n = if toks.len() == 2 && pos == 0 { toks[1].parse::<u32>().unwrap_or(1) } else { 1 };
        *m.entry(toks[pos].to_string()).or_insert(0) += n;
    }
    m
}

/// One genome's CheckM2 result.
#[derive(Clone, Debug)]
pub struct Checkm2Result {
    pub name: String,
    pub prediction: Prediction,
    /// `None` when proteins were supplied (`--genes`), as CheckM2 omits
    /// the gene-calling columns then.
    pub stats: Option<GenomeStats>,
    /// Distinct KO features with a count > 0.
    pub kos_annotated: usize,
    pub ko_source: String,
}

/// `quality_report.tsv` (CheckM2's default `auto` mode columns).
pub fn quality_report_tsv(results: &[Checkm2Result]) -> String {
    let with_stats = results.iter().all(|r| r.stats.is_some());
    let mut out = String::from("Name\tCompleteness\tContamination\tCompleteness_Model_Used");
    if with_stats {
        out.push_str("\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length");
    }
    out.push_str("\tAdditional_Notes\n");
    let mut sorted: Vec<&Checkm2Result> = results.iter().collect();
    sorted.sort_by(|a, b| a.name.cmp(&b.name));
    for r in sorted {
        let p = &r.prediction;
        out.push_str(&format!(
            "{}\t{}\t{}\t{}",
            r.name,
            py_float(np_round(p.completeness, 2)),
            py_float(np_round(p.contamination, 2)),
            p.model
        ));
        if let (true, Some(s)) = (with_stats, &r.stats) {
            out.push_str(&format!(
                "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                s.ttable,
                py_float(np_round(s.coding_density, 3)),
                s.contig_n50 as i64,
                py_float(s.average_gene_length),
                s.genome_size,
                py_float(np_round(s.gc, 2)),
                s.total_coding_sequences,
                s.total_contigs,
                s.max_contig_length
            ));
        }
        out.push('\t');
        out.push_str(&p.notes);
        out.push('\n');
    }
    out
}

/// Unified-report record (MIMAG tier from completeness and contamination).
pub fn to_record(r: &Checkm2Result, length: usize, domain: Domain) -> QcRecord {
    let p = &r.prediction;
    let mut rec = QcRecord::new(&r.name, length, domain, "checkm2");
    rec.gene_count = r.stats.as_ref().map(|s| s.total_coding_sequences).unwrap_or(0);
    rec.method = p.model.into();
    if r.kos_annotated == 0 {
        // metadata-only features: the models are outside their training domain
        rec.completeness = None;
        rec.contamination = None;
        rec.quality = Quality::NotDetermined;
        rec.warnings.push(format!("no KO annotations ({}); completeness/contamination not reported", r.ko_source));
        rec = rec
            .with_extra("metadata_only_completeness", py_float(np_round(p.completeness, 2)))
            .with_extra("metadata_only_contamination", py_float(np_round(p.contamination, 2)));
    } else {
        rec.completeness = Some(np_round(p.completeness, 2));
        rec.contamination = Some(np_round(p.contamination, 2));
        rec.quality = Quality::from_mimag(rec.completeness, rec.contamination);
        if p.contamination > 10.0 {
            rec.warnings.push(format!("CheckM2 contamination {}%", py_float(np_round(p.contamination, 2))));
        }
    }
    if p.notes != "None" {
        rec.warnings.push(p.notes.clone());
    }
    rec = rec
        .with_extra("completeness_general", py_float(np_round(p.general, 2)))
        .with_extra("ko_source", r.ko_source.clone())
        .with_extra("kos_annotated", r.kos_annotated.to_string());
    if let Some(s) = p.specific {
        rec = rec.with_extra("completeness_specific", py_float(np_round(s as f64, 2)));
    }
    if let Some(c) = p.cosine {
        rec = rec.with_extra("cosine_similarity", py_float(np_round(c, 2)));
    }
    if let Some(s) = &r.stats {
        rec = rec
            .with_extra("translation_table", s.ttable.to_string())
            .with_extra("coding_density", py_float(np_round(s.coding_density, 3)))
            .with_extra("contig_n50", (s.contig_n50 as i64).to_string())
            .with_extra("gc_content", py_float(np_round(s.gc, 2)))
            .with_extra("max_contig_length", s.max_contig_length.to_string());
    }
    rec
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_space_is_21241() {
        let fs = FeatureSpace::load().unwrap();
        assert_eq!(fs.metadata_order.len(), 22);
        assert_eq!(fs.ko_genes.len(), 19999);
        assert_eq!(fs.pathways.len(), 416);
        assert_eq!(fs.modules.len(), 757);
        assert_eq!(fs.categories.len(), 47);
        assert_eq!(fs.len(), 21241);
    }

    #[test]
    fn metadata_counts_stops_in_length() {
        let fs = FeatureSpace::load().unwrap();
        let md = fs.metadata(&[b"AAAC*".to_vec(), b"GGX".to_vec()]);
        assert_eq!((md[0], md[1], md[5]), (3.0, 1.0, 2.0));
        assert_eq!(md[20], 8.0); // AALength includes '*' and X
        assert_eq!(md[21], 2.0);
    }

    #[test]
    fn ko_counts_and_modules_use_counts() {
        let fs = FeatureSpace::load().unwrap();
        // pick a module and give each of its KOs two copies
        let (kos, n) = fs.module_kos.iter().find(|(k, n)| !k.is_empty() && k.len() == *n).cloned().unwrap();
        let mut ko = vec![0u32; fs.ko_genes.len()];
        for &i in &kos {
            ko[i] = 2;
        }
        let md = vec![0.0; 22];
        let v = fs.build_vector(&md, &ko);
        let mi = fs.module_kos.iter().position(|(k, _)| *k == kos).unwrap();
        let col = 22 + fs.ko_genes.len() + fs.pathways.len() + mi;
        assert!((v[col] - 2.0 * kos.len() as f64 / n as f64).abs() < 1e-12, "{}", v[col]);
        assert_eq!(v[22 + kos[0]], 2.0);
    }

    #[test]
    fn model_choice_rules() {
        // highly reduced, low completeness -> general
        assert_eq!(choose_model(40.0, 42.0, 0.9, 1000.0 * 41.0).1, GENERAL);
        // high completeness, familiar -> specific
        assert_eq!(choose_model(98.0, 99.0, 0.95, 1e9).1, SPECIFIC);
        // high completeness but novel (cos small) -> general
        assert_eq!(choose_model(98.0, 99.0, 0.5, 1e9).1, GENERAL);
        assert_eq!(
            divergence_note(90.0, 60.0),
            "Low confidence prediction - substantial (30%) disagreement between completeness prediction models"
        );
        assert_eq!(divergence_note(90.0, 40.0), "None");
    }

    #[test]
    fn n50_and_ko_file() {
        // length-weighted median: [1,1,1,2,2] -> 1; [1,2,2,3,3,3] -> 2.5
        assert_eq!(checkm2_n50(&[2, 1, 1, 1]), 1.0);
        assert_eq!(checkm2_n50(&[3, 2, 1]), 2.5);
        assert_eq!(checkm2_n50(&[3, 2]), 3.0);
        let m = parse_ko_file("K00001\nK00001\ngene1\tK00002\nko:K00003\t5\nK00004 3\n#x\n");
        assert_eq!(m["K00001"], 2);
        assert_eq!(m["K00002"], 1);
        assert_eq!(m["K00003"], 5);
        assert_eq!(m["K00004"], 3);
    }

    #[test]
    fn coding_density_drops_the_last_base() {
        use crate::common::Gene;
        let c = Contig { id: "c".into(), header: "c".into(), seq: vec![b'A'; 100] };
        let g = |b: usize, e: usize| Gene {
            contig: "c".into(),
            contig_index: 0,
            index: 1,
            begin: b,
            end: e,
            strand: 1,
            aa: vec![],
            partial_left: false,
            partial_right: false,
            stop: true,
        };
        let gs = GeneSet { genes: vec![g(1, 30), g(21, 60)], trans_table: 11 };
        // union 1..60 = 60 bases, minus the last coding base -> 59
        assert!((coding_density(std::slice::from_ref(&c), &gs) - 0.59).abs() < 1e-12);
    }
}
