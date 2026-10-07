//! # sentinel-checkm2
//!
//! Bacterial / archaeal completeness and contamination in the spirit of
//! **CheckM2**. The full 21,241-feature vector is assembled exactly as CheckM2
//! orders it —
//!
//! ```text
//! [ Metadata(22) | KO_Genes(19999) | KO_Pathways(416) | KO_Modules(757) | KO_Categories(47) ]
//! ```
//!
//! — where *Metadata* is the proteome's amino-acid composition plus total
//! length and CDS count (all database-free), and the KO blocks come from KEGG
//! KOfam annotation of the predicted proteins. Pathways, modules and categories
//! are aggregated from the KO presence/absence vector using the shipped KEGG
//! hierarchy. The assembled vector is scored by the **real** CheckM2 LightGBM
//! models through [`sentinel_gbm`].
//!
//! The KO annotations themselves require the (large) CheckM2 KOfam database;
//! when it is absent the KO blocks are zero and only the metadata contributes,
//! so callers should treat predictions as valid only when KO annotations were
//! supplied.

use std::collections::{HashMap, HashSet};

use anyhow::{Context, Result};
use sentinel_core::{Domain, Gene, QcRecord, Quality};
use sentinel_gbm::Booster;

const FEATURE_ORDERING: &str = include_str!("../../../assets/checkm2/feature_ordering.json");
const MODULE_DEFS: &str = include_str!("../../../assets/checkm2/module_definitions.json");
const PATH_CAT_MAP: &str = include_str!("../../../assets/checkm2/kegg_path_category_mapping.json");

/// The 20 standard amino acids, in the order CheckM2 stores their counts.
pub const AMINO_ACIDS: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

/// The KEGG feature hierarchy and column ordering loaded from CheckM2's data.
pub struct FeatureSpace {
    /// Metadata column order (22): the 20 AA codes then `AALength`, `CDS`.
    pub metadata_order: Vec<String>,
    /// KO gene column order (19999) and a KO → column-index map.
    pub ko_genes: Vec<String>,
    pub ko_index: HashMap<String, usize>,
    /// Pathway column order and, per pathway, the list of KEGG IDs it covers.
    pub pathways: Vec<String>,
    pathway_kos: Vec<Vec<String>>,
    /// Category column order and its KEGG-ID lists.
    pub categories: Vec<String>,
    category_kos: Vec<Vec<String>>,
    /// Module column order and its KEGG-ID lists.
    pub modules: Vec<String>,
    module_kos: Vec<Vec<String>>,
}

impl FeatureSpace {
    /// Total feature-vector length (should be 21,241).
    pub fn len(&self) -> usize {
        self.metadata_order.len()
            + self.ko_genes.len()
            + self.pathways.len()
            + self.modules.len()
            + self.categories.len()
    }
    /// True if no features are defined.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Load the embedded CheckM2 feature definitions.
    pub fn load() -> Result<Self> {
        let fo: serde_json::Value =
            serde_json::from_str(FEATURE_ORDERING).context("parsing feature_ordering.json")?;
        let md: serde_json::Value =
            serde_json::from_str(MODULE_DEFS).context("parsing module_definitions.json")?;
        let pcm: serde_json::Value =
            serde_json::from_str(PATH_CAT_MAP).context("parsing kegg_path_category_mapping.json")?;

        let as_str_list = |v: &serde_json::Value| -> Vec<String> {
            v.as_array()
                .map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect())
                .unwrap_or_default()
        };

        let metadata_order = as_str_list(&fo["Metadata"]);
        let ko_genes = as_str_list(&fo["KO_Genes"]);
        let pathways = as_str_list(&fo["KO_Pathways"]);
        let categories = as_str_list(&fo["KO_Categories"]);
        let modules = as_str_list(&fo["KO_Modules"]);

        let ko_index: HashMap<String, usize> = ko_genes
            .iter()
            .enumerate()
            .map(|(i, k)| (k.clone(), i))
            .collect();

        // path_category_mapping: parallel index-keyed maps Kegg_ID / KO_Pathways
        // / KO_Categories. Reconstruct the aligned rows, then bucket KEGG IDs by
        // pathway and by category (preserving multiplicity, as CheckM2 does).
        let kegg_id = &pcm["Kegg_ID"];
        let row_pathway = &pcm["KO_Pathways"];
        let row_cat = &pcm["KO_Categories"];

        let mut pathway_to_kos: HashMap<String, Vec<String>> = HashMap::new();
        let mut category_to_kos: HashMap<String, Vec<String>> = HashMap::new();
        if let Some(obj) = kegg_id.as_object() {
            for (idx, kid) in obj {
                let kid = kid.as_str().unwrap_or("").to_string();
                if let Some(p) = row_pathway.get(idx).and_then(|v| v.as_str()) {
                    pathway_to_kos.entry(p.to_string()).or_default().push(kid.clone());
                }
                if let Some(c) = row_cat.get(idx).and_then(|v| v.as_str()) {
                    category_to_kos.entry(c.to_string()).or_default().push(kid.clone());
                }
            }
        }
        let pathway_kos = pathways
            .iter()
            .map(|p| pathway_to_kos.get(p).cloned().unwrap_or_default())
            .collect();
        let category_kos = categories
            .iter()
            .map(|c| category_to_kos.get(c).cloned().unwrap_or_default())
            .collect();

        // module_definitions: module -> [KO,...]
        let module_kos = modules
            .iter()
            .map(|m| {
                md.get(m)
                    .and_then(|v| v.as_array())
                    .map(|a| a.iter().filter_map(|x| x.as_str().map(String::from)).collect())
                    .unwrap_or_default()
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

    /// Compute the 22 metadata features from a genome's predicted proteins:
    /// raw counts of each amino acid across the concatenated proteome, then
    /// `AALength` (total residues) and `CDS` (gene count).
    pub fn metadata(&self, genes: &[Gene]) -> Vec<f64> {
        let mut aa_counts = [0u64; 20];
        let mut aa_total = 0u64;
        for g in genes {
            for &b in &g.aa {
                aa_total += 1;
                if let Some(pos) = AMINO_ACIDS.iter().position(|&a| a == b) {
                    aa_counts[pos] += 1;
                }
            }
        }
        let mut out = Vec::with_capacity(self.metadata_order.len());
        for col in &self.metadata_order {
            let v = match col.as_str() {
                "AALength" => aa_total as f64,
                "CDS" => genes.len() as f64,
                aa if aa.len() == 1 => {
                    let b = aa.as_bytes()[0];
                    AMINO_ACIDS
                        .iter()
                        .position(|&a| a == b)
                        .map(|p| aa_counts[p] as f64)
                        .unwrap_or(0.0)
                }
                _ => 0.0,
            };
            out.push(v);
        }
        out
    }

    /// Assemble the full ordered feature vector from the metadata and the set of
    /// KOs annotated as present in the genome.
    pub fn build_vector(&self, genes: &[Gene], ko_present: &HashSet<String>) -> Vec<f64> {
        let mut v = Vec::with_capacity(self.len());

        // 1. Metadata.
        v.extend(self.metadata(genes));

        // 2. KO_Genes presence/absence (1 if present).
        let mut ko_bits = vec![0.0f64; self.ko_genes.len()];
        for ko in ko_present {
            if let Some(&i) = self.ko_index.get(ko) {
                ko_bits[i] = 1.0;
            }
        }
        v.extend_from_slice(&ko_bits);

        // helper: fraction of a KEGG-ID list that is present.
        let frac_present = |kos: &[String]| -> f64 {
            if kos.is_empty() {
                return 0.0;
            }
            let n = kos.iter().filter(|k| ko_present.contains(*k)).count();
            n as f64 / kos.len() as f64
        };

        // 3. KO_Pathways.
        for kos in &self.pathway_kos {
            v.push(frac_present(kos));
        }
        // 4. KO_Modules.
        for kos in &self.module_kos {
            v.push(frac_present(kos));
        }
        // 5. KO_Categories.
        for kos in &self.category_kos {
            v.push(frac_present(kos));
        }

        v
    }
}

/// The two real CheckM2 LightGBM models.
pub struct Models {
    pub comp: Booster,
    pub cont: Booster,
}

impl Models {
    /// Load completeness and contamination models from their `.gbm` paths.
    pub fn load(comp_path: &str, cont_path: &str) -> Result<Self> {
        Ok(Models {
            comp: Booster::load(comp_path).context("loading completeness model")?,
            cont: Booster::load(cont_path).context("loading contamination model")?,
        })
    }

    /// Predict `(completeness, contamination)` with CheckM2's output clamps
    /// (completeness to `[0, 100]`, contamination to `≥ 0`).
    pub fn predict(&self, features: &[f64]) -> (f64, f64) {
        let comp = self.comp.predict(features).clamp(0.0, 100.0);
        let cont = self.cont.predict(features).max(0.0);
        (comp, cont)
    }
}

/// Build one genome's QC record from an assembled feature vector and models.
pub fn assess(
    genome_id: &str,
    total_len: usize,
    domain: Domain,
    n_genes: usize,
    features: &[f64],
    models: &Models,
    ko_annotated: bool,
) -> QcRecord {
    let (comp, cont) = models.predict(features);
    let mut rec = QcRecord::new(genome_id, total_len, domain, "checkm2");
    rec.gene_count = n_genes;
    rec.completeness = Some(comp);
    rec.contamination = Some(cont);
    rec.method = "ml:general-lightgbm".into();
    rec.quality = Quality::from_completeness(Some(comp));
    if cont > 10.0 {
        rec.warnings
            .push(format!("elevated contamination {cont:.1}%"));
    }
    if !ko_annotated {
        rec.warnings.push(
            "no KO annotations supplied (KOfam DB absent) — prediction is metadata-only".into(),
        );
        rec.quality = Quality::NotDetermined;
    }
    rec
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gene_with_protein(seq: &[u8]) -> Gene {
        Gene {
            contig: "c".into(),
            index: 1,
            begin: 1,
            end: seq.len() * 3,
            strand: 1,
            aa: seq.to_vec(),
            partial: false,
        }
    }

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
    fn metadata_counts_amino_acids() {
        let fs = FeatureSpace::load().unwrap();
        // proteome "AAAC" + "GG": A=3,C=1,G=2; AALength=6; CDS=2
        let genes = vec![gene_with_protein(b"AAAC"), gene_with_protein(b"GG")];
        let md = fs.metadata(&genes);
        // metadata_order starts A,C,D,E,F,G,...
        assert_eq!(md[0], 3.0); // A
        assert_eq!(md[1], 1.0); // C
        assert_eq!(md[5], 2.0); // G (index 5 in ACDEFG...)
        // AALength, CDS are the last two
        let n = md.len();
        assert_eq!(md[n - 2], 6.0); // AALength
        assert_eq!(md[n - 1], 2.0); // CDS
    }

    #[test]
    fn vector_length_and_ko_presence() {
        let fs = FeatureSpace::load().unwrap();
        let genes = vec![gene_with_protein(b"ACDEFG")];
        let mut present = HashSet::new();
        // Use a KO that exists in the ordering.
        let some_ko = fs.ko_genes[0].clone();
        present.insert(some_ko.clone());
        let v = fs.build_vector(&genes, &present);
        assert_eq!(v.len(), 21241);
        // The KO bit for ko_genes[0] sits right after the 22 metadata features.
        assert_eq!(v[22], 1.0);
    }
}
