// Build several diverse 21,241-dim vectors (varying KO load) and predict each.
use sentinel_checkm2::{FeatureSpace, Models};
use sentinel_core::{call_genes, read_fasta, Domain};
use std::collections::HashSet;
fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).unwrap();
    let comp_p = std::env::args().nth(2).unwrap();
    let cont_p = std::env::args().nth(3).unwrap();
    let contigs = read_fasta(&path)?;
    let gs = call_genes(&contigs, Domain::Bacteria, 1)?;
    let fs = FeatureSpace::load()?;
    let models = Models::load(&comp_p, &cont_p)?;
    let ko_counts = [0usize, 50, 150, 500, 1500, 5000];
    let mut vec_lines = Vec::new();
    let mut preds = Vec::new();
    for &n in &ko_counts {
        let mut ko = HashSet::new();
        // spread the selection across the KO space for variety
        let step = (fs.ko_genes.len() / n.max(1)).max(1);
        for k in fs.ko_genes.iter().step_by(step).take(n) { ko.insert(k.clone()); }
        let v = fs.build_vector(&gs.genes, &ko);
        let (comp, cont) = models.predict(&v);
        vec_lines.push(v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"));
        preds.push((n, comp, cont));
        println!("SENT\t{}\t{:.10}\t{:.10}", n, comp, cont);
    }
    std::fs::write("/tmp/multi_vectors.tsv", vec_lines.join("\n"))?;
    let js: Vec<_> = preds.iter().map(|(n,c,k)| serde_json::json!({"ko":n,"comp":c,"cont":k})).collect();
    std::fs::write("/tmp/sentinel_multi_pred.json", serde_json::to_string(&js)?)?;
    Ok(())
}
