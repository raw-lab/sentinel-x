use sentinel_checkm2::{FeatureSpace, Models};
use sentinel_core::{call_genes, read_fasta, Domain};
use std::collections::HashSet;
fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: full_vector <fasta> <comp.gbm> <cont.gbm>");
    let comp_p = std::env::args().nth(2).unwrap();
    let cont_p = std::env::args().nth(3).unwrap();
    let contigs = read_fasta(&path)?;
    let gs = call_genes(&contigs, Domain::Bacteria, 1)?;
    let fs = FeatureSpace::load()?;
    // make ~150 KOs present to exercise KO + pathway/module/category blocks
    let mut ko = HashSet::new();
    for k in fs.ko_genes.iter().step_by(130).take(150) { ko.insert(k.clone()); }
    let v = fs.build_vector(&gs.genes, &ko);
    // dump the vector for the Python side
    std::fs::write("/tmp/full_vector.tsv", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t"))?;
    // Sentinel prediction
    let models = Models::load(&comp_p, &cont_p)?;
    let (comp, cont) = models.predict(&v);
    println!("SENTINEL_PRED\t{:.10}\t{:.10}", comp, cont);
    std::fs::write("/tmp/sentinel_pred.json", serde_json::json!({"comp":comp,"cont":cont}).to_string())?;
    Ok(())
}
