use sentinel_checkm2::FeatureSpace;
use sentinel_core::{call_genes, read_fasta, Domain};
use std::collections::HashSet;
fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).unwrap();
    let ko_file = std::env::args().nth(2).unwrap();
    let contigs = read_fasta(&path)?;
    let gs = call_genes(&contigs, Domain::Bacteria, 1)?;
    let fs = FeatureSpace::load()?;
    let ko: HashSet<String> = std::fs::read_to_string(&ko_file)?
        .lines().map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect();
    let v = fs.build_vector(&gs.genes, &ko);
    std::fs::write("/tmp/sentinel_full_vector.tsv",
        v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\n"))?;
    eprintln!("Sentinel vector length {}, nonzero {}", v.len(), v.iter().filter(|&&x| x!=0.0).count());
    Ok(())
}
