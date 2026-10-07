use sentinel_checkm2::FeatureSpace;
use sentinel_core::{call_genes, read_fasta, Domain};
fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: metadata_faa <fasta>");
    let contigs = read_fasta(&path)?;
    let gs = call_genes(&contigs, Domain::Bacteria, 1)?;
    // write proteins to .faa exactly as called (no stop codon)
    let mut faa = String::new();
    for g in &gs.genes {
        faa.push_str(&format!(">{}\n{}\n", g.id(), String::from_utf8_lossy(&g.aa)));
    }
    std::fs::write("/tmp/sentinel_proteins.faa", &faa)?;
    // Sentinel metadata (22-dim: 20 AA counts + AALength + CDS)
    let fs = FeatureSpace::load()?;
    let md = fs.metadata(&gs.genes);
    let names = &fs.metadata_order;
    let mut out = serde_json::Map::new();
    for (n, v) in names.iter().zip(md.iter()) {
        out.insert(n.clone(), serde_json::json!(*v));
    }
    std::fs::write("/tmp/sentinel_metadata.json", serde_json::to_string(&out)?)?;
    eprintln!("wrote {} proteins; metadata: AALength={} CDS={}",
        gs.genes.len(), md[md.len()-2], md[md.len()-1]);
    Ok(())
}
