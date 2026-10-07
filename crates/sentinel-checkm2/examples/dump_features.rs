use sentinel_checkm2::FeatureSpace;
use std::collections::HashSet;
fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let kos: HashSet<String> = std::fs::read_to_string(&args[1])?
        .lines().map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect();
    let fs = FeatureSpace::load()?;
    let v = fs.build_vector(&[], &kos);
    let p0 = 22 + fs.ko_genes.len();
    let m0 = p0 + fs.pathways.len();
    let c0 = m0 + fs.modules.len();
    for (i, name) in fs.pathways.iter().enumerate() {
        let val = v[p0 + i];
        if val > 0.0 { println!("P\t{}\t{:.10}", name, val); }
    }
    for (i, name) in fs.modules.iter().enumerate() {
        let val = v[m0 + i];
        if val > 0.0 { println!("M\t{}\t{:.10}", name, val); }
    }
    for (i, name) in fs.categories.iter().enumerate() {
        let val = v[c0 + i];
        if val > 0.0 { println!("C\t{}\t{:.10}", name, val); }
    }
    Ok(())
}
