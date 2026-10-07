use sentinel_checkv::{detect_terminal_repeat, TrKind};
use sentinel_core::read_fasta;
fn main() -> anyhow::Result<()> {
    let path = std::env::args().nth(1).expect("usage: tr_compare <fasta>");
    let contigs = read_fasta(&path)?;
    let mut out = serde_json::Map::new();
    let mut n_tr = 0;
    for c in &contigs {
        let tr = detect_terminal_repeat(&c.seq);
        let (t, l) = match tr.kind {
            TrKind::Dtr => ("DTR", tr.length),
            TrKind::Itr => ("ITR", tr.length),
            TrKind::None => ("None", 0),
        };
        if t != "None" { n_tr += 1; println!("{:<42}{:>8}{:>6}{:>8}", c.id, c.len(), t, l); }
        out.insert(c.id.clone(), serde_json::json!([c.len(), t, l]));
    }
    eprintln!("Sentinel: {}/{} contigs have a terminal repeat", n_tr, contigs.len());
    std::fs::write("/tmp/sentinel_tr.json", serde_json::to_string(&out)?)?;
    Ok(())
}
