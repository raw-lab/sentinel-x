use sentinel_core::{read_fasta, genes::call_genes, seqio::Domain};
fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let contigs = read_fasta(&args[1])?;
    println!("read {} contigs", contigs.len());
    let gs = call_genes(&contigs, Domain::Virus, 1)?;
    println!("called {} genes (trans_table {})", gs.genes.len(), gs.trans_table);
    for g in gs.genes.iter().take(3) {
        println!("  {} {}..{} strand {} protein {} aa: {}...",
            g.id(), g.begin, g.end, g.strand, g.protein_len(),
            std::str::from_utf8(&g.aa[..g.aa.len().min(30)]).unwrap_or("?"));
    }
    Ok(())
}
