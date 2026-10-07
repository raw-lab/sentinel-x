//! Validation harness: load a LightGBM model + a TSV of feature vectors, print
//! predictions. Used to cross-check against Python LightGBM.
//! Usage: validate <model.gbm> <vectors.tsv>
use sentinel_gbm::Booster;
fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let booster = Booster::load(&args[1])?;
    let text = std::fs::read_to_string(&args[2])?;
    for line in text.lines() {
        let row: Vec<f64> = line.split('\t').map(|t| t.trim().parse().unwrap_or(0.0)).collect();
        println!("{:.8}", booster.predict(&row));
    }
    Ok(())
}
