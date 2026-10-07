//! Quickstart: the Rust analogue of the Python README's `slow_square` demo.
//!
//! Run it locally:
//! ```text
//! cargo run --release --example quickstart
//! ```
//!
//! Turn it into a cluster node without touching the code — HydraMPP reads its
//! own `--hydra-*` flags from the command line:
//! ```text
//! # on the head node
//! cargo run --release --example quickstart -- --hydra-host
//! # on each worker (pointing at the head's IP)
//! cargo run --release --example quickstart -- --hydra-client 10.0.0.1
//! ```

use std::thread;
use std::time::{Duration, Instant};

use hydra_mpp_core::prelude::*;

fn main() -> anyhow::Result<()> {
    // A deliberately slow task so the parallel speed-up is visible.
    fn slow_square(x: u64) -> u64 {
        thread::sleep(Duration::from_millis(200));
        x * x
    }

    let hydra = Hydra::new();
    hydra.register("slow_square", slow_square);

    // `Config::from_args()` consumes any `--hydra-*` flags and leaves the rest
    // for your own parser; with no flags this is plain local mode.
    let cfg = Config::from_args();
    hydra.init(cfg)?;

    let started = Instant::now();

    // How many CPUs the local node actually has (so the demo never asks for
    // more than exist — a job reserving more CPUs than any node provides would
    // wait until a big enough node joins).
    let local_cpus = hydra.nodes()?[0].total_cpus.max(1);
    let heavy_cpus = local_cpus.min(4);

    // Fire 20 ordinary jobs ... (note the explicit `u64`: bincode is not
    // self-describing, so the argument type must match the task signature).
    let mut ids: Vec<JobId> = (0..20u64)
        .map(|i| hydra.remote("slow_square", &i))
        .collect::<Result<_, _>>()?;
    // ... plus one "heavy" job that reserves several CPU slots.
    ids.push(hydra.task("slow_square").cpus(heavy_cpus).submit(&99u64)?);

    println!("submitted {} jobs; draining as they finish...\n", ids.len());

    // Drain: wait for one result at a time, then get() it (which pops it).
    let mut pending = ids.clone();
    while !pending.is_empty() {
        let (ready, still) = hydra.wait(&pending, Some(Duration::from_secs(30)), 1)?;
        for id in ready {
            let r = hydra.get(id)?;
            match r.value::<u64>() {
                Ok(v) => println!(
                    "  job {id:>2}: {} -> {v:<5} ({:.3}s on {})",
                    r.func_name,
                    r.runtime_secs,
                    r.hostname.as_deref().unwrap_or("?")
                ),
                Err(e) => println!("  job {id:>2}: FAILED: {e}"),
            }
        }
        pending = still;
    }

    println!(
        "\nall {} jobs done in {:.2}s (≈{:.2}s of work serialized would be {:.1}s)",
        ids.len(),
        started.elapsed().as_secs_f64(),
        0.2,
        ids.len() as f64 * 0.2
    );

    for n in hydra.nodes()? {
        println!(
            "  node {}{}: {}/{} CPUs free, {}/{} GPUs free",
            n.hostname,
            if n.is_local { " (local)" } else { "" },
            n.free_cpus,
            n.total_cpus,
            n.free_gpus,
            n.total_gpus
        );
    }

    hydra.shutdown();
    Ok(())
}
