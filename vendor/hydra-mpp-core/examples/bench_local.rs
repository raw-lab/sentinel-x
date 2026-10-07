//! Local micro-benchmark: measures HydraMPP's scheduling overhead and the
//! speed-up on CPU-bound work, using **real** wall-clock timings (nothing here
//! is fabricated — run it and read the numbers off your own machine).
//!
//! ```text
//! cargo run --release --example bench_local
//! ```
//!
//! It reports three things:
//!   1. dispatch overhead: how long ~50k trivial jobs take end-to-end, i.e. the
//!      per-job cost of submit → schedule → run → collect;
//!   2. parallel speed-up: a fixed amount of CPU-bound work done serially vs.
//!      through HydraMPP across all cores;
//!   3. effective throughput in jobs/second.

use std::time::Instant;

use hydra_mpp_core::prelude::*;

/// Pure CPU work: a small integer mixing loop (no allocation, no I/O).
fn burn(iters: u64) -> u64 {
    let mut acc = 0u64;
    for i in 0..iters {
        acc = acc.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(i | 1);
    }
    acc
}

fn drain_all(hydra: &Hydra, ids: &[JobId]) -> anyhow::Result<()> {
    // Collect results in bounded batches: each wait() returns as soon as a
    // batch is ready and get() pops them, so neither the queue nor the pending
    // set ever grows unbounded.
    let mut pending = ids.to_vec();
    while !pending.is_empty() {
        let batch = pending.len().min(256);
        let (ready, still) = hydra.wait(&pending, None, batch)?;
        for id in ready {
            let _ = hydra.get(id)?; // pops the finished job
        }
        pending = still;
    }
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let hydra = Hydra::new();
    hydra.register("noop", |_: u64| 0u64);
    hydra.register("burn", burn);
    hydra.init(Config::local().quiet(true))?;

    let cpus = hydra.nodes()?[0].total_cpus;
    println!("HydraMPP local benchmark — {cpus} CPUs\n");

    // ---- 1. Dispatch overhead (trivial jobs) ------------------------------
    let n_overhead = 50_000u64;
    let t = Instant::now();
    let ids: Vec<JobId> = (0..n_overhead)
        .map(|i| hydra.remote("noop", &i))
        .collect::<Result<_, _>>()?;
    drain_all(&hydra, &ids)?;
    let elapsed = t.elapsed().as_secs_f64();
    println!("dispatch overhead:");
    println!("  {n_overhead} no-op jobs in {elapsed:.3}s");
    println!("  = {:.1} jobs/s, {:.2} µs/job\n", n_overhead as f64 / elapsed, elapsed * 1e6 / n_overhead as f64);

    // ---- 2. Parallel speed-up on CPU-bound work ---------------------------
    let n_jobs = (cpus as u64) * 8;
    let iters = 5_000_000u64;

    // Serial baseline (call the function directly).
    let t = Instant::now();
    let mut sink = 0u64;
    for _ in 0..n_jobs {
        sink = sink.wrapping_add(burn(iters));
    }
    let serial = t.elapsed().as_secs_f64();
    std::hint::black_box(sink);

    // Through HydraMPP across all cores.
    let t = Instant::now();
    let ids: Vec<JobId> = (0..n_jobs)
        .map(|_| hydra.remote("burn", &iters))
        .collect::<Result<_, _>>()?;
    drain_all(&hydra, &ids)?;
    let parallel = t.elapsed().as_secs_f64();

    println!("parallel speed-up on CPU-bound work:");
    println!("  {n_jobs} jobs × {iters} iters");
    println!("  serial   : {serial:.3}s");
    println!("  HydraMPP : {parallel:.3}s");
    println!("  speed-up : {:.2}× on {cpus} CPUs\n", serial / parallel);

    hydra.shutdown();
    Ok(())
}
