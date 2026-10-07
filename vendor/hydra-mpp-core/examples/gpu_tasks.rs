//! GPU scheduling demo.
//!
//! HydraMPP treats GPUs as a reservable resource and pins concrete device ids to
//! each job. This example registers a task that reads its assigned devices
//! (race-free, via a thread-local) and reports the `CUDA_VISIBLE_DEVICES` string
//! it would hand to a GPU library or subprocess.
//!
//! Run it on any machine — even one with **no GPUs** — by simulating a GPU box:
//! ```text
//! cargo run --release --example gpu_tasks            # autodetect real GPUs
//! cargo run --release --example gpu_tasks -- --hydra-gpus 4   # simulate 4 GPUs
//! ```
//! On a real cluster the very same binary fans GPU jobs across nodes:
//! ```text
//! cargo run --release --example gpu_tasks -- --hydra-slurm "$SLURM_JOB_NODELIST"
//! ```

use std::thread;
use std::time::Duration;

use hydra_mpp_core::prelude::*;

/// A stand-in for real GPU work. It reads the devices HydraMPP pinned to this
/// task and returns the `CUDA_VISIBLE_DEVICES` value a kernel/subprocess would use.
fn gpu_work(idx: u64) -> String {
    let devices = hydra_mpp_core::current_gpus();
    let cvd = hydra_mpp_core::cuda_visible_devices();

    // In real code you would now launch your kernel or a subprocess, e.g.:
    //   std::process::Command::new("my_aligner")
    //       .env("CUDA_VISIBLE_DEVICES", hydra_mpp_core::cuda_visible_devices())
    //       .status()?;
    thread::sleep(Duration::from_millis(200)); // pretend to compute

    format!("job {idx:>2} → CUDA_VISIBLE_DEVICES=\"{cvd}\" (devices {devices:?})")
}

fn main() -> anyhow::Result<()> {
    let hydra = Hydra::new();
    hydra.register("gpu_work", gpu_work);

    // --hydra-gpus N simulates a GPU box on CPU-only hardware so the demo runs
    // anywhere; with no flag, GPUs are autodetected (nvidia-smi / CUDA_VISIBLE_DEVICES).
    hydra.init(Config::from_args())?;

    let node = &hydra.nodes()?[0];
    let n_gpus = node.total_gpus;
    println!(
        "node {} — {} CPUs, {} GPUs {:?}\n",
        node.hostname, node.total_cpus, n_gpus, node.gpu_ids
    );

    if n_gpus == 0 {
        println!("No GPUs detected. Re-run with `-- --hydra-gpus 4` to simulate a GPU node.");
        hydra.shutdown();
        return Ok(());
    }

    // Fire more GPU jobs than there are GPUs: at most `n_gpus` run at once, each
    // pinned to a distinct device. `.cpus(0)` makes them GPU-gated, not CPU-gated,
    // so device-level concurrency shows even on a single-CPU box. One heavier job
    // reserves two GPUs.
    let mut ids: Vec<JobId> = (0..(n_gpus as u64 * 2))
        .map(|i| hydra.task("gpu_work").cpus(0).gpus(1).submit(&i))
        .collect::<Result<_, _>>()?;
    if n_gpus >= 2 {
        ids.push(hydra.task("gpu_work").cpus(0).gpus(2).submit(&999u64)?);
    }

    println!("submitted {} GPU jobs (≤{} run concurrently):\n", ids.len(), n_gpus);

    let mut pending = ids;
    while !pending.is_empty() {
        let (ready, still) = hydra.wait(&pending, Some(Duration::from_secs(60)), 1)?;
        for id in ready {
            let r = hydra.get(id)?;
            // The devices the scheduler reserved (r.gpu_ids) match what the task
            // saw via the thread-local — proof the pinning is consistent.
            match r.value::<String>() {
                Ok(msg) => println!("  {msg}   [scheduler pinned {:?}]", r.gpu_ids),
                Err(e) => println!("  job failed: {e}"),
            }
        }
        pending = still;
    }

    // Demonstrate fail-fast: asking for more GPUs than exist errors immediately
    // instead of hanging the queue.
    let greedy = hydra.task("gpu_work").gpus(n_gpus + 100).submit(&0u64)?;
    let (_ready, _) = hydra.wait(&[greedy], Some(Duration::from_secs(5)), 1)?;
    if let Some(err) = hydra.get(greedy)?.error() {
        println!("\nfail-fast check: over-subscribed GPU job rejected → {err}");
    }

    hydra.shutdown();
    Ok(())
}
