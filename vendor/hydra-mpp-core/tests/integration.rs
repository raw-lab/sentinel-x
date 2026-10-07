//! In-process integration tests (local mode — no network required).
//!
//! These lock down the behaviours most prone to regression, including the bugs
//! fixed during the audit: fail-fast admission (no queue wedge), panic
//! isolation, `cpus(0)` gating, and the bounded UDP status snapshot.

use std::time::Duration;

use hydra_mpp_core::prelude::*;

/// Each test uses a distinct status port so the per-instance UDP status servers
/// never collide when cargo runs tests in parallel.
fn local(port: u16) -> Config {
    Config::local().port(port).quiet(true)
}

#[test]
fn basic_compute_is_correct() {
    let h = Hydra::new();
    h.register("square", |x: u64| x * x);
    h.init(local(24611)).unwrap();
    let ids: Vec<JobId> = (0..16u64).map(|i| h.remote("square", &i).unwrap()).collect();
    let (ready, pending) = h.wait(&ids, Some(Duration::from_secs(30)), ids.len()).unwrap();
    assert!(pending.is_empty(), "all jobs should finish");
    assert_eq!(ready.len(), 16);

    let mut got: Vec<u64> = ids.iter().map(|&id| h.get(id).unwrap().value::<u64>().unwrap()).collect();
    got.sort_unstable();
    let want: Vec<u64> = (0..16u64).map(|i| i * i).collect();
    assert_eq!(got, want);
    h.shutdown();
}

#[test]
fn fail_fast_when_reserving_more_cpus_than_exist() {
    let h = Hydra::new();
    h.register("noop", |x: u64| x);
    h.init(local(24612)).unwrap();

    // No node has a million CPUs → must fail fast, not hang the queue.
    let id = h.task("noop").cpus(1_000_000).submit(&1u64).unwrap();
    let (ready, _) = h.wait(&[id], Some(Duration::from_secs(10)), 1).unwrap();
    assert_eq!(ready.len(), 1, "unplaceable job must resolve, not hang");

    let r = h.get(id).unwrap();
    assert!(!r.is_ok());
    assert!(r.error().unwrap().contains("no node provides"));
    h.shutdown();
}

#[test]
fn fail_fast_when_reserving_more_gpus_than_exist() {
    let h = Hydra::new();
    h.register("noop", |x: u64| x);
    // Simulate a 2-GPU box, then ask for 9.
    h.init(local(24613).num_gpus(2)).unwrap();

    let id = h.task("noop").cpus(0).gpus(9).submit(&1u64).unwrap();
    let (ready, _) = h.wait(&[id], Some(Duration::from_secs(10)), 1).unwrap();
    assert_eq!(ready.len(), 1);
    assert!(h.get(id).unwrap().error().unwrap().contains("no node provides"));
    h.shutdown();
}

#[test]
fn panicking_task_is_isolated_and_runtime_survives() {
    let h = Hydra::new();
    h.register("boom", |x: u64| -> u64 {
        if x == 0 {
            panic!("intentional test panic");
        }
        x * 2
    });
    h.init(local(24614)).unwrap();

    let bad = h.task("boom").submit(&0u64).unwrap();
    let (ready, _) = h.wait(&[bad], Some(Duration::from_secs(10)), 1).unwrap();
    assert_eq!(ready.len(), 1);
    let r = h.get(bad).unwrap();
    assert!(!r.is_ok());
    assert!(r.error().unwrap().to_lowercase().contains("panic"));

    // The runtime must still work after a task panic.
    let good = h.remote("boom", &21u64).unwrap();
    assert_eq!(h.get_blocking::<u64>(good).unwrap(), 42);
    h.shutdown();
}

#[test]
fn cpus_zero_jobs_are_not_throttled() {
    let h = Hydra::new();
    h.register("id", |x: u64| x);
    h.init(local(24615)).unwrap();

    // Even on a 1-CPU machine, cpus(0) jobs are gated only by (absent) GPUs.
    let ids: Vec<JobId> = (0..32u64)
        .map(|i| h.task("id").cpus(0).submit(&i).unwrap())
        .collect();
    let (ready, pending) = h.wait(&ids, Some(Duration::from_secs(30)), ids.len()).unwrap();
    assert!(pending.is_empty());
    assert_eq!(ready.len(), 32);
    h.shutdown();
}

#[test]
fn put_injects_a_finished_result() {
    let h = Hydra::new();
    h.init(local(24616)).unwrap();
    let id = h.put("constant", &vec![1u8, 2, 3]).unwrap();
    let r = h.get(id).unwrap();
    assert!(r.finished && r.is_ok());
    assert_eq!(r.value::<Vec<u8>>().unwrap(), vec![1, 2, 3]);
    h.shutdown();
}

#[test]
fn status_snapshot_counts_are_bounded_and_correct() {
    let port = 24617;
    let h = Hydra::new();
    // A task that blocks briefly so we can observe non-trivial counts.
    h.register("slow", |x: u64| {
        std::thread::sleep(Duration::from_millis(50));
        x
    });
    h.init(local(port)).unwrap();

    // Submit far more than the sample cap to exercise truncation + counts.
    let n = 500u64;
    let ids: Vec<JobId> = (0..n).map(|i| h.task("slow").cpus(0).submit(&i).unwrap()).collect();

    // Query the live UDP status server.
    let snap = hydra_mpp_core::query_status("127.0.0.1", port, Duration::from_secs(2)).unwrap();
    assert_eq!(snap.total, n as usize, "counts cover the whole queue");
    assert!(snap.queue.len() <= 100, "sample is bounded");
    assert!(snap.queue_truncated, "sample omits entries for a large queue");
    assert_eq!(snap.pending + snap.running + snap.done, n as usize);

    // Drain so shutdown is clean.
    let (_r, pending) = h.wait(&ids, Some(Duration::from_secs(60)), ids.len()).unwrap();
    assert!(pending.is_empty());
    h.shutdown();
}

#[test]
fn map_returns_results_in_input_order() {
    let h = Hydra::new();
    h.register("dbl", |x: u64| x * 2);
    h.init(local(24618)).unwrap();

    // Out-of-order completion must not scramble results.
    let out: Vec<u64> = h.map("dbl", &[5u64, 1, 4, 2, 3]).unwrap();
    assert_eq!(out, vec![10, 2, 8, 4, 6]);

    // Empty input is a no-op, not an error.
    let empty: Vec<u64> = h.map("dbl", &[] as &[u64]).unwrap();
    assert!(empty.is_empty());
    h.shutdown();
}

#[test]
fn map_propagates_a_failing_task() {
    let h = Hydra::new();
    h.register("inv", |x: u64| -> u64 {
        if x == 0 {
            panic!("boom");
        }
        100 / x
    });
    h.init(local(24619)).unwrap();

    // One panicking element must surface as an error for the whole batch.
    let res = h.map::<u64, u64>("inv", &[4u64, 2, 0, 1]);
    assert!(res.is_err(), "a panicking task should abort the batch");
    h.shutdown();
}

#[test]
fn warm_pool_computes_correctly_when_forced() {
    // Auto won't enable the pool on a single core, so force it on to exercise
    // the warm-pool code path (enqueue → worker → finish) for CPU-bound jobs.
    let h = Hydra::new();
    h.register("square", |x: u64| x * x);
    h.init(local(24620).pool_mode(PoolMode::Always)).unwrap();

    let inputs: Vec<u64> = (0..32u64).collect();
    let out: Vec<u64> = h.map("square", &inputs).unwrap();
    assert_eq!(out, inputs.iter().map(|i| i * i).collect::<Vec<_>>());
    h.shutdown();
}

#[test]
fn warm_pool_still_pins_gpus() {
    // A CPU-bound job that *also* reserves a GPU is routed through the pool;
    // the pool worker must still pin a concrete device for the task to read.
    let h = Hydra::new();
    h.register("which_gpu", |_x: u64| hydra_mpp_core::current_gpus());
    h.init(local(24621).num_gpus(2).pool_mode(PoolMode::Always)).unwrap();

    let id = h.task("which_gpu").cpus(1).gpus(1).submit(&0u64).unwrap();
    let gpus: Vec<usize> = h.get_blocking(id).unwrap();
    assert_eq!(gpus.len(), 1, "pooled GPU job must be pinned to one device");
    h.shutdown();
}

#[test]
fn remote_many_matches_per_task_and_keeps_order() {
    let h = Hydra::new();
    h.register("inc", |x: u64| x + 1);
    h.init(local(24622)).unwrap();

    let args: Vec<u64> = (0..50).collect();
    let ids = h.remote_many("inc", &args, 1).unwrap();
    assert_eq!(ids.len(), 50);
    let (ready, pending) = h.wait(&ids, Some(Duration::from_secs(30)), ids.len()).unwrap();
    assert!(pending.is_empty());
    assert_eq!(ready.len(), 50);

    // ids are in input order, so result[k] corresponds to args[k].
    let got: Vec<u64> = ids
        .iter()
        .map(|&id| h.get(id).unwrap().value::<u64>().unwrap())
        .collect();
    assert_eq!(got, (1..=50u64).collect::<Vec<_>>());
    h.shutdown();
}

#[test]
fn map_owned_zero_copy_local_matches_map() {
    let h = Hydra::new();
    h.register("double", |x: u64| x * 2);
    h.init(local(24623)).unwrap();
    let args: Vec<u64> = (0..40).collect();
    // local runtime → exercises the zero-copy typed payload path
    let got: Vec<u64> = h.map_owned("double", args.clone()).unwrap();
    let expect: Vec<u64> = args.iter().map(|x| x * 2).collect();
    assert_eq!(got, expect);
    h.shutdown();
}

#[test]
fn map_owned_handles_owned_heap_values() {
    let h = Hydra::new();
    h.register("sum", |v: Vec<f64>| v.iter().sum::<f64>());
    h.init(local(24624)).unwrap();
    let args: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64; 100]).collect();
    let got: Vec<f64> = h.map_owned("sum", args).unwrap();
    assert_eq!(got.len(), 10);
    assert_eq!(got[3], 3.0 * 100.0);
    h.shutdown();
}

#[test]
fn for_each_streams_in_order_and_isolates_failures() {
    let h = Hydra::new();
    h.register("checked", |x: i64| {
        if x == 3 {
            panic!("boom on 3")
        } else {
            x * 10
        }
    });
    h.init(local(24625)).unwrap();
    let args: Vec<i64> = (0..6).collect();
    let mut seen: Vec<(usize, bool, i64)> = Vec::new();
    h.for_each::<i64, i64, _>("checked", &args, 2, |i, r| match r {
        Ok(v) => seen.push((i, true, v)),
        Err(_) => seen.push((i, false, 0)),
    })
    .unwrap();
    assert_eq!(seen.len(), 6);
    // delivered in input order
    for (k, (i, _, _)) in seen.iter().enumerate() {
        assert_eq!(*i, k);
    }
    // index 3 panicked → delivered as Err; the batch was NOT aborted
    assert!(!seen[3].1, "task 3 should have failed");
    assert!(seen[0].1 && seen[5].1, "others should succeed");
    assert_eq!(seen[0].2, 0);
    assert_eq!(seen[5].2, 50);
    h.shutdown();
}

#[test]
fn map_owned_zero_copy_result_is_typed_roundtrip() {
    let h = Hydra::new();
    // big owned arg AND big owned result: both sides are zero-copy locally
    h.register("scale", |v: Vec<f64>| v.iter().map(|x| x * 2.0).collect::<Vec<f64>>());
    h.init(local(24626)).unwrap();
    let args: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64; 64]).collect();
    let out: Vec<Vec<f64>> = h.map_owned("scale", args).unwrap();
    assert_eq!(out.len(), 5);
    assert_eq!(out[2], vec![4.0; 64]); // 2*2
    assert_eq!(out[4][0], 8.0); // 4*2
    h.shutdown();
}

#[test]
fn concurrent_waiters_wake_on_completion_under_submit_storm() {
    use std::sync::Arc as StdArc;
    use std::thread;
    use std::time::Duration;

    let h = StdArc::new(Hydra::new());
    h.register("work", |x: u64| {
        thread::sleep(Duration::from_millis(2)); // ensure waiters actually block
        x * 2
    });
    h.init(local(24630)).unwrap();

    // 16 threads each submit a job and block on ITS result. With the condvar
    // split, a submit wakes only the scheduler — these waiters must wake on
    // their own completions, not on the submit storm below.
    let mut handles = Vec::new();
    for i in 0..16u64 {
        let h = h.clone();
        handles.push(thread::spawn(move || -> u64 {
            let id = h.remote("work", &i).unwrap();
            h.get_blocking::<u64>(id).unwrap()
        }));
    }
    // Submit storm from another thread (wakes work_cv only).
    let h2 = h.clone();
    let storm = thread::spawn(move || {
        for j in 1000..1300u64 {
            let _ = h2.remote("work", &j).unwrap();
        }
    });

    let sum: u64 = handles.into_iter().map(|h| h.join().unwrap()).sum();
    storm.join().unwrap();
    assert_eq!(sum, (0..16u64).map(|i| i * 2).sum::<u64>());
    h.shutdown();
}

#[test]
fn blocking_submit_does_not_hang() {
    use std::time::Duration;
    let h = Hydra::new();
    h.register("slow", |x: u64| {
        std::thread::sleep(Duration::from_millis(15));
        x + 1
    });
    h.init(local(24631)).unwrap();
    // .blocking(true) -> await_dispatch waits on done_cv, woken by dispatch_pass
    // on Pending->Running. Must return (not hang), then the job completes.
    let id = h.task("slow").blocking(true).submit(&41u64).unwrap();
    let v: u64 = h.get_blocking(id).unwrap();
    assert_eq!(v, 42);
    h.shutdown();
}
