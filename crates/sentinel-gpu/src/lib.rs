//! # sentinel-gpu
//!
//! Compute acceleration for Sentinel's embarrassingly-parallel kernels, with an
//! **always-available CPU fallback**. The heavy GPU backend (`wgpu`) is gated
//! behind the `gpu` cargo feature, so the default build stays light and fast to
//! compile on constrained machines; the public API is identical either way.
//!
//! The accelerated kernels are the ones that dominate large runs:
//!
//! * [`pairwise_distances`] — all-vs-all (1 − cosine) over composition vectors
//!   (used by CHIMERA-SCAN and SAGE on large genomes);
//! * [`batch_gc`] — GC content of many sequences at once.
//!
//! Select a backend with [`Backend::best`] (GPU if compiled in and available,
//! else CPU) or force one explicitly.

use rayon::prelude::*;

/// Which compute backend to use.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Backend {
    /// Multi-threaded CPU (always available).
    Cpu,
    /// GPU via wgpu (only when built with the `gpu` feature **and** a device is
    /// present). Falls back to CPU automatically if unavailable.
    Gpu,
}

impl Backend {
    /// Pick the best available backend.
    pub fn best() -> Backend {
        #[cfg(feature = "gpu")]
        {
            if gpu::device_available() {
                return Backend::Gpu;
            }
        }
        Backend::Cpu
    }

    /// Human label including whether the GPU path is compiled in.
    pub fn describe(&self) -> String {
        match self {
            Backend::Cpu => "CPU (rayon)".to_string(),
            Backend::Gpu => {
                #[cfg(feature = "gpu")]
                {
                    "GPU (wgpu)".to_string()
                }
                #[cfg(not(feature = "gpu"))]
                {
                    "GPU requested but not compiled (rebuild with --features gpu); using CPU"
                        .to_string()
                }
            }
        }
    }
}

/// All-vs-all (1 − cosine) distance matrix over row vectors, flattened
/// row-major (`n * n`). Symmetric, zero diagonal.
pub fn pairwise_distances(vectors: &[Vec<f32>], backend: Backend) -> Vec<f32> {
    match backend {
        Backend::Gpu => {
            #[cfg(feature = "gpu")]
            {
                if let Some(out) = gpu::pairwise_distances(vectors) {
                    return out;
                }
            }
            pairwise_distances_cpu(vectors)
        }
        Backend::Cpu => pairwise_distances_cpu(vectors),
    }
}

fn pairwise_distances_cpu(vectors: &[Vec<f32>]) -> Vec<f32> {
    let n = vectors.len();
    // Precompute norms.
    let norms: Vec<f32> = vectors
        .iter()
        .map(|v| v.iter().map(|x| x * x).sum::<f32>().sqrt())
        .collect();
    let mut out = vec![0.0f32; n * n];
    // Parallelise over rows; write each row's upper triangle, mirror later.
    let rows: Vec<(usize, Vec<f32>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut row = vec![0.0f32; n];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dot: f32 = vectors[i]
                    .iter()
                    .zip(vectors[j].iter())
                    .map(|(a, b)| a * b)
                    .sum();
                let denom = norms[i] * norms[j];
                let cos = if denom > 1e-12 { dot / denom } else { 0.0 };
                row[j] = 1.0 - cos;
            }
            (i, row)
        })
        .collect();
    for (i, row) in rows {
        out[i * n..(i + 1) * n].copy_from_slice(&row);
    }
    out
}

/// GC content (fraction over A/C/G/T) for each sequence.
pub fn batch_gc(seqs: &[&[u8]], backend: Backend) -> Vec<f32> {
    // GC counting is memory-bound; the CPU path is optimal. The GPU branch is
    // kept for API symmetry and falls through to CPU.
    let _ = backend;
    seqs.par_iter()
        .map(|s| {
            let mut gc = 0u32;
            let mut at = 0u32;
            for &b in *s {
                match b {
                    b'G' | b'g' | b'C' | b'c' => gc += 1,
                    b'A' | b'a' | b'T' | b't' => at += 1,
                    _ => {}
                }
            }
            let tot = gc + at;
            if tot == 0 {
                0.0
            } else {
                gc as f32 / tot as f32
            }
        })
        .collect()
}

#[cfg(feature = "gpu")]
mod gpu {
    //! wgpu compute backend. Compiled only with `--features gpu`.
    use bytemuck::{Pod, Zeroable};

    #[repr(C)]
    #[derive(Clone, Copy, Pod, Zeroable)]
    struct Meta {
        n: u32,
        dim: u32,
    }

    /// Probe whether a GPU adapter can be acquired.
    pub fn device_available() -> bool {
        pollster::block_on(async {
            let instance = wgpu::Instance::default();
            instance
                .request_adapter(&wgpu::RequestAdapterOptions::default())
                .await
                .is_some()
        })
    }

    /// GPU all-vs-all (1 − cosine). Returns `None` on any device error so the
    /// caller can fall back to CPU.
    pub fn pairwise_distances(vectors: &[Vec<f32>]) -> Option<Vec<f32>> {
        if vectors.is_empty() {
            return Some(Vec::new());
        }
        let n = vectors.len();
        let dim = vectors[0].len();
        let mut flat = Vec::with_capacity(n * dim);
        for v in vectors {
            if v.len() != dim {
                return None;
            }
            flat.extend_from_slice(v);
        }
        pollster::block_on(run(&flat, n, dim))
    }

    const SHADER: &str = r#"
struct Meta { n: u32, dim: u32 };
@group(0) @binding(0) var<storage, read> data: array<f32>;
@group(0) @binding(1) var<storage, read_write> out: array<f32>;
@group(0) @binding(2) var<uniform> meta: Meta;

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x; let j = gid.y;
    if (i >= meta.n || j >= meta.n) { return; }
    var dot = 0.0; var na = 0.0; var nb = 0.0;
    for (var k = 0u; k < meta.dim; k = k + 1u) {
        let a = data[i * meta.dim + k];
        let b = data[j * meta.dim + k];
        dot = dot + a * b; na = na + a * a; nb = nb + b * b;
    }
    var cosv = 0.0;
    let denom = sqrt(na) * sqrt(nb);
    if (denom > 1e-12) { cosv = dot / denom; }
    out[i * meta.n + j] = 1.0 - cosv;
}
"#;

    async fn run(flat: &[f32], n: usize, dim: usize) -> Option<Vec<f32>> {
        use wgpu::util::DeviceExt;
        let instance = wgpu::Instance::default();
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions::default())
            .await?;
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default(), None)
            .await
            .ok()?;

        let data_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("data"),
            contents: bytemuck::cast_slice(flat),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let out_size = (n * n * std::mem::size_of::<f32>()) as u64;
        let out_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("out"),
            size: out_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let staging = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("staging"),
            size: out_size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let meta = Meta { n: n as u32, dim: dim as u32 };
        let meta_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("meta"),
            contents: bytemuck::bytes_of(&meta),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("cosine"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("pipeline"),
            layout: None,
            module: &module,
            entry_point: "main",
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: data_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: out_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: meta_buf.as_entire_binding() },
            ],
        });

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor::default());
            pass.set_pipeline(&pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            let groups = ((n as u32) + 7) / 8;
            pass.dispatch_workgroups(groups, groups, 1);
        }
        encoder.copy_buffer_to_buffer(&out_buf, 0, &staging, 0, out_size);
        queue.submit(Some(encoder.finish()));

        let slice = staging.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx.send(r); });
        device.poll(wgpu::Maintain::Wait);
        rx.recv().ok()?.ok()?;
        let data = slice.get_mapped_range();
        let result: Vec<f32> = bytemuck::cast_slice(&data).to_vec();
        drop(data);
        staging.unmap();
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpu_pairwise_is_symmetric_zero_diag() {
        let v = vec![
            vec![1.0f32, 0.0, 0.0],
            vec![0.0f32, 1.0, 0.0],
            vec![1.0f32, 1.0, 0.0],
        ];
        let d = pairwise_distances(&v, Backend::Cpu);
        let n = v.len();
        for i in 0..n {
            assert!((d[i * n + i]).abs() < 1e-6); // zero diagonal
            for j in 0..n {
                assert!((d[i * n + j] - d[j * n + i]).abs() < 1e-6); // symmetric
            }
        }
        // orthogonal vectors -> distance 1
        assert!((d[0 * n + 1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn batch_gc_correct() {
        let seqs: Vec<&[u8]> = vec![b"GGCC", b"ATAT", b"ACGT"];
        let g = batch_gc(&seqs, Backend::Cpu);
        assert!((g[0] - 1.0).abs() < 1e-6);
        assert!((g[1] - 0.0).abs() < 1e-6);
        assert!((g[2] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn best_backend_without_gpu_feature_is_cpu() {
        // With default features (no `gpu`), best() must be CPU.
        #[cfg(not(feature = "gpu"))]
        assert_eq!(Backend::best(), Backend::Cpu);
    }
}
