use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_friedman_test::{Matrix, friedman};
use std::hint::black_box;

fn bench_friedman(c: &mut Criterion) {
    // 200k blocks × 5 treatments, small integer range to force within-block ties.
    let k = 5;
    let n = 200_000usize;
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut rows = Vec::with_capacity(n);
    for _ in 0..n {
        let mut row = Vec::with_capacity(k);
        for _ in 0..k {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            row.push((state % 7) as f64);
        }
        rows.push(row);
    }
    let m = Matrix { rows, k };

    c.bench_function("friedman_200k_5trt", |b| {
        b.iter(|| {
            let r = friedman(black_box(&m)).unwrap();
            black_box(r.q)
        });
    });
}

criterion_group!(benches, bench_friedman);
criterion_main!(benches);
