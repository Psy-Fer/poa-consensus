/// Speed comparison harness — run with:
///   cargo test --release --test speed -- --nocapture
///
/// Each workload is a function that builds a PoaGraph from N reads and calls
/// consensus(). Wall time is measured across REPS repetitions. No external deps.
use poa_consensus::{PoaConfig, PoaGraph};
use std::time::Instant;

const REPS: u32 = 5;

fn rng(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn random_base(state: &mut u64) -> u8 {
    b"ACGT"[(rng(state) % 4) as usize]
}

fn make_seq(len: usize, seed: u64) -> Vec<u8> {
    let mut s = seed;
    (0..len).map(|_| random_base(&mut s)).collect()
}

fn mutate(seq: &[u8], error_rate: f64, seed: u64) -> Vec<u8> {
    let mut s = seed;
    seq.iter()
        .map(|&b| {
            if (rng(&mut s) % 10000) < (error_rate * 10000.0) as u64 {
                random_base(&mut s)
            } else {
                b
            }
        })
        .collect()
}

fn repeat_seq(unit: &[u8], count: usize) -> Vec<u8> {
    unit.repeat(count)
}

fn time_workload<F: Fn() -> ()>(label: &str, f: F) {
    // warmup
    f();
    let t0 = Instant::now();
    for _ in 0..REPS {
        f();
    }
    let elapsed = t0.elapsed();
    let per_rep = elapsed / REPS;
    println!("  {label:<45} {per_rep:>8.2?} / rep");
}

// ─── Workloads ────────────────────────────────────────────────────────────────

/// 100 reads × 500 bp non-repetitive: diagonal skip should dominate.
fn workload_clean_reads() {
    let seq = make_seq(500, 42);
    let cfg = PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..100).map(|i| mutate(&seq, 0.05, 1000 + i)).collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 50 reads × CAG×40 (120 bp) repetitive: slide-and-lock / lookahead active.
fn workload_cag_repeat() {
    let seq = repeat_seq(b"CAG", 40);
    let cfg = PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..50).map(|i| mutate(&seq, 0.05, 2000 + i)).collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 50 reads × AAGGG×30 (150 bp): RFC1-like pentanucleotide repeat.
fn workload_aaggg_repeat() {
    let seq = repeat_seq(b"AAGGG", 30);
    let cfg = PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..50).map(|i| mutate(&seq, 0.05, 3000 + i)).collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 200 reads × 200 bp: stale-spine recompute benefit at high depth.
fn workload_high_depth() {
    let seq = make_seq(200, 99);
    let cfg = PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..200).map(|i| mutate(&seq, 0.05, 4000 + i)).collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 30 reads × 600 bp: longer reads, adaptive band, moderate depth.
fn workload_long_reads() {
    let seq = make_seq(600, 77);
    let cfg = PoaConfig {
        min_reads: 3,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..30).map(|i| mutate(&seq, 0.05, 5000 + i)).collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 20 reads × CAG×20 (60 bp), two alleles balanced: multi-allele path.
fn workload_multi_allele() {
    let short: Vec<u8> = [b"ACGTACGT".as_slice(), &repeat_seq(b"CAG", 10), b"TTTTGGGG"].concat();
    let long: Vec<u8> = [b"ACGTACGT".as_slice(), &repeat_seq(b"CAG", 20), b"TTTTGGGG"].concat();
    let cfg = PoaConfig {
        min_reads: 3,
        min_allele_freq: 0.2,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let mut reads: Vec<Vec<u8>> = (0..10).map(|i| mutate(&short, 0.05, 6000 + i)).collect();
    reads.extend((0..10).map(|i| mutate(&long, 0.05, 7000 + i)));
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus_multi(&refs, 0, &cfg).unwrap();
}

// ─── Entry point ─────────────────────────────────────────────────────────────

#[test]
fn speed_comparison() {
    println!("\n=== POA speed workloads (release, {REPS} reps each) ===\n");
    time_workload("clean reads   (100r × 500bp, 5% err)", workload_clean_reads);
    time_workload("CAG×40        ( 50r × 120bp, 5% err)", workload_cag_repeat);
    time_workload(
        "AAGGG×30      ( 50r × 150bp, 5% err)",
        workload_aaggg_repeat,
    );
    time_workload("high depth    (200r × 200bp, 5% err)", workload_high_depth);
    time_workload("long reads    ( 30r × 600bp, 5% err)", workload_long_reads);
    time_workload(
        "multi-allele  ( 20r × ~90bp, 5% err)",
        workload_multi_allele,
    );
    println!();
}
