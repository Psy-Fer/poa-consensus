/// Speed comparison harness — run with:
///   cargo test --release --test speed -- --nocapture
///
/// Each workload is a function that builds a PoaGraph from N reads and calls
/// consensus(). Wall time and peak heap allocation are measured across REPS
/// repetitions. No external deps.
use poa_consensus::{PoaConfig, PoaGraph};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

// ─── Tracking allocator ───────────────────────────────────────────────────────

struct TrackingAllocator;

/// Live bytes currently on the heap.
static LIVE: AtomicUsize = AtomicUsize::new(0);
/// Baseline live bytes at the most recent reset_peak() call.
static BASE: AtomicUsize = AtomicUsize::new(0);
/// Maximum live bytes seen since the most recent reset_peak() call.
static PEAK: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for TrackingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc(layout) };
        if !ptr.is_null() {
            let live = LIVE.fetch_add(layout.size(), Ordering::Relaxed) + layout.size();
            let mut peak = PEAK.load(Ordering::Relaxed);
            while live > peak {
                match PEAK.compare_exchange_weak(peak, live, Ordering::Relaxed, Ordering::Relaxed) {
                    Ok(_) => break,
                    Err(p) => peak = p,
                }
            }
        }
        ptr
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) };
        LIVE.fetch_sub(layout.size(), Ordering::Relaxed);
    }
}

#[global_allocator]
static ALLOCATOR: TrackingAllocator = TrackingAllocator;

/// Reset the peak counter to the current live baseline.
fn reset_peak() {
    let live = LIVE.load(Ordering::SeqCst);
    BASE.store(live, Ordering::SeqCst);
    PEAK.store(live, Ordering::SeqCst);
}

/// Peak heap increase above baseline since the last reset_peak() call (bytes).
fn peak_alloc() -> usize {
    PEAK.load(Ordering::SeqCst)
        .saturating_sub(BASE.load(Ordering::SeqCst))
}

// ─── Harness ─────────────────────────────────────────────────────────────────

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

/// Build a random sequence of `len` bp that contains no occurrence of `unit`.
/// Used for flanks so that unit-counting in the consensus is unambiguous.
fn make_flank(len: usize, unit: &[u8], seed: u64) -> Vec<u8> {
    let mut s = seed;
    let mut seq: Vec<u8> = (0..len).map(|_| random_base(&mut s)).collect();
    let k = unit.len();
    let mut i = 0;
    while i + k <= seq.len() {
        if &seq[i..i + k] == unit {
            let j = i + k - 1;
            let orig = seq[j];
            seq[j] = b"ACGT".iter().copied().find(|&b| b != orig).unwrap();
            i += k;
        } else {
            i += 1;
        }
    }
    seq
}

/// Assemble left_flank + repeat×count + right_flank into a single template.
fn flanked_template(unit: &[u8], count: usize, flank: usize, lseed: u64, rseed: u64) -> Vec<u8> {
    let left = make_flank(flank, unit, lseed);
    let right = make_flank(flank, unit, rseed);
    [left, unit.repeat(count), right].concat()
}

fn time_workload<F: Fn()>(label: &str, f: F) {
    f(); // warmup (also seeds allocator baseline)
    reset_peak();
    let t0 = Instant::now();
    for _ in 0..REPS {
        f();
    }
    let elapsed = t0.elapsed();
    let per_rep = elapsed / REPS;
    let peak_mb = peak_alloc() as f64 / (1024.0 * 1024.0);
    println!("  {label:<45} {per_rep:>8.2?} / rep   peak {peak_mb:>6.1} MB");
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

/// 15 reads × 5 kb: medium-long reads, adaptive band (w≈60).
fn workload_5k_reads() {
    let seq = make_seq(5000, 111);
    let cfg = PoaConfig {
        min_reads: 3,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..15).map(|i| mutate(&seq, 0.05, 8000 + i)).collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 10 reads × 10 kb: long reads, adaptive band (w≈110).
fn workload_10k_reads() {
    let seq = make_seq(10000, 222);
    let cfg = PoaConfig {
        min_reads: 3,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..10).map(|i| mutate(&seq, 0.05, 9000 + i)).collect();
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

// ─── Flanked repeat workloads ────────────────────────────────────────────────
//
// Same repeat as the bare variants above, but each read includes non-repetitive
// flanking sequence on both sides.  The flanks give the minimizer anchor engine
// unique k-mers to lock onto at the repeat boundaries.
//
// With 5% substitution error and k=15, ~46% of k-mers survive unmodified, giving
// ≈ total_flank_bp × 0.46 / 10 anchors.  The density gate threshold is 15:
//   100 bp/side (200 bp total) → ~9 anchors → gate rejects, no anchor benefit.
//   250 bp/side (500 bp total) → ~23 anchors → gate passes, anchors fire.

/// 50 reads × (100bp flank + CAG×40 + 100bp flank) — anchors below threshold.
fn workload_cag_flanked_100() {
    let template = flanked_template(b"CAG", 40, 100, 30_001, 40_001);
    let cfg = PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..50)
        .map(|i| mutate(&template, 0.05, 50_000 + i))
        .collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 50 reads × (250bp flank + CAG×40 + 250bp flank) — anchors above threshold.
fn workload_cag_flanked_250() {
    let template = flanked_template(b"CAG", 40, 250, 30_002, 40_002);
    let cfg = PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..50)
        .map(|i| mutate(&template, 0.05, 51_000 + i))
        .collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 50 reads × (100bp flank + AAGGG×30 + 100bp flank) — anchors below threshold.
fn workload_aaggg_flanked_100() {
    let template = flanked_template(b"AAGGG", 30, 100, 30_003, 40_003);
    let cfg = PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..50)
        .map(|i| mutate(&template, 0.05, 52_000 + i))
        .collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

/// 50 reads × (250bp flank + AAGGG×30 + 250bp flank) — anchors above threshold.
fn workload_aaggg_flanked_250() {
    let template = flanked_template(b"AAGGG", 30, 250, 30_004, 40_004);
    let cfg = PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    };
    let reads: Vec<Vec<u8>> = (0..50)
        .map(|i| mutate(&template, 0.05, 53_000 + i))
        .collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let _ = poa_consensus::consensus(&refs, 0, &cfg).unwrap();
}

// ─── Entry point ─────────────────────────────────────────────────────────────

#[test]
fn speed_comparison() {
    println!("\n=== POA speed workloads (release, {REPS} reps each) ===\n");
    println!(
        "  {:<45} {:>14}   {:>12}",
        "workload", "time/rep", "peak alloc"
    );
    println!("  {}", "-".repeat(77));
    time_workload("clean reads   (100r × 500bp, 5% err)", workload_clean_reads);
    time_workload("CAG×40        ( 50r × 120bp, 5% err)", workload_cag_repeat);
    time_workload(
        "AAGGG×30      ( 50r × 150bp, 5% err)",
        workload_aaggg_repeat,
    );
    time_workload("high depth    (200r × 200bp, 5% err)", workload_high_depth);
    time_workload("long reads    ( 30r × 600bp, 5% err)", workload_long_reads);
    time_workload("5k reads      ( 15r ×  5kbp, 5% err)", workload_5k_reads);
    time_workload("10k reads     ( 10r × 10kbp, 5% err)", workload_10k_reads);
    time_workload(
        "multi-allele  ( 20r × ~90bp, 5% err)",
        workload_multi_allele,
    );
    println!("  {}", "-".repeat(77));
    println!("  Flanked repeat variants  (anchors inactive at 100bp/side, active at 250bp/side)");
    println!("  {}", "-".repeat(77));
    time_workload(
        "CAG×40 +100bp flank (50r × 320bp, 5% err)",
        workload_cag_flanked_100,
    );
    time_workload(
        "CAG×40 +250bp flank (50r × 620bp, 5% err)",
        workload_cag_flanked_250,
    );
    time_workload(
        "AAGGG×30 +100bp flank (50r × 350bp, 5% err)",
        workload_aaggg_flanked_100,
    );
    time_workload(
        "AAGGG×30 +250bp flank (50r × 650bp, 5% err)",
        workload_aaggg_flanked_250,
    );
    println!();
}

// ═══ Performance matrix: isolate cost drivers (length × depth × band) ══════════
//
// The workloads above are fixed scenarios. This second harness varies ONE axis
// at a time so the scaling in each dimension is visible, reports median + spread
// (not a single noisy sample), and attributes cost to a phase (DP alignment vs
// graph construction vs consensus extraction = heaviest_path + coverage trim +
// interior filter). Run with:
//   cargo test --release --test speed perf_matrix -- --nocapture

/// depth reads of `len` bp, 5% substitution error, deterministic.
fn gen_reads(len: usize, depth: usize, seed: u64) -> Vec<Vec<u8>> {
    let base = make_seq(len, seed);
    (0..depth)
        .map(|i| mutate(&base, 0.05, seed.wrapping_add(1_000 + i as u64)))
        .collect()
}

fn cfg_adaptive() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    }
}

fn cfg_unbanded() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        band_width: 0,
        adaptive_band: false,
        warn_on_long_unbanded: false,
        ..Default::default()
    }
}

fn cfg_fixed(w: usize) -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        band_width: w,
        adaptive_band: false,
        ..Default::default()
    }
}

/// Run `f` `reps` times (after one warmup), returning (median, min, max) wall
/// time and the peak heap increase (MB) observed across the timed reps.
fn bench<F: Fn()>(reps: usize, f: F) -> (Duration, Duration, Duration, f64) {
    f(); // warmup + allocator baseline
    reset_peak();
    let mut times = Vec::with_capacity(reps);
    for _ in 0..reps {
        let t0 = Instant::now();
        f();
        times.push(t0.elapsed());
    }
    let peak_mb = peak_alloc() as f64 / (1024.0 * 1024.0);
    times.sort();
    (
        times[times.len() / 2],
        times[0],
        times[times.len() - 1],
        peak_mb,
    )
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}

/// One end-to-end consensus() measurement row (median ms, spread ms, peak MB).
fn row_consensus(label: &str, reads: &[Vec<u8>], cfg: &PoaConfig, reps: usize) {
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let (med, lo, hi, peak) = bench(reps, || {
        let _ = poa_consensus::consensus(&refs, 0, cfg).unwrap();
    });
    println!(
        "  {label:<40} {:>9.3} ms  [{:>7.3}–{:<7.3}]  {peak:>7.1} MB",
        ms(med),
        ms(lo),
        ms(hi),
    );
}

/// Phase attribution: total consensus() vs build (new+add_read) vs extract
/// (consensus() on a prebuilt graph) vs pure DP alignment (align_read_ops of
/// every read against the full graph). Build ≈ DP + graph mutation; extract =
/// heaviest_path + coverage trim + interior filter.
fn phase_breakdown(label: &str, reads: &[Vec<u8>], cfg: &PoaConfig, reps: usize) {
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();

    let (total, _, _, total_pk) = bench(reps, || {
        let _ = poa_consensus::consensus(&refs, 0, cfg).unwrap();
    });
    let (build, _, _, build_pk) = bench(reps, || {
        let mut g = PoaGraph::new(&reads[0], cfg.clone()).unwrap();
        for r in &reads[1..] {
            g.add_read(r).unwrap();
        }
        std::hint::black_box(&g);
    });

    // Build once (untimed) for extract + pure-align measurements.
    let mut g = PoaGraph::new(&reads[0], cfg.clone()).unwrap();
    for r in &reads[1..] {
        g.add_read(r).unwrap();
    }
    let (extract, _, _, extract_pk) = bench(reps, || {
        let _ = g.consensus().unwrap();
    });
    let (align, _, _, _) = bench(reps, || {
        for r in &refs {
            let _ = g.align_read_ops(r);
        }
    });

    println!("  {label}");
    println!(
        "     total consensus()          {:>9.3} ms   peak {:>7.1} MB",
        ms(total),
        total_pk
    );
    println!(
        "     build (new + add_read)     {:>9.3} ms   peak {:>7.1} MB   ({:.0}% of total)",
        ms(build),
        build_pk,
        100.0 * ms(build) / ms(total).max(1e-9)
    );
    println!(
        "       └ DP align (full graph)  {:>9.3} ms                    ({:.0}% of build)",
        ms(align),
        100.0 * ms(align) / ms(build).max(1e-9)
    );
    println!(
        "     extract (HB + trim + filt) {:>9.3} ms   peak {:>7.1} MB   ({:.0}% of total)",
        ms(extract),
        extract_pk,
        100.0 * ms(extract) / ms(total).max(1e-9)
    );
}

#[test]
fn perf_matrix() {
    println!("\n=== POA performance matrix (release; median over reps, [min–max]) ===");

    // ── Read-length scaling (depth 20, adaptive band) ─────────────────────────
    println!("\n  READ-LENGTH scaling  (depth 20, adaptive band; 5% error)");
    println!("  {}", "-".repeat(80));
    for &(len, reps) in &[(150usize, 15usize), (600, 9), (2000, 5), (5000, 3)] {
        let reads = gen_reads(len, 20, 0x100 + len as u64);
        row_consensus(
            &format!("len {len:>5}bp × 20 reads"),
            &reads,
            &cfg_adaptive(),
            reps,
        );
    }

    // ── Depth scaling (600bp, adaptive band) ──────────────────────────────────
    println!("\n  DEPTH scaling  (600bp reads, adaptive band; 5% error)");
    println!("  {}", "-".repeat(80));
    for &(depth, reps) in &[(10usize, 11usize), (20, 9), (50, 5)] {
        let reads = gen_reads(600, depth, 0x200 + depth as u64);
        row_consensus(
            &format!("600bp × {depth:>3} reads"),
            &reads,
            &cfg_adaptive(),
            reps,
        );
    }

    // ── Band-config comparison (600bp, depth 20) ──────────────────────────────
    println!("\n  BAND-CONFIG  (600bp × 20 reads; 5% error)");
    println!("  {}", "-".repeat(80));
    let r600 = gen_reads(600, 20, 0x300);
    row_consensus("adaptive (w=10+0.01L, floor 50)", &r600, &cfg_adaptive(), 9);
    row_consensus("fixed band w=200", &r600, &cfg_fixed(200), 9);
    row_consensus("unbanded (band_width=0)", &r600, &cfg_unbanded(), 9);

    // ── Band-config at larger size (2kb, depth 15) — where banding pays off ────
    println!("\n  BAND-CONFIG at 2kb  (2000bp × 15 reads; 5% error)");
    println!("  {}", "-".repeat(80));
    let r2k = gen_reads(2000, 15, 0x400);
    row_consensus("adaptive (w≈30)", &r2k, &cfg_adaptive(), 5);
    row_consensus("fixed band w=200", &r2k, &cfg_fixed(200), 5);
    row_consensus("unbanded (band_width=0)", &r2k, &cfg_unbanded(), 3);

    // ── Phase attribution ─────────────────────────────────────────────────────
    println!("\n  PHASE ATTRIBUTION  (adaptive band; where does time/memory go?)");
    println!("  {}", "-".repeat(80));
    phase_breakdown(
        "600bp × 20 reads",
        &gen_reads(600, 20, 0x500),
        &cfg_adaptive(),
        9,
    );
    phase_breakdown(
        "2000bp × 15 reads",
        &gen_reads(2000, 15, 0x600),
        &cfg_adaptive(),
        5,
    );
    phase_breakdown(
        "150bp × 50 reads (sub-kb STR regime)",
        &gen_reads(150, 50, 0x700),
        &cfg_adaptive(),
        11,
    );

    // Overall process peak RSS (monotonic high-water mark; heap peaks above are
    // per-workload and resettable, this is the whole-process ceiling).
    if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
        for line in status.lines() {
            if line.starts_with("VmHWM") {
                println!(
                    "\n  process peak RSS (VmHWM): {}",
                    line.trim_start_matches("VmHWM:").trim()
                );
            }
        }
    }
    println!();
}
