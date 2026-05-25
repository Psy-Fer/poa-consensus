/// Flank-size sweep: measure how non-repetitive flanking context affects
/// consensus accuracy and alignment time across repeat sizes and repeat-count
/// heterogeneity (unstable repeats where the number of units varies per read).
///
/// Key questions:
///   - At what flank length do minimizer anchors fire (density gate ≥ 15)?
///   - Does flank context improve accuracy for heterogeneous repeats?
///   - Where does performance improve vs just paying for longer reads?
///
/// Run with:
///   cargo test --release --test flank_sweep -- --nocapture --include-ignored
///
/// Error model: 4% substitution + 1% insertion + 1% deletion (approx ONT R10).
/// Heterogeneity sigma: per-read repeat count drawn from Normal(modal, σ).
/// Accuracy is measured over ACC_REPS seeds and reported as the median Δunits.
use poa_consensus::PoaConfig;
use std::time::Instant;

const TIME_REPS: u32 = 3; // timing reps (sigma=0 reads)
const ACC_REPS: u32 = 5;  // accuracy reps (het sigma reads); median is reported
const N_READS: usize = 50;
const SUB: f64 = 0.04;
const INS: f64 = 0.01;
const DEL: f64 = 0.01;

// ─── RNG ─────────────────────────────────────────────────────────────────────

fn rng(s: &mut u64) -> u64 {
    *s ^= *s << 13;
    *s ^= *s >> 7;
    *s ^= *s << 17;
    *s
}

fn rand_f64(s: &mut u64) -> f64 {
    rng(s) as f64 / u64::MAX as f64
}

fn random_base(s: &mut u64) -> u8 {
    b"ACGT"[(rng(s) % 4) as usize]
}

/// Normal(0, 1) approximation via sum of 12 Uniform[0,1) minus 6 (CLT).
fn normal01(s: &mut u64) -> f64 {
    let sum: f64 = (0..12)
        .map(|_| (rng(s) % 1_000_000) as f64 / 1_000_000.0)
        .sum();
    sum - 6.0
}

// ─── Error simulator ─────────────────────────────────────────────────────────

fn simulate(template: &[u8], sub: f64, ins: f64, del: f64, seed: u64) -> Vec<u8> {
    let mut s = seed;
    let mut out = Vec::with_capacity(template.len());
    for &base in template {
        if rand_f64(&mut s) < ins {
            out.push(random_base(&mut s));
        }
        if rand_f64(&mut s) < del {
            continue;
        }
        if rand_f64(&mut s) < sub {
            out.push(random_base(&mut s));
        } else {
            out.push(base);
        }
    }
    out
}

// ─── Sequence helpers ─────────────────────────────────────────────────────────

/// Random `len`-bp sequence guaranteed to contain no occurrence of `unit`.
/// Flanks built this way let `count_units` measure the consensus unambiguously.
fn make_flank(len: usize, unit: &[u8], seed: u64) -> Vec<u8> {
    if len == 0 {
        return vec![];
    }
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

/// Count non-overlapping occurrences of `unit` in `seq`.
fn count_units(seq: &[u8], unit: &[u8]) -> usize {
    if unit.is_empty() {
        return 0;
    }
    let k = unit.len();
    let mut n = 0;
    let mut i = 0;
    while i + k <= seq.len() {
        if &seq[i..i + k] == unit {
            n += 1;
            i += k;
        } else {
            i += 1;
        }
    }
    n
}

// ─── Read generator ───────────────────────────────────────────────────────────

/// Build `n_reads` simulated reads: left_flank + repeat×count + right_flank.
///
/// When `het_sigma > 0` each read's repeat count is drawn independently from
/// Normal(modal_count, het_sigma), rounded and clamped to [1, 3×modal].  This
/// models a locus with somatic instability or genuine allele-length variance.
///
/// Flanks are fixed per call (same reference sequence for all reads).  Error
/// simulation is applied per-read with distinct seeds.
fn make_reads(
    unit: &[u8],
    modal_count: usize,
    flank_len: usize,
    sub: f64,
    ins: f64,
    del: f64,
    het_sigma: f64,
    base_seed: u64,
) -> Vec<Vec<u8>> {
    let left = make_flank(flank_len, unit, base_seed.wrapping_add(1_000_000));
    let right = make_flank(flank_len, unit, base_seed.wrapping_add(2_000_000));
    (0..N_READS)
        .map(|i| {
            let mut rs = base_seed.wrapping_add(i as u64 * 37 + 13);
            let count = if het_sigma > 0.0 {
                let offset = (normal01(&mut rs) * het_sigma).round() as isize;
                (modal_count as isize + offset).max(1).min(modal_count as isize * 3) as usize
            } else {
                modal_count
            };
            let repeat = unit.repeat(count);
            let template: Vec<u8> = left
                .iter()
                .chain(repeat.iter())
                .chain(right.iter())
                .cloned()
                .collect();
            simulate(&template, sub, ins, del, base_seed.wrapping_add(i as u64 * 7 + 31))
        })
        .collect()
}

// ─── Cell runners ─────────────────────────────────────────────────────────────

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 3,
        band_width: 50,
        adaptive_band: true,
        ..Default::default()
    }
}

/// Mean alignment time using uniform reads (sigma=0) over `TIME_REPS` seeds.
fn time_cell(unit: &[u8], modal: usize, flank: usize, base_seed: u64) -> std::time::Duration {
    let c = cfg();
    let mut total = std::time::Duration::ZERO;
    for rep in 0..TIME_REPS {
        let reads =
            make_reads(unit, modal, flank, SUB, INS, DEL, 0.0, base_seed + rep as u64 * 1_000);
        let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
        let t0 = Instant::now();
        let _ = poa_consensus::consensus(&refs, 0, &c).unwrap();
        total += t0.elapsed();
    }
    total / TIME_REPS
}

/// Median Δunits (consensus_count - modal) over `ACC_REPS` seeds.
fn accuracy_cell(
    unit: &[u8],
    modal: usize,
    flank: usize,
    het_sigma: f64,
    base_seed: u64,
) -> isize {
    let c = cfg();
    let mut deltas: Vec<isize> = (0..ACC_REPS)
        .map(|rep| {
            let reads = make_reads(
                unit, modal, flank, SUB, INS, DEL, het_sigma,
                base_seed + rep as u64 * 1_000,
            );
            let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
            let result = poa_consensus::consensus(&refs, 0, &c).unwrap();
            count_units(&result.sequence, unit) as isize - modal as isize
        })
        .collect();
    deltas.sort();
    deltas[ACC_REPS as usize / 2] // median
}

// ─── Sweep ───────────────────────────────────────────────────────────────────

struct Case {
    unit: &'static [u8],
    modal: usize,
    label: &'static str,
}

const CASES: &[Case] = &[
    Case { unit: b"CAG",   modal: 20,  label: "CAG×20  ( 60 bp repeat)" },
    Case { unit: b"CAG",   modal: 40,  label: "CAG×40  (120 bp repeat)" },
    Case { unit: b"CAG",   modal: 100, label: "CAG×100 (300 bp repeat)" },
    Case { unit: b"AAGGG", modal: 30,  label: "AAGGG×30 (150 bp repeat)" },
];

const FLANK_SIZES: &[usize] = &[0, 50, 100, 150, 200, 300];
const HET_SIGMAS: &[f64] = &[0.0, 2.0, 5.0];

/// True when the density gate passes: total flank bp × 0.46 / MINI_W ≥ 15.
/// Approximately: both flanks together need ≥ 326 bp of non-repetitive context.
fn anchor_fires(flank_per_side: usize) -> bool {
    (2 * flank_per_side) as f64 * 0.46 / 10.0 >= 15.0
}

#[test]
#[ignore] // run explicitly: cargo test --release --test flank_sweep -- --nocapture --include-ignored
fn flank_sweep() {
    println!(
        "\n=== Flank size sweep ({N_READS} reads/rep, {:.0}% sub {:.0}% ins {:.0}% del) ===",
        SUB * 100.0,
        INS * 100.0,
        DEL * 100.0,
    );
    println!(
        "  timing: {TIME_REPS} reps (σ=0)  |  accuracy: {ACC_REPS} reps, median Δunits reported"
    );
    println!("  ⚓ = minimizer anchors active (density gate ≥ 15 anchors)\n");

    for (ci, case) in CASES.iter().enumerate() {
        let repeat_bp = case.modal * case.unit.len();
        println!("  ── {} ──────────────────────────────", case.label);
        // Header
        let sigma_headers: String = HET_SIGMAS
            .iter()
            .map(|&s| format!("  {:>4}", format!("σ={:.0}", s)))
            .collect();
        println!("  {:>10}  {:>8}{}", "flank/side", "time", sigma_headers);
        println!("  {:>10}  {:>8}{}", "──────────", "────────",
            HET_SIGMAS.iter().map(|_| "  ────").collect::<String>());

        for (fi, &flank) in FLANK_SIZES.iter().enumerate() {
            let anchor_marker = if anchor_fires(flank) { "⚓" } else { "  " };
            let read_bp = 2 * flank + repeat_bp;
            // Use a seed that varies per (case, flank) to avoid systematic bad seeds.
            let cell_seed = 100_000 + ci as u64 * 10_000 + fi as u64 * 1_000;

            let time = time_cell(case.unit, case.modal, flank, cell_seed);
            let deltas: Vec<isize> = HET_SIGMAS
                .iter()
                .enumerate()
                .map(|(si, &sig)| {
                    accuracy_cell(
                        case.unit, case.modal, flank, sig,
                        cell_seed + si as u64 * 100_000,
                    )
                })
                .collect();

            let delta_cols: String = deltas
                .iter()
                .map(|&d| format!("  {:>+4}", d))
                .collect();
            println!(
                "  {:>7} bp{}  {:>6.1?}{}   ({read_bp} bp/read)",
                flank, anchor_marker, time, delta_cols,
            );
        }
        println!();
    }
}
