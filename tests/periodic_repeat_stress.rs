//! Stress tests for periodic-repeat consensus — the regime the single-allele
//! diagonal-skip gating targets (see PoaConfig::multi_allele and
//! design/vntr_overcall_delete_edge_visibility.md).
//!
//! Reads are generated deterministically in-test (a fixed-seed LCG + a simple
//! ONT-like sub/ins/del error model) so there are no external-tool or data-file
//! dependencies and results are reproducible across platforms. Each test builds
//! a repeat of known unit count and asserts the consensus length lands within a
//! sub-unit tolerance of the truth — the property periodic over-calling breaks.
//!
//! These are safety tests: the single-allele default (diagonal-skip off) must
//! stay accurate across unit size, copy number, depth, and phase, and the
//! multi-allele path must still separate two length alleles. The real-data
//! over-call guard lives in tests/hg002_muc1_real_data.rs.

use poa_consensus::{PoaConfig, consensus, consensus_multi};

/// Tiny deterministic PRNG (SplitMix64) — reproducible, no external deps.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn frac(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn base(&mut self, exclude: u8) -> u8 {
        loop {
            let b = b"ACGT"[(self.next() % 4) as usize];
            if b != exclude {
                return b;
            }
        }
    }
}

const LEFT: &[u8] = b"AAAGACCCAACCCTATGACTTTAACTTCTTACAGCTACCACAGCCCCTAAACCCG";
const RIGHT: &[u8] = b"TCCATTCTCAATTCCCAGCCACCACTCTGATACTCCTACCACCCTTGCCAGCCAT";

fn truth(unit: &[u8], n: usize) -> Vec<u8> {
    let mut v = LEFT.to_vec();
    for _ in 0..n {
        v.extend_from_slice(unit);
    }
    v.extend_from_slice(RIGHT);
    v
}

/// Apply a sub/ins/del error model at rate `err`.
fn mutate(seq: &[u8], err: f64, rng: &mut Rng) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &c in seq {
        let r = rng.frac();
        if r < err * 0.45 {
            out.push(rng.base(c)); // substitution
        } else if r < err * 0.72 {
            // deletion: skip
        } else if r < err {
            out.push(rng.base(0)); // insertion
            out.push(c);
        } else {
            out.push(c);
        }
    }
    out
}

fn reads_for(unit: &[u8], n_units: usize, depth: usize, err: f64, seed: u64) -> Vec<Vec<u8>> {
    let t = truth(unit, n_units);
    let mut rng = Rng(seed);
    (0..depth).map(|_| mutate(&t, err, &mut rng)).collect()
}

fn consensus_len(reads: &[Vec<u8>]) -> usize {
    let slices: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
    // seed 0 is a real read; default config = single-allele (diag-skip off).
    consensus(&slices, 0, &PoaConfig::default())
        .expect("consensus")
        .sequence
        .len()
}

// A divergent 60 bp MUC1-like unit and a 3 bp unit, to cover both regimes.
const U60: &[u8] = b"CGGCTCCTGGTGCCCAGGCCCACCCTCCGCACCCTCACGTCCTGCTCCAGGCACCCACGT";
const U3: &[u8] = b"CAG";

/// Single-allele periodic consensus must track the true length across unit
/// count and depth (including the high-depth accumulation regime) — no
/// systematic over-call from phantom repeat units.
#[test]
fn periodic_single_allele_length_accurate() {
    // (unit, n_units, depth, err, seed, tolerance_bp)
    let cases: &[(&[u8], usize, usize, f64, u64, usize)] = &[
        (U60, 20, 20, 0.05, 1, 60),
        (U60, 65, 30, 0.05, 2, 90),
        (U60, 65, 80, 0.06, 3, 90), // high depth: accumulation regime
        (U60, 100, 40, 0.05, 4, 120),
        (U3, 40, 30, 0.05, 5, 6),
        (U3, 100, 30, 0.06, 6, 9),
    ];
    for &(unit, n, depth, err, seed, tol) in cases {
        let reads = reads_for(unit, n, depth, err, seed);
        let truth_bp = truth(unit, n).len();
        let len = consensus_len(&reads);
        assert!(
            len.abs_diff(truth_bp) <= tol,
            "unit={}bp n={n} depth={depth} err={err}: consensus {len} bp not within \
             {tol} bp of truth {truth_bp} bp (Δ={})",
            unit.len(),
            len as isize - truth_bp as isize,
        );
    }
}

/// Phase-shifted reads (each starting at a different offset within the repeat,
/// as happens without a flanking anchor) must not inflate the length. This is
/// the class the diagonal-skip greedy match used to over-call.
#[test]
fn periodic_phase_shifted_reads_not_overcalled() {
    let n_units = 60usize;
    let t = truth(U60, n_units);
    let mut rng = Rng(42);
    let mut reads = Vec::new();
    // Rotate the repeat body by a random number of units per read so reads
    // enter at different phases, then apply error.
    for _ in 0..30 {
        let shift = (rng.next() as usize % n_units) * U60.len();
        // keep flanks, rotate only the repeat body
        let body_start = LEFT.len();
        let body_end = t.len() - RIGHT.len();
        let body = &t[body_start..body_end];
        let rotated: Vec<u8> = body[shift..]
            .iter()
            .chain(&body[..shift])
            .copied()
            .collect();
        let mut r = LEFT.to_vec();
        r.extend_from_slice(&rotated);
        r.extend_from_slice(RIGHT);
        reads.push(mutate(&r, 0.05, &mut rng));
    }
    let truth_bp = t.len();
    let len = consensus_len(&reads);
    assert!(
        len <= truth_bp + 2 * U60.len(),
        "phase-shifted periodic reads over-called: {len} bp vs truth {truth_bp} bp",
    );
}

/// The multi-allele path must still SEPARATE two length alleles — the diag-skip
/// gating keeps diagonal-skip ON in multi-allele mode (see PoaConfig::multi_allele)
/// precisely so short-allele reads don't drift onto the long allele's extra
/// units. This guards that the single-allele gating did not collapse or shatter
/// allele separation.
///
/// Scope note: this asserts *separation* (two clearly length-distinct alleles),
/// NOT per-allele length precision. Exact periodic multi-allele length calling
/// is a separate, harder problem (cross-bubble phasing on a shared periodic
/// prefix; the short allele here is called a few units short) that single_off
/// does not touch — its accuracy is covered by bench/validate.py's `multi_*`
/// scenarios, not this unit test.
#[test]
fn multi_allele_two_lengths_separated() {
    let (short_n, long_n) = (30usize, 50usize);
    let mut rng = Rng(7);
    let mut reads = Vec::new();
    for _ in 0..20 {
        reads.push(mutate(&truth(U60, short_n), 0.05, &mut rng));
    }
    for _ in 0..20 {
        reads.push(mutate(&truth(U60, long_n), 0.05, &mut rng));
    }
    let slices: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
    let alleles = consensus_multi(&slices, 0, &PoaConfig::default()).expect("consensus_multi");
    assert_eq!(
        alleles.len(),
        2,
        "expected exactly two length alleles (not collapsed, not over-split), got {}",
        alleles.len()
    );
    let mut lens: Vec<usize> = alleles.iter().map(|a| a.sequence.len()).collect();
    lens.sort_unstable();
    // The two alleles must be clearly distinct in length (well over one unit
    // apart) — the separation property. (Truth gap is 20 units ≈ 1200 bp.)
    assert!(
        lens[1] - lens[0] >= 10 * U60.len(),
        "alleles not clearly length-distinct: {} vs {} bp",
        lens[0],
        lens[1],
    );
}
