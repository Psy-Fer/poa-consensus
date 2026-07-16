//! Regression family for seed-length sensitivity in periodic/homogeneous
//! tandem repeats: an auto-selected seed that is atypically short relative
//! to the true read population can make `heaviest_path` systematically
//! under-call a homogeneous repeat's true length, because the majority's
//! extra content (relative to the short seed) has nowhere unambiguous to
//! insert -- any position in a homogeneous run is an equally valid place,
//! so different reads scatter their inserts across different positions and
//! no single insertion position accumulates enough coverage to survive the
//! majority/coverage floor on its own.
//!
//! Confirmed real-data root cause and confirmed limitation of the existing
//! remedy: see CHANGELOG's "Seed-length sensitivity in periodic/homogeneous
//! tandem repeats" entry and CLAUDE.md's "Seed-Sensitivity Retry" section.
//! `bench/validate.py`'s tracked `cag50_d20` scenario is the clearest named
//! example -- its own best-scoring candidate among
//! `consensus_fit_scored`'s four remedies ties with pass-1, both wrong.
//!
//! This file exists because that single named scenario, and the two ad hoc
//! seed-sweep draws referenced in the CHANGELOG (`cag60_d20_s42`,
//! `cag50_d20_s2`), are one-off anecdotes: a batch generalization check this
//! session (25 fresh random draws per scenario class, `CAG` at 50 and 60
//! units and `GAA` at 55 units, depth 20, ONT R10 error model) found the
//! underlying baseline failure rate for this class is very high (~55-57%
//! of fresh draws), and the currently-shipped `consensus_fit_scored` remedy
//! fixes only ~29-31% of those failures. A single hardcoded scenario cannot
//! represent that -- these tests assert on an aggregate pass count across
//! multiple deterministic constructions per sub-case, documenting the
//! CURRENT shipped behavior honestly (including where it currently fails)
//! rather than a single pass/fail gate. No pre-clustering or other new
//! mechanism is exercised here; that is still under design review (see
//! `design/seed_sensitivity_preclustering.md`) and not wired into
//! `consensus_fit_scored`/`consensus_adaptive`.
//!
//! Three sub-cases, deliberately distinguished because they stress
//! meaningfully different points along the two confirmed-relevant axes
//! (how much of the population carries scattered noise, and whether that
//! noise clusters into sub-populations or spreads continuously):
//!
//! - **Few outliers**: only a small minority of reads carry any
//!   length-altering noise; the rest are byte-identical copies of the true
//!   sequence. Confirmed (see doc comment on the helper) to be the most
//!   benign shape -- a single outlier's own deletion/insertion gives the
//!   majority a single, unambiguous position to agree on.
//! - **Many scattered outliers**: most/all reads carry independent,
//!   differently-positioned noise (the `cag50_d20`-like continuous-spread
//!   shape) -- the shape that actually triggers vote fragmentation.
//! - **Bimodal outliers**: noise clusters into two sub-populations with a
//!   fixed length offset from each other (not a continuous spread), while
//!   still being simulated from one single true allele -- a distinct shape
//!   from continuous scatter, worth testing separately since a real
//!   bimodal split is also the shape a *genuine* second allele would
//!   produce, and this file is specifically about the single-allele path.
//!
//! No real patient data; deterministic (no external RNG), same discipline
//! as `tests/band_width_zero_unbanded.rs` and
//! `tests/multi_allele_periodic_phasing.rs`.

use poa_consensus::{AlignmentMode, PoaConfig, SeedSelection, consensus_fit_scored, select_seed};

const LEFT: &[u8] = b"ACGTACGTCGATCGATTAGCTAGCGCTAGCTA"; // 32bp unique flank
const RIGHT: &[u8] = b"ATCGATCGCGATCGATTAGCTAGCTGCATGCA"; // 32bp unique flank

fn make_true_read(true_units: usize, unit: &[u8]) -> Vec<u8> {
    let mut r = LEFT.to_vec();
    r.extend_from_slice(&unit.repeat(true_units));
    r.extend_from_slice(RIGHT);
    r
}

/// Deterministic scattered indel noise, same shape as
/// `tests/multi_allele_periodic_phasing.rs`'s `with_scattered_noise`:
/// inserts/deletes single bases at positions that vary by `read_idx` and a
/// cycling stride, so different reads disrupt different columns. `density`
/// controls what fraction of reads get any noise at all (1-in-`density`);
/// `density = 1` means every read gets noise.
fn with_scatter(
    seq: &[u8],
    read_idx: usize,
    region: std::ops::Range<usize>,
    density: usize,
) -> Vec<u8> {
    if read_idx % density != 0 {
        return seq.to_vec();
    }
    let span = region.end - region.start;
    let mut out = seq.to_vec();
    let pos = region.start + (read_idx * 7) % span;
    if read_idx % 2 == 0 {
        if pos < out.len() {
            out.remove(pos);
        }
    } else {
        out.insert(pos.min(out.len()), b'A');
    }
    let pos2 = region.start + (read_idx * 13 + 5) % span;
    if read_idx % 2 == 0 {
        out.insert(pos2.min(out.len()), b'A');
    } else if pos2 < out.len() {
        out.remove(pos2);
    }
    out
}

/// Fixed, deterministic bimodal offset: even-indexed reads get the true
/// sequence, odd-indexed reads are shifted `offset_units` shorter (a
/// systematic, not scattered, length difference) -- simulates two
/// consistent error sub-modes within what is still one true allele, as
/// distinct from `with_scatter`'s continuous per-read spread.
fn with_bimodal_shift(
    true_units: usize,
    unit: &[u8],
    read_idx: usize,
    offset_units: usize,
) -> Vec<u8> {
    let units = if read_idx % 2 == 0 {
        true_units
    } else {
        true_units - offset_units
    };
    make_true_read(units, unit)
}

fn cfg() -> PoaConfig {
    PoaConfig {
        band_width: 50,
        adaptive_band: true,
        alignment_mode: AlignmentMode::SemiGlobal,
        min_reads: 2,
        ..PoaConfig::default()
    }
}

fn count_repeat(seq: &[u8], unit: &[u8]) -> usize {
    let mut n = 0;
    let mut i = 0;
    while i + unit.len() <= seq.len() {
        if &seq[i..i + unit.len()] == unit {
            n += 1;
            i += unit.len();
        } else {
            i += 1;
        }
    }
    n
}

/// Runs `consensus_fit_scored` (the crate's own shipped seed-sensitivity
/// remedy, the same entry point this crate's CLI uses) and returns whether
/// the recovered unit count exactly matches `true_units`.
fn recovers_true_units(reads: &[Vec<u8>], true_units: usize, unit: &[u8]) -> bool {
    let slices: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let seed_idx = select_seed(&slices, &SeedSelection::Auto).unwrap();
    let result = consensus_fit_scored(&slices, seed_idx, &cfg()).unwrap();
    let units = count_repeat(&result.consensuses[0].sequence, unit);
    units == true_units
}

/// Sub-case 1: few outliers (~1-in-7 reads carry any noise). Swept across
/// several true unit counts and read-count/noise-stride variations to
/// avoid the "one lucky/unlucky draw" trap this session found repeatedly.
#[test]
fn seed_sensitivity_few_outliers_family() {
    let unit: &[u8] = b"CAG";
    let cases: [(usize, usize); 6] = [(50, 20), (55, 20), (60, 24), (50, 24), (45, 18), (58, 22)];
    let mut passed = 0;
    let mut results = Vec::new();
    for &(true_units, n_reads) in &cases {
        let true_read = make_true_read(true_units, unit);
        let region = LEFT.len()..(LEFT.len() + true_units * unit.len());
        let reads: Vec<Vec<u8>> = (0..n_reads)
            .map(|i| with_scatter(&true_read, i, region.clone(), 7))
            .collect();
        let ok = recovers_true_units(&reads, true_units, unit);
        results.push((true_units, n_reads, ok));
        if ok {
            passed += 1;
        }
    }
    // Documenting current behavior: the "few outliers" shape is the most
    // benign of the three sub-cases (a single outlier's deletion/insertion
    // gives the majority one unambiguous position to agree on, rather than
    // scattering across many). Confirmed empirically to pass all 6 cases at
    // the time this test was written -- if this regresses below all 6, that
    // is worth investigating as a real behavior change, not just noise.
    assert_eq!(
        passed,
        cases.len(),
        "few-outliers sub-case expected to be robust (all {} cases), got {}/{}: {:?}",
        cases.len(),
        passed,
        cases.len(),
        results
    );
}

/// Sub-case 2: many scattered outliers (every read carries independent
/// noise) -- the `cag50_d20`-like continuous-spread shape that actually
/// triggers vote fragmentation. Documenting current shipped behavior
/// honestly, not aspirationally: a 74-draw batch generalization check this
/// session (across `CAG` at 50/60 units and `GAA` at 55 units, depth 20,
/// ONT R10) found the baseline fails ~55-57% of fresh draws in this class,
/// and `consensus_fit_scored`'s existing four-candidate remedy fixes only
/// ~29-31% of those failures. This synthetic, deterministic reconstruction
/// of the same shape is, if anything, a harder stress test than that batch
/// (every read carries noise here, not just a random fraction) --
/// confirmed to fail all 8 constructed cases with the currently-shipped
/// remedy. `#[ignore]`d as a known, real, confirmed-fragile class rather
/// than weakened into a trivial always-passing assertion; the assertion
/// below documents the bar this sub-case SHOULD clear once a real fix (see
/// `design/seed_sensitivity_preclustering.md`, still under review) is
/// wired in, not the bar it currently clears (0/8).
#[test]
#[ignore = "known limitation: seed-length sensitivity under continuous \
            scattered-outlier noise is not fixed by the currently-shipped \
            consensus_fit_scored remedy (confirmed 0/8 on this synthetic \
            family, consistent with a 74-draw real-pipeline batch check \
            finding ~29-31% fix rate on failures of this general shape); \
            see design/seed_sensitivity_preclustering.md for the proposed, \
            not-yet-implemented fix"]
fn seed_sensitivity_many_scattered_outliers_family() {
    let unit: &[u8] = b"CAG";
    let cases: [(usize, usize); 8] = [
        (50, 20),
        (55, 20),
        (60, 24),
        (50, 24),
        (45, 18),
        (58, 22),
        (52, 21),
        (48, 19),
    ];
    let mut passed = 0;
    let mut results = Vec::new();
    for &(true_units, n_reads) in &cases {
        let true_read = make_true_read(true_units, unit);
        let region = LEFT.len()..(LEFT.len() + true_units * unit.len());
        let reads: Vec<Vec<u8>> = (0..n_reads)
            .map(|i| with_scatter(&true_read, i, region.clone(), 1))
            .collect();
        let ok = recovers_true_units(&reads, true_units, unit);
        results.push((true_units, n_reads, ok));
        if ok {
            passed += 1;
        }
    }
    // Aspirational bar (a real fix should clear at least half): currently
    // 0/8, hence #[ignore]d above rather than weakened to a trivial
    // always-true assertion. Un-ignore this once a real fix for the
    // seed-length-sensitivity class lands and re-confirm the actual rate.
    assert!(
        passed >= 4,
        "many-scattered-outliers sub-case: expected at least 4/{} to pass once a \
         real fix for seed-length sensitivity lands; got {passed}/{}: {:?} \
         (confirmed 0/8 with the currently-shipped remedy alone, which is why \
         this test is #[ignore]d)",
        cases.len(),
        cases.len(),
        results
    );
}

/// Sub-case 3: bimodal outliers -- a systematic length offset between two
/// read sub-populations (not continuous scatter), while still nominally one
/// true allele. Distinct failure shape from sub-case 2: worth testing
/// separately since a real bimodal split in read lengths is also what a
/// *genuine* second allele looks like, and this file is specifically about
/// the single-allele path's behavior on noise that merely resembles that.
#[test]
fn seed_sensitivity_bimodal_outliers_family() {
    let unit: &[u8] = b"CAG";
    // (true_units, n_reads, offset_units between the two sub-populations)
    let cases: [(usize, usize, usize); 6] = [
        (50, 20, 2),
        (55, 20, 3),
        (60, 24, 2),
        (50, 24, 4),
        (45, 18, 2),
        (58, 22, 3),
    ];
    let mut passed = 0;
    let mut results = Vec::new();
    for &(true_units, n_reads, offset) in &cases {
        let reads: Vec<Vec<u8>> = (0..n_reads)
            .map(|i| with_bimodal_shift(true_units, unit, i, offset))
            .collect();
        let ok = recovers_true_units(&reads, true_units, unit);
        results.push((true_units, n_reads, offset, ok));
        if ok {
            passed += 1;
        }
    }
    // Confirmed empirically to pass all 6 cases at the time this test was
    // written -- a systematic (not scattered) length offset between two
    // sub-populations turns out to be as benign as the few-outliers
    // sub-case, not as fragile as continuous scatter (sub-case 2). If this
    // regresses below all 6, that is worth investigating as a real behavior
    // change, not just noise -- see the aggregate-assertion rationale in
    // this file's module doc comment for why every sub-case here checks a
    // count, not a single hardcoded scenario.
    assert_eq!(
        passed,
        cases.len(),
        "bimodal-outliers sub-case expected to be robust (all {} cases), got {}/{}: {:?}",
        cases.len(),
        passed,
        cases.len(),
        results
    );
}
