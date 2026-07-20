//! Correctness/discrimination coverage for the diagnostics & confidence layer
//! (the crate's key differentiator vs abPOA/SPOA, per the architecture
//! comparison). Complements tests/sv_analysis.rs (which covers
//! consensus_confidence, allele_fractions, credible_interval, max_accuracy,
//! has_competing_allele, should_call_multiallele, bridged_consensus, and many
//! SV bubble cases). This file targets the GAPS:
//!
//!   * GraphStats: assert each field *discriminates* clean vs two-allele vs
//!     noisy input, not merely that it is populated.
//!   * diagnose()/ConsensusWarnings: every signal fires when it should AND does
//!     NOT false-fire on clean input (both directions), plus is_clean().
//!   * consensus_fit: a correct candidate scores strictly better (lower) than a
//!     deliberately-wrong one against the same reads.
//!   * count_credible_interval / max_achievable_accuracy: numeric edge cases
//!     (empty, single obs, zero variance, wide variance, n=0, exact obs) with
//!     hand-computed expected values.
//!
//! All inputs are deterministic synthetic constructions (no real data).

use poa_consensus::analysis::{
    DiagnoseConfig, consensus_fit, count_credible_interval, diagnose, max_achievable_accuracy,
};
use poa_consensus::{
    BubbleSite, Consensus, CoverageGap, GapKind, GraphStats, PoaConfig, consensus,
};

// ── Deterministic helpers ────────────────────────────────────────────────────

fn xs(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

/// Copy `truth`, substituting each base with probability ~`pct`% (deterministic).
fn mutate_subs(truth: &[u8], pct: u64, seed: u64) -> Vec<u8> {
    let mut state = seed | 1;
    let bases = b"ACGT";
    truth
        .iter()
        .map(|&b| {
            if xs(&mut state) % 100 < pct {
                let mut nb = bases[(xs(&mut state) % 4) as usize];
                if nb == b {
                    nb = bases[(xs(&mut state) % 4) as usize];
                }
                nb
            } else {
                b
            }
        })
        .collect()
}

fn cfg() -> PoaConfig {
    PoaConfig {
        min_reads: 2,
        ..PoaConfig::default()
    }
}

/// A clean single-allele Consensus literal: every position fully supported,
/// no gaps, no bubbles. `median_input_read_len` == sequence length.
fn clean_consensus(seq: &[u8], n_reads: usize) -> Consensus {
    let len = seq.len();
    Consensus {
        sequence: seq.to_vec(),
        coverage: vec![n_reads as u32; len],
        path_weights: vec![n_reads as i32; len],
        n_reads,
        graph_stats: GraphStats {
            median_input_read_len: len,
            ..GraphStats::default()
        },
        gaps: vec![],
        bubble_sites: vec![],
        read_indices: vec![],
    }
}

// ── GraphStats discrimination ─────────────────────────────────────────────────

#[test]
fn graphstats_clean_unanimous_all_fields_pristine() {
    // 12 identical reads → linear chain, every stat at its "no disagreement" value.
    let truth: &[u8] = b"ACGTACGTGGCCAATTACGTACGT";
    let reads: Vec<&[u8]> = vec![truth; 12];
    let g = consensus(&reads, 0, &cfg()).unwrap().graph_stats;

    assert_eq!(g.bubble_count, 0, "no bubbles on unanimous input");
    assert_eq!(g.max_bubble_depth, 0);
    assert_eq!(g.longest_bubble_span, 0);
    assert_eq!(
        g.single_support_fraction, 0.0,
        "no singleton-supported nodes on unanimous input"
    );
    assert_eq!(g.coverage_variance, 0.0, "uniform coverage → zero variance");
    assert!(
        g.mean_column_entropy < 1e-9,
        "one base per column → ~zero entropy; got {}",
        g.mean_column_entropy
    );
    assert!(
        g.edge_weight_gini < 1e-9,
        "uniform edge weights → gini ~0; got {}",
        g.edge_weight_gini
    );
    assert!(
        (g.coverage_mean - 12.0).abs() < 1e-6,
        "every node covered by all 12 reads; got {}",
        g.coverage_mean
    );
}

#[test]
fn graphstats_two_allele_snp_raises_bubble_and_entropy() {
    // 6 reads allele A, 6 reads allele B differing by a single SNP (pos 10).
    let a: &[u8] = b"ACGTACGTACATACGTACGT";
    let b: &[u8] = b"ACGTACGTACGTACGTACGT"; // differs at index 10: A vs G
    let reads: Vec<&[u8]> = vec![a, a, a, a, a, a, b, b, b, b, b, b];
    let two = consensus(&reads, 0, &cfg()).unwrap().graph_stats;
    let clean = consensus(&[b; 12], 0, &cfg()).unwrap().graph_stats;

    assert!(
        two.bubble_count >= 1,
        "a real 6-vs-6 SNP split must register a bubble; got {}",
        two.bubble_count
    );
    assert!(
        two.max_bubble_depth >= 6,
        "minority arm carries 6 reads; got max_bubble_depth={}",
        two.max_bubble_depth
    );
    assert!(
        two.edge_weight_gini > clean.edge_weight_gini,
        "a bubble makes edge weights less uniform (gini up): {} vs {}",
        two.edge_weight_gini,
        clean.edge_weight_gini
    );
}

#[test]
fn graphstats_delete_disagreement_raises_column_entropy() {
    // mean_column_entropy measures per-node Match-vs-Delete disagreement (not
    // base substitution). Construct a single-base indel: 8 reads carry an
    // interior base, 4 reads delete it. That node has coverage 8 / delete 4,
    // binary entropy H(8/12) ≈ 0.918 bits, so the mean rises above clean.
    let left = b"ACGTGCATTGCAAGTC";
    let right = b"GTCATGCAAGTCATGC";
    let x = b"A"; // flanked by C (left end) and G (right start): unambiguous indel
    let present: Vec<u8> = [left.as_slice(), x, right.as_slice()].concat();
    let absent: Vec<u8> = [left.as_slice(), right.as_slice()].concat();

    let reads: Vec<&[u8]> = vec![
        &present, &present, &present, &present, &present, &present, &present, &present, &absent,
        &absent, &absent, &absent,
    ];
    let g = consensus(&reads, 0, &cfg()).unwrap().graph_stats;
    let clean = consensus(&[present.as_slice(); 12], 0, &cfg())
        .unwrap()
        .graph_stats;

    assert_eq!(
        clean.mean_column_entropy, 0.0,
        "unanimous → zero column entropy"
    );
    assert!(
        g.mean_column_entropy > 0.01,
        "an 8-vs-4 Match/Delete node must raise mean_column_entropy; got {}",
        g.mean_column_entropy
    );
}

#[test]
fn graphstats_length_allele_sets_longest_bubble_span() {
    // A length-differing (structural) two-allele set: a 9bp insertion arm.
    let short: &[u8] = b"ACGTACGTACGTACGTACGT";
    let long: &[u8] = b"ACGTACGTACAAAAAAAAAGTACGTACGTACGT";
    let reads: Vec<&[u8]> = vec![
        short, short, short, short, short, short, long, long, long, long, long, long,
    ];
    let g = consensus(&reads, 0, &cfg()).unwrap().graph_stats;
    let clean = consensus(&[short; 12], 0, &cfg()).unwrap().graph_stats;

    assert_eq!(
        clean.longest_bubble_span, 0,
        "no structural bubble when unanimous"
    );
    assert!(
        g.longest_bubble_span > 0,
        "a length-changing allele must set longest_bubble_span > 0; got {}",
        g.longest_bubble_span
    );
}

#[test]
fn graphstats_noisy_input_raises_single_support_and_variance() {
    // 15 reads = truth with ~6% per-base substitution noise. Scattered singleton
    // error nodes should push single_support_fraction and coverage_variance up.
    let truth: &[u8] = b"ACGTACGTGGCCAATTACGTACGTTTAACCGGACGTACGT";
    let mut reads_owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..15 {
        reads_owned.push(mutate_subs(truth, 6, 0x1234_5678 + i * 97));
    }
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let noisy = consensus(&reads, 0, &cfg()).unwrap().graph_stats;
    let clean = consensus(&vec![truth; 15], 0, &cfg()).unwrap().graph_stats;

    assert_eq!(clean.single_support_fraction, 0.0);
    assert!(
        noisy.single_support_fraction > 0.05,
        "noisy input should have many singleton-supported nodes; got {}",
        noisy.single_support_fraction
    );
    assert!(
        noisy.coverage_variance > clean.coverage_variance,
        "noisy coverage is uneven: variance {} should exceed clean {}",
        noisy.coverage_variance,
        clean.coverage_variance
    );
}

// ── diagnose() / ConsensusWarnings: both directions per signal ─────────────────

#[test]
fn diagnose_clean_consensus_is_clean() {
    let w = diagnose(
        &clean_consensus(b"ACGTACGTACGTACGT", 12),
        &DiagnoseConfig::default(),
    );
    assert!(
        w.is_clean(),
        "a fully-supported deep consensus must be clean: {:?}",
        w
    );
    assert!(w.low_depth.is_none());
    assert!(!w.has_coverage_gaps);
    assert!(w.interior_low_support.is_none());
    assert!(w.structural_competing.is_none());
    assert!(w.truncation_suspected.is_none());
}

#[test]
fn diagnose_low_depth_fires_below_threshold_only() {
    // Below the default warn threshold (10) and critical threshold (5).
    let shallow = clean_consensus(b"ACGTACGTACGTACGT", 3);
    let w = diagnose(&shallow, &DiagnoseConfig::default());
    let ld = w
        .low_depth
        .as_ref()
        .expect("3 reads is below the depth_warn_threshold of 10");
    assert_eq!(ld.n_reads, 3);
    assert!(
        ld.is_critical,
        "3 < depth_critical_threshold (5) → critical"
    );
    assert!(!w.is_clean());

    // At/above the threshold: no depth warning.
    let deep = clean_consensus(b"ACGTACGTACGTACGT", 12);
    assert!(
        diagnose(&deep, &DiagnoseConfig::default())
            .low_depth
            .is_none(),
        "12 reads is above threshold → no low_depth warning"
    );
}

#[test]
fn diagnose_coverage_gaps_fires_only_with_gaps() {
    let mut c = clean_consensus(b"ACGTACGTACGTACGT", 12);
    assert!(!diagnose(&c, &DiagnoseConfig::default()).has_coverage_gaps);

    c.gaps.push(CoverageGap {
        start: 8,
        end: 8,
        kind: GapKind::Unknown,
    });
    let w = diagnose(&c, &DiagnoseConfig::default());
    assert!(
        w.has_coverage_gaps,
        "a non-empty gaps vec must set has_coverage_gaps"
    );
    assert!(!w.is_clean());
}

#[test]
fn diagnose_interior_low_support_fires_on_weak_interior_only() {
    // 10 reads; one interior position supported by a single read (fraction 0.1
    // < default 0.15), the rest fully supported.
    let n = 10;
    let mut c = clean_consensus(b"ACGTACGTACGTACGT", n); // len 16
    let mid = 8;
    c.path_weights[mid] = 1; // fraction 1/10 = 0.10

    let w = diagnose(&c, &DiagnoseConfig::default());
    let iw = w
        .interior_low_support
        .expect("an interior position at fraction 0.10 must flag");
    assert_eq!(iw.position, mid);
    assert!(
        (iw.fraction - 0.10).abs() < 1e-6,
        "got fraction {}",
        iw.fraction
    );

    // The same weak position sitting in the boundary margin (default 20%) must
    // NOT flag: only the middle 60% is checked.
    let mut edge = clean_consensus(b"ACGTACGTACGTACGT", n);
    edge.path_weights[1] = 1; // within the leading 20% boundary margin
    assert!(
        diagnose(&edge, &DiagnoseConfig::default())
            .interior_low_support
            .is_none(),
        "a weak position inside the boundary margin should not flag as interior"
    );
}

#[test]
fn diagnose_structural_competing_fires_in_single_allele_mode_only() {
    let mut c = clean_consensus(b"ACGTACGTACGTACGT", 12);
    c.bubble_sites.push(BubbleSite {
        consensus_pos: 8,
        arm_read_counts: vec![8, 4],
        arm_sequences: vec![b"AAAA".to_vec(), b"".to_vec()],
        is_structural: true,
    });

    // Single-allele mode: the structural bubble is flagged.
    let single = diagnose(&c, &DiagnoseConfig::default());
    let sc = single
        .structural_competing
        .as_ref()
        .expect("a structural bubble must be flagged in single-allele mode");
    assert_eq!(sc.site_count, 1);
    assert_eq!(sc.min_arm_reads, 4, "weakest arm carries 4 reads");
    assert!(!single.is_clean());

    // Allele-partition mode suppresses it (already in multi mode).
    let allele_cfg = DiagnoseConfig {
        is_allele_partition: true,
        ..DiagnoseConfig::default()
    };
    assert!(
        diagnose(&c, &allele_cfg).structural_competing.is_none(),
        "structural_competing must be suppressed in allele-partition mode"
    );
}

#[test]
fn diagnose_truncation_fires_when_consensus_far_shorter_than_reads() {
    // Consensus 30bp, median input read 100bp → ratio 0.30 < default 0.60.
    let mut c = clean_consensus(&[b'A'; 30], 12);
    c.graph_stats.median_input_read_len = 100;
    let w = diagnose(&c, &DiagnoseConfig::default());
    let tw = w
        .truncation_suspected
        .expect("30/100 = 0.30 < 0.60 must flag truncation");
    assert_eq!(tw.consensus_len, 30);
    assert_eq!(tw.median_read_len, 100);
    assert!((tw.ratio - 0.30).abs() < 1e-6);

    // Consensus length ~ median read length → no truncation flag.
    let mut ok = clean_consensus(&[b'A'; 98], 12);
    ok.graph_stats.median_input_read_len = 100;
    assert!(
        diagnose(&ok, &DiagnoseConfig::default())
            .truncation_suspected
            .is_none(),
        "0.98 ratio should not flag truncation"
    );
}

#[test]
fn diagnose_real_shallow_consensus_flags_low_depth() {
    // End-to-end (not a literal): 3 real reads → low_depth fires via the real pipeline.
    let reads: Vec<&[u8]> = vec![b"ACGTACGTACGTACGT"; 3];
    let c = consensus(&reads, 0, &cfg()).unwrap();
    let w = diagnose(&c, &DiagnoseConfig::default());
    assert!(w.low_depth.is_some(), "3-read consensus must warn on depth");
}

// ── consensus_fit: correct candidate beats a wrong one ────────────────────────

#[test]
fn consensus_fit_ranks_correct_candidate_below_wrong() {
    let truth: &[u8] = b"ACGTACGTGGCCAATTACGTACGTTTAACCGG";
    let mut reads_owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..12 {
        reads_owned.push(mutate_subs(truth, 3, 0x00AB_CD00 + i * 131));
    }
    let reads: Vec<&[u8]> = reads_owned.iter().map(|r| r.as_slice()).collect();
    let c = cfg();

    // A wrong candidate: truth with a chunk deleted (length very different).
    let wrong: Vec<u8> = truth[..16].to_vec();
    // A different wrong candidate: an unrelated sequence of similar length.
    let unrelated: Vec<u8> = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA".to_vec();

    let fit_correct = consensus_fit(truth, &reads, &c);
    let fit_wrong = consensus_fit(&wrong, &reads, &c);
    let fit_unrelated = consensus_fit(&unrelated, &reads, &c);

    assert!(
        fit_correct < fit_wrong,
        "correct candidate should fit better (lower) than a truncated one: {} vs {}",
        fit_correct,
        fit_wrong
    );
    assert!(
        fit_correct < fit_unrelated,
        "correct candidate should fit better than an unrelated one: {} vs {}",
        fit_correct,
        fit_unrelated
    );
}

// ── count_credible_interval: numeric edge cases (hand-computed) ────────────────

#[test]
fn credible_interval_empty_is_nan() {
    let (lo, hi) = count_credible_interval(&[], 0.95);
    assert!(lo.is_nan() && hi.is_nan(), "empty input → (NaN, NaN)");
}

#[test]
fn credible_interval_single_observation_is_point() {
    let (lo, hi) = count_credible_interval(&[42.0], 0.95);
    assert_eq!(
        (lo, hi),
        (42.0, 42.0),
        "one observation → degenerate point interval"
    );
}

#[test]
fn credible_interval_zero_variance_is_point() {
    let (lo, hi) = count_credible_interval(&[5.0, 5.0, 5.0, 5.0], 0.95);
    assert_eq!(
        (lo, hi),
        (5.0, 5.0),
        "identical observations → zero-width interval at the mean"
    );
}

#[test]
fn credible_interval_hand_computed() {
    // values [10,12,14]: mean 12, sample var = (4+0+4)/2 = 4, std = 2.
    // 95% → z = 1.959964; half-width = z*std/sqrt(n) = 1.959964*2/sqrt(3) ≈ 2.2632.
    let (lo, hi) = count_credible_interval(&[10.0, 12.0, 14.0], 0.95);
    let mid = (lo + hi) / 2.0;
    assert!(
        (mid - 12.0).abs() < 1e-9,
        "interval must be centered on the mean; got {mid}"
    );
    let hw = (hi - lo) / 2.0;
    assert!(
        (hw - 2.2632).abs() < 1e-3,
        "half-width should match z*std/sqrt(n) ≈ 2.2632; got {hw:.4}"
    );
}

#[test]
fn credible_interval_wider_when_variance_is_wider() {
    let tight = count_credible_interval(&[9.0, 10.0, 11.0], 0.95);
    let wide = count_credible_interval(&[0.0, 10.0, 20.0], 0.95);
    let w_tight = tight.1 - tight.0;
    let w_wide = wide.1 - wide.0;
    assert!(
        w_wide > w_tight,
        "higher variance must yield a wider interval: {w_wide:.3} vs {w_tight:.3}"
    );
}

// ── max_achievable_accuracy: numeric edge cases (hand-computed) ────────────────

#[test]
fn max_accuracy_zero_observations_is_zero() {
    assert_eq!(
        max_achievable_accuracy(0, 1.37),
        0.0,
        "no observations → 0 accuracy"
    );
}

#[test]
fn max_accuracy_exact_observations_is_one() {
    assert_eq!(
        max_achievable_accuracy(20, 0.0),
        1.0,
        "sigma=0 (exact) → certainty"
    );
    assert_eq!(
        max_achievable_accuracy(5, -1.0),
        1.0,
        "non-positive sigma → certainty"
    );
}

#[test]
fn max_accuracy_monotonic_in_depth() {
    let a = max_achievable_accuracy(5, 1.37);
    let b = max_achievable_accuracy(20, 1.37);
    let d = max_achievable_accuracy(50, 1.37);
    assert!(
        a < b && b < d,
        "more observations → higher ceiling: {a:.3} < {b:.3} < {d:.3}"
    );
}

#[test]
fn max_accuracy_hand_computed_values() {
    // n=20, sigma=1.37: arg = sqrt(20)/(2*1.37) ≈ 1.6321; 2*Φ(1.6321)-1 ≈ 0.897.
    let acc = max_achievable_accuracy(20, 1.37);
    assert!((acc - 0.897).abs() < 0.01, "expected ≈0.897; got {acc:.4}");

    // Hard locus n=20, sigma=3.0: arg ≈ 0.7454; 2*Φ(0.7454)-1 ≈ 0.544.
    let hard = max_achievable_accuracy(20, 3.0);
    assert!(
        (hard - 0.544).abs() < 0.01,
        "expected ≈0.544; got {hard:.4}"
    );
    assert!(hard < acc, "a wider per-read sigma lowers the ceiling");
}
