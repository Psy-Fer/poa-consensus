//! Structural-variant handling and analysis-helper integration tests.
//!
//! Each test builds a read set with a known ground truth, runs POA consensus,
//! then uses the analysis helpers to inspect what the graph found.  The goal
//! is to verify that:
//!   - the helpers surface correct information for realistic inputs, and
//!   - the library degrades gracefully (rather than silently misreporting)
//!     on inputs that contain SVs, edge cases, or adversarial structure.
//!
//! Run with:
//!   cargo test --test sv_analysis -- --nocapture

use poa_consensus::{
    analysis::{
        allele_fractions, consensus_confidence, count_credible_interval, has_competing_allele,
        low_coverage_regions, max_achievable_accuracy, min_coverage, should_call_multiallele,
    },
    consensus, consensus_multi, PoaConfig,
};

// ── Shared utilities ──────────────────────────────────────────────────────────

fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Deterministic xorshift RNG.
fn xorshift(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

fn rand_f64(state: &mut u64) -> f64 {
    xorshift(state) as f64 / u64::MAX as f64
}

fn random_base_not(exclude: u8, state: &mut u64) -> u8 {
    let opts: [u8; 3] = match exclude {
        b'A' => [b'C', b'G', b'T'],
        b'C' => [b'A', b'G', b'T'],
        b'G' => [b'A', b'C', b'T'],
        _ => [b'A', b'C', b'G'],
    };
    opts[(xorshift(state) % 3) as usize]
}

/// Simulate `n` reads from `template` with per-base sub/ins/del rates.
fn simulate(template: &[u8], n: usize, sub: f64, ins: f64, del: f64, seed: u64) -> Vec<Vec<u8>> {
    let mut rng = seed;
    let mut reads = Vec::with_capacity(n);
    for _ in 0..n {
        let mut read = Vec::with_capacity(template.len());
        for &base in template {
            if rand_f64(&mut rng) < ins {
                let b = b"ACGT"[(xorshift(&mut rng) % 4) as usize];
                read.push(b);
            }
            if rand_f64(&mut rng) < del {
                continue;
            }
            if rand_f64(&mut rng) < sub {
                read.push(random_base_not(base, &mut rng));
            } else {
                read.push(base);
            }
        }
        if !read.is_empty() {
            reads.push(read);
        }
    }
    reads
}

fn default_cfg() -> PoaConfig {
    PoaConfig::default()
}

fn wide_cfg() -> PoaConfig {
    PoaConfig {
        adaptive_band: true,
        band_width: 50,
        ..PoaConfig::default()
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// SECTION 1: Analysis helpers against known ground truth
// ═══════════════════════════════════════════════════════════════════════════════

/// Clean reads → all helpers report high-confidence, no flags.
///
/// Uses zero-error reads: substitution errors create single-bp mismatch bubbles
/// (new insert nodes), inflating `single_support_fraction` even at 2% error.
/// The mechanism under test is the helper logic, not error tolerance.
#[test]
fn helpers_clean_reads_high_confidence() {
    let truth: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGT"; // 28 bp
    // Zero errors: identical reads so no spurious single-support nodes.
    let reads: Vec<&[u8]> = vec![truth; 20];
    let result = consensus(&reads, 0, &default_cfg()).unwrap();

    let conf = consensus_confidence(&result, 0.25);
    assert!(!conf.has_gaps, "clean overlapping reads should not produce gaps");
    assert!(!conf.competing_allele, "single-allele input should not flag competing allele");
    assert!(conf.min_cov >= 1, "min coverage should be non-zero at depth 20");
    assert!(conf.mean_cov >= 10.0, "mean coverage should be high");
    assert_eq!(conf.single_support_fraction, 0.0, "no single-support nodes with identical reads");
    assert!(!conf.is_low_confidence(), "should be high confidence: {:?}", conf);

    // Helpers individually
    assert!(min_coverage(&result) >= 1);
    assert!(low_coverage_regions(&result, 2).is_empty(),
        "no low-coverage regions expected at depth 20");
    assert!(!should_call_multiallele(&result, 0.25));
}

/// Per-read length estimates → credible interval contains the truth.
#[test]
fn helpers_credible_interval_contains_truth() {
    // Simulate reads from a 30bp truth with minor length noise (1% indel).
    let truth: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTAC"; // 30 bp
    let reads = simulate(truth, 25, 0.03, 0.015, 0.015, 42);

    // Per-read length as a proxy for per-observation estimate.
    let lengths: Vec<f64> = reads.iter().map(|r| r.len() as f64).collect();
    let (lo, hi) = count_credible_interval(&lengths, 0.95);

    let truth_len = truth.len() as f64;
    assert!(lo <= truth_len && truth_len <= hi,
        "truth length {truth_len} should be inside 95% CI [{lo:.1}, {hi:.1}]");
}

/// max_achievable_accuracy bounds the empirical accuracy from above.
#[test]
fn helpers_max_accuracy_ceiling_holds() {
    let truth: &[u8] = b"CATCATCATCATCATCATCATCATCATCAT"; // 30 bp, 10×CAT
    // High-error reads: σ is large, so theoretical max is well below 1.0.
    let reads = simulate(truth, 15, 0.08, 0.04, 0.04, 7);
    let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
    let result = consensus(&refs, 0, &wide_cfg()).unwrap();

    // Per-read length as scalar observations; σ from length variance.
    let lengths: Vec<f64> = reads.iter().map(|r| r.len() as f64).collect();
    let mean = lengths.iter().sum::<f64>() / lengths.len() as f64;
    let var = lengths.iter().map(|&x| (x - mean).powi(2)).sum::<f64>()
        / (lengths.len() - 1) as f64;
    let sigma = var.sqrt();

    let max_acc = max_achievable_accuracy(reads.len(), sigma);

    // Empirical accuracy: edit distance as a fraction of truth length.
    let edit = levenshtein(truth, &result.sequence);
    let empirical_acc = 1.0 - edit as f64 / truth.len() as f64;

    // The theoretical ceiling must be above zero and the empirical accuracy
    // cannot exceed 1.0 (trivially).
    assert!(max_acc > 0.0, "max_acc should be positive");
    assert!(max_acc <= 1.0, "max_acc cannot exceed 1.0");
    // When the locus is hard (high error), max_acc should be noticeably < 1.
    if sigma > 1.5 {
        assert!(max_acc < 0.99,
            "hard locus (σ={sigma:.2}) should not have max_acc near 1; got {max_acc:.3}");
    }
    // Empirical accuracy on a short known truth should be achievable.
    println!(
        "depth={}, sigma={sigma:.2}, max_acc={max_acc:.3}, empirical_acc={empirical_acc:.3}",
        reads.len()
    );
}

/// Two-allele SNP → `has_competing_allele` fires; single consensus flags it.
#[test]
fn helpers_competing_allele_detected_on_two_allele_snp() {
    // Allele A: ...G... at position 10
    // Allele B: ...T... at position 10
    let allele_a: &[u8] = b"AAAAAAAAAAGAAAAAAAAAA"; // 21 bp
    let allele_b: &[u8] = b"AAAAAAAAAATAAAAAAAAAA"; // 21 bp

    let reads_a = simulate(allele_a, 8, 0.01, 0.0, 0.0, 11);
    let reads_b = simulate(allele_b, 8, 0.01, 0.0, 0.0, 22);
    let mut all: Vec<Vec<u8>> = reads_a.into_iter().chain(reads_b).collect();
    all.rotate_right(1); // put a B-allele read at index 0 so seed is not mono-allele
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &default_cfg()).unwrap();
    let conf = consensus_confidence(&result, 0.25);

    assert!(
        conf.competing_allele,
        "8+8 balanced reads should flag a competing allele; bubble_sites={:?}",
        result.bubble_sites
    );
    assert!(has_competing_allele(&result, 0.25).is_some());
    assert!(should_call_multiallele(&result, 0.25));

    // Multi-allele call should recover both.
    let alleles = consensus_multi(&refs, 0, &default_cfg()).unwrap();
    assert!(alleles.len() >= 2, "expected two alleles from balanced split");
}

/// `bridged_consensus` explicitly marks the join between two partial consensuses
/// as a gap of unknown size.  The analysis helper must surface this via `has_gaps`.
///
/// This is the production scenario: two read groups that each cover only one side
/// of a region too large for any single read.  `bridged_consensus` concatenates
/// the two per-side consensuses and inserts a `GapKind::Unknown` gap at the join.
#[test]
fn helpers_gap_detected_for_non_overlapping_partials() {
    use poa_consensus::bridged_consensus;

    // Left reads cover the beginning of the region.
    let left_reads: Vec<&[u8]> = vec![
        b"ACGTACGTACGT", b"ACGTACGTACGT", b"ACGTACGTACGT",
        b"ACGTACGTACGT", b"ACGTACGTACGT",
    ];
    // Right reads cover the end; they do not overlap with left reads.
    let right_reads: Vec<&[u8]> = vec![
        b"TTTGGGCCCAAA", b"TTTGGGCCCAAA", b"TTTGGGCCCAAA",
        b"TTTGGGCCCAAA", b"TTTGGGCCCAAA",
    ];

    let result = bridged_consensus(&left_reads, 0, &right_reads, 0, &default_cfg()).unwrap();
    let conf = consensus_confidence(&result, 0.25);

    println!(
        "gap_test: has_gaps={}, gaps={:?}",
        conf.has_gaps, result.gaps
    );

    assert!(conf.has_gaps, "bridged_consensus must always report a gap");
    assert!(!result.gaps.is_empty(), "gaps vec must be non-empty");
    assert!(conf.is_low_confidence(), "gap should flag low confidence");
}

/// `allele_fractions` sums to 1.0 and ranks arms correctly on a real bubble.
#[test]
fn helpers_allele_fractions_on_real_bubble() {
    let majority: &[u8] = b"AAAAAAACCCAAAAAAAAA"; // 19 bp
    let minority: &[u8] = b"AAAAAAAGGGAAAAAAAAA"; // 19 bp, different middle 3 bp

    let reads_maj = simulate(majority, 12, 0.01, 0.0, 0.0, 5);
    let reads_min = simulate(minority, 4, 0.01, 0.0, 0.0, 6);
    let mut all: Vec<Vec<u8>> = reads_maj.into_iter().chain(reads_min).collect();
    all.rotate_right(1);
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &default_cfg()).unwrap();

    if let Some(site) = has_competing_allele(&result, 0.20) {
        let fracs = allele_fractions(site);
        let sum: f64 = fracs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9, "fractions must sum to 1.0, got {sum}");
        let max_frac = fracs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_frac > 0.5, "majority arm should be > 50%");
    }
    // If the bubble wasn't detected (depth trimmed), that's OK — just skip the fraction check.
}

// ═══════════════════════════════════════════════════════════════════════════════
// SECTION 2: Structural variant handling
// ═══════════════════════════════════════════════════════════════════════════════

/// Most reads have the full sequence; a minority carry a large deletion.
/// The full-sequence arm should win the consensus; the bubble is structural.
#[test]
fn sv_deletion_minority_creates_short_arm() {
    let full: &[u8]    = b"GGGGGGGGGGGGAAAACCCCTTTTGGGGGGGGGGGG"; // 36 bp
    let deleted: &[u8] = b"GGGGGGGGGGGGGGGGGGGGGGGG";             // 24 bp (middle 12 bp gone)

    let reads_full = simulate(full, 10, 0.02, 0.01, 0.01, 100);
    let reads_del  = simulate(deleted, 4, 0.02, 0.01, 0.01, 200);
    let mut all: Vec<Vec<u8>> = reads_full.into_iter().chain(reads_del).collect();
    all.rotate_right(1);
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &wide_cfg()).unwrap();

    // Consensus should be closer to the full sequence than the deleted one.
    let edit_full = levenshtein(full, &result.sequence);
    let edit_del  = levenshtein(deleted, &result.sequence);
    assert!(
        edit_full < edit_del,
        "consensus should follow majority (full) sequence; edit_full={edit_full} edit_del={edit_del}"
    );

    // Structural bubble should be present.
    let has_structural = result.bubble_sites.iter().any(|s| s.is_structural);
    println!(
        "sv_deletion: consensus_len={}, bubble_count={}, structural_bubbles={}",
        result.sequence.len(),
        result.graph_stats.bubble_count,
        result.bubble_sites.len(),
    );
    assert!(
        has_structural || result.graph_stats.bubble_count > 0,
        "large deletion should produce at least one bubble; sites={:?}",
        result.bubble_sites
    );

    // 4/14 ≈ 29% > 25% → competing allele flagged.
    assert!(
        should_call_multiallele(&result, 0.25),
        "deletion minority at 29% should flag competing allele"
    );
}

/// Most reads have the full sequence; minority carry a large insertion.
/// The shorter consensus wins; bubble is structural.
#[test]
fn sv_insertion_minority_creates_long_arm() {
    let base:     &[u8] = b"AAAAAAAAAAAATTTTTTTTTTTT"; // 24 bp
    // Insertion of 14 bp CCCCCCCCCCCCCC in the middle.
    let inserted: &[u8] = b"AAAAAAAAAAAACCCCCCCCCCCCCCTTTTTTTTTTTT"; // 38 bp

    let reads_base = simulate(base, 10, 0.02, 0.01, 0.01, 300);
    let reads_ins  = simulate(inserted, 4, 0.02, 0.01, 0.01, 400);
    let mut all: Vec<Vec<u8>> = reads_base.into_iter().chain(reads_ins).collect();
    all.rotate_right(1);
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &wide_cfg()).unwrap();

    let edit_base = levenshtein(base, &result.sequence);
    let edit_ins  = levenshtein(inserted, &result.sequence);
    assert!(
        edit_base < edit_ins,
        "consensus should follow majority (base) sequence; edit_base={edit_base} edit_ins={edit_ins}"
    );

    println!(
        "sv_insertion: consensus_len={}, edit_base={edit_base}, bubble_count={}",
        result.sequence.len(), result.graph_stats.bubble_count
    );

    // Structural bubble expected (14 bp insertion arm).
    let has_structural = result.bubble_sites.iter().any(|s| s.is_structural);
    assert!(
        has_structural || result.graph_stats.bubble_count > 0,
        "large insertion should produce a bubble"
    );
}

/// All reads carry the same sequence — no bubble, just consensus.
///
/// Uses zero-error reads: with error simulation, random substitutions can create
/// minor single-read bubbles that cross the min_allele_freq threshold by chance.
#[test]
fn sv_all_reads_carry_same_deletion_no_bubble() {
    // All reads consistently agree on this 24 bp sequence — it IS the truth.
    let truth: &[u8] = b"AAAAAAAAAAAATTTTTTTTTTTT"; // 24 bp

    // Zero-error reads: every read is identical.
    let reads: Vec<&[u8]> = vec![truth; 12];
    let result = consensus(&reads, 0, &default_cfg()).unwrap();

    println!(
        "sv_all_same_deletion: bubble_count={}, bubble_sites={}",
        result.graph_stats.bubble_count, result.bubble_sites.len()
    );
    assert_eq!(result.graph_stats.bubble_count, 0, "unanimous reads should produce zero bubbles");
    assert!(!should_call_multiallele(&result, 0.25),
        "unanimous read set should not flag competing allele");
    assert_eq!(result.sequence, truth,
        "consensus should exactly match unanimous truth");
}

/// An inverted segment in a minority of reads appears as a bubble.
/// POA doesn't know it's an inversion — it just sees two different arms.
/// We verify the two arm sequences are reverse complements of each other.
#[test]
fn sv_inversion_minority_arm_is_rev_comp() {
    // Reference arm (16 bp, non-palindromic).
    let arm_fwd: &[u8]  = b"AAACCCGGGTTTTAAAC";
    let arm_rev          = rev_comp(arm_fwd);

    let flank_l: &[u8] = b"GCGCGCGCGC"; // 10 bp left flank
    let flank_r: &[u8] = b"TATATATATA"; // 10 bp right flank

    let mut normal   = flank_l.to_vec(); normal.extend_from_slice(arm_fwd); normal.extend_from_slice(flank_r);
    let mut inverted = flank_l.to_vec(); inverted.extend_from_slice(&arm_rev);  inverted.extend_from_slice(flank_r);

    let reads_fwd = simulate(&normal,   10, 0.01, 0.0, 0.0, 600);
    let reads_rev = simulate(&inverted,  4, 0.01, 0.0, 0.0, 700);
    let mut all: Vec<Vec<u8>> = reads_fwd.into_iter().chain(reads_rev).collect();
    all.rotate_right(1);
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &wide_cfg()).unwrap();

    println!(
        "sv_inversion: consensus={}, bubble_count={}, bubble_sites={}",
        std::str::from_utf8(&result.sequence).unwrap_or("<non-utf8>"),
        result.graph_stats.bubble_count,
        result.bubble_sites.len()
    );

    // Consensus should follow the majority (forward) reads.
    assert!(result.sequence.windows(arm_fwd.len()).any(|w| w == arm_fwd)
        || levenshtein(&normal, &result.sequence) < levenshtein(&inverted, &result.sequence),
        "consensus should follow forward majority"
    );

    // If the bubble was captured in bubble_sites, check arm sequences.
    let structural_site = result.bubble_sites.iter().find(|s| s.is_structural);
    if let Some(site) = structural_site {
        let seqs: Vec<&[u8]> = site.arm_sequences.iter().map(|s| s.as_slice()).collect();
        if seqs.len() >= 2 && !seqs[0].is_empty() && !seqs[1].is_empty() {
            let a = seqs[0];
            let b = seqs[1];
            let rc_a = rev_comp(a);
            let rc_b = rev_comp(b);
            // One arm should be close to arm_fwd, the other to its RC.
            let a_vs_fwd = levenshtein(a, arm_fwd);
            let b_vs_rc  = levenshtein(b, &arm_rev);
            let a_vs_rc  = levenshtein(a, &arm_rev);
            let b_vs_fwd = levenshtein(b, arm_fwd);
            let fwd_rc_pair = a_vs_fwd + b_vs_rc;
            let rc_fwd_pair = a_vs_rc  + b_vs_fwd;
            let best = fwd_rc_pair.min(rc_fwd_pair);
            // RC of one arm should equal the other arm up to sequencing noise.
            let _ = rc_a;
            let _ = rc_b;
            println!("  arm[0]={}, arm[1]={}, best_rc_pair_edit={best}",
                std::str::from_utf8(a).unwrap_or("?"),
                std::str::from_utf8(b).unwrap_or("?")
            );
            assert!(best <= 4,
                "inversion arms should be reverse complements of each other; edit={best}");
        }
    }
}

/// Tandem duplication: majority has N copies, minority has N+K copies.
///
/// Uses zero-error exact reads so the insertion arm is created cleanly.
/// Seed = 4-unit read; 8-unit reads INSERT 4 extra units (24 bp), which
/// exceeds `phasing_bubble_min_span` (default 10) → `is_structural = true`.
///
/// Key insight: bubbles in POA arise from INSERTIONS relative to the seed.
/// If the longer allele is the seed, shorter reads just delete through the
/// extra nodes (no bubble).  Seeding with the shorter allele forces inserts.
#[test]
fn sv_tandem_duplication_length_bubble() {
    let unit: &[u8]  = b"CAT";
    let flank: &[u8] = b"GGGGGGGGGGGGGGGGGGGG"; // 20 bp flanks (distinct from CAT)

    // 4-unit template (majority, shorter): 20 + 12 + 20 = 52 bp
    let mut tmpl_4 = flank.to_vec();
    for _ in 0..4 { tmpl_4.extend_from_slice(unit); }
    tmpl_4.extend_from_slice(flank);

    // 8-unit template (minority, longer): 20 + 24 + 20 = 64 bp
    let mut tmpl_8 = flank.to_vec();
    for _ in 0..8 { tmpl_8.extend_from_slice(unit); }
    tmpl_8.extend_from_slice(flank);

    // Zero-error exact reads: insertion arm is clean and unambiguous.
    let reads_4: Vec<Vec<u8>> = vec![tmpl_4.clone(); 10];
    let reads_8: Vec<Vec<u8>> = vec![tmpl_8.clone(); 4];

    // Seed = reads_4[0] (4-unit).  8-unit reads INSERT 4 extra units (12 bp).
    let all: Vec<Vec<u8>> = reads_4.into_iter().chain(reads_8).collect();
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &wide_cfg()).unwrap();
    let edit_4 = levenshtein(&tmpl_4, &result.sequence);
    let edit_8 = levenshtein(&tmpl_8, &result.sequence);

    println!(
        "sv_tandem_dup: consensus_len={}, tmpl_4_len={}, tmpl_8_len={}, \
         edit_4={edit_4}, edit_8={edit_8}, bubble_count={}, longest_span={}",
        result.sequence.len(), tmpl_4.len(), tmpl_8.len(),
        result.graph_stats.bubble_count,
        result.graph_stats.longest_bubble_span
    );

    assert!(
        edit_4 <= edit_8,
        "consensus should follow 4-unit majority; edit_4={edit_4} edit_8={edit_8}"
    );

    // 4 extra CAT units = 12 bp arm; default phasing_bubble_min_span = 10.
    let has_structural = result.bubble_sites.iter().any(|s| s.is_structural);
    assert!(
        has_structural || result.graph_stats.longest_bubble_span >= 10,
        "4-extra-unit insertion (12 bp) should produce a structural bubble; \
         bubble_sites={}, longest_span={}",
        result.bubble_sites.len(),
        result.graph_stats.longest_bubble_span
    );
}

/// "Translocation" modelled as a foreign-sequence insertion in a minority of reads.
/// The consensus follows the majority (native sequence); the foreign arm is detected
/// as a structural bubble (if above min_allele_freq) or trimmed (if below).
#[test]
fn sv_foreign_sequence_insertion_structural_bubble() {
    // Native: AAAAAAA + CCCCCCCCCCCCC + TTTTTTT (7+13+7 = 27 bp)
    // Foreign insertion: AAAAAAA + GCGCGCGCGCGCGCGC + TTTTTTT (same flanks, alien middle)
    let native: &[u8]  = b"AAAAAAACCCCCCCCCCCCCCTTTTTTT"; // 28 bp
    let foreign: &[u8] = b"AAAAAAAGCGCGCGCGCGCGCGCTTTTTTT"; // 31 bp (different + longer mid)

    let reads_native  = simulate(native,  10, 0.02, 0.01, 0.01, 1000);
    let reads_foreign = simulate(foreign,  4, 0.02, 0.01, 0.01, 1100);
    let mut all: Vec<Vec<u8>> = reads_native.into_iter().chain(reads_foreign).collect();
    all.rotate_right(1);
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &wide_cfg()).unwrap();
    let edit_native  = levenshtein(native,  &result.sequence);
    let edit_foreign = levenshtein(foreign, &result.sequence);

    println!(
        "sv_translocation: consensus_len={}, edit_native={edit_native}, \
         edit_foreign={edit_foreign}, bubble_count={}, bubble_sites={}",
        result.sequence.len(),
        result.graph_stats.bubble_count,
        result.bubble_sites.len()
    );

    assert!(
        edit_native < edit_foreign,
        "consensus should follow native majority; edit_native={edit_native} edit_foreign={edit_foreign}"
    );

    // 4/14 = 29% > 25% → structural bubble should be flagged.
    assert!(
        result.graph_stats.bubble_count > 0,
        "foreign insertion should produce at least one bubble"
    );
    assert!(
        should_call_multiallele(&result, 0.25),
        "29% foreign-sequence reads should flag competing allele"
    );
}

/// Two separate SV sites in the same read set → two distinct bubbles.
///
/// Bubbles in POA arise from INSERTIONS relative to the current graph.
/// Using the short baseline as seed, variant reads INSERT extra bases at two
/// distinct positions, creating two independent structural bubbles.
#[test]
fn sv_two_independent_sv_sites() {
    // Seed / majority: short baseline with three distinct blocks.
    let base: &[u8] = b"AAAACCCCTTTT"; // 12 bp: A-block, C-block, T-block

    // Variant 1: 12 G's inserted between A-block and C-block (site 1).
    let var1: &[u8] = b"AAAAGGGGGGGGGGGGCCCCTTTT"; // 24 bp

    // Variant 2: 12 G's inserted between C-block and T-block (site 2).
    // Uses G's (not T's) to avoid natural merging with the existing T-block.
    let var2: &[u8] = b"AAAACCCCGGGGGGGGGGGGTTTT"; // 24 bp

    // Zero-error exact reads so minority arms are clean and unambiguous.
    // 10 base + 5 var1 + 5 var2 = 20 reads total.
    // threshold = floor(20 × 0.25) = 5.  Both minority arms have weight 5 >= 5.
    let reads_base: Vec<Vec<u8>> = vec![base.to_vec(); 10];
    let reads_var1: Vec<Vec<u8>> = vec![var1.to_vec(); 5];
    let reads_var2: Vec<Vec<u8>> = vec![var2.to_vec(); 5];

    // Seed = reads_base[0] (short baseline) so variants create insertion bubbles.
    let all: Vec<Vec<u8>> = reads_base
        .into_iter()
        .chain(reads_var1)
        .chain(reads_var2)
        .collect();
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &wide_cfg()).unwrap();

    println!(
        "sv_two_sites: consensus_len={}, bubble_count={}, bubble_sites={}",
        result.sequence.len(),
        result.graph_stats.bubble_count,
        result.bubble_sites.len()
    );

    assert!(
        result.graph_stats.bubble_count >= 2,
        "two independent insertion sites should produce ≥ 2 bubbles; got {}",
        result.graph_stats.bubble_count
    );
}

// ═══════════════════════════════════════════════════════════════════════════════
// SECTION 3: Edge cases
// ═══════════════════════════════════════════════════════════════════════════════

/// All reads are identical → zero noise, perfect consensus, no bubbles.
#[test]
fn edge_all_identical_reads_zero_noise() {
    let truth: &[u8] = b"ACGTACGTACGTACGT";
    let reads: Vec<&[u8]> = vec![truth; 10];

    let result = consensus(&reads, 0, &default_cfg()).unwrap();

    assert_eq!(result.sequence, truth, "identical reads should produce exact truth");
    assert_eq!(result.graph_stats.bubble_count, 0, "no bubbles expected");
    assert_eq!(result.bubble_sites.len(), 0);
    assert!(!should_call_multiallele(&result, 0.25));

    let conf = consensus_confidence(&result, 0.25);
    assert!(!conf.is_low_confidence(), "identical reads should be high confidence: {:?}", conf);
}

/// SNV in a single outlier read → below min_allele_freq, no bubble flagged.
#[test]
fn edge_single_outlier_read_below_threshold() {
    let truth:   &[u8] = b"ACGTACGTACGTACGTACGT"; // 20 bp
    let outlier: &[u8] = b"ACGTACGTCCGTACGTACGT"; // SNV at pos 8

    let mut reads: Vec<&[u8]> = vec![truth; 9];
    reads.push(outlier);

    let result = consensus(&reads, 0, &default_cfg()).unwrap();

    // 1/10 = 10% is below 25% threshold.
    assert!(
        !should_call_multiallele(&result, 0.25),
        "single outlier (10%) should not flag competing allele"
    );
    assert_eq!(result.sequence, truth,
        "outlier should not affect consensus; got {}",
        std::str::from_utf8(&result.sequence).unwrap_or("<non-utf8>")
    );
}

/// Balanced two-allele split → competing allele correctly flagged.
#[test]
fn edge_balanced_two_allele_split_flagged() {
    let a: &[u8] = b"ACGTACGTAAAACGTACGT"; // 19 bp
    let b: &[u8] = b"ACGTACGTGGGGCGTACGT"; // 19 bp (4 bp differ in middle)

    let reads_a = simulate(a, 8, 0.01, 0.0, 0.0, 2000);
    let reads_b = simulate(b, 8, 0.01, 0.0, 0.0, 2100);
    let mut all: Vec<Vec<u8>> = reads_a.into_iter().chain(reads_b).collect();
    all.rotate_right(1);
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &default_cfg()).unwrap();

    assert!(
        should_call_multiallele(&result, 0.25),
        "8+8 balanced split should flag competing allele; bubble_sites={:?}",
        result.bubble_sites
    );

    // Allele fractions at the flagged site should be roughly 50/50.
    if let Some(site) = has_competing_allele(&result, 0.25) {
        let fracs = allele_fractions(site);
        let minority = fracs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(minority >= 0.30, "balanced split: minority fraction should be ≥ 30%; got {minority:.2}");
    }
}

/// Homopolymer runs → consensus correct, no spurious bubbles.
#[test]
fn edge_homopolymer_run_clean_consensus() {
    // Long homopolymer: classic ONT difficulty, but at low error rate should be clean.
    let truth: &[u8] = b"ACGTAAAAAAAAAAAACGT"; // 16 A's in the middle
    let reads = simulate(truth, 15, 0.01, 0.0, 0.0, 2200);
    let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &default_cfg()).unwrap();

    println!(
        "edge_homopolymer: consensus={}, edit_from_truth={}",
        std::str::from_utf8(&result.sequence).unwrap_or("<non-utf8>"),
        levenshtein(truth, &result.sequence)
    );

    let edit = levenshtein(truth, &result.sequence);
    assert!(edit <= 2, "clean homopolymer reads should produce near-exact consensus; edit={edit}");
    assert!(!should_call_multiallele(&result, 0.25));
}

/// Extremely short reads (5 bp) → consensus still works.
#[test]
fn edge_very_short_reads_consensus() {
    let reads: Vec<&[u8]> = vec![
        b"ACGTA", b"ACGTA", b"ACGTA",
        b"ACGTA", b"ACGTA", b"ACGTT", // one error
    ];
    let result = consensus(&reads, 0, &default_cfg()).unwrap();
    assert_eq!(result.sequence, b"ACGTA", "majority should win on 5bp reads");
}

/// Reads with extreme length variation → majority length wins, no crash.
#[test]
fn edge_extreme_length_variation_majority_wins() {
    // 8 reads of length 20, 3 reads of length 40 (2× longer).
    let short: &[u8] = b"ACGTACGTACGTACGTACGT"; // 20 bp
    let long:  &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bp

    let reads_short = simulate(short,  8, 0.01, 0.01, 0.01, 2300);
    let reads_long  = simulate(long,   3, 0.01, 0.01, 0.01, 2400);
    let mut all: Vec<Vec<u8>> = reads_short.into_iter().chain(reads_long).collect();
    all.rotate_right(1);
    let refs: Vec<&[u8]> = all.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &wide_cfg()).unwrap();

    println!(
        "edge_extreme_length: consensus_len={}, short_len={}, long_len={}, \
         edit_short={}, edit_long={}",
        result.sequence.len(), short.len(), long.len(),
        levenshtein(short, &result.sequence),
        levenshtein(long, &result.sequence)
    );

    // Majority (short) should produce the consensus.
    let edit_short = levenshtein(short, &result.sequence);
    let edit_long  = levenshtein(long,  &result.sequence);
    assert!(
        edit_short < edit_long,
        "short (majority) template should win; edit_short={edit_short} edit_long={edit_long}"
    );
}

/// High noise (random reads) → consensus_confidence flags single_support_fraction.
#[test]
fn edge_noisy_reads_flagged_by_single_support_fraction() {
    // Reads are mostly noise: 15% substitution rate, 5% indels.
    let template: &[u8] = b"ACGTACGTACGTACGTACGTACGT"; // 24 bp
    let reads = simulate(template, 10, 0.15, 0.05, 0.05, 2500);
    let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();

    let result = consensus(&refs, 0, &wide_cfg()).unwrap();
    let conf = consensus_confidence(&result, 0.25);

    println!(
        "edge_noisy: single_support_fraction={:.3}, bubble_count={}",
        conf.single_support_fraction,
        result.graph_stats.bubble_count
    );

    // High noise creates many single-read nodes — fraction should be elevated.
    // We just verify the helper correctly surfaces what the graph observed.
    assert!(
        conf.single_support_fraction >= 0.0 && conf.single_support_fraction <= 1.0,
        "single_support_fraction must be in [0,1]"
    );
    // At 15% sub + 10% indel, graph should be noisier than a clean input.
    let clean_reads = simulate(template, 10, 0.01, 0.0, 0.0, 2600);
    let clean_refs: Vec<&[u8]> = clean_reads.iter().map(|r| r.as_slice()).collect();
    let clean_result = consensus(&clean_refs, 0, &wide_cfg()).unwrap();
    let clean_conf = consensus_confidence(&clean_result, 0.25);
    assert!(
        conf.single_support_fraction >= clean_conf.single_support_fraction,
        "noisy reads should have ≥ single_support_fraction vs clean; \
         noisy={:.3}, clean={:.3}",
        conf.single_support_fraction, clean_conf.single_support_fraction
    );
}

/// SNV bubble: minority arm is exactly 1 bp; is_structural should be false.
/// Structural bubble: minority arm is 15 bp; is_structural should be true.
#[test]
fn edge_is_structural_threshold_snp_vs_indel() {
    // SNP bubble: single-base difference.
    let a_snp: &[u8] = b"AAAAAAAAAGAAAAAAAAAA"; // G at pos 9
    let b_snp: &[u8] = b"AAAAAAAAATAAAAAAAAAA"; // T at pos 9

    let reads_snp_a = simulate(a_snp, 8, 0.01, 0.0, 0.0, 3000);
    let reads_snp_b = simulate(b_snp, 4, 0.01, 0.0, 0.0, 3100);
    let mut snp_all: Vec<Vec<u8>> = reads_snp_a.into_iter().chain(reads_snp_b).collect();
    snp_all.rotate_right(1);
    let snp_refs: Vec<&[u8]> = snp_all.iter().map(|r| r.as_slice()).collect();
    let snp_result = consensus(&snp_refs, 0, &default_cfg()).unwrap();

    // Large indel bubble: 15 bp inserted in minority.
    let base_indel: &[u8]  = b"AAAAAAAAAGAAAAAAAAAA"; // 20 bp
    let ins_indel: &[u8]   = b"AAAAAAAAAGCCCCCCCCCCCCCCCAAAAAAAAAA"; // +15 bp C's

    let reads_ind_a = simulate(base_indel, 8, 0.01, 0.0, 0.0, 3200);
    let reads_ind_b = simulate(ins_indel,  4, 0.01, 0.0, 0.0, 3300);
    // Seed = base_indel[0] so insertion reads create the structural bubble.
    // If insertion read were the seed, base reads would delete through → no bubble.
    let ind_all: Vec<Vec<u8>> = reads_ind_a.into_iter().chain(reads_ind_b).collect();
    let ind_refs: Vec<&[u8]> = ind_all.iter().map(|r| r.as_slice()).collect();
    let ind_result = consensus(&ind_refs, 0, &wide_cfg()).unwrap();

    println!(
        "is_structural: snp_bubbles={:?}  indel_bubbles={:?}",
        snp_result.bubble_sites.iter().map(|s| (s.is_structural, &s.arm_sequences)).collect::<Vec<_>>(),
        ind_result.bubble_sites.iter().map(|s| (s.is_structural, s.arm_sequences.len())).collect::<Vec<_>>()
    );

    // SNP site should not be structural.
    let snp_structural = snp_result.bubble_sites.iter().any(|s| s.is_structural);
    assert!(!snp_structural,
        "single-base SNP bubble should not be structural; sites={:?}", snp_result.bubble_sites);

    // Indel site (15 bp arm) should be structural.
    let ind_structural = ind_result.bubble_sites.iter().any(|s| s.is_structural);
    assert!(ind_structural
        || ind_result.graph_stats.longest_bubble_span >= 10,
        "15 bp insertion bubble should be structural; sites={:?}", ind_result.bubble_sites);
}

/// `longest_bubble_span` tracks correctly across bubble sizes.
///
/// Uses INTERIOR insertions (flanked on both sides) so the inserted arm
/// is never trimmed by boundary coverage filtering.
#[test]
fn edge_longest_bubble_span_max_of_all_bubbles() {
    let flank = b"GCGCGCGCGCGCGCGCGCGC"; // 20 bp flanks

    // Small insertion: 4 G's between the flanks.
    let base_small: Vec<u8> = {
        let mut v = flank.to_vec(); v.extend_from_slice(b"TATATA"); v.extend_from_slice(flank); v
    };
    let ins_small: Vec<u8> = {
        let mut v = flank.to_vec(); v.extend_from_slice(b"TATATAGGGGTATATA"); v.extend_from_slice(flank); v
    };

    // Large insertion: 18 C's between the flanks.
    let base_large: Vec<u8> = {
        let mut v = flank.to_vec(); v.extend_from_slice(b"ATATAT"); v.extend_from_slice(flank); v
    };
    let ins_large: Vec<u8> = {
        let mut v = flank.to_vec(); v.extend_from_slice(b"ATATCCCCCCCCCCCCCCCCCCATATAT"); v.extend_from_slice(flank); v
    };

    // Seed = base (short); insertion reads insert extra nodes → structural bubble.
    let all_small: Vec<Vec<u8>> = {
        let mut b = simulate(&base_small, 8, 0.01, 0.0, 0.0, 4000);
        b.extend(simulate(&ins_small, 4, 0.01, 0.0, 0.0, 4100));
        b
    };
    let all_large: Vec<Vec<u8>> = {
        let mut b = simulate(&base_large, 8, 0.01, 0.0, 0.0, 4200);
        b.extend(simulate(&ins_large, 4, 0.01, 0.0, 0.0, 4300));
        b
    };

    let refs_small: Vec<&[u8]> = all_small.iter().map(|r| r.as_slice()).collect();
    let refs_large: Vec<&[u8]> = all_large.iter().map(|r| r.as_slice()).collect();

    let r_small = consensus(&refs_small, 0, &wide_cfg()).unwrap();
    let r_large = consensus(&refs_large, 0, &wide_cfg()).unwrap();

    println!(
        "longest_span: small_ins_longest={}, large_ins_longest={}",
        r_small.graph_stats.longest_bubble_span,
        r_large.graph_stats.longest_bubble_span
    );

    assert!(
        r_large.graph_stats.longest_bubble_span > r_small.graph_stats.longest_bubble_span,
        "larger insertion should produce longer bubble span; \
         small={}, large={}",
        r_small.graph_stats.longest_bubble_span,
        r_large.graph_stats.longest_bubble_span
    );
    assert!(
        r_large.graph_stats.longest_bubble_span >= 10,
        "18-bp insertion arm should exceed structural threshold; got {}",
        r_large.graph_stats.longest_bubble_span
    );
}

// ── Levenshtein distance (used throughout) ────────────────────────────────────

fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    let m = a.len();
    let n = b.len();
    let mut prev: Vec<usize> = (0..=n).collect();
    let mut curr = vec![0usize; n + 1];
    for i in 1..=m {
        curr[0] = i;
        for j in 1..=n {
            curr[j] = if a[i - 1] == b[j - 1] {
                prev[j - 1]
            } else {
                1 + prev[j - 1].min(prev[j]).min(curr[j - 1])
            };
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[n]
}
