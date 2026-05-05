/// Extract the repeat segment from `read` by locating `left_flank` and `right_flank`
/// using semi-global pairwise alignment.
///
/// Returns a slice of `read` that begins immediately after the best-scoring position of
/// `left_flank` and ends immediately before the best-scoring position of `right_flank`.
/// Returns `None` if either flank cannot be located with sufficient confidence, or if the
/// located positions are inconsistent (right starts before left ends).
///
/// The alignment is semi-global: the flank must be fully consumed (no free end gaps on
/// the flank side), but the read may have unaligned prefix and suffix (free end gaps on
/// the read side). Scoring: match +2, mismatch -1, gap -2.  The gap penalty is set to
/// -2 (rather than -1) so that skipping two read bases never outscores a direct
/// mismatch, keeping the anchor tight.  The score threshold requires the best alignment
/// to reach at least `flank.len()`, so noisy-but-present flanks pass while completely
/// absent flanks are rejected.
pub fn extract_flanked_region<'a>(
    read: &'a [u8],
    left_flank: &[u8],
    right_flank: &[u8],
) -> Option<&'a [u8]> {
    let seg_start = find_flank_end(read, left_flank)?;
    let seg_end = find_right_flank_start(read, right_flank)?;
    if seg_start >= seg_end {
        return None;
    }
    Some(&read[seg_start..seg_end])
}

// Gap penalty used by find_flank_end.  Must be more negative than half the match
// bonus (+2) so that two insertions (cost 2*GAP = -4) never beat a mismatch (-1).
// With GAP = -2: a 5-match+2-gap path scores 6 < 7 (4-match+1-mismatch), so the
// DP anchors at the compact correct position rather than spanning extra read chars.
const GAP: i32 = -2;

/// Returns the position in `read` immediately after the best semi-global alignment of
/// `flank` ends (i.e. where the repeat segment should start).
fn find_flank_end(read: &[u8], flank: &[u8]) -> Option<usize> {
    if flank.is_empty() {
        return Some(0);
    }
    if read.len() < flank.len() {
        return None;
    }

    let m = flank.len();
    let n = read.len();
    let cols = n + 1;
    let mut dp = vec![0i32; (m + 1) * cols];

    // Consuming flank chars without advancing in read costs GAP per base.
    for i in 1..=m {
        dp[i * cols] = GAP * i as i32;
    }
    // dp[0][j] = 0 for all j: alignment may start at any read position (free prefix).

    for i in 1..=m {
        for j in 1..=n {
            let sub = if flank[i - 1] == read[j - 1] {
                2i32
            } else {
                -1i32
            };
            let mat = dp[(i - 1) * cols + (j - 1)] + sub;
            let del = dp[(i - 1) * cols + j] + GAP; // consume flank char, skip read base
            let ins = dp[i * cols + (j - 1)] + GAP; // consume read char, skip flank base
            dp[i * cols + j] = mat.max(del).max(ins);
        }
    }

    // Take the FIRST position with the maximum score so the flank anchors as early
    // as possible.  max_by_key breaks ties by returning the last element, which
    // overshoots the anchor when the flank matches at several equivalent positions.
    let mut best_j = 1;
    let mut best_score = dp[m * cols + 1];
    for j in 2..=n {
        let s = dp[m * cols + j];
        if s > best_score {
            best_score = s;
            best_j = j;
        }
    }

    // Require score ≥ flank length: ~75% identity at zero indels, or near-perfect
    // with a few gaps.  Completely absent flanks score near zero or negative.
    if best_score < flank.len() as i32 {
        return None;
    }

    Some(best_j)
}

/// Returns the position in `read` where the right flank begins (i.e. where the repeat
/// segment should end).  Implemented by reversing both sequences and calling
/// `find_flank_end`, then converting back to the original coordinate space.
fn find_right_flank_start(read: &[u8], right_flank: &[u8]) -> Option<usize> {
    if right_flank.is_empty() {
        return Some(read.len());
    }
    let rev_read: Vec<u8> = read.iter().rev().copied().collect();
    let rev_flank: Vec<u8> = right_flank.iter().rev().copied().collect();
    let end_in_rev = find_flank_end(&rev_read, &rev_flank)?;
    Some(read.len() - end_in_rev)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_flanks() {
        let read = b"GGGGGCATCATCATCATAAAAA";
        let seg =
            extract_flanked_region(read, b"GGGGG", b"AAAAA").expect("should find both flanks");
        assert_eq!(seg, b"CATCATCATCAT");
    }

    #[test]
    fn single_base_flanks() {
        let read = b"GCATCATCATA";
        let seg = extract_flanked_region(read, b"G", b"A").expect("single-base flanks");
        // left: G ends at position 1; right: A starts at position 10
        assert_eq!(seg, b"CATCATCAT");
    }

    #[test]
    fn no_left_flank_returns_none() {
        let read = b"GGGGGCATCATCATCATAAAAA";
        assert!(extract_flanked_region(read, b"TTTTT", b"AAAAA").is_none());
    }

    #[test]
    fn no_right_flank_returns_none() {
        let read = b"GGGGGCATCATCATCATAAAAA";
        assert!(extract_flanked_region(read, b"GGGGG", b"TTTTT").is_none());
    }

    #[test]
    fn read_shorter_than_flank_returns_none() {
        assert!(extract_flanked_region(b"ACG", b"GGGGG", b"AAAAA").is_none());
    }

    #[test]
    fn noisy_flanks() {
        // One mismatch in each flank (sequencing error scenario)
        // left:  GGGGG → GGGAG (1 mismatch)
        // right: AAAAA → AAATA (1 mismatch)
        let read = b"GGGAGCATCATCATCATAAATA";
        let seg = extract_flanked_region(read, b"GGGGG", b"AAAAA")
            .expect("noisy flanks should still align");
        assert_eq!(seg, b"CATCATCATCAT");
    }

    #[test]
    fn overlapping_flanks_returns_none() {
        // Flanks overlap so there is no valid segment between them
        let read = b"GGGGGAAAAA";
        assert!(extract_flanked_region(read, b"GGGGGAAAAA", b"GGGGGAAAAA").is_none());
    }

    #[test]
    fn empty_repeat_segment_returns_none() {
        // Left flank ends exactly where right flank begins → zero-length segment
        let read = b"GGGGGAAAAA";
        assert!(extract_flanked_region(read, b"GGGGG", b"AAAAA").is_none());
    }

    #[test]
    fn longer_read_with_flanked_str() {
        // Simulate a realistic read: 10-bp flanks + 15-bp repeat
        let read = b"ATCGATCGATCATCATCATCATCATTAGCTAGCTA";
        //           ^--left 10--^^---repeat 15---^^--right 10--^
        // left_flank = ATCGATCGAT (10 bp), right_flank = TAGCTAGCTA (10 bp)
        let seg = extract_flanked_region(read, b"ATCGATCGAT", b"TAGCTAGCTA")
            .expect("10-bp flanks + 15-bp repeat");
        assert_eq!(seg, b"CATCATCATCATCAT");
    }

    // ── Error tolerance vs flank length ──────────────────────────────────────
    //
    // The score threshold is flank.len().  Perfect score = 2L.  Each mismatch
    // costs 3 relative to perfect, so failure requires > L/3 errors.
    //
    //  5 bp  → fails at ≥ 2 mismatches (40% error)
    // 10 bp  → fails at ≥ 4 mismatches (40%)
    // 20 bp  → fails at ≥ 7 mismatches (35%) — safe through ONT R9.4 (~10%)
    // 30 bp  → fails at ≥ 11 mismatches (37%)
    //
    // At ONT R10 (~8% error), expected errors per flank:
    //  5 bp → 0.4 expected  (fails ~6% of reads — unacceptable)
    // 10 bp → 0.8 expected  (fails ~0.5% of reads — marginal)
    // 20 bp → 1.6 expected  (fails ~0.003% — safe)
    // 30 bp → 2.4 expected  (essentially never fails)

    #[test]
    fn short_5bp_flank_fails_at_two_mismatches() {
        // 5 bp flank, 2 mismatches: score = 2×3 − 2 = 4 < 5 → None.
        // Represents a partial read where sequencing errors corrupt a short flank.
        // Two errors in 5 bp is only 40% error rate, well above ONT R10's 8%
        // average, but easily reached by chance.
        let read = b"GAGAGCATCATCATCATAAAAA"; // "GGGGG" → "GAGAG" (2 errors: pos 1 G→A, pos 3 G→A)
        assert!(
            extract_flanked_region(read, b"GGGGG", b"AAAAA").is_none(),
            "5 bp flank with 2 mismatches should fail the score threshold"
        );
    }

    #[test]
    fn medium_10bp_flank_survives_two_mismatches() {
        // 10 bp flank, 2 mismatches: score = 2×8 − 2 = 14 ≥ 10 → Some.
        // The same absolute error count that kills a 5 bp flank is easily
        // absorbed by a 10 bp flank.
        // Left flank: ATCGATCGAT → ATCTATCGAA (2 errors: pos 3 G→T, pos 9 T→A)
        let read = b"ATCTATCGAACATCATCATCATCATTAGCTAGCTA";
        //            ^---flank with 2 errors---^^--repeat--^^--right flank--^
        let seg = extract_flanked_region(read, b"ATCGATCGAT", b"TAGCTAGCTA")
            .expect("10 bp flank should survive 2 mismatches");
        assert_eq!(seg, b"CATCATCATCATCAT");
    }

    #[test]
    fn medium_10bp_flank_fails_at_four_mismatches() {
        // 10 bp flank, 4 mismatches: score = 2×6 − 4 = 8 < 10 → None.
        // Four errors in 10 bp is 40% — above the failure threshold.
        // Flanks chosen to avoid motif overlap with the CAG repeat (no ATC subsequence),
        // preventing gap-rescued alignments from reaching the threshold.
        // GGTTAACCGG → TGTAAAACGT: errors at pos 0 (G→T), 3 (T→A), 6 (C→A), 9 (G→T)
        // score = 6×2 − 4×1 = 8 < 10 → None (all gapped alignments also score < 10)
        let read = b"TGTAAAACGTCAGCAGCAGCAGCAGCCTTAAGGCC";
        assert!(
            extract_flanked_region(read, b"GGTTAACCGG", b"CCTTAAGGCC").is_none(),
            "10 bp flank with 4 mismatches should fail"
        );
    }

    #[test]
    fn long_20bp_flank_survives_six_mismatches() {
        // 20 bp flank, 6 mismatches: score = 2×14 − 6 = 22 ≥ 20 → Some.
        // Six errors in 20 bp is 30% — beyond ONT R9.4 worst-case, still passes.
        // This is the recommended minimum for ONT data.
        let lf = b"AATTCCGGAATTCCGGAATT";
        let rf = b"TTAAGGCCTTAAGGCCTTAA";
        let repeat = b"CAGCAGCAG";
        let mut lf_noisy = *lf;
        lf_noisy[0] ^= 1; // flip lsb to change base
        lf_noisy[3] ^= 1;
        lf_noisy[7] ^= 1;
        lf_noisy[11] ^= 1;
        lf_noisy[14] ^= 1;
        lf_noisy[18] ^= 1;
        let mut read2: Vec<u8> = Vec::new();
        read2.extend_from_slice(&lf_noisy);
        read2.extend_from_slice(repeat);
        read2.extend_from_slice(rf);
        // Only proceed if we actually have 6 mismatches (sanity check in test)
        let mismatches = lf_noisy
            .iter()
            .zip(lf.iter())
            .filter(|(a, b)| a != b)
            .count();
        assert_eq!(
            mismatches, 6,
            "test setup: expected 6 mismatches in noisy flank"
        );
        let seg = extract_flanked_region(&read2, lf, rf)
            .expect("20 bp flank should survive 6 mismatches (30% error)");
        assert_eq!(seg, repeat.as_ref());
    }

    #[test]
    fn long_20bp_flank_fails_at_seven_mismatches() {
        // 20 bp flank, 7 mismatches: score = 2×13 − 7 = 19 < 20 → None.
        // Just past the failure threshold.
        let lf = b"AATTCCGGAATTCCGGAATT";
        let rf = b"TTAAGGCCTTAAGGCCTTAA";
        let mut lf_noisy = *lf;
        for pos in [0, 3, 7, 11, 14, 17, 19] {
            lf_noisy[pos] ^= 1;
        }
        let mut read: Vec<u8> = Vec::new();
        read.extend_from_slice(&lf_noisy);
        read.extend_from_slice(b"CAGCAGCAG");
        read.extend_from_slice(rf);
        assert!(
            extract_flanked_region(&read, lf, rf).is_none(),
            "20 bp flank with 7 mismatches should fail"
        );
    }

    // ── Partial-read / minimum-captured-flank tests ───────────────────────────
    //
    // A partial read has limited flanking sequence at one or both ends.  The
    // query flank must be no longer than the flank the read actually captured;
    // if the query is longer than the captured region the alignment will fail.
    //
    // Practical rule: set query flank length ≤ (read_len − repeat_len) / 2.

    #[test]
    fn partial_read_5bp_captured_passes_5bp_query() {
        // Read barely spans: only 5 bp of each flank captured.
        // A 5 bp query finds them perfectly.
        let read = b"GGGGGCATCATCATCATAAAAA";
        let seg = extract_flanked_region(read, b"GGGGG", b"AAAAA")
            .expect("5 bp query matches 5 bp captured flank");
        assert_eq!(seg, b"CATCATCATCAT");
    }

    #[test]
    fn partial_read_5bp_captured_fails_20bp_query() {
        // Same read, but the query asks for 20 bp of left flank.
        // The read only has 5 bp, so the alignment can't score ≥ 20.
        // score ≈ 2×5 − 15×2 = −20  (15 gap penalties) → well below threshold.
        let read = b"GGGGGCATCATCATCATAAAAA";
        assert!(
            extract_flanked_region(read, b"AAAAAAAAAAAAAAAAGGGGG", b"AAAAA").is_none(),
            "20 bp query must fail when read only has 5 bp of captured flank"
        );
    }

    #[test]
    fn partial_read_10bp_captured_passes_10bp_query_with_one_error() {
        // Read has 10 bp of each flank captured; 1 sequencing error in left flank.
        // A 10 bp query survives 1 error: score = 2×9 − 1 = 17 ≥ 10.
        let read_1err = b"ATCGTTCGATCATCATCATCATCATTAGCTAGCTA"; // pos 4 A→T (ATCGATCGAT → ATCGTTCGAT)
        let seg = extract_flanked_region(read_1err, b"ATCGATCGAT", b"TAGCTAGCTA")
            .expect("10 bp query should survive 1 mismatch");
        assert_eq!(seg, b"CATCATCATCATCAT");
    }

    // ── Flank containing repeat k-mers ────────────────────────────────────────
    //
    // The flank QUERY must contain enough unique (non-repeat) sequence that the
    // first-maximum anchor lands at the true boundary.  If the flank query is
    // composed entirely of the repeat unit the anchor will be wrong.

    #[test]
    fn flank_ending_in_repeat_unit_anchors_correctly_with_unique_prefix() {
        // Left flank ends in "CAT" (= the repeat unit) but has a unique 5-bp
        // prefix "GGGGG".  The unique prefix scores higher at the true boundary
        // than at any shifted repeat position, so the anchor is correct.
        // Read:  GGGGG + CAT (part of left query) + CATCATCAT (repeat) + GGGGG (right flank)
        // Layout: 0-4=GGGGG, 5-7=CAT, 8-16=CATCATCAT, 17-21=GGGGG  (22 bytes total)
        let read = b"GGGGGCATCATCATCATGGGGG"; // 4 × CAT
        // Left query "GGGGGCAT" (8 bp) anchors at 0-7 → seg_start = 8
        // Right query "GGGGG" anchors at 17-21 → seg_end = 17
        let seg = extract_flanked_region(read, b"GGGGGCAT", b"GGGGG")
            .expect("unique prefix keeps anchor at correct boundary");
        assert_eq!(seg, b"CATCATCAT");
    }

    #[test]
    fn flank_query_identical_to_repeat_unit_anchors_at_first_occurrence() {
        // Degenerate case: the flank query IS the repeat unit.
        // The anchor lands at the first occurrence, which may not be the true
        // boundary.  This test documents the known failure mode so callers know
        // that flanks must contain unique (non-repeat) sequence.
        let read = b"CATCATCATCATGGG";
        //           CATCAT  repeat, GGG = right flank
        // Left query = "CAT" = repeat unit itself.
        // First occurrence of "CAT" is at position 0 → seg_start = 3.
        // seg = "CATCATCAT" (3 units), but a caller who set the boundary at
        // the start of the repeat would expect 4 units.
        let seg = extract_flanked_region(read, b"CAT", b"GGG")
            .expect("degenerate flank returns something");
        // Anchors at first CAT, missing the leading repeat unit.
        assert_eq!(
            seg.len() % 3,
            0,
            "result is still a whole number of CAT units"
        );
        assert!(
            seg.len() < 12,
            "but it is shorter than the full repeat (known anchor drift)"
        );
    }

    // ── Realistic STR scenarios ───────────────────────────────────────────────

    #[test]
    fn cag_repeat_10bp_flanks_one_error_each() {
        // HD / SCA-like locus: CAG repeat, 10 bp flanks, 1 error per flank.
        // Represents a typical HiFi read (rare errors).
        let left_true = b"ATCGATCGAT"; // 10 bp unique
        let right_true = b"TAGCTAGCTA"; // 10 bp unique
        let repeat = b"CAGCAGCAGCAGCAGCAGCAG"; // 7 × CAG
        // 1 error each flank:
        let left_noisy = b"ATCGATCGTT"; // pos 8 A→T
        let right_noisy = b"TAGCTAGCCA"; // pos 8 T→C
        let read: Vec<u8> = [left_noisy.as_ref(), repeat, right_noisy].concat();
        let seg = extract_flanked_region(&read, left_true, right_true)
            .expect("10 bp flanks + 1 error each should pass for HiFi");
        assert_eq!(seg, repeat.as_ref());
    }

    #[test]
    fn gaa_repeat_frda_like_20bp_flanks_with_noise() {
        // FRDA-like locus: GAA repeat, 20 bp flanks, 2 errors in left flank.
        // Left flank ends in "GAA" (= repeat unit) but unique prefix saves it.
        // Represents a safe ONT R10 scenario.
        let left_true = b"TTCCTGCAGTTCCTGCAGAA"; // 20 bp, last 3 = GAA
        let right_true = b"TTCTTGCAGTTCTTGCAGTT"; // 20 bp unique
        let repeat = b"GAAGAAGAAGAAGAAGAAGAA"; // 7 × GAA
        // 2 errors in left flank (10% of 20 bp = realistic ONT rate):
        let mut left_noisy = *left_true;
        left_noisy[2] ^= 2; // flip a bit to change base
        left_noisy[9] ^= 2;
        let read: Vec<u8> = [left_noisy.as_ref(), repeat, right_true.as_ref()].concat();
        let seg = extract_flanked_region(&read, left_true, right_true)
            .expect("20 bp flank ending in repeat unit + 2 errors should pass");
        assert_eq!(seg, repeat.as_ref());
    }

    #[test]
    fn ont_r10_worst_case_20bp_flank_passes() {
        // Worst-case ONT R10 scenario: 6 errors in a 20 bp flank.
        // 6/20 = 30% error rate, just inside the pass threshold (fail at ≥ 7).
        // Demonstrates that 20 bp is the safe minimum for ONT R10.
        let lf = b"AATTCCGGAATTCCGGAATT";
        let rf = b"TTAAGGCCTTAAGGCCTTAA";
        let mut lf_noisy = *lf;
        for pos in [1, 4, 8, 12, 15, 18] {
            // 6 positions
            lf_noisy[pos] ^= 1;
        }
        let mismatches = lf_noisy
            .iter()
            .zip(lf.iter())
            .filter(|(a, b)| a != b)
            .count();
        assert_eq!(mismatches, 6);
        let mut read: Vec<u8> = Vec::new();
        read.extend_from_slice(&lf_noisy);
        read.extend_from_slice(b"CAGCAGCAGCAGCAG");
        read.extend_from_slice(rf);
        let seg = extract_flanked_region(&read, lf, rf)
            .expect("20 bp flank survives 6 mismatches (30% = ONT worst-case)");
        assert_eq!(seg, b"CAGCAGCAGCAGCAG");
    }
}
