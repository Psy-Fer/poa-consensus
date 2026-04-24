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
}
