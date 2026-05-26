use crate::error::PoaError;
use std::collections::{HashMap, HashSet};

/// Strategy for selecting the seed read used to initialise the POA graph.
///
/// The seed read is built into the graph as a linear chain of nodes.  All
/// other reads are aligned into this backbone, so the choice affects bubble
/// visibility: reads *longer* than the seed produce `Insert` operations
/// (explicit bubble arms visible in `arm_sequences`), while reads *shorter*
/// than the seed produce `Delete` operations (skip edges with empty
/// `arm_sequences`).  Picking the shortest fully-spanning read as seed
/// maximises `Insert`-type bubble coverage from all longer reads.
///
/// # Which strategy to use
///
/// | Scenario | Strategy |
/// |---|---|
/// | Caller already knows which read is spanning (e.g., bladerunner) | `Explicit(idx)` |
/// | All reads are full-length, similar lengths | `Shortest` |
/// | Mixed spanning + partial reads, source unknown | `Auto` |
/// | Expansion larger than any read (no spanning reads) | `Auto` → `Err(NoSpanningReads)` → use `bridged_consensus` |
#[derive(Debug, Clone)]
pub enum SeedSelection {
    /// Terminal k-mer frequency heuristic: find the shortest read that has
    /// common k-mers at *both* its left and right terminal regions.
    ///
    /// The first and last [`TERM_LEN`] bases of every read are k-mer-sampled.
    /// A k-mer is a valid anchor only if it is *end-specific*: common at one
    /// end (≥ `max(2, floor(n × 0.3))` reads) but rare at the other
    /// (`< threshold` reads).  K-mers frequent at both ends are interior
    /// repetitive sequence and are discarded.
    ///
    /// **Fallback hierarchy:**
    /// 1. Spanning candidates found → shortest candidate.
    /// 2. Non-overlapping left-only and right-only groups, no bridge →
    ///    `Err(PoaError::NoSpanningReads)` — use [`bridged_consensus`].
    /// 3. No cluster structure (all reads too short or highly repetitive) →
    ///    longest read.
    ///
    /// Complexity: O(n × [`TERM_LEN`]) — no alignment performed.
    ///
    /// [`bridged_consensus`]: crate::bridged_consensus
    Auto,

    /// Use the read at `idx` as the seed.
    ///
    /// Returns `Err(SeedOutOfBounds)` if `idx >= reads.len()`.
    Explicit(usize),

    /// Use the shortest read as the seed without any k-mer analysis.
    ///
    /// All longer reads will produce `Insert` operations (explicit bubble
    /// arms).  Appropriate when all reads are known to be fully spanning and
    /// the caller simply wants to maximise bubble visibility.
    Shortest,
}

// ── Constants ─────────────────────────────────────────────────────────────────

/// k-mer length for terminal anchor sampling.  15-mers are specific enough to
/// be nearly unique in non-repetitive flanking sequence.
const TERM_K: usize = 15;

/// Number of bases sampled from each end of a read.  50 bp covers typical STR
/// flanking anchors while remaining cheap even for short reads.
///
/// Reads shorter than `2 × TERM_LEN` have overlapping terminal windows; for
/// those reads the end-specific filter suppresses the shared k-mers, so they
/// contribute to anchoring only if their flanks are distinct from the interior.
const TERM_LEN: usize = 50;

/// Fraction of reads that must share a terminal k-mer for it to count as a
/// common anchor.  30 % works down to depth ≈7 (threshold clamps to 2).
const ANCHOR_FRAC: f64 = 0.3;

// ── Public API ────────────────────────────────────────────────────────────────

/// Select the index of the seed read from `reads` using `selection`.
///
/// Returns the chosen index on success.  May return:
/// - [`PoaError::EmptyInput`] — `reads` is empty.
/// - [`PoaError::SeedOutOfBounds`] — `Explicit(idx)` with `idx >= reads.len()`.
/// - [`PoaError::NoSpanningReads`] — `Auto` detected two non-overlapping read
///   groups with no bridging read; use [`bridged_consensus`] instead.
///
/// [`bridged_consensus`]: crate::bridged_consensus
pub fn select_seed(reads: &[&[u8]], selection: &SeedSelection) -> Result<usize, PoaError> {
    if reads.is_empty() {
        return Err(PoaError::EmptyInput);
    }
    match selection {
        SeedSelection::Explicit(idx) => {
            if *idx >= reads.len() {
                Err(PoaError::SeedOutOfBounds { index: *idx, len: reads.len() })
            } else {
                Ok(*idx)
            }
        }
        SeedSelection::Shortest => Ok(shortest(reads)),
        SeedSelection::Auto => auto_select(reads),
    }
}

// ── Internal ──────────────────────────────────────────────────────────────────

fn shortest(reads: &[&[u8]]) -> usize {
    reads.iter().enumerate().min_by_key(|(_, r)| r.len()).map(|(i, _)| i).unwrap()
}

fn longest(reads: &[&[u8]]) -> usize {
    reads.iter().enumerate().max_by_key(|(_, r)| r.len()).map(|(i, _)| i).unwrap()
}

fn auto_select(reads: &[&[u8]]) -> Result<usize, PoaError> {
    let n = reads.len();
    // Need at least 2 reads to agree on a terminal k-mer; 30 % gives good
    // sensitivity at depth 10 while suppressing noise from damaged reads.
    let threshold = ((n as f64 * ANCHOR_FRAC) as usize).max(2).min(n);

    // Collect deduplicated terminal k-mer sets per read (one set per end).
    let left_sets: Vec<HashSet<u64>> = reads.iter().map(|r| terminal_kmers(r, true)).collect();
    let right_sets: Vec<HashSet<u64>> = reads.iter().map(|r| terminal_kmers(r, false)).collect();

    // Global frequency tables: how many reads contribute each k-mer at each end.
    let mut left_freq: HashMap<u64, usize> = HashMap::new();
    let mut right_freq: HashMap<u64, usize> = HashMap::new();
    for s in &left_sets  { for &h in s { *left_freq.entry(h).or_insert(0) += 1; } }
    for s in &right_sets { for &h in s { *right_freq.entry(h).or_insert(0) += 1; } }

    // End-specific anchors: k-mers common at one end but rare at the other.
    // K-mers frequent at BOTH ends are interior repetitive sequence (e.g. the
    // repeat unit itself) and must be discarded — otherwise partial reads whose
    // interior happens to resemble the opposite flank look spuriously spanning.
    let is_left_anchor = |h: &u64| -> bool {
        *left_freq.get(h).unwrap_or(&0) >= threshold
            && *right_freq.get(h).unwrap_or(&0) < threshold
    };
    let is_right_anchor = |h: &u64| -> bool {
        *right_freq.get(h).unwrap_or(&0) >= threshold
            && *left_freq.get(h).unwrap_or(&0) < threshold
    };

    let left_ok: Vec<bool> = left_sets.iter()
        .map(|s| s.iter().any(is_left_anchor))
        .collect();
    let right_ok: Vec<bool> = right_sets.iter()
        .map(|s| s.iter().any(is_right_anchor))
        .collect();

    // Spanning candidates: end-specific anchor on both sides.
    let candidates: Vec<usize> = (0..n).filter(|&i| left_ok[i] && right_ok[i]).collect();
    if !candidates.is_empty() {
        // Shortest spanning read: all longer reads will INSERT, maximising
        // explicit bubble arms in arm_sequences.
        return Ok(*candidates.iter().min_by_key(|&&i| reads[i].len()).unwrap());
    }

    // Two-cluster detection: distinct left-only and right-only groups with no
    // bridging read.  Caller should use bridged_consensus.
    let left_only = (0..n).filter(|&i| left_ok[i] && !right_ok[i]).count();
    let right_only = (0..n).filter(|&i| !left_ok[i] && right_ok[i]).count();
    if left_only > 0 && right_only > 0 {
        return Err(PoaError::NoSpanningReads { left_depth: left_only, right_depth: right_only });
    }

    // Fallback: no cluster structure detected (reads too short, highly
    // repetitive flanks, or uniform noise).  Longest read is a safer seed
    // than shortest when spanning status is unknown.
    Ok(longest(reads))
}

/// Deduplicated set of k-mer hashes from the terminal region of `read`.
///
/// `left = true` → first [`TERM_LEN`] bases; `left = false` → last [`TERM_LEN`] bases.
/// Returns an empty set for reads shorter than [`TERM_K`].
fn terminal_kmers(read: &[u8], left: bool) -> HashSet<u64> {
    if read.len() < TERM_K {
        return HashSet::new();
    }
    let term_len = TERM_LEN.min(read.len());
    let slice = if left {
        &read[..term_len]
    } else {
        &read[read.len() - term_len..]
    };
    kmers_of(slice)
}

/// Rolling 2-bit k-mer hash over `seq`, returning a deduplicated set of hashes.
/// Non-ACGT bytes reset the accumulator; k-mers spanning a reset are skipped.
fn kmers_of(seq: &[u8]) -> HashSet<u64> {
    let mut out = HashSet::new();
    if seq.len() < TERM_K {
        return out;
    }
    let mask = (1u64 << (2 * TERM_K)) - 1;
    let mut h: u64 = 0;
    let mut valid: usize = 0;
    for &b in seq {
        let bits = match b {
            b'A' | b'a' => { valid += 1; 0u64 }
            b'C' | b'c' => { valid += 1; 1 }
            b'G' | b'g' => { valid += 1; 2 }
            b'T' | b't' => { valid += 1; 3 }
            _ => { valid = 0; h = 0; continue; }
        };
        h = ((h << 2) | bits) & mask;
        if valid >= TERM_K {
            out.insert(h);
        }
    }
    out
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Explicit / Shortest ───────────────────────────────────────────────────

    #[test]
    fn explicit_valid() {
        let reads: &[&[u8]] = &[b"AAA", b"CCC", b"GGG"];
        assert_eq!(select_seed(reads, &SeedSelection::Explicit(2)).unwrap(), 2);
    }

    #[test]
    fn explicit_out_of_bounds() {
        let reads: &[&[u8]] = &[b"AAA", b"CCC"];
        assert!(matches!(
            select_seed(reads, &SeedSelection::Explicit(5)),
            Err(PoaError::SeedOutOfBounds { index: 5, len: 2 })
        ));
    }

    #[test]
    fn shortest_picks_min_len() {
        let reads: &[&[u8]] = &[b"AAAAAAA", b"AAA", b"AAAAA"];
        assert_eq!(select_seed(reads, &SeedSelection::Shortest).unwrap(), 1);
    }

    #[test]
    fn empty_input_errors() {
        let reads: &[&[u8]] = &[];
        assert!(matches!(select_seed(reads, &SeedSelection::Auto), Err(PoaError::EmptyInput)));
        assert!(matches!(select_seed(reads, &SeedSelection::Shortest), Err(PoaError::EmptyInput)));
        assert!(matches!(
            select_seed(reads, &SeedSelection::Explicit(0)),
            Err(PoaError::EmptyInput)
        ));
    }

    // ── Auto: spanning detection ──────────────────────────────────────────────

    // 30-bp flanks: long enough to yield multiple 15-mers, short enough that
    // reads stay manageable.
    const LEFT_FLANK: &[u8]  = b"ACGTACGTACGTACGTACGTACGTACGTAC"; // 30 bp
    const RIGHT_FLANK: &[u8] = b"TGCATGCATGCATGCATGCATGCATGCATG"; // 30 bp
    const REPEAT: &[u8]      = b"CAT";                              // 3 bp unit

    // Reads must be long enough that the TERM_LEN=50 window at each end does
    // not reach the opposite flank.  Required: read_len > 2 × TERM_LEN = 100 bp.
    // With LEFT(30) + RIGHT(30) = 60 bp of flanks, need ≥ 41 bp of repeat
    // in the middle (≥14 CAT units).

    fn make_spanning(n_repeat: usize) -> Vec<u8> {
        let mut v = LEFT_FLANK.to_vec();
        for _ in 0..n_repeat { v.extend_from_slice(REPEAT); }
        v.extend_from_slice(RIGHT_FLANK);
        v
    }

    fn make_left_only(n_repeat: usize) -> Vec<u8> {
        // Right terminal = pure REPEAT, no RIGHT_FLANK k-mers.
        // Need read_len > TERM_LEN + LEFT_FLANK_len = 50 + 30 = 80 bp, so
        // the left flank does not appear in the right terminal window.
        let mut v = LEFT_FLANK.to_vec();
        for _ in 0..n_repeat { v.extend_from_slice(REPEAT); }
        v
    }

    fn make_right_only(n_repeat: usize) -> Vec<u8> {
        // Left terminal = pure REPEAT, no LEFT_FLANK k-mers.
        let mut v: Vec<u8> = Vec::new();
        for _ in 0..n_repeat { v.extend_from_slice(REPEAT); }
        v.extend_from_slice(RIGHT_FLANK);
        v
    }

    #[test]
    fn auto_picks_shortest_spanning() {
        // Three spanning reads of different lengths; shortest should be selected.
        // All > 100 bp so the TERM_LEN windows do not overlap.
        let r_short  = make_spanning(15); // 30+45+30 = 105 bp
        let r_medium = make_spanning(20); // 120 bp
        let r_long   = make_spanning(25); // 135 bp
        // Place r_short at index 2 to verify index correctness.
        let reads: Vec<&[u8]> = vec![&r_medium, &r_long, &r_short];
        let idx = select_seed(&reads, &SeedSelection::Auto).unwrap();
        assert_eq!(idx, 2, "expected shortest spanning read at index 2");
    }

    #[test]
    fn auto_ignores_partials_picks_spanning() {
        // Mix: left-only and right-only partials alongside two spanning reads.
        // Auto should select the shorter spanning read as seed.
        let s_short = make_spanning(15); // 105 bp — spanning, shortest
        let s_long  = make_spanning(20); // 120 bp — spanning
        // 30 CAT units = 90 bp of repeat → total left_only = 30+90 = 120 bp > 80 ✓
        let left1 = make_left_only(30);
        let left2 = make_left_only(30);
        let left3 = make_left_only(30);
        let right1 = make_right_only(30);
        let right2 = make_right_only(30);

        // n=7, threshold = max(2, floor(7×0.3)) = 2.
        // LEFT_FLANK k-mers: in left terminal of s_short, s_long, left1/2/3 → freq=5 ≥ 2 ✓
        // RIGHT_FLANK k-mers: in right terminal of s_short, s_long, right1/2 → freq=4 ≥ 2 ✓
        // REPEAT k-mers appear at BOTH ends → end-specific filter discards them.
        // Spanning candidates: s_short (idx 0) and s_long (idx 1).
        // Shortest spanning = s_short at index 0.
        let reads: Vec<&[u8]> = vec![&s_short, &s_long, &left1, &left2, &left3, &right1, &right2];
        let idx = select_seed(&reads, &SeedSelection::Auto).unwrap();
        assert_eq!(reads[idx].len(), s_short.len(), "expected the shortest spanning read");
    }

    #[test]
    fn auto_two_cluster_errors() {
        // All reads are either left-only or right-only, no spanning read.
        // Auto should return NoSpanningReads.
        let left1  = make_left_only(30);
        let left2  = make_left_only(30);
        let left3  = make_left_only(30);
        let right1 = make_right_only(30);
        let right2 = make_right_only(30);
        let right3 = make_right_only(30);

        // n=6, threshold = max(2, floor(6×0.3)) = 2.
        let reads: Vec<&[u8]> = vec![&left1, &left2, &left3, &right1, &right2, &right3];
        match select_seed(&reads, &SeedSelection::Auto) {
            Err(PoaError::NoSpanningReads { left_depth, right_depth }) => {
                assert_eq!(left_depth, 3);
                assert_eq!(right_depth, 3);
            }
            other => panic!("expected NoSpanningReads, got {other:?}"),
        }
    }

    #[test]
    fn auto_fallback_to_longest_when_no_cluster() {
        // Reads too short to extract TERM_K=15 k-mers → no anchors, no clusters.
        // Fallback should return the longest read.
        let reads: Vec<&[u8]> = vec![b"ACGTACG", b"ACGT", b"ACGTACGTACG"];
        // All reads < 15 bp → terminal_kmers returns empty → left/right_freq empty.
        // No left_ok, no right_ok → no clusters → fallback to longest (index 2).
        let idx = select_seed(&reads, &SeedSelection::Auto).unwrap();
        assert_eq!(idx, 2, "expected longest read as fallback");
    }
}
