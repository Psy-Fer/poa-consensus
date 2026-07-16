//! Regression fixture for a severe bug discovered while validating the
//! `BandTooNarrow` fix (see `tests/adaptive_band_collapse.rs`): once
//! `align_with_retry` started actually widening the band for reads that
//! needed it, a *different* read set (this fixture) hit an OOM/crash --
//! "memory allocation of 2147483648 bytes failed" (and, with a looser
//! `ulimit`, 4294967296, then 8589934592 -- doubling each time, the
//! signature of an unbounded `Vec` growth loop, not a single oversized
//! allocation).
//!
//! # Root cause, traced precisely
//!
//! The crash was not in `align()`'s own traceback (which already has a
//! defensive bound, see `band_too_narrow_tests` in `src/graph.rs`) but in
//! `heaviest_path`'s predecessor-pointer traceback (`src/graph.rs`), called
//! from `add_read` to refresh the cached spine. A debug build with
//! `RUST_BACKTRACE=full` pinned the crash to
//! `heaviest_path -> RawVec::grow_one`, and a permanent-invariant check (the
//! `debug_assert_eq!` now in `add_read`, comparing `topological_order`'s
//! output length against `self.nodes.len()`) confirmed the graph had a
//! genuine cycle by the time this fixture's 14th read was processed: 310 of
//! 798 nodes were excluded from the topological order, i.e. stuck with
//! in-degree that Kahn's algorithm could never reduce to zero.
//!
//! Dumping every edge `add_to_graph` created for that read (gated on
//! `read_idx`) found the exact back edge: an ordinary `Match` op wired an
//! edge from a node at rank 471 to a node at rank 466 -- backwards. Tracing
//! *why* found that `align()`'s own traceback is provably rank-monotonic in
//! the node identities it names (confirmed empirically: dumping the full,
//! unmodified `ops` sequence showed no rank ever decreasing), but
//! `add_to_graph`'s single-base *substitution reuse* (`try_reuse_arm`, the
//! content-addressed fork-arm cache) can redirect `prev` to a different,
//! already-existing node -- and that substitute's rank is not necessarily
//! less than the rank of whatever node `align()`'s traceback visits *next*.
//! The very next, ordinary `Match`/`Delete` edge is created via
//! `increment_or_add_edge` with no rank check at all (it never needed one
//! before reuse existed), so it silently wired a real back edge once the
//! substitute's rank happened to exceed the next node's rank -- corrupting
//! the graph into a cycle that only surfaced later, when the spine was next
//! refreshed and `heaviest_path` walked straight into it.
//!
//! # Fix implemented
//!
//! `try_reuse_arm` gained two rank checks (both against the entry-snapshot
//! `rank_of` `add_read` already computes once per call, threaded through
//! `add_to_graph`): (1) the existing local check that `fork`'s rank is less
//! than the reused chain's first node's rank, and (2) the new
//! `reuse_would_create_back_edge`, which additionally rejects the reuse if
//! its *last* node's rank would not be less than the rank of the very next
//! `Match`/`Delete` reference in the read's own remaining ops -- exactly the
//! edge that would otherwise go backward. Rejecting falls through to
//! ordinary fresh-node creation, the same safe fallback `try_reuse_arm`
//! already used for its other failure cases. As a second, independent line
//! of defense (not a substitute for the real fix), `add_read` now carries a
//! permanent, zero-cost-in-release `debug_assert_eq!` confirming the graph
//! is still a DAG on every call, so any future regression of this class
//! fails loudly and immediately in debug/test builds rather than hanging or
//! OOMing downstream.
//!
//! # Fixture provenance
//!
//! `fixtures/reuse_back_edge_cag50.fa` (22 reads, CAG×50, HiFi error model)
//! is a `bench/validate.py`-style synthetic pbsim3 simulation (seed 2001,
//! from this session's own batch-severity re-check harness) -- no real
//! patient data, same discipline as every other fixture in this directory.

use poa_consensus::{AlignmentMode, PoaConfig, consensus};

fn parse_fasta(s: &str) -> Vec<Vec<u8>> {
    let mut reads: Vec<Vec<u8>> = Vec::new();
    let mut current: Vec<u8> = Vec::new();
    for line in s.lines() {
        if line.starts_with('>') {
            if !current.is_empty() {
                reads.push(std::mem::take(&mut current));
            }
        } else {
            current.extend_from_slice(line.trim().as_bytes());
        }
    }
    if !current.is_empty() {
        reads.push(current);
    }
    reads
}

#[test]
fn reuse_back_edge_no_longer_corrupts_graph_into_cycle() {
    let reads = parse_fasta(include_str!("fixtures/reuse_back_edge_cag50.fa"));
    let slices: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();

    assert_eq!(
        slices.len(),
        22,
        "fixture drifted -- this test assumes exactly 22 reads"
    );

    // The exact config and seed that reproduced the OOM: the CLI's own
    // default band, seeded on the population's median-length read (357bp)
    // rather than the shortest, which is what originally drove the
    // substitution-reuse's redirected node to a higher rank than the next
    // node in that read's own traceback.
    let median_idx = 16usize;
    assert_eq!(
        slices[median_idx].len(),
        357,
        "fixture drifted -- this test assumes seed index 16 is the known 357bp read"
    );

    let cfg = PoaConfig {
        band_width: 50,
        adaptive_band: true,
        alignment_mode: AlignmentMode::SemiGlobal,
        min_reads: 2,
        ..PoaConfig::default()
    };

    // Before the fix, this call either hung consuming unbounded memory or
    // was killed by the OS (OOM) partway through graph construction, well
    // before ever returning. Completing at all, in bounded time and memory,
    // is the primary assertion here.
    let result = consensus(&slices, median_idx, &cfg).unwrap();

    // Not just "didn't crash": the true CAG×50 window is comfortably within
    // a few tens of bp of the input reads' own lengths (a homogeneous
    // 50-unit CAG repeat plus non-repetitive flanks). A collapsed or
    // garbage consensus from a corrupted graph would look nothing like
    // this.
    let len = result.sequence.len();
    assert!(
        (250..=400).contains(&len),
        "expected a plausible CAG×50-window consensus length, got {len}bp"
    );
}
