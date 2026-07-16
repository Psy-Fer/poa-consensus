//! Regression fixture for a severe, confirmed bug found incidentally during
//! this session's seed-length-sensitivity batch generalization check (Part
//! A of that investigation): under the CLI's own **default** config
//! (`band_width: 50, adaptive_band: true`), a specific read set collapses
//! the entire consensus to a single base, silently -- no error, no warning.
//!
//! # Fixture provenance
//!
//! `fixtures/adaptive_band_collapse_gaa55.fa` (28 reads, GAA×55, depth 20,
//! ONT R10 error model) is the exact read set from one of the 74 fresh
//! random draws in the seed-sensitivity batch check (`bench/validate.py`'s
//! own synthetic pbsim3 simulation pipeline, `pbsim3 seed=1022`) -- no real
//! patient data, same discipline as `tests/band_width_zero_unbanded.rs` and
//! every other synthetic fixture in this directory. 3 of 74 draws in that
//! batch (~4%) hit this exact collapse; a follow-up scope check found it
//! recurs at ~13% (2/15) on longer repeats (`CAG×200`, ONT R10) and not at
//! all (0/15) under the HiFi error model or on a short, "clean" control
//! (`CAG×20`) -- see this file's root-cause trace below for why longer
//! repeats and higher error rates are the aggravating factors.
//!
//! # Root cause, traced precisely (not guessed)
//!
//! `select_seed(Auto)` (correctly, by its own documented design) picks the
//! shortest spanning read as seed -- 292bp here, against a population
//! ranging up to 434bp (a 142bp spread, driven by ordinary ONT R10 indel
//! error accumulated over a 55-unit homogeneous repeat). `align()`'s banded
//! DP, at the default `band_width=50`/`adaptive_band=true` (effective band
//! ≈ 50, since `band_width` floors the adaptive formula's own ≈14), cannot
//! track the diagonal shift needed once a read is more than ~50bp longer
//! than the current graph. Confirmed by direct sweep: **the exact
//! transition point is 73bp** (the precise length difference between the
//! 292bp seed and one specific 365bp read in this fixture) -- any band
//! width ≥ 73 aligns that read correctly (Match≈289, Insert≈76, Delete≈3);
//! any band width < 73, including the adaptive formula's own unfloored
//! value, produces `align_read_ops` returning an **empty `Vec<AlignOp>`**
//! (confirmed directly via the public `PoaGraph::align_read_ops` API, not
//! inferred) -- not an error, a silently empty alignment.
//!
//! Traced to the exact code location: `align()`'s "Find best terminal cell
//! at column l" step searches every graph node's own banded row-window for
//! a cell at query column `l` (the full read length); when the band is too
//! narrow for the diagonal shift, no cell anywhere reaches column `l`, so
//! that search returns `None`. The fallback for `None`,
//! `(best_t, best_state) = (0, State::M)`, is not a valid starting point at
//! all (node 0's own row-window is nowhere near column `l`), so the very
//! first traceback lookup finds an `UNSET` cell and the loop breaks having
//! pushed zero ops. `add_to_graph` on an empty ops list is a complete
//! no-op: no coverage, no new nodes, nothing recorded for that read. Once
//! every non-seed read (all longer than this seed by more than the band)
//! independently hits this, `PoaGraph` ends up holding only the seed's own
//! linear 292-node chain, each node at `coverage=1` -- confirmed directly
//! via `PoaGraph::stats()` (`single_support_fraction: 1.0`,
//! `coverage_mean: 1.0`) and `PoaGraph::node_count()` (292, exactly the
//! seed's own length). `consensus()`'s boundary trim, given every node
//! below the depth-derived floor, then collapses to a 1-base output.
//!
//! **`PoaError::BandTooNarrow` is never actually raised anywhere in
//! `src/graph.rs`** (confirmed by exhaustive `grep`) despite being fully
//! defined in `src/error.rs`, documented in `src/lib.rs`'s crate docs, and
//! having a ready-made, user-facing hint in `src/main.rs`'s
//! `explain_error` ("try `--adaptive-band`, or `--band-width {required}`").
//! `TODO.md`'s "Smart retry (3-pass)" checkbox describes an automatic
//! `BandTooNarrow` → retry-wider → retry-unbanded loop as already shipped;
//! no such retry exists anywhere in `add_read()` or elsewhere in the
//! current code -- this appears to be stale documentation (the same
//! documentation-vs-reality drift already found and flagged elsewhere this
//! session, e.g. Known Bug #6's stale `#[ignore]` and `phasing_groups`'s
//! stale doc comment), not a mechanism that regressed.
//!
//! # Why the shipped safety net (`consensus_fit_scored`) does not catch this
//!
//! `stats().single_support_fraction` on the collapsed graph is `1.0`, well
//! above the `0.3` trigger `consensus_fit_scored`/`consensus_adaptive` both
//! use -- the trigger *does* fire. But `seed_sensitivity_retry` still
//! returns the broken 1-base candidate (`action: PassThrough`), for two
//! independently-confirmed reasons:
//!
//! 1. **`analysis::consensus_fit` itself hits the identical bug when
//!    scoring the degenerate candidate.** It builds a throwaway graph
//!    seeded on the 1-base candidate sequence and aligns every original
//!    read against it -- and aligning a 348bp read against a 1-node graph
//!    is exactly the same "band too narrow for the length gap" condition,
//!    so it *also* returns empty ops (confirmed directly). Zero
//!    Insert+Delete ops against zero reads that "fit" scores as `0.0` --
//!    literally the best possible score `consensus_fit` can produce --
//!    for what is actually the worst candidate. The one candidate that
//!    *is* genuinely correct in this fixture (`AlternateSeedRetry`,
//!    re-seeded on the population's median-length read, which does not
//!    trigger the band-vs-length-gap condition) scores `0.039254` -- a
//!    real, non-zero, *worse*-looking number than the broken candidates'
//!    fake `0.0`, so it loses the comparison.
//! 2. **`consensus_fit_scored` has no truncation-ratio check at all.**
//!    `consensus_adaptive` (the fuller two-pass entry point) does, and
//!    confirmed directly on this exact fixture+seed to catch and fully fix
//!    this: `action: TruncationRetry { recovered: true }`, full 364bp
//!    output. The CLI (`src/main.rs`) deliberately calls
//!    `consensus_fit_scored` instead of `consensus_adaptive`, specifically
//!    to avoid `consensus_adaptive`'s multi-allele bubble detection firing
//!    unexpectedly on a het SNP -- an unrelated, deliberate tradeoff whose
//!    side effect is silently losing this specific protection too.
//!
//! # Fix implemented
//!
//! (a) At the `terminal_best == None` branch in `align()`, return
//! `Err(PoaError::BandTooNarrow { configured, required })` instead of the
//! bogus `(0, State::M)` substitution -- restores the *already-designed*,
//! already-documented, already-user-facing-hinted error path instead of
//! silent corruption. `required` is estimated from `best_j_per_t` (already
//! tracked per row for the diagonal-skip machinery): the row that got
//! closest to query column `l` fell short by `l - max(best_j_per_t)`, so
//! widening the margin by that amount is the natural single-shot estimate.
//! (b) A new `align_with_retry` wrapper (used by `add_read`,
//! `align_read_ops`, and `align_read_ops_unbanded` -- so both graph
//! construction *and* `analysis::consensus_fit`'s own internal
//! throwaway-graph scoring, which calls `align_read_ops`, are covered by
//! the same fix) catches `BandTooNarrow` and retries: pass 2 at the
//! error's own `required` estimate (not trusted as a proof of sufficiency,
//! since it is a heuristic); pass 3, only if pass 2 *also* errors, fully
//! unbanded (`band_width = 0, adaptive_band = false`), which is
//! mathematically guaranteed sufficient. See CHANGELOG for the full trace
//! and validated impact.

use poa_consensus::{AlignmentMode, PoaConfig, SeedSelection, consensus_fit_scored, select_seed};

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
fn adaptive_band_default_config_no_longer_collapses_gaa55() {
    let reads = parse_fasta(include_str!("fixtures/adaptive_band_collapse_gaa55.fa"));
    let slices: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();

    // The CLI's own default config: band_width=50, adaptive_band=true,
    // SemiGlobal (the crate's own default alignment mode since v0.2.0).
    let cfg = PoaConfig {
        band_width: 50,
        adaptive_band: true,
        alignment_mode: AlignmentMode::SemiGlobal,
        min_reads: 2,
        ..PoaConfig::default()
    };

    let seed_idx = select_seed(&slices, &SeedSelection::Auto).unwrap();
    // Confirmed this is the pathological seed: 292bp against a population
    // extending to 434bp, a gap larger than the effective band.
    assert_eq!(
        slices[seed_idx].len(),
        292,
        "fixture drifted -- this test assumes seed_idx picks the known 292bp read"
    );

    let result = consensus_fit_scored(&slices, seed_idx, &cfg).unwrap();

    // Full recovery, not just "not collapsed": the true GAA×55 window is
    // ~364bp (full extracted-read width including flanks, confirmed via
    // this exact fixture's own median_input_read_len). A tolerance of a
    // few bp accommodates ordinary per-read indel noise in the consensus's
    // own boundary, the same way every other accuracy test in this suite
    // does -- this is not a loose "at least not 1bp" check.
    let len = result.consensuses[0].sequence.len();
    assert!(
        len.abs_diff(364) <= 10,
        "expected the BandTooNarrow retry to recover a ~364bp consensus \
         (the true GAA×55 window), got {len}bp"
    );
}
