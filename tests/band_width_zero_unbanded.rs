//! Regression for `PoaConfig::band_width = 0` not actually being unbanded.
//!
//! `band_width` is documented as "0 = unbanded (full NW over DAG)", but
//! `align()`'s spine-margin computation used to silently fall back to a fixed
//! 50-column window whenever `band_width = 0` and `adaptive_band = false` --
//! identical in practice to `band_width = 50` -- so the documented
//! correctness-baseline escape hatch was unreachable through the public API.
//!
//! This was found investigating a real long-repeat locus where the CLI
//! default (banded, adaptive on) silently fabricated an extra copy of a
//! tandem-repeat unit beyond what any input read, abPOA, or SPOA supported --
//! a concrete instance of the documented "silent wrong alignment on narrow
//! band in repetitive regions" failure mode. That investigation used real
//! patient sequencing data that cannot be redistributed with this crate, so
//! this test instead uses a synthetic tandem-repeat expansion: a seed with no
//! expansion, and reads carrying a large expansion the seed doesn't have --
//! the same shape of problem (a majority-supported variant requiring a wide
//! DP window to place correctly), without needing the exact real-world
//! diagonal-ambiguity mechanism to reproduce a clear regression signal.
//!
//! A fixed `band_width = 50` could not place the expansion at all when this
//! test was first written (the graph never recovered a coherent consensus).
//! That is no longer true as of the `BandTooNarrow`/`align_with_retry` fix
//! (see `tests/adaptive_band_collapse.rs`): `align()` now detects exactly
//! this "band too narrow to reach the end of the read" condition and
//! retries with a wider band automatically, so `band_width = 50` now also
//! reconstructs the expansion exactly, the same as genuinely unbanded DP
//! (`band_width = 0`, `adaptive_band = false`). Both are tested below.

use poa_consensus::{AlignmentMode, PoaConfig, PoaGraph};

fn build_graph(band_width: usize, seed: &[u8]) -> PoaGraph {
    let cfg = PoaConfig {
        band_width,
        adaptive_band: false,
        min_reads: 2,
        alignment_mode: AlignmentMode::SemiGlobal,
        ..PoaConfig::default()
    };
    PoaGraph::new(seed, cfg).unwrap()
}

/// A short tandem repeat (60 CAG units) and a majority-supported expansion
/// (+70 more CAG units inserted mid-repeat) not present in the seed.
fn repeat_and_expansion() -> (Vec<u8>, Vec<u8>) {
    let base: Vec<u8> = "CAG".repeat(60).into_bytes(); // 180bp, no expansion
    let extra: Vec<u8> = "CAG".repeat(70).into_bytes(); // +210bp expansion
    let half = base.len() / 2;
    let mut expanded = base[..half].to_vec();
    expanded.extend_from_slice(&extra);
    expanded.extend_from_slice(&base[half..]);
    (base, expanded)
}

#[test]
fn band_width_zero_reconstructs_expansion_exactly() {
    let (base, expanded) = repeat_and_expansion();
    let mut graph = build_graph(0, &base);
    for _ in 0..6 {
        graph.add_read(&expanded).unwrap();
    }
    let result = graph.consensus().unwrap();
    assert_eq!(
        result.sequence,
        expanded,
        "band_width=0 (adaptive_band=false) must reconstruct the full \
         expansion via a genuinely unbanded window; got {} bp, expected {} bp",
        result.sequence.len(),
        expanded.len()
    );
}

#[test]
fn fixed_band_width_50_now_recovers_via_band_too_narrow_retry() {
    // This test originally documented the opposite: a merely-50-wide fixed
    // band could not correctly place a +210bp expansion the seed lacks
    // entirely, and asserted `assert_ne!` against the correct answer with a
    // note that read "if this now passes, band_width=50 alone has become
    // sufficient and this contrast test should be revisited." It now passes
    // -- revisited here, per that note, rather than left stale.
    //
    // The reason is not that band_width=50 alone became wide enough on its
    // own; it is that `align()` now correctly returns
    // `Err(PoaError::BandTooNarrow)` instead of silently returning an
    // empty/degenerate alignment when a read needs a wider window than
    // configured (see `tests/adaptive_band_collapse.rs` for the severe
    // silent-corruption bug this closes), and a new retry wrapper
    // (`align_with_retry`, used by `add_read`/`align_read_ops`/
    // `align_read_ops_unbanded`) catches that error and retries with a
    // wider band -- escalating to fully unbanded if needed -- so every
    // caller gets a correct answer transparently instead of either a silent
    // wrong one or a hard failure.
    let (base, expanded) = repeat_and_expansion();
    let mut graph = build_graph(50, &base);
    for _ in 0..6 {
        graph.add_read(&expanded).unwrap();
    }
    let result = graph.consensus().unwrap();
    assert_eq!(
        result.sequence, expanded,
        "band_width=50 should now reconstruct the expansion exactly via the \
         BandTooNarrow retry, the same as genuinely unbanded (band_width=0)"
    );
}
