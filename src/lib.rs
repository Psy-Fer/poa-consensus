//! # poa-consensus
//!
//! A pure-Rust banded Partial Order Alignment (POA) library for building a
//! consensus sequence from a set of reads.
//!
//! POA builds a directed acyclic graph (DAG) from the reads, aligns each read
//! into the graph using dynamic programming, and extracts a consensus by
//! following the heaviest (most-supported) path through the graph.
//!
//! ---
//!
//! ## Affine gap penalties
//!
//! The alignment scoring uses **affine gap penalties**. This section explains
//! what that means and why it matters, because the standard explanation in most
//! papers is needlessly opaque.
//!
//! ### The problem with linear gaps
//!
//! A linear gap model charges a flat penalty per gap base, say `-1` for every
//! missing or inserted base. This means a single deletion of 4 bases and four
//! separate single-base deletions cost exactly the same:
//!
//! ```text
//! One 4-base deletion:       -1 -1 -1 -1   = -4
//! Four 1-base deletions:     -1 -1 -1 -1   = -4
//! ```
//!
//! In biology these are very different events. A 4-base deletion is one
//! mutation: a single copying mistake that dropped 4 bases together. Four
//! separate 1-base deletions are four independent mutations. The linear model
//! treats them identically, so it cannot prefer the more plausible explanation.
//!
//! ### Affine gaps: one price to start, a smaller price to continue
//!
//! Affine scoring uses two parameters:
//!
//! - **`gap_open`**: a one-time penalty paid at the moment a gap starts.
//!   This represents the cost of the mutation event itself.
//! - **`gap_extend`**: a per-base penalty paid for every base inside the gap.
//!   This represents the cost of each missing or inserted base.
//!
//! A gap of length `k` costs: `gap_open + k * gap_extend`
//!
//! ### Example
//!
//! Using `gap_open = -2` and `gap_extend = -1`:
//!
//! ```text
//! One gap of length 4:
//!     gap_open + 4 * gap_extend = -2 + (-4) = -6
//!
//! Four gaps of length 1:
//!     4 * (gap_open + 1 * gap_extend) = 4 * (-3) = -12
//! ```
//!
//! The single 4-base deletion costs -6. Four separate 1-base deletions cost
//! -12. The aligner now strongly prefers the single-event explanation, which
//! matches biological reality.
//!
//! With linear gaps (gap = -1) both scenarios cost -4 and the aligner cannot
//! distinguish them at all.
//!
//! ### Worked alignment example
//!
//! Say we are aligning a read to a reference with a 4-base deletion:
//!
//! ```text
//! Reference: A C G T T T T T A C G T
//! Read:      A C G T - - - - A C G T
//! ```
//!
//! Scoring (match = +1, mismatch = -1, gap_open = -2, gap_extend = -1):
//!
//! ```text
//! 8 matching bases:  8 * (+1)           =  +8
//! 1 gap of length 4: -2 + 4 * (-1)     =  -6
//! Total:                                =  +2
//! ```
//!
//! Now consider a different alignment that avoids the gap by accepting 4
//! mismatches instead:
//!
//! ```text
//! Reference: A C G T T T T T A C G T
//! Read:      A C G T X X X X A C G T   (X = mismatch)
//! ```
//!
//! ```text
//! 8 matching bases:   8 * (+1)          =  +8
//! 4 mismatches:       4 * (-1)          =  -4
//! Total:                                =  +4
//! ```
//!
//! Here the mismatched alignment scores higher (+4 vs +2), so the aligner
//! would choose it. Whether that is the right call depends on the biology: if
//! the read truly has a 4-base deletion, these parameters would produce the
//! wrong alignment. Tuning `gap_open` and `gap_extend` controls this
//! trade-off. A more negative `gap_open` makes the aligner more willing to
//! open a gap rather than accept mismatches; a less negative `gap_extend`
//! makes longer gaps cheaper relative to short ones.
//!
//! ### Why individual alignment ambiguity matters less in POA
//!
//! In a single pairwise alignment there is often no way to know whether a gap
//! or a mismatch is the correct interpretation. POA largely sidesteps this
//! problem because it aligns many reads into the same graph and then extracts
//! a consensus, rather than making a definitive call on any one alignment in
//! isolation.
//!
//! If a gap in one read reflects a real deletion, most other reads will carry
//! the same gap. The corresponding graph edges accumulate high weight, and the
//! heaviest-path consensus will include the deletion. If the gap is a
//! sequencing error, other reads will traverse the same graph node as a match,
//! and the consensus will correct it.
//!
//! This means parameter choices affect alignment quality at the margins, but
//! the coverage threshold and heaviest-path extraction together resolve most
//! of the ambiguity that would be fatal to a single pairwise alignment. The
//! defaults (`gap_open = -2`, `gap_extend = -1`) are calibrated for Oxford
//! Nanopore and PacBio HiFi reads, where indels are the dominant error type.
//!
//! ### How affine scoring changes the dynamic programming
//!
//! Standard Needleman-Wunsch uses a single DP table. Affine scoring requires
//! three tables, one per "state" the aligner can be in at any position:
//!
//! - **`M[i][j]`**: best score where query position `i` aligns to graph node
//!   `j` as a match or mismatch.
//! - **`I[i][j]`**: best score where query position `i` is an insertion (a
//!   base in the read with no corresponding node in the graph).
//! - **`D[i][j]`**: best score where graph node `j` is a deletion (a node in
//!   the graph skipped by the read).
//!
//! The recurrences are:
//!
//! ```text
//! M[i][j] = max(M[i-1][j-1], I[i-1][j-1], D[i-1][j-1]) + score(query[i], node[j])
//!
//! I[i][j] = max(M[i][j-1] + gap_open + gap_extend,
//!               I[i][j-1] + gap_extend)
//!
//! D[i][j] = max(M[i-1][j] + gap_open + gap_extend,
//!               D[i-1][j] + gap_extend)
//! ```
//!
//! The key insight in `I` and `D`: transitioning from `M` to a gap state pays
//! `gap_open + gap_extend` (opening cost plus first base). Staying in the same
//! gap state pays only `gap_extend` (continuing an existing gap). The aligner
//! tracks which state it is in at every cell, so it knows whether to apply the
//! full open penalty or just the extend penalty.
//!
//! Traceback follows whichever state and predecessor gave the maximum score at
//! each cell, producing a sequence of Match, Insert, and Delete operations.
//!
//! For POA over a graph rather than a linear reference, `j` indexes a graph
//! node in topological order rather than a column in a matrix, but the
//! recurrence is identical.
//!
//! ### Choosing parameter values
//!
//! The defaults used in this crate (`gap_open = -2`, `gap_extend = -1`) work
//! well for Oxford Nanopore and PacBio HiFi reads at short tandem repeat loci.
//! As a rule of thumb:
//!
//! - Increase the magnitude of `gap_open` (e.g. `-4`) to make the aligner
//!   more reluctant to open gaps, preferring mismatches instead.
//! - Decrease the magnitude of `gap_extend` (e.g. `-0.5`, or use integer
//!   scaling) to make long gaps cheaper, useful when reads have systematic
//!   length variation.
//! - Setting `gap_open = 0` recovers linear gap behaviour where only
//!   `gap_extend` matters.

pub mod config;
pub mod error;
pub mod flank;
pub mod graph;
pub mod orient;
pub mod types;

#[cfg(test)]
mod tests;

pub use config::{AlignmentMode, ConsensusMode, PoaConfig};
pub use error::PoaError;
pub use flank::extract_flanked_region;
pub use graph::PoaGraph;
pub use orient::{auto_orient, orient_to_seed, reverse_complement};
pub use types::{Consensus, GraphStats, Strand};

// ── Internal helper ───────────────────────────────────────────────────────────

fn build_graph(reads: &[&[u8]], seed_idx: usize, config: PoaConfig) -> Result<PoaGraph, PoaError> {
    let mut graph = PoaGraph::new(reads[seed_idx], config)?;
    for (i, read) in reads.iter().enumerate() {
        if i != seed_idx {
            graph.add_read(read)?;
        }
    }
    Ok(graph)
}

fn validate(reads: &[&[u8]], seed_idx: usize) -> Result<(), PoaError> {
    if reads.is_empty() {
        return Err(PoaError::EmptyInput);
    }
    if seed_idx >= reads.len() {
        return Err(PoaError::SeedOutOfBounds {
            index: seed_idx,
            len: reads.len(),
        });
    }
    Ok(())
}

// ── Public convenience wrappers ───────────────────────────────────────────────

/// Build a single-allele consensus from `reads`.
///
/// `seed_idx` is the index of the read used to initialise the graph; choose a
/// median-length read for best results.
pub fn consensus(
    reads: &[&[u8]],
    seed_idx: usize,
    config: &PoaConfig,
) -> Result<Consensus, PoaError> {
    validate(reads, seed_idx)?;
    build_graph(reads, seed_idx, config.clone())?.consensus()
}

/// Build a multi-allele consensus from `reads`.
///
/// Returns one [`Consensus`] per detected allele.  If no heterozygous bubble is
/// found the result is a single-element `Vec` equivalent to calling [`consensus`].
pub fn consensus_multi(
    reads: &[&[u8]],
    seed_idx: usize,
    config: &PoaConfig,
) -> Result<Vec<Consensus>, PoaError> {
    validate(reads, seed_idx)?;
    build_graph(reads, seed_idx, config.clone())?.consensus_multi()
}

/// Two-pass adaptive consensus.
///
/// **Pass 1** builds a graph with the supplied config and computes [`GraphStats`].
/// **Pass 2** is selected by inspecting those stats:
///
/// | Condition | Pass-2 action |
/// |---|---|
/// | 1-3 bubbles, minority arm ≥ `min_allele_freq × n` | `consensus_multi` on pass-1 graph |
/// | `single_support_fraction > 0.3` | Tighten `min_coverage_fraction` to ≥ 0.6, rebuild |
/// | Coverage CV > 1.5 and mode is `Global` | Switch to `SemiGlobal`, rebuild |
/// | Otherwise | Return pass-1 single consensus; no rebuild |
///
/// Always returns a `Vec<Consensus>`: one element for single-allele outcomes,
/// two for diploid.
pub fn consensus_adaptive(
    reads: &[&[u8]],
    seed_idx: usize,
    config: &PoaConfig,
) -> Result<Vec<Consensus>, PoaError> {
    validate(reads, seed_idx)?;

    // ── Pass 1 ───────────────────────────────────────────────────────────────
    let graph = build_graph(reads, seed_idx, config.clone())?;
    let stats = graph.stats();

    // ── Decision ─────────────────────────────────────────────────────────────
    let n = reads.len();
    let allele_threshold = (n as f64 * config.min_allele_freq).ceil() as usize;

    // Multi-allele: few bubbles with a well-supported minority arm.
    // Re-use the pass-1 graph; consensus_multi rebuilds per-allele sub-graphs.
    if stats.bubble_count >= 1
        && stats.bubble_count <= 3
        && stats.max_bubble_depth >= allele_threshold
    {
        return graph.consensus_multi();
    }

    // Noisy input: high fraction of singleton-supported nodes.
    // Tighten the coverage threshold and rebuild.
    if stats.single_support_fraction > 0.3 {
        let mut cfg2 = config.clone();
        cfg2.min_coverage_fraction = cfg2.min_coverage_fraction.max(0.6);
        return Ok(vec![build_graph(reads, seed_idx, cfg2)?.consensus()?]);
    }

    // Uneven boundary coverage (high CV) suggests partial reads in global mode.
    // Switch to semi-global alignment and rebuild.
    let cv = if stats.coverage_mean > 0.0 {
        stats.coverage_variance.sqrt() / stats.coverage_mean
    } else {
        0.0
    };
    if cv > 1.5 && config.alignment_mode == AlignmentMode::Global {
        let mut cfg2 = config.clone();
        cfg2.alignment_mode = AlignmentMode::SemiGlobal;
        return Ok(vec![build_graph(reads, seed_idx, cfg2)?.consensus()?]);
    }

    // No second pass needed.
    Ok(vec![graph.consensus()?])
}
