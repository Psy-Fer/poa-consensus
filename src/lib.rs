//! # poa-consensus
//!
//! A pure-Rust banded Partial Order Alignment (POA) library for building a
//! consensus sequence from a set of reads.
//!
//! POA builds a directed acyclic graph (DAG) from the reads, aligns each
//! subsequent read into the graph using affine-gap dynamic programming, and
//! extracts a consensus by following the heaviest (most-supported) path through
//! the graph.  It handles length variation between reads naturally: inserts and
//! deletions create separate branches in the graph, and the heaviest-path
//! extraction resolves them by read support rather than by fixed rules.
//!
//! ## When to use this crate
//!
//! This crate is optimised for **short to medium reads** (50 bp to ~20 kb with
//! banded DP) where the graph stays small enough for in-memory DP — short
//! tandem repeat (STR) loci, amplicon consensus, per-locus nanopore or HiFi
//! read sets.
//!
//! The nearest alternative on crates.io is
//! [`poasta`](https://crates.io/crates/poasta) (Broad Institute, pure Rust,
//! gap-affine A\* alignment).  `poasta` excels at larger graphs such as
//! bacterial genes and HLA loci.  For STR reads where graphs are small and
//! throughput across many loci matters, a well-tuned banded DP is faster and
//! simpler.  Both crates are pure Rust; neither wraps a C library.
//!
//! ## Quick start
//!
//! ### Functional API (single call)
//!
//! ```rust
//! use poa_consensus::{consensus, consensus_multi, PoaConfig};
//!
//! let reads: Vec<&[u8]> = vec![
//!     b"CATCATCAT",
//!     b"CATCATCAT",
//!     b"CATCGTCAT",
//!     b"CATCATCAT",
//! ];
//!
//! // Single-allele consensus; seed_idx=0 seeds the graph with reads[0].
//! let result = consensus(&reads, 0, &PoaConfig::default())?;
//! println!("{}", String::from_utf8_lossy(&result.sequence));
//!
//! // Multi-allele: returns one Consensus per detected allele.
//! let alleles = consensus_multi(&reads, 0, &PoaConfig::default())?;
//! # Ok::<(), poa_consensus::PoaError>(())
//! ```
//!
//! ### Stateful API (inspect graph between reads)
//!
//! ```rust
//! use poa_consensus::{PoaGraph, PoaConfig};
//!
//! let reads: &[&[u8]] = &[b"CATCATCAT", b"CATCATCAT", b"CATCGTCAT"];
//!
//! let mut graph = PoaGraph::new(reads[0], PoaConfig::default())?;
//! for read in &reads[1..] {
//!     graph.add_read(read)?;
//! }
//! let consensus = graph.consensus()?;
//! let stats     = graph.stats();
//! println!("bubbles: {}", stats.bubble_count);
//! # Ok::<(), poa_consensus::PoaError>(())
//! ```
//!
//! ### Two-pass adaptive mode
//!
//! ```rust
//! use poa_consensus::{consensus_adaptive, PoaConfig};
//!
//! let reads: Vec<&[u8]> = vec![
//!     b"CATCATCAT", b"CATCATCAT", b"CATCATCAT",
//!     b"CATCGTCAT", b"CATCGTCAT", b"CATCGTCAT",
//! ];
//! // Pass 1 builds the graph and computes GraphStats.
//! // Pass 2 is selected automatically: multi-allele split, noise tightening,
//! // or semi-global switch, depending on what the stats reveal.
//! let alleles = consensus_adaptive(&reads, 0, &PoaConfig::default())?;
//! # Ok::<(), poa_consensus::PoaError>(())
//! ```
//!
//! ## Seed selection
//!
//! The seed read initialises the graph as a linear chain of nodes.  All other
//! reads are aligned into this initial structure, so a poor seed degrades
//! alignment quality for every subsequent read.
//!
//! **Choose a median-length read.**  The seed acts as the backbone; reads
//! shorter than the seed produce terminal deletions and reads longer than the
//! seed produce extensions.  A median-length seed minimises both.  In
//! high-throughput use (many loci), `reads.iter().enumerate().min_by_key(|(_, r)| r.len().abs_diff(median_len))` is a one-liner.
//!
//! **Avoid outliers.**  The longest and shortest reads are the most likely to
//! be error-prone or to span a different number of repeat units.  Using one as
//! the seed skews the initial graph in a direction that the alignment of
//! subsequent reads must then correct.
//!
//! **Orient reads before selecting the seed.**  Mixed-strand input produces a
//! garbage graph silently.  Call [`auto_orient`] before POA construction; it
//! uses k-mer matching (O(n) per read, no alignment) and returns borrowed
//! slices for reads already on the correct strand.
//!
//! Seed selection is the caller's responsibility.  The API takes an explicit
//! `seed_idx` so that the selection logic can live in the caller and be tuned
//! per application.
//!
//! ## Band width and scale
//!
//! The DP matrix has one row per read base and one column per graph node.
//! Unbanded alignment is O(read_len × graph_nodes) memory and time — manageable
//! for reads up to ~1 kb but prohibitive above that.
//!
//! | Read length | Band width | Memory per read (3 matrices, i32) |
//! |---|---|---|
//! | 600 bp | unbanded | ~1.4 MB |
//! | 600 bp | 100 | ~2.4 KB |
//! | 20 kb | unbanded | ~9.6 GB |
//! | 20 kb | adaptive (w≈210) | ~200 MB |
//!
//! **Rules of thumb:**
//!
//! - Reads ≤ 1 kb: unbanded (`band_width = 0`) is fine.
//! - Reads 1 kb–20 kb: set `band_width` to at least twice the expected
//!   length difference between reads, or enable `adaptive_band = true` (abPOA
//!   formula: `w = adaptive_band_b + adaptive_band_f × read_len`).
//! - Reads > 20 kb: adaptive banding is required; consider `poasta` for graphs
//!   that approach bacterial-gene size.
//!
//! A band that is too narrow returns `Err(PoaError::BandTooNarrow)` when the
//! terminal column is unreachable.  The library never silently produces a wrong
//! alignment — it errors instead.  [`PoaGraph::warnings_emitted`] counts how
//! many times a long-unbanded warning fired during `add_read` calls.
//!
//! ## Coverage and depth
//!
//! **`min_reads`** (default: 1) is the minimum number of reads required to
//! attempt consensus.  `consensus()` returns `Err(InsufficientDepth)` below
//! this threshold.  For reliable results, use at least 5 reads; 10+ is
//! preferable for heterozygous sites.
//!
//! **Boundary trim** removes leading and trailing nodes whose coverage falls
//! below the majority threshold `(n/2 + 1).max(2)` (or
//! `min_coverage_fraction × n` if set explicitly).  This corrects for seed
//! reads that happen to be longer or shorter than the true consensus length.
//!
//! **In multi-allele mode**, `min_reads` is applied per allele group, not to
//! the total read count.  With 10 reads split 5/5, each group must independently
//! satisfy `min_reads`.
//!
//! ## Known limitations
//!
//! **Phase-shift majority trim.** When the majority of reads are phase-shifted
//! relative to the seed (e.g. most reads start one repeat unit later), boundary
//! trim may incorrectly remove the seed's first node because its Match coverage
//! is low.  The long-term fix is to use [`extract_flanked_region`] to anchor
//! reads to a common reference point before POA; this eliminates phase ambiguity
//! entirely.
//!
//! **Rotation-ambiguous repeats without flanking sequence.** Trinucleotide
//! repeats like GAA can appear as GAA, AAG, or AGA depending on where the read
//! starts relative to the period.  Without flanking anchors, POA on a mixed
//! rotational-phase read set produces unreliable output.  The fix is the same:
//! run [`extract_flanked_region`] first so that every read enters POA already
//! anchored at the same phase.
//!
//! **Long expansions and partial reads.** When the repeat is longer than most
//! reads, reads that do not span the full locus create artifactual terminal
//! deletions that do not increment node coverage, making boundary nodes appear
//! low-coverage and vulnerable to trim.  Use `AlignmentMode::SemiGlobal` (free
//! terminal gaps) to prevent partial reads from distorting the boundary, or
//! apply [`extract_flanked_region`] to exclude non-spanning reads before POA.
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
