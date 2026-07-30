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
//! // Single-allele consensus. `seed_idx` is validated for bounds but the
//! // engine self-seeds on the median-length read, so it is not load-bearing.
//! let result = consensus(&reads, 0, &PoaConfig::default())?;
//! println!("{}", String::from_utf8_lossy(&result.sequence));
//!
//! // Multi-allele: returns one Consensus per detected allele.
//! let alleles = consensus_multi(&reads, 0, &PoaConfig::default())?;
//! # Ok::<(), poa_consensus::PoaError>(())
//! ```
//!
//! ### Inspecting graph statistics
//!
//! Every [`Consensus`] carries a [`GraphStats`] describing the graph it was
//! built from (bubble count, coverage, noise level).  Use it to assess read
//! set quality or to decide whether a second, differently-configured call is
//! warranted.
//!
//! ```rust
//! use poa_consensus::{consensus, PoaConfig};
//!
//! let reads: &[&[u8]] = &[b"CATCATCAT", b"CATCATCAT", b"CATCGTCAT"];
//!
//! let result = consensus(reads, 0, &PoaConfig::default())?;
//! println!("bubbles: {}", result.graph_stats.bubble_count);
//! # Ok::<(), poa_consensus::PoaError>(())
//! ```
//!
//! ## Seed selection
//!
//! The seed read initialises the graph as a linear chain of nodes; all other
//! reads are aligned into this initial structure.  A median-length seed is the
//! best backbone (reads shorter than the seed produce terminal deletions,
//! longer reads produce extensions; the median minimises both), so **the
//! single-allele [`consensus`] entry point self-seeds on the median-length read
//! internally** — `seed_idx` is validated for bounds but does not bias the
//! result.  You do not need to pick a seed for the common case.
//!
//! For lower-level or multi-allele use, the [`select_seed`] helper implements
//! the same median/terminal-anchor heuristics if you want to choose explicitly.
//!
//! **Orient reads first, either way.**  Mixed-strand input produces a garbage
//! graph silently.  Call [`auto_orient`] before POA construction; it uses k-mer
//! matching (O(n) per read, no alignment) and returns borrowed slices for reads
//! already on the correct strand.
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
//! The banded DP uses a static-diagonal-union corridor (abPOA-style anti-fold):
//! the graph-geometry diagonal is always kept in-band, so the terminal column is
//! always reachable and a too-narrow band can never collapse the consensus to an
//! empty/all-gap alignment — at worst it yields a *suboptimal* alignment. On
//! pure-repetitive sequence several diagonals score equally and the DP may pick
//! one that differs by a repeat unit; `diagnose()` surfaces gross truncation so
//! the caller can rebuild unbanded and compare.
//!
//! ## Coverage and depth
//!
//! **`min_reads`** (default: 3) is the minimum number of reads required to
//! attempt consensus.  `consensus()` returns `Err(InsufficientDepth)` below
//! this threshold (the `(n/2 + 1).max(2)` coverage floor breaks down at depth
//! 1-2, so a lower floor produces silently unreliable output).  For reliable
//! results, use at least 5 reads; 10+ is preferable for heterozygous sites.
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
//! Using this crate's defaults, `gap_open = -4` and `gap_extend = -3`:
//!
//! ```text
//! One gap of length 4:
//!     gap_open + 4 * gap_extend = -4 + (-12) = -16
//!
//! Four gaps of length 1:
//!     4 * (gap_open + 1 * gap_extend) = 4 * (-7) = -28
//! ```
//!
//! The single 4-base deletion costs -16. Four separate 1-base deletions cost
//! -28. The aligner strongly prefers the single-event explanation, which
//! matches biological reality.
//!
//! With linear gaps (a flat `-1` per base) both scenarios cost -4 and the
//! aligner cannot distinguish them at all.
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
//! Scoring with the defaults (match = +2, mismatch = -4, gap_open = -4,
//! gap_extend = -3):
//!
//! ```text
//! 8 matching bases:  8 * (+2)           = +16
//! 1 gap of length 4: -4 + 4 * (-3)     = -16
//! Total:                                =   0
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
//! 8 matching bases:   8 * (+2)          = +16
//! 4 mismatches:       4 * (-4)          = -16
//! Total:                                =   0
//! ```
//!
//! With these defaults the two alignments tie at this length: a 4-base gap and
//! 4 mismatches score equally. That balance is deliberate — the gap penalties
//! are harsh enough that the aligner does not open cheap spurious gaps in
//! homopolymer or periodic runs (the failure mode that over-called repeat
//! lengths under looser scoring), but not so harsh that a genuine deletion is
//! always rejected in favour of mismatches. Shorter gaps tilt toward
//! mismatches, longer gaps toward the single deletion event. Crucially, POA
//! does not have to get this call right on any one read (next section).
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
//! defaults (`match = +2`, `mismatch = -4`, `gap_open = -4`, `gap_extend = -3`)
//! are calibrated for Oxford Nanopore and PacBio HiFi reads at tandem-repeat
//! loci, where cheap gaps otherwise over-call repeat length.
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
//! The defaults used in this crate (`match = +2`, `mismatch = -4`,
//! `gap_open = -4`, `gap_extend = -3`) work well for Oxford Nanopore and PacBio
//! HiFi reads at short tandem repeat loci. As a rule of thumb:
//!
//! - Increase the magnitude of `gap_open`/`gap_extend` to make the aligner even
//!   more reluctant to open gaps, preferring mismatches instead.
//! - Decrease the magnitude of `gap_extend` to make long gaps cheaper per base,
//!   useful when reads have systematic length variation.
//! - Setting `gap_open = 0` recovers linear gap behaviour where only
//!   `gap_extend` matters.

pub mod analysis;
pub mod config;
pub mod error;
pub mod flank;
pub mod orient;
pub mod seed;
pub mod types;

/// Internal POA engine. **Not part of the public API** and exempt from semver:
/// its types and functions may change or be removed in any release. Use the
/// crate-root functions ([`consensus`], [`consensus_multi`], [`bridged_consensus`])
/// instead. Exposed only so the crate's own benchmark/robustness test harness can
/// switch between engine variants.
#[doc(hidden)]
pub mod poa2;

/// Internal linkage-phasing primitives behind [`consensus_multi`]. **Not part of
/// the public API** and exempt from semver. See [`poa2`].
#[doc(hidden)]
pub mod phasing;

/// Internal multi-allele consensus pipeline (structural-propose + linkage-confirm)
/// on top of the [`poa2`] engine. **Not part of the public API** and exempt from
/// semver; reached via the crate-root [`consensus_multi`].
#[doc(hidden)]
pub mod multi;

pub use analysis::{
    ConsensusWarnings, DiagnoseConfig, InteriorSupportWarning, LowDepthWarning,
    StructuralCompetingSummary, TruncationWarning, diagnose,
};
pub use config::{AlignmentMode, ConsensusMode, PoaConfig};
pub use error::PoaError;
pub use flank::extract_flanked_region;
pub use orient::{auto_orient, orient_to_seed, reverse_complement};
pub use seed::{SeedSelection, select_seed};
pub use types::{BubbleSite, Consensus, CoverageGap, GapKind, GraphStats, Strand};

// ── Internal helpers ─────────────────────────────────────────────────────────

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

/// Emit a one-line stderr warning when the caller has forced fully unbanded
/// alignment (`band_width == 0` and `adaptive_band == false`) on long reads,
/// where the DP allocates O(nodes × read_len) memory. Gated by
/// [`PoaConfig::warn_on_long_unbanded`] (set it `false` to silence).
fn warn_long_unbanded(reads: &[&[u8]], config: &PoaConfig) {
    const LONG_READ_BP: usize = 1000;
    if config.warn_on_long_unbanded && config.band_width == 0 && !config.adaptive_band {
        let max_len = reads.iter().map(|r| r.len()).max().unwrap_or(0);
        if max_len > LONG_READ_BP {
            eprintln!(
                "poa-consensus: warning: unbanded alignment (band_width = 0) on reads up \
                 to {max_len} bp uses O(nodes × read_len) memory; set band_width or enable \
                 adaptive_band for long reads (silence with warn_on_long_unbanded = false)"
            );
        }
    }
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
    if reads.len() < config.min_reads {
        return Err(PoaError::InsufficientDepth {
            got: reads.len(),
            min: config.min_reads,
        });
    }
    warn_long_unbanded(reads, config);
    // Cutover Step 2: single-allele consensus now runs on the clean engine
    // (poa2). The clean engine is seed-robust, so it picks its own backbone --
    // the MEDIAN-length read -- rather than trusting the caller's `seed_idx`
    // hint (a length-outlier seed would bias a periodic repeat toward that
    // outlier's length; median avoids it). `seed_idx` is still validated for
    // bounds above, preserving the error contract.
    let _ = seed_idx;
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by_key(|&i| reads[i].len());
    let med = order[order.len() / 2];
    let mut g = crate::poa2::Poa::new(reads[med], config.clone());
    for (i, r) in reads.iter().enumerate() {
        if i != med {
            g.add_read(r);
        }
    }
    // Consensus extraction honours `consensus_mode`:
    //  - BestFit (default): best-fit — compute both a heaviest-path and a
    //    majority-frequency consensus and keep whichever the reads better
    //    support (majority-frequency recovers homopolymer/length-variable
    //    repeats the heaviest path over-calls).
    //  - MajorityFrequency: force the MSA-column majority consensus (useful for
    //    high-depth amplicons where column majority is trusted outright).
    // Per-position output fields are made consistent with the chosen sequence.
    Ok(match config.consensus_mode {
        crate::config::ConsensusMode::MajorityFrequency => g.consensus_full_majority(),
        crate::config::ConsensusMode::BestFit => g.consensus_full_best_fit(reads),
    })
}

/// Build a multi-allele consensus from `reads`.
///
/// Returns one [`Consensus`] per detected allele.  If no heterozygous bubble is
/// found the result is a single-element `Vec` equivalent to calling [`consensus`].
///
/// An allele below [`PoaConfig::min_allele_freq`] is merged into the majority and
/// does not appear in the output; to detect a low-frequency / mosaic allele, lower
/// `min_allele_freq` and re-run (see that field's docs).
pub fn consensus_multi(
    reads: &[&[u8]],
    seed_idx: usize,
    config: &PoaConfig,
) -> Result<Vec<Consensus>, PoaError> {
    validate(reads, seed_idx)?;
    warn_long_unbanded(reads, config);
    // Multi-allele runs on the clean poa2 HYBRID engine: structural-bubble
    // discovery proposes splits, linkage discovery + consensus-difference /
    // bimodality confirmation refines them and rejects false splits. It picks its
    // own median-length seed and emits external read indices directly, so
    // `seed_idx` is only used for input validation and no remap is needed. On the
    // robustness matrix it beats the legacy engine on read-assignment cleanliness
    // and single-allele safety at near-equal split sensitivity. (Same-length
    // substitution haplotypes are the one weaker case; tracked separately.)
    crate::multi::hybrid_consensus_multi(reads, config)
}

/// Fraction of reads that must fail to span both flanks before the flanked
/// entry points switch to anchored mode (restrict to spanning reads).
const FLANK_ANCHOR_TRIGGER: f64 = 0.1;

/// Indices of the reads that span BOTH flanks (i.e. are not partial). We anchor
/// by restricting to these reads but keep them WHOLE — the flanks provide unique
/// anchoring context the aligner and (especially) the multi-allele length-phasing
/// need; trimming to the bare repeat first makes periodic segments align
/// ambiguously and collapses length alleles together (measured 2026-07-30:
/// keeping reads whole calls 182/216 length-diploids correct vs 103 when trimmed).
/// The output is sliced back to the repeat region afterwards for a consistent
/// flanks-excluded result.
fn spanning_indices(reads: &[&[u8]], left: &[u8], right: &[u8]) -> Vec<usize> {
    (0..reads.len())
        .filter(|&i| crate::flank::flank_span(reads[i], left, right).is_some())
        .collect()
}

/// Anchor only when a meaningful fraction of reads fail to span (so partials are
/// distorting the raw consensus) AND enough reads still span to build from (so we
/// don't trade a partials problem for a depth problem — the raw fallback covers
/// the rest).
fn should_anchor(total: usize, spanning: usize, config: &PoaConfig) -> bool {
    if total == 0 {
        return false;
    }
    let partial_frac = 1.0 - (spanning as f64 / total as f64);
    partial_frac >= FLANK_ANCHOR_TRIGGER && spanning >= config.min_reads
}

/// Slice a whole-region `Consensus` down to the repeat segment between `left` and
/// `right` (located in the consensus's own low-error sequence), so a flanked
/// result always describes the sequence BETWEEN the flanks. Returns the consensus
/// unchanged if the flanks cannot be located in it (rare; the consensus is clean).
fn slice_to_repeat(c: Consensus, left: &[u8], right: &[u8]) -> Consensus {
    let Some((s, e)) = crate::flank::flank_span(&c.sequence, left, right) else {
        return c;
    };
    if e > c.sequence.len() || s >= e {
        return c;
    }
    Consensus {
        sequence: c.sequence[s..e].to_vec(),
        coverage: c
            .coverage
            .get(s..e)
            .map(<[u32]>::to_vec)
            .unwrap_or_default(),
        path_weights: c
            .path_weights
            .get(s..e)
            .map(<[i32]>::to_vec)
            .unwrap_or_default(),
        n_reads: c.n_reads,
        graph_stats: c.graph_stats,
        gaps: c
            .gaps
            .into_iter()
            .filter(|g| g.start >= s && g.end <= e)
            .map(|g| CoverageGap {
                start: g.start - s,
                end: g.end - s,
                kind: g.kind,
            })
            .collect(),
        bubble_sites: c
            .bubble_sites
            .into_iter()
            .filter(|b| b.consensus_pos >= s && b.consensus_pos < e)
            .map(|mut b| {
                b.consensus_pos -= s;
                b
            })
            .collect(),
        read_indices: c.read_indices,
    }
}

/// Single-allele consensus of the repeat region between `left_flank` and
/// `right_flank`, auto-anchoring when the read set is partial-heavy.
///
/// Reads that do not span both flanks (partial / non-spanning reads) distort a
/// raw consensus. When a meaningful fraction of reads fail to span, this builds
/// only from the reads that DO span (kept whole, so the flanks still anchor the
/// alignment); otherwise it builds from all reads. Either way the result is sliced
/// so the returned [`Consensus`] describes the sequence **between** the flanks
/// (flanks excluded). Falls back to the raw all-reads consensus when too few reads
/// span to anchor safely. Flanks of at least ~20 bp of unique sequence are
/// recommended.
pub fn consensus_flanked(
    reads: &[&[u8]],
    left_flank: &[u8],
    right_flank: &[u8],
    config: &PoaConfig,
) -> Result<Consensus, PoaError> {
    validate(reads, 0)?;
    let span = spanning_indices(reads, left_flank, right_flank);
    let full = if should_anchor(reads.len(), span.len(), config) {
        // Partial-heavy: build only from the spanning reads (kept whole).
        let refs: Vec<&[u8]> = span.iter().map(|&i| reads[i]).collect();
        consensus(&refs, 0, config)
    } else {
        // Spanning-heavy or too few spanned to anchor safely: all reads.
        consensus(reads, 0, config)
    }?;
    Ok(slice_to_repeat(full, left_flank, right_flank))
}

/// Multi-allele version of [`consensus_flanked`]: returns one [`Consensus`] per
/// detected allele, each describing the repeat region between the flanks. Uses
/// the same partial-heavy auto-anchoring; `read_indices` are always reported as
/// indices into the input `reads` slice (remapped from the anchored subset).
pub fn consensus_multi_flanked(
    reads: &[&[u8]],
    left_flank: &[u8],
    right_flank: &[u8],
    config: &PoaConfig,
) -> Result<Vec<Consensus>, PoaError> {
    validate(reads, 0)?;
    let span = spanning_indices(reads, left_flank, right_flank);
    if should_anchor(reads.len(), span.len(), config) {
        // Partial-heavy: phase only the spanning reads (kept whole so the flanks
        // still anchor the length-bubble detection), then slice each allele to the
        // repeat and remap read indices from the spanning subset back to originals.
        let refs: Vec<&[u8]> = span.iter().map(|&i| reads[i]).collect();
        if let Ok(mut alleles) = consensus_multi(&refs, 0, config) {
            for a in &mut alleles {
                for idx in &mut a.read_indices {
                    *idx = span[*idx];
                }
            }
            return Ok(alleles
                .into_iter()
                .map(|c| slice_to_repeat(c, left_flank, right_flank))
                .collect());
        }
    }
    let full = consensus_multi(reads, 0, config)?;
    Ok(full
        .into_iter()
        .map(|c| slice_to_repeat(c, left_flank, right_flank))
        .collect())
}

/// Build a consensus from two non-overlapping read groups with a gap of
/// unknown length between them.
///
/// Use this when reads cover only the left end and right end of a template
/// with no read bridging the middle — for example, when a repeat expansion
/// is longer than any individual read.
///
/// Each group is assembled independently using [`consensus`].  The results
/// are concatenated and a single [`CoverageGap`] with
/// [`GapKind::Unknown`] is inserted at the join point.
///
/// The returned `Consensus.n_reads` is the sum of both groups.
/// `Consensus.graph_stats` reflects the left group's graph.
///
/// # Size interpretation
///
/// ```text
/// |=left_consensus=|???unknown???|=right_consensus=|
/// total ≥ left_consensus.len() + right_consensus.len()
/// ```
///
/// `gap.min_size()` returns `None` (size unknown).  The caller can report:
/// ```text
/// format!("≥{} bp with a gap of unknown length",
///     cons.sequence.len())
/// ```
pub fn bridged_consensus(
    left_reads: &[&[u8]],
    left_seed_idx: usize,
    right_reads: &[&[u8]],
    right_seed_idx: usize,
    config: &PoaConfig,
) -> Result<Consensus, PoaError> {
    validate(left_reads, left_seed_idx)?;
    validate(right_reads, right_seed_idx)?;

    let left = consensus(left_reads, left_seed_idx, config)?;
    let right = consensus(right_reads, right_seed_idx, config)?;

    let join = left.sequence.len();

    let mut sequence = left.sequence.clone();
    sequence.extend_from_slice(&right.sequence);

    let mut coverage = left.coverage.clone();
    coverage.extend_from_slice(&right.coverage);

    let mut path_weights = left.path_weights.clone();
    path_weights.extend_from_slice(&right.path_weights);

    // Carry forward any spanning gaps from each segment (offset right gaps by join).
    let mut gaps = left.gaps.clone();
    for gap in &right.gaps {
        gaps.push(CoverageGap {
            start: gap.start + join,
            end: gap.end + join,
            kind: gap.kind,
        });
    }
    // Insert the unknown gap at the join point, then sort by position.
    gaps.push(CoverageGap {
        start: join,
        end: join,
        kind: GapKind::Unknown,
    });
    gaps.sort_by_key(|g| g.start);

    Ok(Consensus {
        sequence,
        coverage,
        path_weights,
        n_reads: left.n_reads + right.n_reads,
        graph_stats: left.graph_stats,
        gaps,
        bubble_sites: vec![],
        read_indices: vec![],
    })
}
