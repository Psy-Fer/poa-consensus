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
