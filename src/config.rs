/// How each read is aligned to the graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentMode {
    /// End-to-end: both the read and the graph must be fully consumed, so
    /// terminal gaps are charged. Only appropriate when every read spans exactly
    /// the same template (e.g. trimmed amplicons).
    Global,
    /// Default. Free terminal gaps: a read may start and end anywhere in the
    /// graph (traceback starts from the best cell in the last row). Required for
    /// extracted STR reads: under `Global`, if the seed carries more leading
    /// identical bases than most reads, shorter reads place their deletions
    /// inside the run and the extra bases stay on the consensus. Semi-global
    /// also keeps partial reads from distorting the boundaries.
    SemiGlobal,
}

/// How the consensus sequence is extracted from the built graph.
///
/// Honoured by the single-allele entry points ([`consensus`](crate::consensus),
/// [`consensus_flanked`](crate::consensus_flanked),
/// [`bridged_consensus`](crate::bridged_consensus)). The multi-allele entry
/// points always emit each allele's heaviest path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsensusMode {
    /// Default. "Best-fit": build both a heaviest-path and a majority-frequency
    /// consensus and realign every read to each. The majority-frequency one is
    /// kept only if it is strictly better on both indel measures: fewer indel
    /// bases, and a lower affine indel penalty (each insert/delete event pays
    /// `|gap_open|`, each indel base `|gap_extend|`). Mismatches are not
    /// counted. Otherwise the heaviest path is kept. The two win on
    /// different inputs (heaviest on clean/short reads, majority on high-error
    /// length-variable repeats), so picking per call usually gets the better of
    /// both; it is a heuristic, not a guarantee.
    BestFit,
    /// Force the majority-frequency (MSA-column) consensus regardless of fit:
    /// each column emits its plurality base, counting read deletions explicitly.
    /// Best when column majority is trusted outright, e.g. high-depth amplicons.
    MajorityFrequency,
}

/// Alignment, banding, consensus and multi-allele parameters.
///
/// `PoaConfig::default()` suits STR / VNTR loci from long reads (band 50 +
/// adaptive, semi-global, scoring 2 / -4 / -4 / -3, best-fit consensus,
/// homopolymer handling and tandem gap costs on). Set `homopolymer_min_run:
/// 0`, `homopolymer_gap: None` and `tandem_gap: None` for the plain engine path
/// (about 2.4x less CPU).
///
/// # Banding
///
/// Banding is on when `band_width > 0 || adaptive_band`. The band half-width
/// (columns each side of the band centre) is
/// `max(band_width, adaptive_band_b + ceil(adaptive_band_f * read_len), 4)`
/// with `adaptive_band` on, and `max(band_width, 4)` with it off.
/// Fully unbanded alignment needs `band_width = 0` **and**
/// `adaptive_band = false`. Reads shorter than 0.8x the graph span (partial
/// reads) are always aligned unbanded.
#[derive(Debug, Clone)]
pub struct PoaConfig {
    /// Minimum band half-width. Default 50. Banding is on when this is > 0 or
    /// `adaptive_band` is set; `0` together with `adaptive_band = false` means
    /// fully unbanded (warns on reads > ~1 kb, see `warn_on_long_unbanded`).
    pub band_width: usize,
    /// abPOA-style adaptive band: widen the half-width to
    /// `adaptive_band_b + ceil(adaptive_band_f * L)` when that exceeds
    /// `band_width` (reads over ~4 kb at the defaults), and band even when
    /// `band_width == 0`. Off: a fixed `band_width` half-width. Default `true`.
    pub adaptive_band: bool,
    /// Base component of the adaptive half-width term (abPOA default: 10).
    pub adaptive_band_b: usize,
    /// Length-proportional component of the adaptive half-width term (abPOA default: 0.01).
    pub adaptive_band_f: f32,
    /// Score for a read base matching a node base. Positive. Default 2.
    pub match_score: i32,
    /// Score for a read base aligned to a different node base. Negative. Default -4.
    pub mismatch_score: i32,
    /// One-time penalty when a gap opens. Negative. Default -4.
    pub gap_open: i32,
    /// Per-base penalty inside a gap. Negative. Default -3.
    pub gap_extend: i32,
    /// Boundary-trim floor as a fraction of reads: leading and trailing
    /// consensus nodes covered by fewer than `ceil(n * min_coverage_fraction)`
    /// reads are trimmed. Interior positions are never removed. `0.0` (default)
    /// uses a strict majority, `(n/2 + 1).max(2)`.
    pub min_coverage_fraction: f64,
    /// Minimum fraction of reads supporting an allele for it to be called (the
    /// minority-arm support floor is `ceil(n * min_allele_freq).max(2)`).
    ///
    /// This is a **hard mosaic-sensitivity boundary**, not just a noise filter: an
    /// allele below this frequency is merged into the majority and does **not**
    /// appear anywhere in the output (its bubble is emitted only at or above this
    /// floor). A single run therefore cannot report a sub-threshold mosaic: to
    /// detect a low-frequency / subclonal allele (e.g. a 10-15% mosaic), lower
    /// `min_allele_freq` and re-run `consensus_multi`. On noisy ONT data raise it
    /// (~0.40) instead, to avoid error-driven false second-allele calls.
    pub min_allele_freq: f64,
    /// Minimum number of reads required to build a consensus (default 3);
    /// below it the entry points return `Err(InsufficientDepth)`. Multi-allele
    /// mode checks the total non-empty read count and every allele group.
    pub min_reads: usize,
    /// Read-to-graph alignment mode. Default [`AlignmentMode::SemiGlobal`].
    pub alignment_mode: AlignmentMode,
    /// Consensus extraction mode. Default [`ConsensusMode::BestFit`].
    pub consensus_mode: ConsensusMode,
    /// Emit a warning to stderr when fully unbanded alignment (`band_width = 0`
    /// and `adaptive_band = false`) is used on reads longer than ~1 kb.
    pub warn_on_long_unbanded: bool,
    /// Minimum arm span (nodes) for a bubble to count as structural: structural
    /// bubbles propose multi-allele splits and set [`BubbleSite::is_structural`](crate::BubbleSite::is_structural).
    /// Shorter bubbles (SNPs, short indels) are left to linkage phasing. Default 10.
    pub phasing_bubble_min_span: usize,
    /// Homopolymer compression. Default `2`; `0` = off. Otherwise reads are rewritten so each run of `>= homopolymer_min_run` identical
    /// bases is one "poly" symbol (`1` = pure homopolymer compression, `2` = the
    /// `A,C,G,T,polyA..polyT` alphabet), a consensus is built in that space, and
    /// each poly symbol is expanded to the (lower) median run length of the
    /// reads in its graph column. If that vote is weak (under half the reads, or
    /// a split vote), the run lengths come from realigning the reads to the
    /// compressed consensus instead.
    ///
    /// The result is only a candidate: it is spliced into the plain consensus
    /// only where the two differ inside homopolymer runs (blocks made of run
    /// bases and single interrupts between runs), and the splice is kept only if
    /// it fits the reads strictly better (fewer mismatch + indel bases per read
    /// base) than the plain consensus. Applies via [`consensus`](crate::consensus),
    /// [`consensus_flanked`](crate::consensus_flanked) and
    /// [`bridged_consensus`](crate::bridged_consensus); the multi-allele entry
    /// points ignore it. With the default `homopolymer_gap` the homopolymer
    /// path costs roughly 2.4x the CPU of the plain engine path (both options
    /// off: `homopolymer_min_run: 0, homopolymer_gap: None`; CLI `--no-homopolymer`).
    pub homopolymer_min_run: usize,
    /// Homopolymer-aware gap costs. Default `Some((-2, -1))`; `None` = off.
    /// `Some((open, extend))` replaces `gap_open`/`gap_extend` for an indel in
    /// a homopolymer context:
    ///
    /// - an inserted read base that equals the read base before or after it;
    /// - a deleted node that has a same-base graph neighbour (predecessor or
    ///   successor) **and** whose base equals the read base on either side of
    ///   the deletion point.
    ///
    /// Requiring a run on both sides keeps a periodic repeat's unit indel
    /// (`CAG|CAG`) at full cost. Cheap enough values (e.g. `(-2, -1)`) make a
    /// base shifted within a run align as indels rather than a mismatch pair.
    /// Guarded like `homopolymer_min_run` (splice inside runs, kept only on a
    /// strictly better fit; splice run threshold 2 when used alone), with the
    /// same entry-point scope.
    pub homopolymer_gap: Option<(i32, i32)>,
    /// Tandem-repeat gap costs. Default `Some((-2, -1))`; `None` = off.
    ///
    /// Gaining or losing one copy of a short repeat unit (polymerase slippage,
    /// sequencing stutter) is the commonest real difference inside a tandem
    /// repeat, but plain affine costs price it like any other indel. Reads
    /// then reach for cheaper-looking alignments: sliding across a repeat
    /// boundary with mismatches and reusing a node another read created, which
    /// can carry a wrong extra unit onto the heaviest path depending on read
    /// order. With this set, an indel that copies its neighbour in a locally
    /// periodic stretch of the read (a 2-6 base unit immediately followed by a
    /// copy of itself) pays `(open, extend)` instead. The test is local to each
    /// alignment cell: no repeat tract or unit count is identified.
    ///
    /// Like the homopolymer options this only builds a candidate, never changes
    /// the plain engine: a consensus is built with these costs, and it replaces
    /// the plain call only if it fits the reads strictly better (fewer
    /// mismatches plus indel bases per read base). Applied wholesale, unguarded, the costs pull
    /// split-length repeats toward the shorter length. The candidate is built
    /// only when the plain call has a base with support quality below 10
    /// (somewhere the reads disagree), so clean loci pay nothing. Applies via
    /// [`consensus`](crate::consensus), [`consensus_flanked`](crate::consensus_flanked)
    /// and [`bridged_consensus`](crate::bridged_consensus); the multi-allele
    /// entry points ignore it. CLI `--tandem-gap OPEN,EXTEND|off`.
    pub tandem_gap: Option<(i32, i32)>,
}

impl Default for PoaConfig {
    fn default() -> Self {
        PoaConfig {
            band_width: 50,
            adaptive_band: true,
            adaptive_band_b: 10,
            adaptive_band_f: 0.01,
            // abPOA-style scoring: gaps and mismatches are harsh relative to
            // match (+2), so the aligner does not open cheap spurious gaps in
            // homopolymer/periodic runs. The old +1/-1/-2/-1 made gaps too cheap
            // (a gap-open cost only 2 matches), scattering homopolymer alignments
            // at high error and over-calling repeats; abPOA-like scoring roughly
            // halves that error with no regression on clean inputs (validated on
            // the robustness matrix + 3-way comparison, 2026-07-28).
            match_score: 2,
            mismatch_score: -4,
            gap_open: -4,
            gap_extend: -3,
            min_coverage_fraction: 0.0,
            min_allele_freq: 0.2,
            min_reads: 3,
            alignment_mode: AlignmentMode::SemiGlobal,
            consensus_mode: ConsensusMode::BestFit,
            warn_on_long_unbanded: true,
            phasing_bubble_min_span: 10,
            // Homopolymer handling on by default (compression k=2 plus
            // homopolymer gap costs): both candidates are guarded by a splice
            // and a strict fit check, so they only change runs the reads
            // support. ~2.4x the plain engine's CPU.
            homopolymer_min_run: 2,
            homopolymer_gap: Some((-2, -1)),
            // Tandem gap costs on by default: a guarded candidate, built only
            // where the plain call has weak support (see the field docs).
            tandem_gap: Some((-2, -1)),
        }
    }
}
