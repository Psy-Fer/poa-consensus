use std::io::{self, BufRead, BufReader, Write};
use std::process;

use clap::{Parser, ValueEnum};

use poa_consensus::{
    AlignmentMode, ConsensusMode, DiagnoseConfig, PoaConfig, PoaError, SeedSelection, auto_orient,
    diagnose, select_seed,
};

/// Consensus extraction mode, mirrored from [`ConsensusMode`] for the CLI.
#[derive(Clone, Copy, ValueEnum)]
enum ConsensusModeArg {
    /// Heaviest-path (edge-weight) traversal. Best default for STR/VNTR data.
    Heaviest,
    /// Most-frequent-base per column (MSA majority). Better for near-equal-length
    /// HiFi read sets; counts Delete traversals explicitly.
    Majority,
}

impl From<ConsensusModeArg> for ConsensusMode {
    fn from(m: ConsensusModeArg) -> Self {
        match m {
            ConsensusModeArg::Heaviest => ConsensusMode::HeaviestPath,
            ConsensusModeArg::Majority => ConsensusMode::MajorityFrequency,
        }
    }
}

/// Build a consensus sequence from FASTA or FASTQ reads using Partial Order
/// Alignment.
///
/// Reads are auto-oriented to the seed strand before alignment.  The default
/// seed is chosen automatically via a terminal k-mer anchor heuristic (see
/// `SeedSelection::Auto`); override with --seed.
#[derive(Parser)]
#[command(name = "poa-consensus", version, allow_negative_numbers = true)]
struct Args {
    /// Input FASTA or FASTQ file. Use '-' for stdin.
    #[arg(default_value = "-")]
    input: String,

    /// Seed read index (0-based). Defaults to an automatically chosen seed
    /// (terminal k-mer anchor heuristic; falls back to the longest read).
    /// Note: the single-allele path re-seeds on the median-length read
    /// internally, so --seed only affects seed validation and multi-allele.
    #[arg(short, long)]
    seed: Option<usize>,

    /// Multi-allele mode: output one FASTA record per detected allele.
    #[arg(short = 'm', long)]
    multi: bool,

    /// Suppress all warnings and notes on stderr.  The consensus sequence is
    /// still written to stdout.  Errors that prevent consensus building are
    /// always printed regardless of this flag.
    #[arg(short = 'q', long)]
    quiet: bool,

    // ── Band ──────────────────────────────────────────────────────────────────
    /// Fixed band width floor for banded DP (0 = unbanded).  Combined with
    /// adaptive band (default on), effective width is max(band_width, adaptive
    /// formula).  Set to 0 with --no-adaptive-band for fully unbanded DP.
    #[arg(short = 'b', long, default_value_t = 50, help_heading = "Band")]
    band_width: usize,

    /// Disable adaptive band width.  Adaptive band (default on) uses the
    /// formula w = b + f × read_len as a floor alongside --band-width.
    #[arg(long, help_heading = "Band")]
    no_adaptive_band: bool,

    /// Adaptive band base component `b` in w = b + f × read_len (abPOA default).
    #[arg(long, default_value_t = 10, help_heading = "Band")]
    adaptive_band_b: usize,

    /// Adaptive band length-proportional component `f` in w = b + f × read_len.
    #[arg(long, default_value_t = 0.01, help_heading = "Band")]
    adaptive_band_f: f32,

    // ── Scoring ───────────────────────────────────────────────────────────────
    /// Match score (positive).
    #[arg(
        long = "match",
        value_name = "SCORE",
        default_value_t = 1,
        help_heading = "Scoring"
    )]
    match_score: i32,

    /// Mismatch score (negative).
    #[arg(long = "mismatch", value_name = "SCORE", default_value_t = -1, help_heading = "Scoring")]
    mismatch_score: i32,

    /// Gap-open penalty (negative; charged once when a gap opens).
    #[arg(long, value_name = "PENALTY", default_value_t = -2, help_heading = "Scoring")]
    gap_open: i32,

    /// Gap-extend penalty (negative; per base inside a gap).
    #[arg(long, value_name = "PENALTY", default_value_t = -1, help_heading = "Scoring")]
    gap_extend: i32,

    // ── Coverage / consensus ──────────────────────────────────────────────────
    /// Minimum reads required to attempt consensus.
    #[arg(long, default_value_t = 3, help_heading = "Coverage / consensus")]
    min_reads: usize,

    /// Fraction of reads that must cover a node for it to appear in the
    /// consensus.  0 uses the strict-majority default (n/2 + 1).max(2).
    #[arg(long, default_value_t = 0.0, help_heading = "Coverage / consensus")]
    min_coverage_fraction: f64,

    /// Consensus extraction mode.
    #[arg(long, value_enum, default_value_t = ConsensusModeArg::Heaviest, help_heading = "Coverage / consensus")]
    consensus_mode: ConsensusModeArg,

    // ── Alignment ─────────────────────────────────────────────────────────────
    /// Global alignment: penalise terminal gaps (use when reads are guaranteed
    /// to span the full locus from identical start/end positions).  The default
    /// is semi-global (free terminal gaps), which is correct for extracted STR
    /// reads that may start or end at slightly different positions.
    #[arg(long, help_heading = "Alignment")]
    global: bool,

    // ── Multi-allele ──────────────────────────────────────────────────────────
    /// Minimum fraction of reads supporting a bubble arm for it to count as an
    /// allele candidate.  Raise to ~0.40 on noisy ONT data to avoid spurious
    /// second-allele calls.
    #[arg(long, default_value_t = 0.2, help_heading = "Multi-allele")]
    min_allele_freq: f64,

    /// Minimum arm span (nodes) for a bubble to count as a structural variant
    /// for cross-bubble phasing; below this uses single-bubble (SNP) partitioning.
    #[arg(long, default_value_t = 10, help_heading = "Multi-allele")]
    phasing_bubble_min_span: usize,

    // ── Diagnostics ───────────────────────────────────────────────────────────
    /// Suppress the long-read unbanded warning (reads > ~1 kb with band 0).
    #[arg(long, help_heading = "Diagnostics")]
    no_long_unbanded_warning: bool,

    /// Read count below which a depth warning is emitted.
    #[arg(long, default_value_t = 10, help_heading = "Diagnostics")]
    depth_warn_threshold: usize,

    /// Read count below which the depth warning is marked critical.
    #[arg(long, default_value_t = 5, help_heading = "Diagnostics")]
    depth_critical_threshold: usize,

    /// Per-allele depth warning threshold (multi-allele mode).
    #[arg(long, default_value_t = 15, help_heading = "Diagnostics")]
    depth_allele_threshold: usize,

    /// Weight-fraction below which an interior-support warning fires.
    #[arg(long, default_value_t = 0.15, help_heading = "Diagnostics")]
    interior_support_threshold: f32,

    /// Fraction of the consensus skipped from each end when checking interior
    /// support (the middle 1 - 2×margin is checked).
    #[arg(long, default_value_t = 0.20, help_heading = "Diagnostics")]
    boundary_margin: f32,

    /// consensus_len / median_read_len ratio below which a truncation warning
    /// fires.  0 disables the check.
    #[arg(long, default_value_t = 0.6, help_heading = "Diagnostics")]
    truncation_ratio_threshold: f32,
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // ── Read all sequences ────────────────────────────────────────────────────
    let reads = load_reads(&args.input)?;

    if reads.is_empty() {
        return Err("no reads in input".into());
    }

    // Warn on non-ACGT bases (POA handles them as mismatches, but callers
    // should know).
    for seq in &reads {
        if seq
            .iter()
            .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
        {
            eprintln!(
                "poa-consensus: warning: read contains non-ACGT bases; \
                 treating as mismatches"
            );
            break;
        }
    }

    // ── Single-record passthrough ─────────────────────────────────────────────
    if reads.len() == 1 {
        let stdout = io::stdout();
        let mut out = stdout.lock();
        writeln!(out, ">consensus reads=1 seed=0 band=unbanded")?;
        out.write_all(&reads[0])?;
        writeln!(out)?;
        return Ok(());
    }

    // ── Seed selection and orientation ────────────────────────────────────────
    let seed_idx = match args.seed {
        Some(idx) => {
            if idx >= reads.len() {
                return Err(format!("seed {} out of range ({} reads)", idx, reads.len()).into());
            }
            idx
        }
        None => {
            let slices: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            match select_seed(&slices, &SeedSelection::Auto) {
                Ok(idx) => idx,
                // `select_seed`'s terminal-k-mer spanning heuristic is unreliable
                // at high error / long repeats and can false-positive
                // `NoSpanningReads` on genuinely-spanning noisy reads. The library
                // `consensus()` re-seeds on the median-length read internally and
                // does not require spanning, so rather than refuse output (abPOA
                // and SPOA both produce a best-effort consensus here), fall back to
                // the longest read as the orientation seed and warn.
                Err(e) => {
                    if !args.quiet {
                        eprintln!(
                            "poa-consensus: warning: automatic seed selection was inconclusive \
                             ({e}); proceeding with the longest read as seed. If reads are truly \
                             split into non-overlapping left/right groups, use bridged_consensus."
                        );
                    }
                    (0..reads.len())
                        .max_by_key(|&i| reads[i].len())
                        .unwrap_or(0)
                }
            }
        }
    };

    let oriented = auto_orient(&reads, seed_idx);
    let slices: Vec<&[u8]> = oriented.iter().map(|r| r.as_ref()).collect();

    // ── Config ────────────────────────────────────────────────────────────────
    let adaptive_band = !args.no_adaptive_band;
    let config = PoaConfig {
        band_width: args.band_width,
        adaptive_band,
        adaptive_band_b: args.adaptive_band_b,
        adaptive_band_f: args.adaptive_band_f,
        match_score: args.match_score,
        mismatch_score: args.mismatch_score,
        gap_open: args.gap_open,
        gap_extend: args.gap_extend,
        min_coverage_fraction: args.min_coverage_fraction,
        min_allele_freq: args.min_allele_freq,
        min_reads: args.min_reads,
        alignment_mode: if args.global {
            AlignmentMode::Global
        } else {
            AlignmentMode::SemiGlobal
        },
        consensus_mode: args.consensus_mode.into(),
        warn_on_long_unbanded: !args.no_long_unbanded_warning,
        phasing_bubble_min_span: args.phasing_bubble_min_span,
        // multi_allele is managed by the library (consensus_multi sets it true);
        // the CLI's --multi selects the code path rather than setting this.
        ..PoaConfig::default()
    };

    // Diagnostic thresholds shared by both the single- and multi-allele paths.
    // `is_allele_partition` is flipped on per allele in the multi branch.
    let diag_cfg = DiagnoseConfig {
        depth_warn_threshold: args.depth_warn_threshold,
        depth_critical_threshold: args.depth_critical_threshold,
        depth_allele_threshold: args.depth_allele_threshold,
        interior_support_threshold: args.interior_support_threshold,
        boundary_margin: args.boundary_margin,
        truncation_ratio_threshold: args.truncation_ratio_threshold,
        is_allele_partition: false,
    };

    let band_desc = if adaptive_band {
        "adaptive".to_string()
    } else if args.band_width == 0 {
        "unbanded".to_string()
    } else {
        args.band_width.to_string()
    };

    // ── Consensus and output ──────────────────────────────────────────────────
    let stdout = io::stdout();
    let mut out = stdout.lock();
    let n = reads.len();

    if args.multi {
        let alleles = poa_consensus::consensus_multi(&slices, seed_idx, &config)
            .inspect_err(|e| explain_error(e, n))?;
        let total = alleles.len();
        let allele_cfg = DiagnoseConfig {
            is_allele_partition: true,
            ..diag_cfg.clone()
        };
        for (i, allele) in alleles.iter().enumerate() {
            if !args.quiet {
                let label = format!("allele {}/{}", i + 1, total);
                emit_warnings(&diagnose(allele, &allele_cfg), &label);
            }
            // allele.n_reads is the per-partition read count, not the total.
            // Callers use this to assess whether a minority allele has adequate depth.
            writeln!(
                out,
                ">allele_{} reads={} total_reads={n} seed={seed_idx} band={band_desc} allele={}/{}",
                i + 1,
                allele.n_reads,
                i + 1,
                total
            )?;
            out.write_all(&allele.sequence)?;
            writeln!(out)?;
        }
    } else {
        // Single-allele: the clean poa2 engine (via `consensus()`). It seeds on
        // the median-length read (seed-robust; avoids the legacy short-seed
        // under-call), and node fusion + the static-diagonal-union band prevent
        // the repeat "fold" at construction time. That construction-level
        // correctness is what retired the legacy workarounds this branch used to
        // carry -- the seed-sensitivity retry (`consensus_fit_scored`), the
        // banded-truncation unbanded-rebuild retry, and the interior filter --
        // none of which the clean engine needs.
        let result = poa_consensus::consensus(&slices, seed_idx, &config)
            .inspect_err(|e| explain_error(e, n))?;
        // `consensus()` picks the median-length read as seed internally; recompute
        // it here purely for the informational FASTA header.
        let mut order: Vec<usize> = (0..slices.len()).collect();
        order.sort_by_key(|&i| slices[i].len());
        let med = order[order.len() / 2];
        if !args.quiet {
            emit_warnings(&diagnose(&result, &diag_cfg), "consensus");
        }
        writeln!(out, ">consensus reads={n} seed={med} band={band_desc}")?;
        out.write_all(&result.sequence)?;
        writeln!(out)?;
    }

    Ok(())
}

// ── Diagnostic helpers ────────────────────────────────────────────────────────

/// Emit actionable hints to stderr when consensus building fails.
fn explain_error(e: &PoaError, n_reads: usize) {
    match e {
        PoaError::InsufficientDepth { got, min } => {
            eprintln!("poa-consensus: error: only {got} read(s) provided, minimum is {min}");
            if *got > 0 {
                eprintln!(
                    "  hint: use --min-reads {got} to lower the floor (accuracy will \
                     suffer at low depth)"
                );
            }
        }
        PoaError::BandTooNarrow {
            configured,
            required,
        } => {
            eprintln!(
                "poa-consensus: error: band width {configured} is too narrow \
                 (need ≥ {required} for this read set)"
            );
            eprintln!("  hint: try --adaptive-band, or --band-width {required}");
        }
        PoaError::NoSpanningReads {
            left_depth,
            right_depth,
        } => {
            eprintln!(
                "poa-consensus: error: no read spans the full locus \
                 ({left_depth} left-only, {right_depth} right-only reads)"
            );
            eprintln!(
                "  hint: if reads are split into two non-overlapping groups, \
                 use bridged_consensus() to assemble each side separately"
            );
        }
        _ => {} // Display impl covers the remaining variants
    }
    let _ = n_reads; // reserved for future per-depth guidance
}

/// Format and print each message from a [`ConsensusWarnings`] to stderr.
fn emit_warnings(warnings: &poa_consensus::ConsensusWarnings, label: &str) {
    for (is_warning, msg) in warnings.messages(label) {
        let level = if is_warning { "warning" } else { "note" };
        eprintln!("poa-consensus: {level}: {msg}");
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Load all sequences from `path` (or stdin when `path == "-"`), auto-detecting
/// FASTA (`>`) vs FASTQ (`@`) by the first non-empty byte.
fn load_reads(path: &str) -> Result<Vec<Vec<u8>>, Box<dyn std::error::Error>> {
    if path == "-" {
        let stdin = io::stdin();
        let mut buf = BufReader::new(stdin.lock());
        parse_reads(&mut buf)
    } else {
        let file = std::fs::File::open(path)?;
        let mut buf = BufReader::new(file);
        parse_reads(&mut buf)
    }
}

fn parse_reads<R: BufRead>(reader: &mut R) -> Result<Vec<Vec<u8>>, Box<dyn std::error::Error>> {
    // Peek at the first non-whitespace byte to detect format.
    let first = first_byte(reader)?;
    match first {
        Some(b'>') => parse_fasta(reader),
        Some(b'@') => parse_fastq(reader),
        Some(b) => Err(format!(
            "unexpected first byte 0x{b:02x}; expected '>' (FASTA) or '@' (FASTQ)"
        )
        .into()),
        None => Ok(vec![]),
    }
}

fn parse_fasta<R: BufRead>(reader: &mut R) -> Result<Vec<Vec<u8>>, Box<dyn std::error::Error>> {
    use noodles::fasta;
    let mut fa_reader = fasta::io::Reader::new(reader);
    let mut reads = Vec::new();
    for result in fa_reader.records() {
        let record = result?;
        reads.push(record.sequence().as_ref().to_vec());
    }
    Ok(reads)
}

fn parse_fastq<R: BufRead>(reader: &mut R) -> Result<Vec<Vec<u8>>, Box<dyn std::error::Error>> {
    use noodles::fastq;
    let mut fq_reader = fastq::io::Reader::new(reader);
    let mut reads = Vec::new();
    for result in fq_reader.records() {
        let record = result?;
        reads.push(record.sequence().to_vec());
    }
    Ok(reads)
}

/// Peek at the first byte of a `BufRead` without consuming it.
fn first_byte<R: BufRead>(reader: &mut R) -> Result<Option<u8>, Box<dyn std::error::Error>> {
    loop {
        let buf = reader.fill_buf()?;
        if buf.is_empty() {
            return Ok(None);
        }
        let b = buf[0];
        if b == b'\n' || b == b'\r' || b == b' ' {
            reader.consume(1);
            continue;
        }
        return Ok(Some(b));
    }
}
