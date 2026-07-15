use std::io::{self, BufRead, BufReader, Write};
use std::process;

use clap::Parser;

use poa_consensus::{
    AdaptiveAction, AlignmentMode, DiagnoseConfig, PoaConfig, PoaError, SeedSelection, auto_orient,
    consensus_fit_scored, diagnose, select_seed,
};

/// Build a consensus sequence from FASTA or FASTQ reads using Partial Order
/// Alignment.
///
/// Reads are auto-oriented to the seed strand before alignment.  The default
/// seed is chosen automatically via a terminal k-mer anchor heuristic (see
/// `SeedSelection::Auto`); override with --seed.
#[derive(Parser)]
#[command(name = "poa-consensus", version)]
struct Args {
    /// Input FASTA or FASTQ file. Use '-' for stdin.
    #[arg(default_value = "-")]
    input: String,

    /// Seed read index (0-based). Defaults to an automatically chosen seed
    /// (terminal k-mer anchor heuristic; falls back to the longest read).
    #[arg(short, long)]
    seed: Option<usize>,

    /// Fixed band width floor for banded DP (0 = unbanded).  Combined with
    /// adaptive band (default on), effective width is max(band_width, adaptive
    /// formula).  Set to 0 with --no-adaptive-band for fully unbanded DP.
    #[arg(short = 'b', long, default_value_t = 50)]
    band_width: usize,

    /// Disable adaptive band width.  Adaptive band (default on) uses the
    /// formula w = 10 + 0.01 × read_len as a floor alongside --band-width.
    #[arg(long)]
    no_adaptive_band: bool,

    /// Minimum reads required to attempt consensus.
    #[arg(long, default_value_t = 3)]
    min_reads: usize,

    /// Multi-allele mode: output one FASTA record per detected allele.
    #[arg(short = 'm', long)]
    multi: bool,

    /// Global alignment: penalise terminal gaps (use when reads are guaranteed
    /// to span the full locus from identical start/end positions).  The default
    /// is semi-global (free terminal gaps), which is correct for extracted STR
    /// reads that may start or end at slightly different positions.
    #[arg(long)]
    global: bool,

    /// Suppress all warnings and notes on stderr.  The consensus sequence is
    /// still written to stdout.  Errors that prevent consensus building are
    /// always printed regardless of this flag.
    #[arg(short = 'q', long)]
    quiet: bool,
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
            select_seed(&slices, &SeedSelection::Auto)
                .inspect_err(|e| explain_error(e, reads.len()))?
        }
    };

    let oriented = auto_orient(&reads, seed_idx);
    let slices: Vec<&[u8]> = oriented.iter().map(|r| r.as_ref()).collect();

    // ── Config ────────────────────────────────────────────────────────────────
    let adaptive_band = !args.no_adaptive_band;
    let config = PoaConfig {
        band_width: args.band_width,
        adaptive_band,
        min_reads: args.min_reads,
        alignment_mode: if args.global {
            AlignmentMode::Global
        } else {
            AlignmentMode::SemiGlobal
        },
        ..PoaConfig::default()
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
            ..DiagnoseConfig::default()
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
        // Single-allele path with seed-sensitivity retry, then truncation retry.
        //
        // Pass 1 goes through `consensus_fit_scored` rather than the plain
        // `consensus()`: when the pass-1 graph looks noisy/fragmented
        // (`single_support_fraction > 0.3` -- common on genuinely repetitive
        // loci, not just broken ones), it builds a few candidate remedies
        // (tightened coverage floor, a different seed near the read
        // population's median length, MajorityFrequency mode) and empirically
        // scores each against the actual reads, keeping whichever fits best
        // (see `consensus_fit_scored`'s doc comment for the full mechanism and
        // validated limitations -- it is not a complete fix). Otherwise it is
        // exactly equivalent to `consensus()`.
        //
        // After that, `diagnose()` checks whether the result is suspiciously
        // short relative to the median input read length (ratio <
        // DiagnoseConfig::truncation_ratio_threshold, default 0.60).  When that
        // fires AND banded DP was used AND median_read_len <= 5000 bp, rebuild
        // with band_width = 0.  Unbanded forces the traceback to the correct
        // graph endpoint, recovering from the wrong-diagonal failure mode. This
        // targets a different, more severe failure mode (catastrophic band
        // truncation) than the seed-sensitivity retry above (a few missing
        // repeat units), so both stay in place independently.
        //
        // Note: consensus_adaptive is intentionally NOT used here.  That function
        // also checks for multi-allele bubbles, which can fire on het SNPs in
        // non-repetitive controls and produce unexpected results in a caller that
        // expects single-allele output.  `consensus_fit_scored` is the standalone
        // entry point that gives the seed-sensitivity retry without that.
        let was_banded = config.band_width > 0 || config.adaptive_band;
        let fit_scored = consensus_fit_scored(&slices, seed_idx, &config)
            .inspect_err(|e| explain_error(e, n))?;
        if !args.quiet {
            match fit_scored.action {
                AdaptiveAction::PassThrough => {}
                ref a => eprintln!(
                    "poa-consensus: note: consensus looked noisy (fragmented graph evidence); \
                     seed-sensitivity retry chose an alternate rebuild ({a:?})"
                ),
            }
        }
        let mut result = fit_scored.consensuses.into_iter().next().expect(
            "consensus_fit_scored always returns exactly one consensus in single-allele mode",
        );
        if was_banded {
            let diag = diagnose(&result, &DiagnoseConfig::default());
            if let Some(ref t) = diag.truncation_suspected {
                if t.median_read_len <= 5_000 {
                    let mut cfg2 = config.clone();
                    cfg2.band_width = 0;
                    cfg2.adaptive_band = false;
                    if let Ok(r2) = poa_consensus::consensus(&slices, seed_idx, &cfg2) {
                        result = r2;
                    }
                } else if !args.quiet {
                    eprintln!(
                        "poa-consensus: warning: consensus ({} bp) is {:.0}% of median \
                         read length ({} bp) — suspected banded DP truncation on a long \
                         read set; retry with --band-width 0 to correct",
                        t.consensus_len,
                        t.ratio * 100.0,
                        t.median_read_len,
                    );
                }
            }
        }
        if !args.quiet {
            emit_warnings(&diagnose(&result, &DiagnoseConfig::default()), "consensus");
        }
        writeln!(out, ">consensus reads={n} seed={seed_idx} band={band_desc}")?;
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
