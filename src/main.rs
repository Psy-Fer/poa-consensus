use std::io::{self, BufRead, BufReader, Write};
use std::process;

use clap::Parser;

use poa_consensus::{AlignmentMode, PoaConfig, auto_orient};

/// Build a consensus sequence from FASTA or FASTQ reads using Partial Order
/// Alignment.
///
/// Reads are auto-oriented to the seed strand before alignment.  The default
/// seed is the median-length read; override with --seed.
#[derive(Parser)]
#[command(name = "poa-consensus", version)]
struct Args {
    /// Input FASTA or FASTQ file. Use '-' for stdin.
    #[arg(default_value = "-")]
    input: String,

    /// Seed read index (0-based). Defaults to the median-length read.
    #[arg(short, long)]
    seed: Option<usize>,

    /// Fixed band width for banded DP (0 = unbanded; set to ≥ expected length
    /// variation between reads).
    #[arg(short = 'b', long, default_value_t = 0)]
    band_width: usize,

    /// Enable adaptive band width: w = 10 + 0.01 × read_len (recommended
    /// above 1 kb).
    #[arg(long)]
    adaptive_band: bool,

    /// Minimum reads required to attempt consensus.
    #[arg(long, default_value_t = 3)]
    min_reads: usize,

    /// Multi-allele mode: output one FASTA record per detected allele.
    #[arg(short = 'm', long)]
    multi: bool,

    /// Semi-global alignment: free terminal gaps (recommended when reads do
    /// not all span the full locus).
    #[arg(long)]
    semi_global: bool,
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
        None => median_seed(&reads),
    };

    let oriented = auto_orient(&reads, seed_idx);
    let slices: Vec<&[u8]> = oriented.iter().map(|r| r.as_ref()).collect();

    // ── Config ────────────────────────────────────────────────────────────────
    let config = PoaConfig {
        band_width: args.band_width,
        adaptive_band: args.adaptive_band,
        min_reads: args.min_reads,
        alignment_mode: if args.semi_global {
            AlignmentMode::SemiGlobal
        } else {
            AlignmentMode::Global
        },
        ..PoaConfig::default()
    };

    let band_desc = if args.adaptive_band {
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
        let alleles = poa_consensus::consensus_multi(&slices, seed_idx, &config)?;
        let total = alleles.len();
        for (i, allele) in alleles.iter().enumerate() {
            writeln!(
                out,
                ">allele_{} reads={n} seed={seed_idx} band={band_desc} allele={}/{}",
                i + 1,
                i + 1,
                total
            )?;
            out.write_all(&allele.sequence)?;
            writeln!(out)?;
        }
    } else {
        let result = poa_consensus::consensus(&slices, seed_idx, &config)?;
        writeln!(out, ">consensus reads={n} seed={seed_idx} band={band_desc}")?;
        out.write_all(&result.sequence)?;
        writeln!(out)?;
    }

    Ok(())
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

/// Return the index of the median-length read.
fn median_seed(reads: &[Vec<u8>]) -> usize {
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_unstable_by_key(|&i| reads[i].len());
    order[order.len() / 2]
}
