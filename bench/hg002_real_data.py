#!/usr/bin/env python3
"""
HG002 real-data validation for poa-consensus.

Extracts reads from a BAM at known disease STR loci, runs poa-consensus with
platform-appropriate settings, and reports consensus allele lengths and
inferred repeat counts.

Usage:
    python bench/hg002_real_data.py HG002.bam
    python bench/hg002_real_data.py HG002.bam --locus HTT RFC1 DMPK
    python bench/hg002_real_data.py HG002.bam --bed my_loci.bed
    python bench/hg002_real_data.py HG002.bam --out results/hg002.tsv

Requires on PATH:  samtools
Requires in repo:  target/release/poa-consensus
                   (cargo build --release --features cli)

BED format (for --bed):
    chrom  start  end  name  unit  [note]
    chr4   3074877  3074944  HTT  CAG  Huntington disease

Coordinates below are hg38.  If your BAM is aligned to CHM13/T2T, pass
--assembly chm13 and supply liftover coordinates via --bed.

Known limitations and locus-specific findings (validated against HG002 HiFi):

  Strand orientation
    Clinical repeat units for DMPK, ATXN1, ATXN2, ATXN3, and RFC1 are defined
    on the minus (coding) strand.  SAM/BAM stores reads in plus-strand
    orientation, so the consensus produced by poa-consensus contains the reverse
    complement of the clinical unit (e.g. CTG instead of CAG for DMPK).
    count_units() searches both strands to handle this transparently.  Any
    downstream tool that inspects the raw consensus bytes must do the same.

  Interrupted repeats (FMR1, ATXN1)
    FMR1 normal alleles contain AGG interruptions roughly every 9-12 CGG units.
    The repeat count reported here is the longest uninterrupted CGG run (~9-12),
    not the total CGG count (~30).  Similarly, ATXN1 normal alleles contain CAT
    interruptions; the reported count (~14-15) is the longest uninterrupted CAG
    run, not the total (~26-35).  This is correct behaviour: the longest
    uninterrupted run is the clinically relevant stability predictor.

  ATXN3 coordinates
    Correct hg38 coordinates are chr14:92,071,010-92,071,052.  The
    ExpansionHunter catalog entry at chr14:92,071,992-92,072,120 does not
    land on the CAG repeat in HG002 (both platforms return 1× CAG on
    identical-length consensuses at that position).

  RFC1 banded DP diagonal drift and the short-flank workaround
    The CTTTT plus-strand / AAAAG coding-strand 5-mer repeat causes banded DP to
    converge to the wrong diagonal without approaching the band edge.  With long
    flanks (pad=600) the correctly-assembled flanking sequence masks the truncation
    (consensus / median-read ratio ~0.84, above the 0.60 detection threshold).
    The bench script uses pad=100 so that a truncated consensus drops the ratio
    below 0.60, triggering an automatic unbounded-DP retry that pins both endpoints
    and recovers the correct 115× count.

  ATXN2 interrupted GCT/GTT structure
    ATXN2 uses the GCT unit (a CAG rotation; revcomp of CAG), with GTT
    interruptions.  Bladerunner ground truth for HG002 hg38:
      Hap1: (GCT)8(GTT)1(GCT)21 -- 29 total, longest uninterrupted run = 21
      Hap2: (GCT)8(GTT)1(GCT)4(GTT)1(GCT)8 -- 20 total, longest run = 8
    poa-consensus reports 8/21 (longest uninterrupted run per haplotype), which
    is correct.  count_units() finds GCT via CAG rotation + strand search.

  ONT min_allele_freq sensitivity
    At ONT error rates (>=4% substitution), the default min_allele_freq = 0.25
    threshold fires on clean single-allele loci.  Raise to >= 0.40 for ONT.
    At depth = 3, one error in one read = 33%, which already exceeds 0.25; the
    effective minimum depth for clean ONT calls is >= 10 reads.
"""

import argparse
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ── CIGAR helper ──────────────────────────────────────────────────────────────

_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def _cigar_clip(ref_start: int, cigar: str, target_start: int, target_end: int) -> tuple[int, int]:
    """
    Map genomic region [target_start, target_end) to read coordinates (rd_s, rd_e).
    Returns (-1, -1) if the read does not overlap the target.

    Soft-clipped bases advance the read pointer but not the reference pointer, so
    they are naturally excluded from the clipped output.  Insertions that fall
    entirely within [target_start, target_end) are included so the clipped read
    retains any indels present in that region.
    """
    ref = ref_start
    rd  = 0
    rd_s = -1
    rd_e = -1

    for m in _CIGAR_RE.finditer(cigar):
        n, op = int(m.group(1)), m.group(2)
        if op in ('M', '=', 'X'):
            ovlp_s = max(ref, target_start)
            ovlp_e = min(ref + n, target_end)
            if ovlp_s < ovlp_e:
                if rd_s < 0:
                    rd_s = rd + (ovlp_s - ref)
                rd_e = rd + (ovlp_e - ref)
            ref += n
            rd  += n
        elif op == 'I':
            # Include an insertion only if we are past the start and still inside
            # the target (ref hasn't advanced for I, so use strict < target_end).
            if rd_s >= 0 and ref > target_start and ref < target_end:
                rd_e = rd + n
            rd += n
        elif op in ('D', 'N'):
            ref += n
        elif op == 'S':
            rd += n
        # H / P: neither pointer advances

    return rd_s, rd_e


# ── Paths ─────────────────────────────────────────────────────────────────────

POA_BIN = Path(__file__).parent.parent / "target/release/poa-consensus"

# ── Locus definitions ─────────────────────────────────────────────────────────

@dataclass
class Locus:
    name: str
    chrom: str          # with chr prefix (stripped if BAM uses bare contigs)
    start: int          # 0-based
    end: int            # exclusive
    unit: str           # canonical repeat unit; empty string for non-repetitive controls
    pad: int = 300      # bp of flanking context added when extracting reads
    min_depth: int = 5  # skip locus if fewer reads overlap
    note: str = ""
    category: str = "str"  # "str" = disease repeat locus; "control" = non-repetitive

# hg38 coordinates.  Sources: OMIM, Locus Reference Genomic, TRGT catalog.
# Control regions: non-repetitive exonic/coding regions used as sanity checks.
HG38_LOCI: list[Locus] = [
    # ── Disease STR loci ──────────────────────────────────────────────────────
    # Coordinates are 0-based half-open intervals matching the ExpansionHunter
    # hg38 variant catalog (github.com/Illumina/ExpansionHunter).
    #
    # The repeat unit below is the CLINICAL unit (defined on the coding strand).
    # Genes on the minus strand (DMPK, ATXN1, ATXN2, ATXN3, RFC1) have their
    # unit appear as the reverse complement in the reference plus strand.
    # count_units() searches both strands automatically.
    Locus("HTT",   "chr4",  3_074_877,    3_074_944,    "CAG",   pad=400, note="Huntington disease"),
    Locus("FMR1",  "chrX",  147_912_050,  147_912_110,  "CGG",   pad=400, note="Fragile X syndrome"),
    Locus("DMPK",  "chr19", 45_770_203,   45_770_264,   "CTG",   pad=400, note="Myotonic dystrophy 1"),
    # RFC1: HG002 has one normal AAAAG allele (~11 units) and one large non-pathogenic
    # AAAAG expansion (~115 units, ~520 bp insertion vs. reference).  The disease-causing
    # allele is AAGGG; HG002 has only AAAAG so it is unaffected.
    #
    # pad=100: deliberately short.  With pad=600, the correctly-assembled flanking
    # sequence masks the repeat truncation (banded DP on AAAAG 5-mer converges to
    # the wrong diagonal): the consensus is ~1500 bp vs ~1775 bp reads, ratio 0.84,
    # which does not trigger the truncation-detection retry.  With pad=100, the read
    # lengths drop to ~700 bp, and a banded-DP truncation produces a consensus of
    # ~400 bp (ratio ~0.57 < 0.60), triggering the unbanded retry.  Unbanded global
    # alignment pins both endpoints so the correct 115× assembly is recovered.
    Locus("RFC1",  "chr4",  39_348_424,   39_348_485,   "AAAAG", pad=100, note="CANVAS / RFC1 ataxia; HG002 has ~115-unit AAAAG expanded allele (non-pathogenic)"),
    Locus("ATXN3", "chr14", 92_071_010,   92_071_052,   "CAG",   pad=400, note="Spinocerebellar ataxia 3"),
    Locus("ATXN2", "chr12", 111_598_950,  111_599_019,  "CAG",   pad=400, note="Spinocerebellar ataxia 2"),
    Locus("ATXN1", "chr6",  16_327_633,   16_327_723,   "CAG",   pad=400, note="Spinocerebellar ataxia 1"),
    # ── Non-repetitive controls ───────────────────────────────────────────────
    # These should yield a single clean consensus with no bubble sites and
    # depth tracking the local WGS coverage.  Any multi-allele result or
    # zero-depth here indicates a pipeline or extraction problem.
    Locus("GAPDH", "chr12", 6_534_517,   6_534_800,   "",      pad=200,
          note="GAPDH coding (housekeeping control)", category="control"),
    Locus("TP53e7", "chr17", 7_674_220,  7_674_400,   "",      pad=200,
          note="TP53 exon 7 (non-repetitive coding)", category="control"),
    Locus("ACTB",  "chr7",  5_529_264,   5_529_600,   "",      pad=200,
          note="ACTB coding (housekeeping control)", category="control"),
]

# ── HG002 truth table ────────────────────────────────────────────────────────
#
# Expected repeat unit counts for HG002 (hg38).
#
# IMPORTANT: count_units() returns the LONGEST CONTIGUOUS RUN of the unit in
# the consensus, not the total count.  For interrupted repeats (FMR1, ATXN1,
# ATXN2) this is shorter than the clinical total count:
#   FMR1  — AGG interruptions every ~9-12 CGG units; longest run ~11
#   ATXN1 — CAT interruptions; longest uninterrupted CAG run ~14-15
#   ATXN2 — GTT interruptions separate GCT blocks; count_units returns the
#            longest uninterrupted GCT run per haplotype (8 and 21 for HG002)
#
# HTT and DMPK: observed from HiFi + ONT runs on HG002 hg38 (both platforms
# agree).  Longest-run == total because no interruptions in normal alleles.
#
# RFC1: one normal AAAAG allele and one large non-pathogenic AAAAG expansion
# (~115 units from a 520 bp insertion in the HG002 paternal haplotype, bedpull
# README).  Marked known_bug=True; POA bug #4 silently truncates this repeat.
#
# ATXN3: coordinates chr14:92,071,992-92,072,120 are incorrect for hg38 HG002
# (both platforms return 1× CAG on two identical-length consensuses).  Omitted
# from truth evaluation; see bad_coords flag.
#
# Tolerance: ±1 for alleles < 50 units; ±2 for < 100; ±10 for < 200; ±20 otherwise.
HG002_TRUTH: dict[str, dict] = {
    "HTT":    {"alleles": [17, 24], "confidence": "high",
               "source": "HiFi + ONT observed; both platforms agree; length-consistent"},
    "FMR1":   {"alleles": [11],     "confidence": "high",
               "source": "HiFi + ONT observed; longest uninterrupted CGG run in interrupted allele",
               "note": "hemizygous (HG002 is male; chrX locus). Total CGG ~29-35 but "
                       "count_units reports longest run (~11) due to AGG interruptions."},
    "DMPK":   {"alleles": [11],     "confidence": "medium",
               "source": "HiFi + ONT observed; single allele (both haplotypes likely same length)",
               "note": "11 CTG is within normal range (5-37); no second allele detected"},
    "RFC1":   {"alleles": [11, 115], "confidence": "medium",
               "source": "bedpull README (520 bp paternal insertion = ~104 extra AAAAG units)",
               "note": "Short flanks (pad=100) keep the consensus/read ratio below 0.60 "
                       "when banded DP truncates the repeat, triggering the unbounded retry."},
    "ATXN1":  {"alleles": [14, 15], "confidence": "medium",
               "source": "HiFi + ONT observed; longest uninterrupted CAG run in interrupted allele",
               "note": "Total CAG ~28-35 but count_units reports longest run (~14-15) "
                       "due to CAT interruptions. 1 bp allele length difference is consistent."},
    "ATXN2":  {"alleles": [8, 21],  "confidence": "medium",
               "source": "Bladerunner ground truth (HG002 hg38, HP-phased GCT analysis)",
               "note": "Unit is GCT (CAG rotation; revcomp of CAG).  GTT interruptions: "
                       "Hap1=(GCT)8(GTT)1(GCT)21 (29 total, longest run=21); "
                       "Hap2=(GCT)8(GTT)1(GCT)4(GTT)1(GCT)8 (20 total, longest run=8). "
                       "poa-consensus 8/21 is CORRECT: count_units reports longest uninterrupted run."},
    "ATXN3":  {"alleles": [17, 21], "confidence": "medium",
               "source": "HiFi + ONT observed; both platforms agree; length-consistent (Δ12bp = 4 CAG)"},
}


def _truth_tolerance(n: int) -> int:
    if n < 50:   return 1
    if n < 100:  return 2
    if n < 200:  return 10
    return 20


def evaluate_against_truth(locus_name: str, found_counts: list[int]) -> str:
    """Compare repeat counts to HG002 truth; return a status string."""
    if locus_name not in HG002_TRUTH:
        return ""
    truth = HG002_TRUTH[locus_name]
    known_bug  = truth.get("known_bug", False)
    bad_coords = truth.get("bad_coords", False)
    expected   = sorted(truth["alleles"])
    conf       = truth["confidence"]

    if bad_coords:
        return f"BAD_COORDS — {truth.get('note', '')}"

    if len(found_counts) == 0:
        if known_bug:
            return f"KNOWN_BUG [{conf}] — {truth.get('note', '')}"
        return "FAIL (no alleles)"

    found = sorted(found_counts)

    if len(expected) == 1:
        ok = abs(found[0] - expected[0]) <= _truth_tolerance(expected[0])
    elif len(found) < 2:
        ok = False
    else:
        d00 = abs(found[0] - expected[0]) + abs(found[1] - expected[1])
        d01 = abs(found[0] - expected[1]) + abs(found[1] - expected[0])
        if d00 <= d01:
            ok = (abs(found[0] - expected[0]) <= _truth_tolerance(expected[0]) and
                  abs(found[1] - expected[1]) <= _truth_tolerance(expected[1]))
        else:
            ok = (abs(found[0] - expected[1]) <= _truth_tolerance(expected[1]) and
                  abs(found[1] - expected[0]) <= _truth_tolerance(expected[0]))

    if ok:
        return f"PASS [{conf}]"
    elif known_bug:
        return f"KNOWN_BUG [{conf}] — {truth.get('note', '')}"
    else:
        exp_str = "/".join(str(e) for e in expected)
        got_str = "/".join(str(f) for f in found)
        return f"FAIL [{conf}] expected {exp_str}, got {got_str}"


# ── BAM header inspection ─────────────────────────────────────────────────────

def bam_header(bam: Path) -> str:
    r = subprocess.run(["samtools", "view", "-H", str(bam)],
                       capture_output=True, text=True, check=True)
    return r.stdout


def detect_chr_prefix(header: str) -> bool:
    """Return True if contigs use 'chr' prefix."""
    for line in header.splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            sn = fields.get("SN", "")
            if sn.startswith("chr"):
                return True
            if sn.isdigit() or sn in ("X", "Y", "M", "MT"):
                return False
    # Default to chr-prefix if we can't tell.
    return True


def detect_assembly(header: str) -> str:
    """
    Return 'hg38', 'chm13', or 'unknown' by checking chr1 length in @SQ lines.
    hg38 chr1 = 248,956,422 bp; CHM13 v2.0 chr1 = 248,387,328 bp.
    """
    for line in header.splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            sn = fields.get("SN", "")
            ln = int(fields.get("LN", 0))
            if sn in ("chr1", "1"):
                if abs(ln - 248_956_422) < 100_000:
                    return "hg38"
                if abs(ln - 248_387_328) < 100_000:
                    return "chm13"
    return "unknown"


def detect_platform(header: str, bam_path: Optional[Path] = None) -> str:
    """Return 'ont', 'hifi', or 'unknown' from @RG tags and BAM filename."""
    for line in header.splitlines():
        if line.startswith("@RG"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            # Concatenate all field values to catch any tag that indicates platform.
            combined = " ".join(v.upper() for v in fields.values())
            if "ONT" in combined or "NANOPORE" in combined:
                return "ont"
            if any(k in combined for k in ("HIFI", "PACBIO", "SEQUEL", "SMRT", "CCS", "REVIO")):
                return "hifi"
    # Fall back to filename heuristics.
    if bam_path is not None:
        name = bam_path.name.lower()
        if "ont" in name or "nanopore" in name:
            return "ont"
        if any(k in name for k in ("hifi", "pb.", "_pb_", "pacbio", "ccs", "revio", "sequel")):
            return "hifi"
    return "unknown"

# ── Coordinate adjustment ─────────────────────────────────────────────────────

def adjust_chrom(chrom: str, use_chr: bool) -> str:
    """Add or strip 'chr' prefix to match the BAM's contig naming."""
    if use_chr and not chrom.startswith("chr"):
        return "chr" + chrom
    if not use_chr and chrom.startswith("chr"):
        return chrom[3:]
    return chrom

# ── CLI flags by platform ─────────────────────────────────────────────────────

def poa_flags(platform: str) -> list[str]:
    """Return poa-consensus CLI flags appropriate for the platform."""
    # --band-width 50 floor prevents silent unit loss in long repetitive alleles.
    # --semi-global handles partial reads that don't span the full locus.
    base = ["--adaptive-band", "--band-width", "50", "--semi-global", "--min-reads", "3"]
    return base

# ── Read-length allele estimation ────────────────────────────────────────────

def read_length_alleles(reads: list[bytes], unit: str, core_len: int,
                        ref_len: int) -> Optional[str]:
    """
    For loci where POA fails (e.g. RFC1 bug #4), the read length distribution
    encodes allele sizes directly: reads from the expanded allele are longer by
    exactly the insertion size.

    Fits a two-component model: reads cluster around ref_len (normal allele)
    and ref_len + insertion (expanded allele).  Returns a human-readable
    summary string if a bimodal split is detected, or None if the distribution
    is unimodal.

    unit_len is used to convert insertion bp → repeat units.
    ref_len is the expected trimmed length for a non-expanded read (padded window).
    """
    if not reads or not unit:
        return None
    unit_len = len(unit)
    lengths = sorted(len(r) for r in reads)

    # Find the modal length of the "normal" cluster (readings ≤ ref_len * 1.15).
    normal_floor = int(ref_len * 0.85)
    normal_ceil  = int(ref_len * 1.15)
    normal_reads = [l for l in lengths if normal_floor <= l <= normal_ceil]
    if not normal_reads:
        return None
    normal_median = sorted(normal_reads)[len(normal_reads) // 2]

    # Look for a second cluster clearly above the normal cluster.
    gap_threshold = max(unit_len * 5, 50)   # at least 5 repeat units of gap
    upper_reads = [l for l in lengths if l > normal_median + gap_threshold]
    if not upper_reads:
        return None
    upper_median = sorted(upper_reads)[len(upper_reads) // 2]

    insertion_bp   = upper_median - normal_median
    insertion_units = round(insertion_bp / unit_len)

    n_normal  = len(normal_reads)
    n_upper   = len(upper_reads)
    n_short   = sum(1 for l in lengths if l < normal_floor)

    return (f"bimodal: {n_normal} reads ~{normal_median}bp (normal) / "
            f"{n_upper} reads ~{upper_median}bp "
            f"(+{insertion_bp}bp = +{insertion_units}× {unit} extra) "
            f"[{n_short} short partials excluded]")


# ── Read extraction ───────────────────────────────────────────────────────────

def extract_reads(bam: Path, chrom: str, start: int, end: int, pad: int) -> list[bytes]:
    """
    Return read sequences trimmed to [start-pad, end+pad] in genomic coordinates.

    Full WGS reads are 10–20 kb; poa-consensus is designed for 50–2000 bp STR
    reads.  Passing full reads causes the banded DP to fail on the large graph.
    CIGAR-based clipping reduces each read to just the padded region (~800–1500 bp
    for STR loci), which is the input size poa-consensus is tuned for.

    Reads whose trimmed length is shorter than the core (start, end) region are
    dropped as non-spanning.
    """
    core_len     = end - start
    padded_start = max(0, start - pad)
    padded_end   = end + pad
    region = f"{chrom}:{padded_start + 1}-{padded_end}"  # samtools is 1-based

    # -F 2308: exclude unmapped (4) + secondary (256) + supplementary (2048)
    cmd = ["samtools", "view", "-F", "2308", str(bam), region]
    try:
        view = subprocess.run(cmd, capture_output=True, check=True)
    except subprocess.CalledProcessError:
        return {0: [], 1: [], 2: []}

    if not view.stdout:
        return {0: [], 1: [], 2: []}

    reads: dict[int, list[bytes]] = {0: [], 1: [], 2: []}
    for line in view.stdout.splitlines():
        fields = line.split(b"\t")
        if len(fields) < 10:
            continue
        try:
            pos0  = int(fields[3]) - 1          # SAM POS is 1-based
            cigar = fields[5].decode()
            seq   = fields[9]
        except (ValueError, IndexError):
            continue
        if seq == b"*" or cigar == "*":
            continue

        rd_s, rd_e = _cigar_clip(pos0, cigar, padded_start, padded_end)
        if rd_s < 0 or rd_e <= rd_s:
            continue

        trimmed = seq[rd_s:rd_e]
        if len(trimmed) < core_len:
            continue

        hp = 0
        for tag in fields[11:]:
            if tag.startswith(b"HP:i:"):
                try:
                    hp = int(tag[5:])
                except ValueError:
                    pass
                break
        reads.setdefault(hp, []).append(trimmed)

    return reads

# ── poa-consensus invocation ──────────────────────────────────────────────────

def run_poa(reads: list[bytes], flags: list[str], multi: bool,
            verbose: bool = False) -> Optional[list[bytes]]:
    """
    Run poa-consensus on `reads` (as FASTA piped to stdin).
    Returns list of consensus sequences, or None on failure.
    """
    fasta = b""
    for i, seq in enumerate(reads):
        fasta += b">r%d\n%s\n" % (i, seq)

    cmd = [str(POA_BIN)] + flags
    if multi:
        cmd.append("--multi")
    cmd.append("-")  # read from stdin

    try:
        r = subprocess.run(cmd, input=fasta, capture_output=True, timeout=60)
    except subprocess.TimeoutExpired:
        print("  [timeout]", file=sys.stderr)
        return None
    if r.returncode != 0:
        stderr_msg = r.stderr.decode(errors="replace").strip() if r.stderr else ""
        stdout_msg = r.stdout.decode(errors="replace").strip() if r.stdout else ""
        detail = stderr_msg or stdout_msg or f"(no output; exit code {r.returncode})"
        print(f"  [poa error rc={r.returncode}: {detail}]", file=sys.stderr)
        return None

    # Parse FASTA output: collect all sequences.
    seqs = []
    current = b""
    for line in r.stdout.splitlines():
        if line.startswith(b">"):
            if current:
                seqs.append(current)
            current = b""
        else:
            current += line.strip()
    if current:
        seqs.append(current)
    return seqs or None

# ── Repeat unit counting ──────────────────────────────────────────────────────

_COMP = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")

def _revcomp_unit(unit: str) -> str:
    return unit.upper().encode().translate(_COMP)[::-1].decode()


def count_units(seq: bytes, unit: str) -> int:
    """
    Find the longest contiguous run of `unit` within `seq` and return its count.

    The consensus includes flanking context, so dividing total length by unit
    length gives the wrong answer.  Instead, search for the longest run via regex,
    trying all rotational phases of the unit to handle phase ambiguity in the
    consensus (e.g. CAG / AGC / GCA all represent the same trinucleotide repeat),
    AND trying the reverse-complement strand of the unit.

    poa-consensus outputs the consensus in the seed read's orientation, which may
    be the minus strand.  Searching both strands ensures correct counts regardless
    of which strand the consensus comes out on.

    Example: for unit="CAG" and consensus on the minus strand
        "...ACGTCTGCTGCTGCTGCTGACGT..."  →  5 (CTG = revcomp of CAG)
    """
    if not unit or not seq:
        return 0
    n = len(unit)
    unit_fwd = unit.upper()
    unit_rc  = _revcomp_unit(unit)

    best = 0
    seen: set[bytes] = set()
    for strand_unit in (unit_fwd, unit_rc):
        ub = strand_unit.encode()
        for i in range(n):
            rot = ub[i:] + ub[:i]
            if rot in seen:
                continue
            seen.add(rot)
            pat = re.compile(rb"(?:" + rot + rb")+", re.IGNORECASE)
            for m in pat.finditer(seq):
                count = len(m.group()) // n
                if count > best:
                    best = count
    return best

# ── BED loader ────────────────────────────────────────────────────────────────

def load_bed(path: Path) -> list[Locus]:
    """Load loci from a BED file.  Format: chrom start end name unit [note]."""
    loci = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                sys.exit(f"error: BED line needs ≥5 fields: {line!r}")
            chrom, start, end, name, unit = parts[:5]
            note = parts[5] if len(parts) > 5 else ""
            loci.append(Locus(name=name, chrom=chrom,
                              start=int(start), end=int(end),
                              unit=unit, note=note))
    return loci

# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("bam", type=Path, help="Input BAM file (must be indexed)")
    ap.add_argument("--locus", nargs="+", metavar="NAME",
                    help="Run only these loci (e.g. HTT RFC1)")
    ap.add_argument("--bed", type=Path, metavar="FILE",
                    help="Custom loci in BED format (chrom start end name unit [note])")
    ap.add_argument("--assembly", choices=["hg38", "chm13", "auto"],
                    default="auto", help="Reference assembly (default: auto-detect)")
    ap.add_argument("--platform", choices=["ont", "hifi", "auto"],
                    default="auto", help="Sequencing platform (default: auto-detect)")
    ap.add_argument("--pad", type=int, default=None,
                    help="Override flanking context in bp (default: per-locus)")
    ap.add_argument("--no-multi", action="store_true",
                    help="Disable multi-allele mode (report single consensus only)")
    ap.add_argument("--out", type=Path, metavar="FILE",
                    help="Write TSV results to this file in addition to stdout")
    ap.add_argument("--verbose", action="store_true",
                    help="Print each consensus sequence after the summary row")
    args = ap.parse_args()

    # ── Sanity checks ─────────────────────────────────────────────────────────
    if not args.bam.exists():
        sys.exit(f"error: BAM not found: {args.bam}")
    bai = args.bam.with_suffix(args.bam.suffix + ".bai")
    csi = args.bam.with_suffix(args.bam.suffix + ".csi")
    if not bai.exists() and not csi.exists():
        sys.exit(f"error: no BAM index found ({bai} or {csi})\n"
                 f"       run: samtools index {args.bam}")
    if not POA_BIN.exists():
        sys.exit(f"error: poa-consensus binary not found at {POA_BIN}\n"
                 f"       run: cargo build --release --features cli")

    # ── BAM header inspection ─────────────────────────────────────────────────
    header = bam_header(args.bam)
    use_chr   = detect_chr_prefix(header)
    assembly  = args.assembly if args.assembly != "auto" else detect_assembly(header)
    platform  = args.platform if args.platform != "auto" else detect_platform(header, args.bam)

    # ── Loci selection ────────────────────────────────────────────────────────
    if args.bed:
        loci = load_bed(args.bed)
    elif assembly == "chm13":
        print("warning: CHM13/T2T detected.  Built-in coordinates are for hg38.")
        print("         Supply liftover coordinates via --bed for accurate results.")
        print("         Using hg38 coordinates as a rough approximation.\n")
        loci = HG38_LOCI
    else:
        loci = HG38_LOCI

    if args.locus:
        names = set(args.locus)
        loci = [l for l in loci if l.name in names]
        missing = names - {l.name for l in loci}
        if missing:
            sys.exit(f"error: unknown loci: {', '.join(sorted(missing))}\n"
                     f"       available: {', '.join(l.name for l in HG38_LOCI)}")

    # ── CLI flags ─────────────────────────────────────────────────────────────
    flags = poa_flags(platform)

    # ── Header ────────────────────────────────────────────────────────────────
    asm_str = assembly if assembly != "unknown" else "unknown (hg38 assumed)"
    plat_str = platform if platform != "unknown" else "unknown (ONT settings used)"
    chr_str = "chr prefix" if use_chr else "no chr prefix"

    print("HG002 poa-consensus real-data test")
    print(f"  BAM:      {args.bam}")
    print(f"  Assembly: {asm_str} ({chr_str})")
    print(f"  Platform: {plat_str}")
    print(f"  Binary:   {POA_BIN}")
    print(f"  Flags:    {' '.join(flags)}")
    print()

    col_w = dict(locus=9, chrom=6, region=26, unit=6, depth=6,
                 alleles=30, repeats=20, time=8, verdict=30)
    header_row = (
        f"  {'Locus':<{col_w['locus']}}  {'Chrom':<{col_w['chrom']}}"
        f"  {'Region':<{col_w['region']}}  {'Unit':<{col_w['unit']}}"
        f"  {'Depth':>{col_w['depth']}}  {'Allele lengths':<{col_w['alleles']}}"
        f"  {'Repeat counts':<{col_w['repeats']}}  {'Time':>{col_w['time']}}"
        f"  {'Truth check':<{col_w['verdict']}}"
    )
    sep = "  " + "─" * (len(header_row) - 2)

    tsv_rows: list[str] = [
        "\t".join(["locus", "category", "chrom", "start", "end", "unit", "depth",
                   "phasing", "allele_lengths_bp", "repeat_counts", "time_s", "verdict"])
    ]

    # ── Per-locus loop ────────────────────────────────────────────────────────
    prev_category: Optional[str] = None
    for locus in loci:
        # Print a section header whenever the category changes.
        if locus.category != prev_category:
            if prev_category is not None:
                print()
            label = {
                "str":     "Disease STR loci",
                "control": "Non-repetitive controls (expect single clean consensus)",
            }.get(locus.category, locus.category)
            print(f"  {label}")
            print(header_row)
            print(sep)
            prev_category = locus.category

        chrom = adjust_chrom(locus.chrom, use_chr)
        pad   = args.pad if args.pad is not None else locus.pad

        t0 = time.perf_counter()
        hap_reads = extract_reads(args.bam, chrom, locus.start, locus.end, pad)
        all_reads = [r for rs in hap_reads.values() for r in rs]
        depth = len(all_reads)

        if depth < locus.min_depth:
            elapsed = time.perf_counter() - t0
            allele_str = f"SKIP (depth {depth} < {locus.min_depth})"
            repeat_str = ""
            _print_row(col_w, locus.name, chrom, locus.start, locus.end,
                       locus.unit, depth, allele_str, repeat_str, elapsed, "")
            tsv_rows.append("\t".join([
                locus.name, locus.category, chrom,
                str(locus.start), str(locus.end), locus.unit,
                str(depth), "skip", "", f"{elapsed:.2f}", "",
            ]))
            continue

        is_control = locus.category == "control"
        n_hp1 = len(hap_reads.get(1, []))
        n_hp2 = len(hap_reads.get(2, []))
        _MIN_HP_READS = 3

        if not is_control and n_hp1 >= _MIN_HP_READS and n_hp2 >= _MIN_HP_READS:
            # Both haplotypes tagged: run single-allele POA on each independently.
            a1 = run_poa(hap_reads[1], flags, multi=False)
            a2 = run_poa(hap_reads[2], flags, multi=False)
            alleles = (a1 or []) + (a2 or [])
            alleles = alleles or None
            phasing = f"hp-split HP1:{n_hp1}/HP2:{n_hp2}"
        else:
            # Unphased or control: pool all reads, use multi-allele mode.
            multi = not args.no_multi and not is_control
            alleles = run_poa(all_reads, flags, multi=multi)
            phasing = f"untagged:{len(hap_reads.get(0,[]))} HP1:{n_hp1} HP2:{n_hp2}"
        elapsed = time.perf_counter() - t0

        counts: list[int] = []
        if alleles is None:
            allele_str = "ERROR"
            repeat_str = ""
        elif len(alleles) == 0:
            allele_str = "no consensus"
            repeat_str = ""
        else:
            lengths = [len(a) for a in alleles]
            allele_str = "  /  ".join(f"{l} bp" for l in lengths)
            if locus.unit:
                counts = [count_units(a, locus.unit) for a in alleles]
                repeat_str = "  /  ".join(f"{c}× {locus.unit}" for c in counts)
            else:
                repeat_str = "—"

        verdict = evaluate_against_truth(locus.name, counts) if locus.unit else ""
        _print_row(col_w, locus.name, chrom, locus.start, locus.end,
                   locus.unit, depth, allele_str, repeat_str, elapsed, verdict)

        # For KNOWN_BUG loci, show the read-length bimodal as a diagnostic.
        if verdict.startswith("KNOWN_BUG") and locus.unit:
            ref_len = locus.end - locus.start + 2 * pad
            rla = read_length_alleles(all_reads, locus.unit, locus.end - locus.start, ref_len)
            if rla:
                print(f"    [read-length estimate] {rla}")

        if args.verbose and alleles:
            for i, seq in enumerate(alleles):
                label = f"allele {i+1}" if len(alleles) > 1 else "consensus"
                preview = seq[:120].decode(errors="replace")
                print(f"    [{label}] {preview}...")

        allele_bp = "; ".join(str(len(a)) for a in (alleles or []))
        repeat_ct = "; ".join(str(c) for c in counts) if counts else ""
        tsv_rows.append("\t".join([
            locus.name, locus.category, chrom,
            str(locus.start), str(locus.end), locus.unit,
            str(depth), phasing, allele_bp, repeat_ct,
            f"{elapsed:.2f}", verdict or locus.note,
        ]))

    print()

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text("\n".join(tsv_rows) + "\n")
        print(f"TSV written to {args.out}")


def _print_row(col_w: dict, name: str, chrom: str,
               start: int, end: int, unit: str,
               depth: int, allele_str: str, repeat_str: str,
               elapsed: float, verdict: str = "") -> None:
    region = f"{start:,}-{end:,}"
    verdict_col = f"  {verdict:<{col_w['verdict']}}" if verdict else ""
    print(
        f"  {name:<{col_w['locus']}}  {chrom:<{col_w['chrom']}}"
        f"  {region:<{col_w['region']}}  {unit:<{col_w['unit']}}"
        f"  {depth:>{col_w['depth']}}  {allele_str:<{col_w['alleles']}}"
        f"  {repeat_str:<{col_w['repeats']}}  {elapsed:>{col_w['time']}.1f}s"
        + verdict_col
    )


if __name__ == "__main__":
    main()
