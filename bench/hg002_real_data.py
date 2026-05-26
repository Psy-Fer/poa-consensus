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
"""

import argparse
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

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
    Locus("HTT",   "chr4",  3_074_877,   3_074_944,   "CAG",   pad=400, note="Huntington disease"),
    Locus("FMR1",  "chrX",  147_912_048, 147_912_126, "CGG",   pad=400, note="Fragile X syndrome"),
    Locus("DMPK",  "chr19", 46_272_986,  46_273_024,  "CTG",   pad=400, note="Myotonic dystrophy 1"),
    Locus("RFC1",  "chr4",  39_348_425,  39_348_494,  "AAGGG", pad=400, note="CANVAS / RFC1 ataxia"),
    Locus("ATXN3", "chr14", 92_071_968,  92_072_115,  "CAG",   pad=400, note="Spinocerebellar ataxia 3"),
    Locus("ATXN2", "chr12", 112_036_756, 112_036_823, "CAG",   pad=400, note="Spinocerebellar ataxia 2"),
    Locus("ATXN1", "chr6",  16_327_867,  16_327_954,  "CAG",   pad=400, note="Spinocerebellar ataxia 1"),
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


def detect_platform(header: str) -> str:
    """Return 'ont', 'hifi', or 'unknown' from @RG PL tags."""
    for line in header.splitlines():
        if line.startswith("@RG"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            pl = fields.get("PL", "").upper()
            pm = fields.get("PM", "").upper()
            ds = fields.get("DS", "").upper()
            combined = " ".join([pl, pm, ds])
            if "ONT" in combined or "NANOPORE" in combined:
                return "ont"
            if "HIFI" in combined or "PACBIO" in combined or "SEQUEL" in combined:
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

# ── Read extraction ───────────────────────────────────────────────────────────

def extract_reads(bam: Path, chrom: str, start: int, end: int, pad: int) -> list[bytes]:
    """
    Return sequences of all non-supplementary reads overlapping the padded region.
    Reads shorter than (end - start) are discarded as non-spanning.
    """
    region_len = end - start
    padded_start = max(0, start - pad)
    padded_end   = end + pad
    region = f"{chrom}:{padded_start + 1}-{padded_end}"  # samtools 1-based

    # samtools view: exclude supplementary (-F 2048) and unmapped (-F 4)
    cmd_view = [
        "samtools", "view", "-F", "2052",
        str(bam), region,
    ]
    try:
        view = subprocess.run(cmd_view, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        return []

    if not view.stdout:
        return []

    # Parse the SAM output and filter by read length.
    reads = []
    for line in view.stdout.splitlines():
        fields = line.split(b"\t")
        if len(fields) < 10:
            continue
        seq = fields[9]
        if seq == b"*":
            continue
        # Keep reads that are at least as long as the core region.
        if len(seq) >= region_len:
            reads.append(seq)

    return reads

# ── poa-consensus invocation ──────────────────────────────────────────────────

def run_poa(reads: list[bytes], flags: list[str], multi: bool) -> Optional[list[bytes]]:
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
        return None
    if r.returncode != 0:
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

def count_units(seq: bytes, unit: str) -> int:
    """
    Find the longest contiguous run of `unit` within `seq` and return its count.

    The consensus includes flanking context, so dividing total length by unit
    length gives the wrong answer.  Instead, search for the longest run via regex,
    trying all rotational phases of the unit to handle phase ambiguity in the
    consensus (e.g. CAG / AGC / GCA all represent the same trinucleotide repeat).

    Example: for unit="CAG" and consensus
        "...ACGTCAGCAGCAGCAGCAGACGT..."  →  5 (the run of 5 CAG units)
    not 8 (which is what dividing total length 25 by 3 would give).
    """
    if not unit or not seq:
        return 0
    unit_b = unit.upper().encode()
    n = len(unit_b)
    best = 0
    # Try every rotational phase: CAG, AGC, GCA
    for i in range(n):
        rot = unit_b[i:] + unit_b[:i]
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
    platform  = args.platform if args.platform != "auto" else detect_platform(header)

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
                 alleles=30, repeats=20, time=8)
    header_row = (
        f"  {'Locus':<{col_w['locus']}}  {'Chrom':<{col_w['chrom']}}"
        f"  {'Region':<{col_w['region']}}  {'Unit':<{col_w['unit']}}"
        f"  {'Depth':>{col_w['depth']}}  {'Allele lengths':<{col_w['alleles']}}"
        f"  {'Repeat counts':<{col_w['repeats']}}  {'Time':>{col_w['time']}}"
    )
    sep = "  " + "─" * (len(header_row) - 2)

    tsv_rows: list[str] = [
        "\t".join(["locus", "category", "chrom", "start", "end", "unit", "depth",
                   "allele_lengths_bp", "repeat_counts", "time_s", "note"])
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
        reads = extract_reads(args.bam, chrom, locus.start, locus.end, pad)
        depth = len(reads)

        if depth < locus.min_depth:
            elapsed = time.perf_counter() - t0
            allele_str = f"SKIP (depth {depth} < {locus.min_depth})"
            repeat_str = ""
            _print_row(col_w, locus.name, chrom, locus.start, locus.end,
                       locus.unit, depth, allele_str, repeat_str, elapsed)
            tsv_rows.append("\t".join([
                locus.name, locus.category, chrom,
                str(locus.start), str(locus.end), locus.unit,
                str(depth), "skip", "", f"{elapsed:.2f}", locus.note,
            ]))
            continue

        # Controls run in single-allele mode: they're not repeat loci, and
        # any multi-allele result from a control would flag a pipeline issue.
        is_control = locus.category == "control"
        multi = not args.no_multi and not is_control
        alleles = run_poa(reads, flags, multi=multi)
        elapsed = time.perf_counter() - t0

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

        _print_row(col_w, locus.name, chrom, locus.start, locus.end,
                   locus.unit, depth, allele_str, repeat_str, elapsed)

        allele_bp = "; ".join(str(len(a)) for a in (alleles or []))
        repeat_ct = "; ".join(
            str(count_units(a, locus.unit)) for a in (alleles or [])
        ) if locus.unit else ""
        tsv_rows.append("\t".join([
            locus.name, locus.category, chrom,
            str(locus.start), str(locus.end), locus.unit,
            str(depth), allele_bp, repeat_ct,
            f"{elapsed:.2f}", locus.note,
        ]))

    print()

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text("\n".join(tsv_rows) + "\n")
        print(f"TSV written to {args.out}")


def _print_row(col_w: dict, name: str, chrom: str,
               start: int, end: int, unit: str,
               depth: int, allele_str: str, repeat_str: str,
               elapsed: float) -> None:
    region = f"{start:,}-{end:,}"
    print(
        f"  {name:<{col_w['locus']}}  {chrom:<{col_w['chrom']}}"
        f"  {region:<{col_w['region']}}  {unit:<{col_w['unit']}}"
        f"  {depth:>{col_w['depth']}}  {allele_str:<{col_w['alleles']}}"
        f"  {repeat_str:<{col_w['repeats']}}  {elapsed:>{col_w['time']}.1f}s"
    )


if __name__ == "__main__":
    main()
