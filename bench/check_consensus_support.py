#!/usr/bin/env python3
"""Read-support QC for a POA consensus.

Aligns the input reads back to the consensus with minimap2 and looks for
positions the reads do NOT support:

  * FABRICATED bases  -- a consensus position where a majority of covering
    reads carry a deletion (the consensus has a base the reads lack).  This is
    exactly the failure mode that produced private data 2 bp flank
    fabrication: greedy path selection threaded through phantom repeat-unit
    nodes whose graph coverage was inflated by repeat-node folding, so the
    graph's own statistics could not see the error, but a read pileup can.

  * MISSING bases     -- a consensus position where a majority of covering
    reads carry an insertion (the reads agree on a base the consensus dropped).

Both are computed from a straight `samtools mpileup` of reads-vs-consensus, so
they are independent of the POA graph's internal accounting -- which is the
whole point: this is an external check that catches errors the builder cannot
see in its own coverage numbers.

Usage:
    check_consensus_support.py CONSENSUS.fa READS.{fa,fq} [--preset map-ont]
        [--min-depth 10] [--frac 0.5] [--quiet]

Exit status is the number of flagged positions (0 == clean), capped at 250 so
it stays a usable shell exit code.  Intended both as a standalone QC tool and
as the objective acceptance test for consensus-accuracy work.

Requires `minimap2` and `samtools` on PATH.
"""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# mpileup encodes a read insertion as e.g. "+2AC" at the column *before* the
# inserted bases; a deletion placeholder is "*".
_INS_RE = re.compile(r"[+](\d+)")


def _run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def align(consensus: Path, reads: Path, preset: str, workdir: Path) -> Path:
    """minimap2 reads -> consensus, sorted+indexed BAM. Returns the BAM path."""
    bam = workdir / "aln.bam"
    mm2 = subprocess.Popen(
        ["minimap2", "-a", f"-x{preset}", "--secondary=no", str(consensus), str(reads)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    _run(["samtools", "sort", "-o", str(bam), "-"], stdin=mm2.stdout, stderr=subprocess.DEVNULL)
    mm2.wait()
    if mm2.returncode not in (0, None):
        raise RuntimeError("minimap2 failed")
    _run(["samtools", "index", str(bam)], stderr=subprocess.DEVNULL)
    return bam


def scan(consensus: Path, bam: Path, min_depth: int, frac: float):
    """Return (fabricated, missing) lists of (pos, base, depth, count, fraction)."""
    proc = subprocess.run(
        ["samtools", "mpileup", "-f", str(consensus), str(bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    fabricated = []
    missing = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        pos = int(f[1])
        base = f[2]
        depth = int(f[3])
        col = f[4]
        if depth < min_depth:
            continue
        dele = col.count("*")
        # count reads carrying an insertion immediately after this column
        ins = sum(int(m.group(1)) > 0 for m in _INS_RE.finditer(col))
        if dele / depth > frac:
            fabricated.append((pos, base, depth, dele, dele / depth))
        if ins / depth > frac:
            missing.append((pos, base, depth, ins, ins / depth))
    return fabricated, missing


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("consensus", type=Path)
    ap.add_argument("reads", type=Path)
    ap.add_argument("--preset", default="map-ont", help="minimap2 -x preset (default map-ont; use map-hifi for HiFi)")
    ap.add_argument("--min-depth", type=int, default=10, help="ignore columns below this depth")
    ap.add_argument("--frac", type=float, default=0.5, help="flag when >this fraction of reads disagree")
    ap.add_argument("--quiet", action="store_true", help="only print the summary line")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as td:
        bam = align(args.consensus, args.reads, args.preset, Path(td))
        fabricated, missing = scan(args.consensus, bam, args.min_depth, args.frac)

    if not args.quiet:
        for pos, base, depth, n, fr in fabricated:
            print(f"FABRICATED pos={pos} base={base} depth={depth} reads_deleting={n} ({fr:.0%})")
        for pos, base, depth, n, fr in missing:
            print(f"MISSING    after pos={pos} base={base} depth={depth} reads_inserting={n} ({fr:.0%})")

    total = len(fabricated) + len(missing)
    status = "CLEAN" if total == 0 else f"{total} unsupported position(s)"
    print(f"{args.consensus.name}: {status} "
          f"({len(fabricated)} fabricated, {len(missing)} missing)")
    return min(total, 250)


if __name__ == "__main__":
    sys.exit(main())
