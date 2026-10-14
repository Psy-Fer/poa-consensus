#!/usr/bin/env python3
"""
Real-read comparison on HG002 against per-haplotype sequence truth.

Truth: bladerunner's `giab_truth.py` table (GIAB v5.0q / T2T-HG002 Q100 v1.1),
per haplotype over each STRchive locus window (BED +-100 bp), with the true repeat
interval from TRF (`rep_ext_hstart/hend`). See ~/Data/refs/hg002/truth/README.md.

Reads: a whatshap-haplotagged BAM. Per locus, the reads over the truth window are
extracted with bedpull (which labels each read h1/h2 from its HP tag) and split by
haplotype; every single-allele caller runs on each haplotype's reads. If either
haplotype has fewer than 3 tagged reads, all reads go into one call (scored against
both truth haplotypes). Calls are paired to truth haplotypes by best total match
(HP numbering is arbitrary), as the truth README does for unphased tools.

TRGT and LongTR run once on the whole BAM with all loci, as they would be used,
and their alleles are spliced into the reference window.

Score per haplotype: fitting edit distance of the truth repeat (TRF extended
interval +-20 bp) against the call. Cells are exact/N and mean edit.

Usage:
    python3 bench/hg002_ont_truth.py --bam HG002_ont...bam [--loci NAME ...] [--jobs 4]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import compare_callers as cc
import compare_homopolymer as ch
import other_callers as oc
import validate as v

TRUTH = Path.home() / "Data/refs/hg002/truth/truth_repeats.grch38.tsv"
REF = Path.home() / "Data/refs/hg38.analysisSet.fa.gz"
PAD = 20
CALLERS = ch.CALLERS
GENOTYPERS = ["TRGT", "LongTR"]


def load_truth(path: Path) -> dict:
    loci: dict[str, dict] = {}
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            if r["route"] != "vcf":
                continue
            L = loci.setdefault(r["name"], {
                "chrom": r["chrom"], "ws": int(r["win_start"]), "we": int(r["win_end"]),
                "in_benchmark": r["in_benchmark"] == "1", "haps": {}})
            seq = r["sequence"]
            rs, re_ = int(r["rep_ext_hstart"]), int(r["rep_ext_hend"])
            L["haps"][r["hap"]] = {
                "seq": seq,
                "target": seq[max(0, rs - PAD):min(len(seq), re_ + PAD)],
                "rep_len": re_ - rs,
                "motif": r["rep_motifs"].split(",")[0],
                "locus": (int(r["locus_hstart"]), int(r["locus_hend"])),
            }
            L["motif"] = L["haps"][r["hap"]]["motif"]
    return loci


def fetch_ref(chrom: str, start: int, end: int) -> str:
    p = subprocess.run(["samtools", "faidx", str(REF), f"{chrom}:{start + 1}-{end}"],
                       capture_output=True, text=True)
    return "".join(p.stdout.splitlines()[1:]).upper()


def call_all(seqs: list[bytes], work: Path) -> dict:
    fa = work / "reads.fa"
    ch.write_fasta(fa, seqs)
    return ch.run_all(fa, seqs, work, "map-ont")


def score_pair(calls: list, truths: list[bytes]) -> list:
    """calls: 1 or 2 sequences (None = no call); truths: 1 or 2 targets. Returns one
    edit distance per truth, pairing calls to truths by best total."""
    def d(c, t):
        return None if not c else cc.fitting_edit_distance(t, c.upper())
    if len(truths) == 1:
        ds = [d(c, truths[0]) for c in calls]
        ds = [x for x in ds if x is not None]
        return [min(ds) if ds else None]
    if len(calls) == 1:
        return [d(calls[0], t) for t in truths]
    a = [d(calls[0], truths[0]), d(calls[1], truths[1])]
    b = [d(calls[1], truths[0]), d(calls[0], truths[1])]
    big = 10 ** 9
    tot = lambda xs: sum(big if x is None else x for x in xs)
    return a if tot(a) <= tot(b) else b


def locus_job(args) -> tuple[str, dict, dict]:
    name, L, bam, workdir = args
    work = Path(workdir) / name
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    info = {"reads": 0, "split": False}
    out: dict[str, list] = {}
    try:
        region = f"{L['chrom']}:{L['ws'] + 1}-{L['we']}"
        sub = work / "sub.bam"
        subprocess.run(["samtools", "view", "-b", "-o", str(sub), bam, region], capture_output=True)
        subprocess.run(["samtools", "index", str(sub)], capture_output=True)
        bed = work / "w.bed"
        bed.write_text(f"{L['chrom']}\t{L['ws']}\t{L['we']}\n")
        fa = work / "ext.fa"
        subprocess.run(["bedpull", "-b", str(sub), "-r", str(bed), "-o", str(fa)], capture_output=True)
        groups: dict[str, list] = {}
        if fa.exists():
            for h, s in v.parse_fasta(fa):
                groups.setdefault(h.split("|")[-1], []).append(s)
        g1, g2 = groups.get("h1", []), groups.get("h2", [])
        allr = [s for g in groups.values() for s in g]
        info["reads"] = len(allr)
        haps = sorted(L["haps"])
        truths = [L["haps"][h]["target"].encode() for h in haps]
        if len(truths) == 2 and len(g1) >= 3 and len(g2) >= 3:
            info["split"] = True
            sets = [g1, g2]
        elif len(allr) >= 3:
            sets = [allr]
        else:
            return name, info, {}
        res = []
        for i, s in enumerate(sets):
            w = work / f"set{i}"
            w.mkdir()
            res.append(call_all(s, w))
        for c in CALLERS:
            out[c] = score_pair([r.get(c) for r in res], truths)
    finally:
        if not os.environ.get("KEEP_WORK"):
            shutil.rmtree(work, ignore_errors=True)
    return name, info, out


# ── TRGT / LongTR on the whole BAM ───────────────────────────────────────────

def genotypers(loci: dict, bam: str, work: Path) -> dict:
    """{locus: {tool: [edit per truth hap]}}"""
    work.mkdir(parents=True, exist_ok=True)
    out: dict[str, dict] = {n: {} for n in loci}
    tbed, lbed = work / "trgt.bed", work / "longtr.bed"
    with open(tbed, "w") as ft, open(lbed, "w") as fl:
        for n, L in sorted(loci.items(), key=lambda x: (x[1]["chrom"], x[1]["ws"])):
            s, e = L["bed"]
            m = re.sub("[^ACGT]", "", L["motif"]) or "A"
            ft.write(f"{L['chrom']}\t{s}\t{e}\tID={n};MOTIFS={m};STRUC=({m})n\n")
            fl.write(f"{L['chrom']}\t{s + 1}\t{e}\t{m}\t{n}\n")
    vcfs = {}
    if oc.TRGT_BIN:
        p = subprocess.run([oc.TRGT_BIN, "genotype", "-g", str(REF), "-r", bam, "-b", str(tbed),
                            "-o", str(work / "trgt"), "--sample-name", "HG002", "-t", "4",
                            "--karyotype", "XY"], capture_output=True)
        (work / "trgt.log").write_bytes(p.stdout + p.stderr)
        vcfs["TRGT"] = work / "trgt.vcf.gz"
    if oc.LONGTR_BIN:
        p = subprocess.run([oc.LONGTR_BIN, "--bams", bam, "--fasta", str(REF), "--regions", str(lbed),
                            "--tr-vcf", str(work / "longtr.vcf.gz"), "--lib-from-samp",
                            "--min-reads", "3", "--haploid-chrs", "chrX,chrY"],
                           capture_output=True, stdin=subprocess.DEVNULL)
        (work / "longtr.log").write_bytes(p.stdout + p.stderr)
        vcfs["LongTR"] = work / "longtr.vcf.gz"
    for tool, vcf in vcfs.items():
        if not vcf.exists():
            continue
        txt = subprocess.run(["bash", "-c", f"gzip -dc {vcf}"], capture_output=True, text=True).stdout
        for line in txt.splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            chrom, pos, ref, alts = f[0], int(f[1]), f[3], [a for a in f[4].split(",") if a != "."]
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = re.split(r"[/|]", fmt.get("GT", "."))
            name = next((n for n, L in loci.items()
                         if L["chrom"] == chrom and L["ws"] <= pos - 1 < L["we"]), None)
            if name is None or not all(g.isdigit() for g in gt):
                continue
            L = loci[name]
            win = fetch_ref(chrom, L["ws"], L["we"])
            off = pos - 1 - L["ws"]
            if off + len(ref) > len(win):
                continue
            calls = []
            for g in gt:
                k = int(g)
                allele = ref if k == 0 else alts[k - 1]
                calls.append((win[:off] + allele + win[off + len(ref):]).encode())
            truths = [L["haps"][h]["target"].encode() for h in sorted(L["haps"])]
            out[name][tool] = score_pair(calls if len(set(calls)) > 1 or len(truths) == 1 else calls[:1],
                                         truths)
    return out


# ── Report ───────────────────────────────────────────────────────────────────

def summarize(label: str, rows: list, cols: list[str]):
    print(f"\n{label}: {len(rows)} haplotypes; exact / mean edit (no call counted separately)")
    for c in cols:
        vals = [r[c] for r in rows if r.get(c) is not None]
        miss = len(rows) - len(vals)
        if not vals:
            print(f"  {c:<12} no calls")
            continue
        ex = sum(1 for x in vals if x == 0)
        print(f"  {c:<12} {ex:>3}/{len(rows)} exact  mean edit {sum(vals) / len(vals):6.2f}"
              + (f"  ({miss} no call)" if miss else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--truth", default=str(TRUTH))
    ap.add_argument("--loci", nargs="*")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--no-genotypers", action="store_true")
    ap.add_argument("--workdir", default="bench/work_hg002_ont")
    args = ap.parse_args()
    loci = load_truth(Path(args.truth))
    # Locus (repeat) coordinates on GRCh38 for the genotypers: the catalog BED.
    for line in (Path(args.truth).parent / "loci.grch38.bed").read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 4 and f[3] in loci:
            loci[f[3]]["bed"] = (int(f[1]), int(f[2]))
    loci = {n: L for n, L in loci.items() if "bed" in L}
    if args.loci:
        loci = {n: L for n, L in loci.items() if n in args.loci}
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    gt = {} if args.no_genotypers else genotypers(loci, args.bam, workdir / "genotypers")
    cols = CALLERS + ([] if args.no_genotypers else GENOTYPERS)
    print(f"{'locus':<22}{'hap':<4}{'motif':<9}{'rep':>5}{'rd':>4}{'sp':>3}  "
          + "".join(f"{c[:9]:<10}" for c in cols))
    rows_at, rows_other = [], []
    jobs = [(n, L, args.bam, str(workdir)) for n, L in sorted(loci.items())]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for name, info, out in ex.map(locus_job, jobs):
            L = loci[name]
            out.update(gt.get(name, {}))
            for i, h in enumerate(sorted(L["haps"])):
                H = L["haps"][h]
                row = {c: (out[c][i] if c in out and i < len(out[c]) else None) for c in cols}
                cells = "".join(f"{'-' if row[c] is None else row[c]:<10}" for c in cols)
                print(f"{name[:21]:<22}{h:<4}{H['motif'][:8]:<9}{H['rep_len']:>5}{info['reads']:>4}"
                      f"{'y' if info['split'] else 'n':>3}  {cells}", flush=True)
                m = H["motif"]
                at_penta = len(m) == 5 and set(m) <= {"A", "T", "C", "G"} and sum(b in "AT" for b in m) >= 4 \
                    and max(m.count("A"), m.count("T")) >= 3
                (rows_at if at_penta else rows_other).append(row)
    summarize("A/T-rich pentanucleotide loci (AAAAT/AAAAG/TTTCA-type)", rows_at, cols)
    summarize("All other loci", rows_other, cols)


if __name__ == "__main__":
    main()
