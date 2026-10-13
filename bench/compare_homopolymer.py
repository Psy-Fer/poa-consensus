#!/usr/bin/env python3
"""
Homopolymer-interrupt loci: poa-consensus modes vs abPOA / SPOA / POASTA.

Loci where a base sits between length-variable homopolymers: a single G between
two A runs, and AAAAT-type repeats whose A runs are genuinely 3/4/5 long (shared
by every read) plus per-read noise. Callers:

  - poa        poa-consensus defaults
  - poa:gap    --homopolymer-gap=-2,-1
  - poa:hpc2g  --homopolymer-min-run 2 --homopolymer-gap=-2,-1
  - abPOA, SPOA (semi-global), POASTA (global, harness majority-vote consensus)
  - non-POA (bench/other_callers.py): MAFFT / Kalign MSA + majority, medoid
    read, medoid + minimap2 + samtools-consensus pileup; with --pbsim/--vntr
    also hifiasm, and TRGT / LongTR on a catalog reference with a different
    repeat count than the truth

Two read sources:

  --pbsim   pbsim3 -> minimap2 -> bedpull, as validate.py (realistic error
            profiles, but run-length noise is whatever the error model gives).
            Several pbsim seeds per locus.
  --jitter  controlled model (same as poa2::tests_hp_methods): each A run gets
            per-read length jitter (P(0)=p0, else +-1 mostly, +-2), plus uniform
            substitutions at rate e. Reads are exactly flank+locus+flank.

Each cell is `exact/N (mean edit)`, scored by fitting edit distance of the truth
locus against the consensus, so flank/anchor padding is free.

Usage:
    python3 bench/compare_homopolymer.py --jitter [--trials 50]
    python3 bench/compare_homopolymer.py --pbsim [--seeds 5]
    python3 bench/compare_homopolymer.py --vntr [--seeds 5]   # MUC1-like control
"""

import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import compare_callers as cc
import other_callers as oc
import validate as v

POA_MODES = [
    ("poa", []),
    ("poa:gap", ["--homopolymer-gap=-2,-1"]),
    ("poa:hpc2g", ["--homopolymer-min-run", "2", "--homopolymer-gap=-2,-1"]),
]
CALLERS = ([m for m, _ in POA_MODES] + ["abPOA", "SPOA", "POASTA"]
           + ["MAFFT", "Kalign", "medoid", "pileup"])
# pbsim modes only: need a BAM (TRGT, LongTR) or real-length reads (hifiasm).
BAM_CALLERS = ["hifiasm", "TRGT", "LongTR"]


# ── Loci ─────────────────────────────────────────────────────────────────────

def runs_locus(runs: list[int], base: str = "A", sep: str = "T") -> str:
    return "".join(base * n + sep for n in runs)


def mixed_runs(rng: random.Random, n_units: int) -> list[int]:
    """~15% A3, ~15% A5, rest A4: genuine variation shared by every read."""
    out = []
    for _ in range(n_units):
        u = rng.random()
        out.append(3 if u < 0.15 else 5 if u < 0.3 else 4)
    return out


def loci() -> list[tuple[str, list[int], str, str, str]]:
    """(name, runs, base, sep, truth). Truth = concatenated base^n + sep."""
    rng = random.Random(20261013)
    out = []
    for name, runs, base, sep in [
        ("aaaat15_pure", [4] * 15, "A", "T"),
        ("aaaat15_mixed", mixed_runs(rng, 15), "A", "T"),
        ("aaaat40_mixed", mixed_runs(rng, 40), "A", "T"),
        ("aaaag11_rfc1ref", [4] * 11, "A", "G"),       # RFC1 reference allele
        ("ttttc20_mixed", mixed_runs(rng, 20), "T", "C"),
    ]:
        out.append((name, runs, base, sep, runs_locus(runs, base, sep)))
    # Single interrupt: A12 G A11 (last run has no trailing separator).
    out.append(("a12_g_a11", [12, 11], "A", "G", "A" * 12 + "G" + "A" * 11))
    return out


# ── Callers ──────────────────────────────────────────────────────────────────

def write_fasta(path: Path, seqs: list[bytes]):
    with open(path, "wb") as fh:
        for i, s in enumerate(seqs):
            fh.write(b">r%d\n%s\n" % (i, s))


def first_fasta_seq(out: bytes) -> bytes:
    seq, started = [], False
    for line in out.split(b"\n"):
        if line.startswith(b">"):
            if started:
                break
            started = True
        else:
            seq.append(line.strip())
    return b"".join(seq)


def run_all(fa: Path, seqs: list[bytes], work: Path, preset: str = "map-ont") -> dict:
    """Consensus from every caller (None if it failed)."""
    res = {}
    for name, flags in POA_MODES:
        cmd = [str(v.POA_BIN), "--band-width", "50", "--min-reads", "2", *flags, str(fa)]
        p = subprocess.run(cmd, capture_output=True)
        res[name] = first_fasta_seq(p.stdout) if p.returncode == 0 else None
    res["abPOA"] = cc.run_abpoa(seqs)
    res["SPOA"] = cc.run_spoa(seqs)
    p = subprocess.run(["poasta", "align", "-O", "fasta", str(fa)], capture_output=True)
    res["POASTA"] = cc._msa_majority_consensus(p.stdout) if p.returncode == 0 else None
    res["medoid"] = oc.medoid(seqs)
    res["pileup"] = oc.pileup(seqs, work / "pileup", preset)
    res["MAFFT"] = oc.mafft(seqs, work / "msa")
    res["Kalign"] = oc.kalign(seqs, work / "msa")
    return res


def score(cons, truth: bytes) -> int | None:
    if not cons:
        return None
    return cc.fitting_edit_distance(truth, cons.upper())


# ── Jitter model ─────────────────────────────────────────────────────────────

def jitter(rng: random.Random, p0: float) -> int:
    if rng.random() < p0:
        return 0
    mag = 1 if rng.random() < 0.8 else 2
    return -mag if rng.random() < 0.5 else mag


def flank(rng: random.Random) -> str:
    # No A/T in the flank so it cannot extend the runs.
    return "C" + "".join(rng.choice("CG") for _ in range(28)) + "C"


def jitter_trial(args) -> dict:
    name, runs, base, sep, truth, p0, e, depth, seed = args
    rng = random.Random(seed)
    lf, rf = flank(rng), flank(rng)
    single = name.startswith("a12_g")
    reads = []
    for _ in range(depth):
        ns = [max(1, n + jitter(rng, p0)) for n in runs]
        core = (base * ns[0] + sep + base * ns[1]) if single else runs_locus(ns, base, sep)
        r = list(lf + core + rf)
        for i, c in enumerate(r):
            if rng.random() < e:
                r[i] = rng.choice([b for b in "ACGT" if b != c])
        reads.append("".join(r).encode())
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "reads.fa"
        write_fasta(fa, reads)
        res = run_all(fa, reads, Path(td))
    t = truth.encode()
    return {c: score(res[c], t) for c in CALLERS}


def run_jitter(args):
    cells = [(p0, e, d) for e in (0.0, 0.03) for p0 in (0.7, 0.4) for d in (10, 20, 40)]
    sel = [l for l in loci() if not args.loci or l[0] in args.loci]
    print(f"jitter model, {args.trials} trials per cell; cells = exact/N (mean edit)\n")
    hdr = f"{'locus':<17}{'e':>5}{'p0':>5}{'d':>4}  " + "".join(f"{c:<16}" for c in CALLERS)
    print(hdr)
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for name, runs, base, sep, truth in sel:
            for p0, e, d in cells:
                jobs = [(name, runs, base, sep, truth, p0, e, d, 1000 * t + 7)
                        for t in range(args.trials)]
                rows = list(ex.map(jitter_trial, jobs))
                print(f"{name:<17}{e:>5.2f}{p0:>5}{d:>4}  " + fmt_row(rows), flush=True)


def fmt_row(rows: list[dict], callers: list[str] = CALLERS) -> str:
    out = []
    for c in callers:
        vals = [r.get(c) for r in rows if r.get(c) is not None]
        if not vals:
            out.append(f"{'N/A':<16}")
            continue
        ex = sum(1 for x in vals if x == 0)
        fail = len(rows) - len(vals)
        s = f"{ex}/{len(rows)} ({sum(vals) / len(vals):.2f})" + (f"!{fail}" if fail else "")
        out.append(f"{s:<16}")
    return "".join(out)


# ── pbsim pipeline ───────────────────────────────────────────────────────────

def _fastq_seqs(path: Path) -> list[bytes]:
    lines = path.read_bytes().split(b"\n")
    return [lines[i] for i in range(1, len(lines), 4) if lines[i]]


def motif_of(name: str, truth: str) -> str:
    if name.startswith("vntr"):
        return VNTR_UNIT
    if name.startswith("a12_g"):
        return "A"
    # aaaat15_mixed -> AAAAT, ttttc20_mixed -> TTTTC, aaaag11_rfc1ref -> AAAAG
    return "".join(ch for ch in name.split("_")[0] if ch.isalpha()).upper()


def use_random_flanks(truth: str, seed: int):
    """Swap validate's periodic flank pool (a 32 bp unit x 300) for random
    sequence in this worker. The periodic pool makes the whole reference one
    tandem repeat: hifiasm assembles nothing (and peaks at 16 GB), and
    reference-based callers lose their unique flank anchors. Flank ends are
    chosen so they cannot extend the locus's first/last run."""
    rng = random.Random(seed)

    def mk(n, avoid):
        f = [rng.choice("ACGT") for _ in range(n)]
        if f[-1] == avoid[0]:
            f[-1] = next(b for b in "ACGT" if b != avoid[0])
        return "".join(f)

    left = mk(len(v._LEFT_FLANK_FULL), truth[0])
    right = mk(len(v._RIGHT_FLANK_FULL), truth[-1])[::-1]
    v._LEFT_FLANK_FULL, v._RIGHT_FLANK_FULL = left, right
    v.LEFT_FLANK, v.RIGHT_FLANK = left[:v.FLANK_LEN], right[:v.FLANK_LEN]


def pbsim_trial(args) -> dict | None:
    name, truth, model_key, depth, seed, workdir = args
    use_random_flanks(truth, 7919 * seed + len(truth))
    # Reads must span long loci: length_mean ~ locus + 2.5 kb (< reference length).
    ref_len = 2 * v.FLANK_LEN + len(truth)
    length_mean = min(max(500, int(ref_len * 0.75), len(truth) + 2500), ref_len - 100)
    sc = v.Scenario(f"{name}_{model_key}_s{seed}",
                    [v.Allele(truth, 1, depth, desc=name)], model_key, seed=seed,
                    length_mean=length_mean)
    work = Path(workdir) / sc.name
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    model = v.ERRMODELS[model_key]
    try:
        ref = v.build_reference(sc, work)
        reads = v.simulate_reads(sc, ref, work, model)
        bam = v.align_and_index(reads, ref, work, model["mm2_preset"])
        fa = v.extract_reads(sc, bam, work)
        seqs = [s for _, s in v.parse_fasta(fa)]
        if len(seqs) < 3:
            return None
        preset = model["mm2_preset"]
        # map-ont for the pileup regardless of model: the bench "hifi" model is
        # 89% accurate and map-hifi leaves most reads unmapped to a read draft.
        res = run_all(fa, seqs, work, "map-ont")
        motif = motif_of(name, truth)
        # No hifiasm here: at the pbsim models' 85-92% accuracy it finds no
        # strong overlaps and assembles nothing (see --replicant for HiFi).
        res.update(oc.bam_callers(reads, v.LEFT_FLANK, v.RIGHT_FLANK,
                                  oc.catalog_allele(motif, len(truth)), motif,
                                  work / "bam_callers", preset))
    finally:
        if not os.environ.get("KEEP_WORK"):
            shutil.rmtree(work, ignore_errors=True)
    t = truth.encode()
    return {c: score(res.get(c), t) for c in CALLERS + BAM_CALLERS}


def run_pbsim(args):
    models = ["ont_r10", "ont_r9", "hifi"]
    sel = [l for l in loci() if not args.loci or l[0] in args.loci]
    print(f"pbsim3 pipeline, depth {args.depth}, {args.seeds} seeds per cell; "
          f"cells = exact/N (mean edit)\n")
    print(f"{'locus':<17}{'model':<9}  " + "".join(f"{c:<16}" for c in CALLERS + BAM_CALLERS))
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for name, _runs, _b, _s, truth in sel:
            for m in models:
                jobs = [(name, truth, m, args.depth, 42 + 101 * s, str(workdir))
                        for s in range(args.seeds)]
                rows = [r for r in ex.map(pbsim_trial, jobs) if r is not None]
                if not rows:
                    print(f"{name:<17}{m:<9}  no usable read sets")
                    continue
                print(f"{name:<17}{m:<9}  " + fmt_row(rows, CALLERS + BAM_CALLERS), flush=True)


# ── MUC1-like VNTR (non-homopolymer control) ─────────────────────────────────

# Made-up GC-rich 60 bp unit (bench/make_synthetic_vntr.py), not a real MUC1 motif.
VNTR_UNIT = "GGCTCAGGTGGCAGCGGAGCCTGGACCTGCAGTGGCAGTCCAGGCTCCACGTGGCAGTGG"


def vntr_locus(n_units: int, seed: int) -> str:
    """n_units copies; ~1 in 4 carries one or two substitutions (variant repeats)."""
    rng = random.Random(seed)
    out = []
    for _ in range(n_units):
        u = list(VNTR_UNIT)
        if rng.random() < 0.25:
            for _ in range(rng.choice((1, 2))):
                i = rng.randrange(len(u))
                u[i] = rng.choice([b for b in "ACGT" if b != u[i]])
        out.append("".join(u))
    return "".join(out)


def run_vntr(args):
    models = ["ont_r10", "ont_r9", "hifi"]
    print(f"MUC1-like VNTR (60 bp GC-rich unit), pbsim3, depth {args.depth}, "
          f"{args.seeds} seeds per cell; cells = exact/N (mean edit)\n")
    print(f"{'locus':<17}{'model':<9}  " + "".join(f"{c:<16}" for c in CALLERS + BAM_CALLERS))
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        for n in (10, 25, 50):
            truth = vntr_locus(n, 900 + n)
            name = f"vntr60x{n}"
            for m in models:
                jobs = [(name, truth, m, args.depth, 42 + 101 * s, str(workdir))
                        for s in range(args.seeds)]
                rows = [r for r in ex.map(pbsim_trial, jobs) if r is not None]
                if not rows:
                    print(f"{name:<17}{m:<9}  no usable read sets")
                    continue
                print(f"{name:<17}{m:<9}  " + fmt_row(rows, CALLERS + BAM_CALLERS), flush=True)


# ── str-replicant HiFi ───────────────────────────────────────────────────────

STR_REPLICANT = os.environ.get("STR_REPLICANT") or shutil.which("str-replicant") or str(
    Path.home() / "repos/str-replicant/target/release/str-replicant")
# Flanks longer than the reads so reads start at staggered positions: identical
# read spans are all "contained" for hifiasm and PCR duplicates for LongTR.
SR_FLANK = 3000


def replicant_loci() -> list[tuple[str, str, list[str]]]:
    """(id, truth, str-replicant states as TOML inline tables)."""
    out = []
    for name, runs, base, sep, truth in loci():
        if name.startswith("a12_g"):
            st = [f'{{ homopolymer = "{base}", units = {runs[0]} }}',
                  f'{{ interruption = "{sep}" }}',
                  f'{{ homopolymer = "{base}", units = {runs[1]} }}']
        else:
            st = []
            for n in runs:
                st += [f'{{ homopolymer = "{base}", units = {n} }}', f'{{ interruption = "{sep}" }}']
        out.append((name, truth, st))
    for n in (10, 25, 50):
        truth = vntr_locus(n, 900 + n)
        units = [truth[i:i + len(VNTR_UNIT)] for i in range(0, len(truth), len(VNTR_UNIT))]
        out.append((f"vntr60x{n}", truth, [f'{{ vntr = "{u}", units = 1 }}' for u in units]))
    return out


def replicant_config(depth: int, seed: int) -> str:
    parts = [f"seed = {seed}\n"]
    sep = ",\n  "
    for lid, truth, states in replicant_loci():
        parts.append(f"""
[[locus]]
id = "{lid}"
depth = {depth}
flanks = {{ synthetic = {{ len = {SR_FLANK}, gc = 0.45 }} }}

[[locus.haplotype]]
states = [
  {sep.join(states)},
]

[[locus.haplotype]]
same_as = 0

[locus.reads]
length = {len(truth) + 2500}
placement = "spanning"

[locus.error]
preset = "pacbio-revio-hg002"
""")
    return "".join(parts)


def replicant_trial(args) -> dict | None:
    lid, truth, ref_seq, reads, workdir = args
    work = Path(workdir) / lid
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        ref = work / "reference.fa"
        ref.write_text(f">{lid}\n{ref_seq}\n")
        fq = work / "reads.fq"
        with open(fq, "wb") as fh:
            for i, (sq, ql) in enumerate(reads):
                fh.write(b"@r%d\n%s\n+\n%s\n" % (i, sq, ql))
        # Extract the locus window (+-ANCHOR_PAD) as validate does, from an
        # alignment to str-replicant's own reference (used only for extraction).
        mm = subprocess.run(["minimap2", "-ax", "map-hifi", str(ref), str(fq)], capture_output=True)
        bam = work / "aligned.bam"
        subprocess.run(["samtools", "sort", "-o", str(bam), "-"], input=mm.stdout, capture_output=True)
        subprocess.run(["samtools", "index", str(bam)], capture_output=True)
        bed = work / "targets.bed"
        lo, hi = SR_FLANK, SR_FLANK + len(truth)
        bed.write_text(f"{lid}\t{max(0, lo - v.ANCHOR_PAD)}\t{hi + v.ANCHOR_PAD}\n")
        fa = work / "extracted.fa"
        subprocess.run(["bedpull", "-b", str(bam), "-r", str(bed), "-o", str(fa)], capture_output=True)
        seqs = [s for _, s in v.parse_fasta(fa)] if fa.exists() else []
        if len(seqs) < 3:
            return None
        res = run_all(fa, seqs, work, "map-hifi")
        res["hifiasm"] = oc.hifiasm([sq for sq, _ in reads], work / "hifiasm", ont=False)
        motif = motif_of(lid, truth)
        res.update(oc.bam_callers(fq, ref_seq[:lo], ref_seq[hi:],
                                  oc.catalog_allele(motif, len(truth)), motif,
                                  work / "bam_callers", "map-hifi"))
    finally:
        if not os.environ.get("KEEP_WORK"):
            shutil.rmtree(work, ignore_errors=True)
    t = truth.encode()
    return {c: score(res.get(c), t) for c in CALLERS + BAM_CALLERS}


def run_replicant(args):
    sr_loci = [l for l in replicant_loci() if not args.loci or l[0] in args.loci]
    keep = {l[0] for l in sr_loci}
    workdir = Path(args.workdir) / "replicant"
    workdir.mkdir(parents=True, exist_ok=True)
    print(f"str-replicant HiFi (pacbio-revio-hg002 preset), depth {args.depth}, "
          f"{args.seeds} seeds per locus; cells = exact/N (mean edit)\n")
    print(f"{'locus':<17}{'':<9}  " + "".join(f"{c:<16}" for c in CALLERS + BAM_CALLERS))
    rows: dict[str, list] = {l[0]: [] for l in sr_loci}
    for s in range(args.seeds):
        sd = workdir / f"seed{s}"
        cfg = workdir / f"seed{s}.toml"
        cfg.write_text(replicant_config(args.depth, 42 + 101 * s))
        p = subprocess.run([STR_REPLICANT, "run", str(cfg), "-o", str(sd), "--format", "fastq",
                            "--force", "-t", "4"], capture_output=True)
        if p.returncode != 0:
            sys.exit(f"str-replicant failed: {p.stderr.decode()[-500:]}")
        refs = {h.split()[0]: sq.decode() for h, sq in v.parse_fasta(sd / "reference.fa")}
        lines = (sd / "reads.fq").read_bytes().split(b"\n")
        reads: dict[str, list] = {}
        for i in range(0, len(lines) - 3, 4):
            lid = lines[i][1:].split(b":")[0].decode()
            reads.setdefault(lid, []).append((lines[i + 1], lines[i + 3]))
        jobs = [(lid, truth, refs[lid], reads.get(lid, []), str(sd))
                for lid, truth, _ in sr_loci if lid in keep]
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            for (lid, *_), r in zip(jobs, ex.map(replicant_trial, jobs)):
                if r is not None:
                    rows[lid].append(r)
    for lid, rs in rows.items():
        if rs:
            print(f"{lid:<17}{'hifi':<9}  " + fmt_row(rs, CALLERS + BAM_CALLERS), flush=True)
        else:
            print(f"{lid:<17}{'hifi':<9}  no usable read sets")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jitter", action="store_true")
    ap.add_argument("--pbsim", action="store_true")
    ap.add_argument("--vntr", action="store_true", help="MUC1-like 60 bp VNTR control")
    ap.add_argument("--replicant", action="store_true",
                    help="HiFi reads from str-replicant's pacbio-revio-hg002 preset")
    ap.add_argument("--trials", type=int, default=50)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--depth", type=int, default=20)
    ap.add_argument("--loci", nargs="*")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--workdir", default="bench/work_homopolymer")
    args = ap.parse_args()
    if not v.POA_BIN.exists():
        sys.exit(f"build the CLI first: cargo build --release --features cli ({v.POA_BIN})")
    print("loci:")
    for name, _r, _b, _s, truth in loci():
        print(f"  {name:<17} {truth}")
    print()
    if args.jitter:
        run_jitter(args)
    if args.pbsim:
        run_pbsim(args)
    if args.vntr:
        run_vntr(args)
    if args.replicant:
        run_replicant(args)
    if not (args.jitter or args.pbsim or args.vntr or args.replicant):
        ap.error("pass --jitter, --pbsim, --vntr and/or --replicant")


if __name__ == "__main__":
    main()
