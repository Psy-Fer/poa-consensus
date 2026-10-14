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
locus plus CTX (20) bp of flank on each side against the consensus, so further
flank/anchor padding is free but a repeat over-call is not.

Usage:
    python3 bench/compare_homopolymer.py --jitter [--trials 50]
    python3 bench/compare_homopolymer.py --pbsim [--seeds 5]
    python3 bench/compare_homopolymer.py --vntr [--seeds 5]   # MUC1-like control
    python3 bench/compare_homopolymer.py --replicant [--depth 20]
    python3 bench/compare_homopolymer.py --mosaic [--conditions sd2 sub0.40]
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


# Flank context scored with the truth. Without it, fitting alignment absorbs an
# over-call of a pure repeat (AAAAT x16 contains AAAAT x15 exactly) or of a run
# at the locus edge, and scores it as exact.
CTX = 20


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
    t = (lf[-CTX:] + truth + rf[:CTX]).encode()
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
    t = (v.LEFT_FLANK[-CTX:] + truth + v.RIGHT_FLANK[:CTX]).encode()
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
    t = ref_seq[lo - CTX:hi + CTX].encode()
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


# ── Somatic mosaicism (str-replicant, HiFi) ──────────────────────────────────
#
# High-accuracy reads where the length variation is biological: each molecule
# draws a skew-normal delta in units (str-replicant [locus.mosaic]), or a
# second allele (an expanded subclone) is present at a fraction. The germline
# allele is the target for the single-consensus callers; per-read molecule
# lengths come from str-replicant's `jp` tag.

def mosaic_loci() -> list[str]:
    return ["aaaat15", "aaaat40", "aaaag11", "aaaat15_mixed", "polyA25", "cag40", "vntr60x12"]


def mosaic_states(lid: str, d: int) -> list[str]:
    """str-replicant states for locus `lid` with `d` units added (or removed).
    For the mixed AAAAT locus a unit is one A4+T pair at the end (the existing
    A3/A4/A5 runs are kept), as str-replicant's own mosaic model does."""
    simple = {"aaaat15": ("motif", "AAAAT", 15), "aaaat40": ("motif", "AAAAT", 40),
              "aaaag11": ("motif", "AAAAG", 11), "cag40": ("motif", "CAG", 40),
              "polyA25": ("homopolymer", "A", 25), "vntr60x12": ("vntr", VNTR_UNIT, 12)}
    if lid in simple:
        kind, unit, n = simple[lid]
        return [f'{{ {kind} = "{unit}", units = {n + d} }}']
    runs = next(r for name, r, *_ in loci() if name == lid)
    runs = runs + [4] * d if d >= 0 else runs[:d]
    st = []
    for n in runs:
        st += [f'{{ homopolymer = "A", units = {n} }}', '{ interruption = "T" }']
    return st


def mosaic_fractions(spread: float, skew: float) -> list[tuple[int, float]]:
    """Discretised skew-normal (location 0) over whole-unit deltas: the
    fraction of molecules carrying each delta. Germline (0) first."""
    import math
    if spread == 0:
        return [(0, 1.0)]
    phi = lambda z: math.exp(-z * z / 2) / math.sqrt(2 * math.pi)
    Phi = lambda z: 0.5 * (1 + math.erf(z / math.sqrt(2)))
    k = int(math.ceil(4 * spread))
    w = {d: 2 / spread * phi(d / spread) * Phi(skew * d / spread) for d in range(-k, k + 1)}
    tot = sum(w.values())
    w = {d: x / tot for d, x in w.items() if x / tot >= 0.005}
    tot = sum(w.values())
    return sorted(((d, x / tot) for d, x in w.items()), key=lambda t: (t[0] != 0, t[0]))


# (label, mosaic (spread, skew) in units, subclone (fraction, +units) or None)
MOSAIC_CONDITIONS = [
    ("sd1", (1, 0), None),
    ("sd2", (2, 0), None),
    ("sd4", (4, 0), None),
    ("sd2+skew", (2, 4), None),
    ("sd4+skew", (4, 4), None),
    ("sub0.25", (0, 0), (0.25, 6)),
    ("sub0.40", (0, 0), (0.40, 6)),
    ("sub0.40+sd2", (2, 0), (0.40, 6)),
]


def mosaic_config(cond, depth: int, seed: int, only: set) -> str:
    """One locus per mosaic locus; molecules are a mixture of haplotypes at
    germline +- d units (str-replicant's own [locus.mosaic] only targets
    period-2..12 motif states, which excludes poly-A, the VNTR and the mixed
    A-run locus). With a subclone, the same spread applies around it."""
    _label, (spread, skew), sub = cond
    parts = [f"seed = {seed}\n"]
    sep = ",\n  "
    fr = mosaic_fractions(spread, skew)
    for lid in mosaic_loci():
        if lid not in only:
            continue
        haps = [(d, f * (1 - sub[0] if sub else 1)) for d, f in fr]
        if sub:
            haps += [(sub[1] + d, f * sub[0]) for d, f in fr]
        parts.append(f"""
[[locus]]
id = "{lid}"
depth = {depth}
flanks = {{ synthetic = {{ len = {SR_FLANK}, gc = 0.45 }} }}
""")
        for d, f in haps:
            parts.append(f"""
[[locus.haplotype]]
fraction = {f:.6f}
states = [
  {sep.join(mosaic_states(lid, d))},
]
""")
        parts.append(f"""
[locus.reads]
length = {len(VNTR_UNIT) * 20 + 2500}
placement = "spanning"

[locus.error]
preset = "pacbio-revio-hg002"
""")
    return "".join(parts)


def _extract(lid: str, ref_seq: str, rep_len: int, reads, work: Path):
    """Write reads, align to str-replicant's reference, bedpull the locus
    window (+-ANCHOR_PAD). Returns (fastq path, extracted fasta, seqs)."""
    ref = work / "reference.fa"
    ref.write_text(f">{lid}\n{ref_seq}\n")
    fq = work / "reads.fq"
    with open(fq, "wb") as fh:
        for i, (sq, ql) in enumerate(reads):
            fh.write(b"@r%d\n%s\n+\n%s\n" % (i, sq, ql))
    mm = subprocess.run(["minimap2", "-ax", "map-hifi", str(ref), str(fq)], capture_output=True)
    bam = work / "aligned.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), "-"], input=mm.stdout, capture_output=True)
    subprocess.run(["samtools", "index", str(bam)], capture_output=True)
    bed = work / "targets.bed"
    lo, hi = SR_FLANK, SR_FLANK + rep_len
    bed.write_text(f"{lid}\t{max(0, lo - v.ANCHOR_PAD)}\t{hi + v.ANCHOR_PAD}\n")
    fa = work / "extracted.fa"
    subprocess.run(["bedpull", "-b", str(bam), "-r", str(bed), "-o", str(fa)], capture_output=True)
    seqs = [s for _, s in v.parse_fasta(fa)] if fa.exists() else []
    return fq, fa, seqs


def _fasta_all(out: bytes) -> list[bytes]:
    recs, cur = [], None
    for line in out.split(b"\n"):
        if line.startswith(b">"):
            if cur is not None:
                recs.append(b"".join(cur))
            cur = []
        elif cur is not None:
            cur.append(line.strip())
    if cur is not None:
        recs.append(b"".join(cur))
    return [r for r in recs if r]


MOSAIC_CALLERS = CALLERS + ["TRGT", "LongTR"]


def mosaic_trial(args) -> dict | None:
    lid, germ, minor, modal, ref_seq, reads, workdir = args
    work = Path(workdir) / lid
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        fq, fa, seqs = _extract(lid, ref_seq, len(germ), reads, work)
        if len(seqs) < 3:
            return None
        res = run_all(fa, seqs, work, "map-hifi")
        motif = VNTR_UNIT if lid.startswith("vntr") else motif_of(lid, germ)
        res.update(oc.bam_callers(fq, ref_seq[:SR_FLANK], ref_seq[SR_FLANK + len(germ):],
                                  oc.catalog_allele(motif, len(germ)), motif,
                                  work / "bam_callers", "map-hifi"))
        # Multi-allele: number of alleles, and whether germline (and the
        # subclone, if any) are each matched exactly by some allele.
        p = subprocess.run([str(v.POA_BIN), "--band-width", "50", "--min-reads", "2",
                            "--multi", str(fa)], capture_output=True)
        alleles = _fasta_all(p.stdout) if p.returncode == 0 else []
    finally:
        if not os.environ.get("KEEP_WORK"):
            shutil.rmtree(work, ignore_errors=True)
    lf, rf = ref_seq[SR_FLANK - CTX:SR_FLANK], ref_seq[SR_FLANK + len(germ):][:CTX]
    g = (lf + germ + rf).encode()
    out = {c: score(res.get(c), g) for c in MOSAIC_CALLERS}
    m = (lf + modal + rf).encode()
    out["_modal"] = {c: score(res.get(c), m) == 0 for c in MOSAIC_CALLERS}
    # Signed length error vs the germline, on the truth-fitted span of the call.
    out["_dlen"] = {c: _repeat_len(res[c], lf, rf) - len(germ) for c in MOSAIC_CALLERS
                    if res.get(c) and out[c] is not None}
    out["_multi"] = (len(alleles),
                     any(score(a, g) == 0 for a in alleles),
                     minor is not None and any(score(a, (lf + minor + rf).encode()) == 0
                                               for a in alleles))
    return out


def mosaic_allele(lid: str, d: int) -> str:
    out = []
    for st in mosaic_states(lid, d):
        kv = dict(x.strip().split(" = ") for x in st.strip("{} ").split(", "))
        if "interruption" in kv:
            out.append(kv["interruption"].strip('"'))
        else:
            unit = next(kv[k] for k in ("motif", "homopolymer", "vntr") if k in kv).strip('"')
            out.append(unit * int(kv["units"]))
    return "".join(out)


def _repeat_len(cons: bytes, lf: str, rf: str) -> int:
    """Length of `cons` between the two flank anchors (the scoring context).
    Fitting the whole flanked truth does not give this: for long units it
    prefers ~CTX mismatches in the flank to a whole-unit indel."""
    import edlib
    c = cons.upper().decode()
    left = edlib.align(lf, c, mode="HW", task="locations")["locations"]
    right = edlib.align(rf, c, mode="HW", task="locations")["locations"]
    return right[-1][0] - (left[0][1] + 1)


def _molecule_deltas(header_lines: list[bytes], germ_len: int) -> list[int]:
    """Per-read molecule repeat length minus germline (bp), from the jp tag."""
    out = []
    for h in header_lines:
        for tag in h.split(b"\t"):
            if tag.startswith(b"jp:B:I,"):
                out.append(int(tag.split(b",")[2]) - germ_len)
    return out


def run_mosaic(args):
    import statistics
    names = [l for l in mosaic_loci() if not args.loci or l in args.loci]
    conds = [c for c in MOSAIC_CONDITIONS if not args.conditions or c[0] in args.conditions]
    workdir = Path(args.workdir) / "mosaic"
    workdir.mkdir(parents=True, exist_ok=True)
    print(f"str-replicant HiFi + somatic mosaicism, depth {args.depth}, {args.seeds} seeds; "
          f"cells = exact germline|exact modal molecule/N (mean signed repeat length error vs germline, bp)\n"
          f"reads: % of molecules at germline length, median/modal molecule delta (bp)\n"
          f"multi: mean #alleles, germline found, subclone found\n")
    print(f"{'locus':<15}{'cond':<13}{'reads':<15}{'multi':<16}"
          + "".join(f"{c:<17}" for c in MOSAIC_CALLERS))
    for cond in conds:
        rows: dict[str, list] = {n: [] for n in names}
        mol: dict[str, list] = {n: [] for n in names}
        for s in range(args.seeds):
            sd = workdir / f"{cond[0]}_seed{s}"
            cfg = workdir / f"{cond[0]}_seed{s}.toml"
            cfg.write_text(mosaic_config(cond, args.depth, 42 + 101 * s, set(names)))
            p = subprocess.run([STR_REPLICANT, "run", str(cfg), "-o", str(sd), "--format", "fastq",
                                "--force", "-t", "2"], capture_output=True)
            if p.returncode != 0:
                sys.exit(f"str-replicant failed: {p.stderr.decode()[-800:]}")
            refs = {h.split()[0]: sq.decode() for h, sq in v.parse_fasta(sd / "reference.fa")}
            lines = (sd / "reads.fq").read_bytes().split(b"\n")
            reads: dict[str, list] = {}
            heads: dict[str, list] = {}
            for i in range(0, len(lines) - 3, 4):
                lid = lines[i][1:].split(b":")[0].decode()
                reads.setdefault(lid, []).append((lines[i + 1], lines[i + 3]))
                heads.setdefault(lid, []).append(lines[i])
            jobs = []
            for lid in names:
                if lid not in refs:
                    continue
                germ = refs[lid][SR_FLANK:-SR_FLANK]
                assert germ == mosaic_allele(lid, 0), lid
                minor = None
                if cond[2]:
                    # The subclone allele: str-replicant's reference holds haplotype 0 only.
                    minor = "".join(mosaic_allele(lid, cond[2][1]))
                mol[lid] += _molecule_deltas(heads.get(lid, []), len(germ))
                # The most common molecule (differs from germline when skewed).
                fr = mosaic_fractions(*cond[1])
                modal = mosaic_allele(lid, max(fr, key=lambda t: t[1])[0])
                jobs.append((lid, germ, minor, modal, refs[lid], reads.get(lid, []), str(sd)))
            with ProcessPoolExecutor(max_workers=args.jobs) as ex:
                for (lid, *_), r in zip(jobs, ex.map(mosaic_trial, jobs)):
                    if r is not None:
                        rows[lid].append(r)
            if not os.environ.get("KEEP_WORK"):
                shutil.rmtree(sd, ignore_errors=True)
        for lid in names:
            rs = rows[lid]
            if not rs:
                continue
            d = mol[lid]
            rd = (f"{100 * sum(x == 0 for x in d) / len(d):.0f}% "
                  f"{statistics.median(d):+.0f}/{statistics.mode(d):+d}")
            m = [r["_multi"] for r in rs]
            ms = (f"{sum(x[0] for x in m) / len(m):.1f} {sum(x[1] for x in m)}/{len(m)}"
                  + (f" {sum(x[2] for x in m)}/{len(m)}" if cond[2] else ""))
            cells = []
            for c in MOSAIC_CALLERS:
                vals = [(r[c], r["_dlen"][c]) for r in rs if r.get(c) is not None]
                if not vals:
                    cells.append(f"{'N/A':<17}")
                    continue
                ex = sum(1 for e, _ in vals if e == 0)
                mo = sum(1 for r in rs if r["_modal"][c])
                dl = sum(x for _, x in vals) / len(vals)
                fail = f"!{len(rs) - len(vals)}" if len(vals) < len(rs) else ""
                cells.append(f"{ex}|{mo}/{len(rs)} ({dl:+.1f}){fail}".ljust(17))
            print(f"{lid:<15}{cond[0]:<13}{rd:<15}{ms:<16}" + "".join(cells), flush=True)


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
    ap.add_argument("--mosaic", action="store_true",
                    help="HiFi (str-replicant) with somatic mosaicism / an expanded subclone")
    ap.add_argument("--conditions", nargs="*", help="--mosaic: subset of condition labels")
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
    if args.mosaic:
        run_mosaic(args)
    if not (args.jitter or args.pbsim or args.vntr or args.replicant or args.mosaic):
        ap.error("pass --jitter, --pbsim, --vntr, --replicant and/or --mosaic")


if __name__ == "__main__":
    main()
