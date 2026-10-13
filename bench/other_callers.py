"""
Non-POA consensus callers for the comparison harnesses, to separate what is a
POA-family strength/weakness from what other approaches do on the same reads.

  medoid   the read with the smallest summed edit distance to the others (edlib);
           a floor: no consensus at all.
  pileup   medoid as draft, reads aligned with minimap2, `samtools consensus`
           (simple majority), two rounds (round 2 uses round 1 as the draft).
  MAFFT, Kalign
           progressive MSA + per-column majority (same majority rule as the
           POASTA cell in compare_callers).
  hifiasm  de novo assembly of the read set (HiFi mode, or --ont for ONT
           models); longest primary contig, oriented to the medoid read.
  TRGT     PacBio TRGT genotype, on reads aligned to a catalog reference.
  LongTR   LongTR, same BAM.

TRGT and LongTR need a BAM against a reference that holds the locus. The bench
references (validate.build_reference) contain the TRUTH allele, which would let
a reference-aware caller echo it, so these get their own catalog reference with
a different repeat count (`catalog_allele`).

Binaries: TRGT_BIN / LONGTR_BIN / HIFIASM_BIN env vars, else `trgt` / `LongTR` /
`hifiasm` on PATH (hifiasm also checked at ~/repos/hifiasm/hifiasm).
"""

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import edlib

TRGT_BIN = os.environ.get("TRGT_BIN") or shutil.which("trgt")
LONGTR_BIN = os.environ.get("LONGTR_BIN") or shutil.which("LongTR")
HIFIASM_BIN = (os.environ.get("HIFIASM_BIN") or shutil.which("hifiasm")
               or next((str(p) for p in [Path.home() / "repos/hifiasm/hifiasm"] if p.exists()), None))


def _run(cmd, timeout: int = 600, **kw) -> subprocess.CompletedProcess:
    # stdin=/dev/null unless input is given: kalign reads stdin when it is not a
    # TTY, which hangs forever under a process pool.
    if "input" not in kw:
        kw["stdin"] = subprocess.DEVNULL
    try:
        return subprocess.run([str(c) for c in cmd], capture_output=True, timeout=timeout, **kw)
    except subprocess.TimeoutExpired:
        return subprocess.CompletedProcess(cmd, 124, b"", b"timeout")


def medoid(seqs: list[bytes]) -> Optional[bytes]:
    if not seqs:
        return None
    s = [x.decode() for x in seqs]
    best, best_i = None, 0
    for i, a in enumerate(s):
        tot = 0
        for j, b in enumerate(s):
            if i != j:
                tot += edlib.align(a, b, mode="NW", task="distance")["editDistance"]
                if best is not None and tot >= best:
                    break
        if best is None or tot < best:
            best, best_i = tot, i
    return seqs[best_i]


def _fasta_first(path: Path) -> Optional[bytes]:
    if not path.exists():
        return None
    seq = []
    started = False
    for line in path.read_bytes().split(b"\n"):
        if line.startswith(b">"):
            if started:
                break
            started = True
        else:
            seq.append(line.strip())
    out = b"".join(seq).replace(b"N", b"")
    return out or None


def pileup(seqs: list[bytes], work: Path, preset: str, rounds: int = 2) -> Optional[bytes]:
    """Align-to-draft majority consensus, starting from the medoid read."""
    draft = medoid(seqs)
    if not draft:
        return None
    work.mkdir(parents=True, exist_ok=True)
    reads = work / "pu_reads.fa"
    with open(reads, "wb") as fh:
        for i, s in enumerate(seqs):
            fh.write(b">r%d\n%s\n" % (i, s))
    for r in range(rounds):
        ref = work / f"pu_draft{r}.fa"
        ref.write_bytes(b">draft\n" + draft + b"\n")
        bam = work / f"pu{r}.bam"
        mm = _run(["minimap2", "-ax", preset, "--secondary=no", ref, reads])
        if mm.returncode != 0:
            return None
        st = _run(["samtools", "sort", "-o", bam, "-"], input=mm.stdout)
        if st.returncode != 0:
            return None
        out = work / f"pu_cons{r}.fa"
        c = _run(["samtools", "consensus", "-m", "simple", "-c", "0.5", "-d", "2",
                  "--show-ins", "yes", "--show-del", "no", "-f", "fasta", "-o", out, bam])
        nxt = _fasta_first(out) if c.returncode == 0 else None
        if not nxt:
            return draft if r > 0 else None
        draft = nxt
    return draft


def _write_reads(path: Path, seqs: list[bytes]):
    with open(path, "wb") as fh:
        for i, s in enumerate(seqs):
            fh.write(b">r%d\n%s\n" % (i, s))


def _majority(aln: bytes) -> Optional[bytes]:
    from compare_callers import _msa_majority_consensus
    return _msa_majority_consensus(aln.upper())


def mafft(seqs: list[bytes], work: Path) -> Optional[bytes]:
    if not shutil.which("mafft"):
        return None
    work.mkdir(parents=True, exist_ok=True)
    fa = work / "msa_in.fa"
    _write_reads(fa, seqs)
    p = _run(["mafft", "--auto", "--quiet", "--thread", "1", fa])
    return _majority(p.stdout) if p.returncode == 0 else None


def kalign(seqs: list[bytes], work: Path) -> Optional[bytes]:
    if not shutil.which("kalign"):
        return None
    work.mkdir(parents=True, exist_ok=True)
    fa, out = work / "msa_in.fa", work / "kalign_out.fa"
    _write_reads(fa, seqs)
    p = _run(["kalign", "-i", fa, "-o", out, "-f", "fasta"])
    return _majority(out.read_bytes()) if p.returncode == 0 and out.exists() else None


def _limit_mem(gb: int):
    def f():
        import resource
        resource.setrlimit(resource.RLIMIT_AS, (gb << 30, gb << 30))
    return f


def _revcomp(s: bytes) -> bytes:
    return s.translate(bytes.maketrans(b"ACGTacgt", b"TGCAtgca"))[::-1]


def hifiasm(seqs: list[bytes], work: Path, ont: bool) -> Optional[bytes]:
    if not HIFIASM_BIN or not seqs:
        return None
    work.mkdir(parents=True, exist_ok=True)
    fa = work / "asm_reads.fq"  # --ont mode rejects FASTA
    with open(fa, "wb") as fh:
        for i, s in enumerate(seqs):
            fh.write(b"@r%d\n%s\n+\n%s\n" % (i, s, b"5" * len(s)))
    prefix = work / "asm"
    # -f0: no k-mer bloom filter (the default -f37 allocates 16 GB even for a
    # 4 kb target; hifiasm recommends -f0 for small genomes).
    cmd = [HIFIASM_BIN, "-o", prefix, "-t", "1", "-f0", "-l0", "--n-hap", "1"]
    if ont:
        cmd.append("--ont")
    if _run(cmd + [fa], preexec_fn=_limit_mem(8)).returncode != 0:
        return None
    gfa = Path(f"{prefix}.bp.p_ctg.gfa")
    if not gfa.exists():
        return None
    ctgs = [l.split("\t")[2].encode() for l in gfa.read_text().splitlines()
            if l.startswith("S\t")]
    if not ctgs:
        return None
    ctg = max(ctgs, key=len)
    ref = medoid(seqs).decode()
    fwd = edlib.align(ref, ctg.decode(), mode="HW", task="distance")["editDistance"]
    rev = edlib.align(ref, _revcomp(ctg).decode(), mode="HW", task="distance")["editDistance"]
    return ctg if fwd <= rev else _revcomp(ctg)


# ── BAM-based STR genotypers ─────────────────────────────────────────────────

def catalog_allele(motif: str, truth_len: int) -> str:
    """Reference allele for the catalog: pure motif, ~75% of the truth length
    (so it differs from the sample and holds no shared interruptions)."""
    n = max(2, round(0.75 * truth_len / len(motif)))
    return motif * n


def bam_callers(reads_fq: Path, left: str, right: str, ref_allele: str, motif: str,
                work: Path, preset: str) -> dict:
    """Align `reads_fq` to left+ref_allele+right and run TRGT / LongTR on it."""
    res = {"TRGT": None, "LongTR": None}
    if not (TRGT_BIN or LONGTR_BIN):
        return res
    work.mkdir(parents=True, exist_ok=True)
    ref = work / "catalog_ref.fa"
    ref.write_text(f">locus\n{left}{ref_allele}{right}\n")
    _run(["samtools", "faidx", ref])
    start, end = len(left), len(left) + len(ref_allele)
    mm = _run(["minimap2", "-ax", preset, "-y", "-R", r"@RG\tID:sim\tSM:sim\tLB:sim", ref, reads_fq])
    bam = work / "catalog.bam"
    if mm.returncode != 0 or _run(["samtools", "sort", "-o", bam, "-"], input=mm.stdout).returncode != 0:
        return res
    _run(["samtools", "index", bam])

    if TRGT_BIN:
        bed = work / "trgt.bed"
        bed.write_text(f"locus\t{start}\t{end}\tID=L;MOTIFS={motif};STRUC=({motif})n\n")
        prefix = work / "trgt"
        p = _run([TRGT_BIN, "genotype", "-g", ref, "-r", bam, "-b", bed, "-o", prefix,
                  "--sample-name", "sim"])
        (work / "trgt.log").write_bytes(p.stdout + p.stderr)
        vcf = Path(f"{prefix}.vcf.gz")
        if p.returncode == 0 and vcf.exists():
            res["TRGT"] = _vcf_allele(vcf)

    if LONGTR_BIN:
        bed = work / "longtr.bed"
        bed.write_text(f"locus\t{start + 1}\t{end}\t{motif}\tL\n")
        vcf = work / "longtr.vcf.gz"
        p = _run([LONGTR_BIN, "--bams", bam, "--fasta", ref, "--regions", bed,
                  "--tr-vcf", vcf, "--haploid-chrs", "locus",
                  # pbsim qualities are low; don't let the default Q30 filter drop every read.
                  "--min-mean-qual", "0", "--min-reads", "3"])
        (work / "longtr.log").write_bytes(p.stdout + p.stderr)
        if p.returncode == 0 and vcf.exists():
            res["LongTR"] = _vcf_allele(vcf)
    return res


def _vcf_allele(vcf: Path) -> Optional[bytes]:
    """First called allele of the first record (REF/ALT chosen by GT)."""
    p = _run(["bash", "-c", f"gzip -dc {vcf} 2>/dev/null || cat {vcf}"])
    for line in p.stdout.decode().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 10:
            return None
        ref, alts = f[3], [a for a in f[4].split(",") if a != "."]
        fmt = f[8].split(":")
        gt = dict(zip(fmt, f[9].split(":"))).get("GT", ".")
        idx = gt.replace("|", "/").split("/")[0]
        if not idx.isdigit():
            return None
        k = int(idx)
        allele = ref if k == 0 else alts[k - 1] if k - 1 < len(alts) else None
        return allele.encode() if allele else None
    return None
