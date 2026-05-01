#!/usr/bin/env python3
"""
poa-consensus validation harness.

Full pipeline per scenario:
  1. Build synthetic reference FASTA: left_flank + allele(s) + right_flank
  2. pbsim3  → long reads (--strategy wgs, --method errhmm)
  3. minimap2 → align reads to reference; samtools sort + index
  4. bedpull  → extract reads that span the STR locus
  5. poa-consensus CLI → single or multi-allele consensus
  6. Evaluate: allele length error, repeat-unit count, edit distance

Usage:
    python bench/validate.py [--scenarios NAME ...] [--workdir PATH] [--keep]

Requirements:
    cargo build --release --features cli   (before first run)
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from itertools import permutations
from pathlib import Path
from typing import Optional

# ── Tool paths ──────────────────────────────────────────────────────────────

PBSIM    = Path("~/install/pbsim3/src/pbsim").expanduser()
PBMODELS = Path("~/install/pbsim3/data").expanduser()
POA_BIN  = Path(__file__).parent.parent / "target/release/poa-consensus"

# All other tools (minimap2, samtools, bedpull) are expected on PATH.

# ── Flanks (2048 bp, no tandem repeats) ─────────────────────────────────────

_U = "ACGTACGTCGATCGATTAGCTAGCGCTAGCTA"   # 32 bp unit
LEFT_FLANK  = (_U * 64)[:2048]
RIGHT_FLANK = (_U[::-1] * 64)[:2048]
FLANK_LEN   = 2048
ANCHOR_PAD  = 100   # bp of flanking context included on each side of the STR locus

# ── Error models ────────────────────────────────────────────────────────────
# Constraints baked in:
#   floor(accuracy_mean * 100 * 1.05) must be ≤ model's max accuracy
#   kappa = (length_mean / length_sd)^2 < ~140  (avoids pow() overflow)
#   length_mean must be < reference length

ERRMODELS: dict[str, dict] = {
    "ont_r9": dict(
        errhmm="ERRHMM-ONT.model",
        accuracy_mean=0.85,          # floor(85*1.05)=89 ≤ 99 ✓
        difference_ratio="39:24:36",
        mm2_preset="map-ont",
        label="ONT R9 (85% acc)",
    ),
    "ont_r10": dict(
        errhmm="ERRHMM-ONT-HQ.model",
        accuracy_mean=0.92,          # floor(92*1.05)=96 ≤ 99 ✓
        difference_ratio="39:24:36",
        mm2_preset="map-ont",
        label="ONT R10 HQ (92% acc)",
    ),
    "hifi": dict(
        errhmm="ERRHMM-SEQUEL.model",
        accuracy_mean=0.89,          # floor(89*1.05)=93 ≤ 94 ✓
        difference_ratio="22:45:33",
        mm2_preset="map-hifi",
        label="HiFi/Sequel (89% acc)",
    ),
}

# ── Scenario definition ──────────────────────────────────────────────────────

@dataclass
class Allele:
    unit:   str          # e.g. "CAG"
    count:  int          # repeat count
    depth:  int          # target read depth for this allele

    @property
    def seq(self) -> str:
        return self.unit * self.count

    @property
    def id(self) -> str:
        return f"{self.unit.lower()}{self.count}"


@dataclass
class Scenario:
    name:       str
    alleles:    list[Allele]           # one allele = single; two+ = multi-allele
    model:      str                    # key into ERRMODELS
    multi:      bool  = False          # pass --multi to poa-consensus
    poa_flags:  list  = field(default_factory=list)
    seed:       int   = 42

    @property
    def is_multi(self) -> bool:
        # Only treat as multi-allele if the user explicitly set multi=True.
        # Scenarios with multiple alleles but multi=False (e.g. SV outlier tests)
        # are evaluated as single-allele against alleles[0].
        return self.multi


# ── Scenario catalogue ───────────────────────────────────────────────────────

SCENARIOS: list[Scenario] = [

    # ── Depth sweep: CAG×20, ONT R10 ───────────────────────────────────────
    Scenario("cag20_d05_r10",  [Allele("CAG", 20,  5)],  "ont_r10"),
    Scenario("cag20_d10_r10",  [Allele("CAG", 20, 10)],  "ont_r10"),
    Scenario("cag20_d20_r10",  [Allele("CAG", 20, 20)],  "ont_r10"),
    Scenario("cag20_d30_r10",  [Allele("CAG", 20, 30)],  "ont_r10"),

    # ── Error model comparison at depth 20 ─────────────────────────────────
    Scenario("cag20_d20_r9",   [Allele("CAG", 20, 20)],  "ont_r9"),
    Scenario("cag20_d20_hifi", [Allele("CAG", 20, 20)],  "hifi"),

    # ── Allele length sweep: ONT R10, depth 20 ─────────────────────────────
    Scenario("cag5_d20",       [Allele("CAG",   5, 20)], "ont_r10"),
    Scenario("cag10_d20",      [Allele("CAG",  10, 20)], "ont_r10"),
    Scenario("cag50_d20",      [Allele("CAG",  50, 20)], "ont_r10"),
    Scenario("cag100_d20",     [Allele("CAG", 100, 20)], "ont_r10"),
    Scenario("cag200_d20",     [Allele("CAG", 200, 20)], "ont_r10"),

    # ── FRDA-like GAA, ONT R10 ──────────────────────────────────────────────
    Scenario("gaa50_d20",      [Allele("GAA",  50, 20)], "ont_r10"),
    Scenario("gaa100_d20",     [Allele("GAA", 100, 20)], "ont_r10"),
    Scenario("gaa200_d10",     [Allele("GAA", 200, 10)], "ont_r10"),

    # ── Multi-allele: two alleles each at depth 20 ──────────────────────────
    Scenario("multi_cag15_25",
             [Allele("CAG", 15, 20), Allele("CAG", 25, 20)],
             "ont_r10", multi=True),

    Scenario("multi_cag20_50",
             [Allele("CAG", 20, 20), Allele("CAG", 50, 20)],
             "ont_r10", multi=True),

    Scenario("multi_gaa30_100",
             [Allele("GAA", 30, 20), Allele("GAA", 100, 20)],
             "ont_r10", multi=True),

    # ── SV outlier: majority allele + a few reads from a very different one ─
    # Expect: consensus = majority allele; the 2-read outlier shouldn't win.
    Scenario("sv_cag20_out60",
             [Allele("CAG", 20, 20), Allele("CAG", 60, 2)],
             "ont_r10", multi=False),  # single-allele mode; outlier should be ignored

    Scenario("sv_gaa50_out200",
             [Allele("GAA", 50, 20), Allele("GAA", 200, 2)],
             "ont_r10", multi=False),

    # ── Heterogeneous depth: skewed multi-allele ratio ──────────────────────
    Scenario("multi_skew_cag20_40",
             [Allele("CAG", 20, 24), Allele("CAG", 40, 8)],
             "ont_r10", multi=True),
]


# ── Helpers ──────────────────────────────────────────────────────────────────

def run(cmd: list, *, cwd=None, capture=False, check=True) -> subprocess.CompletedProcess:
    kwargs = dict(cwd=cwd, check=check)
    if capture:
        kwargs["capture_output"] = True
    else:
        kwargs["stdout"] = subprocess.DEVNULL
        kwargs["stderr"] = subprocess.DEVNULL
    return subprocess.run([str(c) for c in cmd], **kwargs)


def levenshtein(a: bytes, b: bytes, max_len: int = 2000) -> int:
    """Edit distance; approximate with |len(a) - len(b)| for sequences > max_len."""
    if len(a) > max_len or len(b) > max_len:
        return abs(len(a) - len(b))
    dp = list(range(len(b) + 1))
    for ca in a:
        dp2 = [dp[0] + 1] + [0] * len(b)
        for j, cb in enumerate(b):
            dp2[j + 1] = min(dp[j] + (0 if ca == cb else 1), dp2[j] + 1, dp[j + 1] + 1)
        dp = dp2
    return dp[-1]


def count_repeat(seq: bytes, unit: str) -> int:
    """Count non-overlapping occurrences of unit in seq (both uppercase)."""
    u = unit.upper().encode()
    s = seq.upper()
    n, i = 0, 0
    while True:
        i = s.find(u, i)
        if i == -1:
            break
        n += 1
        i += len(u)
    return n


def extract_allele_by_unit(cons: bytes, unit: str) -> Optional[bytes]:
    """Extract the allele region from a consensus by locating the repeat unit.

    The repeat unit (e.g. 'CAG', 'GAA') does not appear in the flanking
    sequences (_U or its reverse).  Isolated spurious occurrences in the flank
    (from sequencing errors) are excluded by finding the largest CLUSTER of
    unit positions, where two positions belong to the same cluster if their
    gap is <= 2*len(unit) (tolerates single-base indel errors in the repeat).
    """
    u = unit.upper().encode()
    s = cons.upper()
    n = len(u)

    positions: list[int] = []
    i = 0
    while True:
        i = s.find(u, i)
        if i == -1:
            break
        positions.append(i)
        i += n

    if not positions:
        return None
    if len(positions) == 1:
        return cons[positions[0]: positions[0] + n]

    # Group into clusters: two consecutive positions are in the same cluster
    # when their gap is <= 2*n (allows one extra/missing base per unit).
    clusters: list[list[int]] = [[positions[0]]]
    for p in positions[1:]:
        if p - clusters[-1][-1] <= 2 * n:
            clusters[-1].append(p)
        else:
            clusters.append([p])

    best = max(clusters, key=len)
    return cons[best[0]: best[-1] + n]


def find_allele_in_consensus(cons: bytes, anchor_k: int = 20) -> Optional[bytes]:
    """Extract the allele region from a consensus using the known flanking anchors.

    The flanking sequences are repetitive (_U repeated), so the anchor k-mer
    appears multiple times in the flank.  Use rfind for the left anchor (last
    occurrence before the allele) and find for the right anchor (first
    occurrence after the allele).
    """
    la = LEFT_FLANK[-anchor_k:].encode()
    ra = RIGHT_FLANK[:anchor_k].encode()
    upper = cons.upper()
    # Right anchor: first occurrence anywhere in the consensus
    ri = upper.find(ra.upper())
    if ri == -1:
        return None
    # Left anchor: last occurrence strictly before the right anchor
    li = upper[:ri].rfind(la.upper())
    if li == -1 or ri <= li + anchor_k:
        return None
    return cons[li + anchor_k: ri]


def parse_fasta(path: Path) -> list[tuple[str, bytes]]:
    """Parse a FASTA file into (header, sequence) pairs."""
    records, header, seq = [], None, []
    with open(path, "rb") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(b">"):
                if header is not None:
                    records.append((header, b"".join(seq)))
                header = line[1:].decode()
                seq = []
            else:
                seq.append(line)
    if header is not None:
        records.append((header, b"".join(seq)))
    return records


# ── Pipeline steps ───────────────────────────────────────────────────────────

def build_reference(scenario: Scenario, work: Path) -> Path:
    """Write one FASTA record per allele: LEFT_FLANK + allele + RIGHT_FLANK."""
    ref = work / "reference.fa"
    with open(ref, "w") as fh:
        for i, allele in enumerate(scenario.alleles):
            seq = LEFT_FLANK + allele.seq + RIGHT_FLANK
            fh.write(f">allele_{i}  unit={allele.unit} count={allele.count}\n{seq}\n")
    return ref


def simulate_reads(scenario: Scenario, ref: Path, work: Path, model: dict) -> Path:
    """Run pbsim on each allele and merge into a single FASTQ."""
    allele_len_max = max(len(a.seq) for a in scenario.alleles)
    ref_len = FLANK_LEN * 2 + allele_len_max

    # Length parameters: mean = ~75 % of reference length; sd = mean / 5.
    length_mean = max(500, int(ref_len * 0.75))
    length_mean = min(length_mean, ref_len - 100)   # must be < ref_len
    length_sd   = max(50, length_mean // 5)          # kappa ≤ 25

    # pbsim outputs one .fq.gz per FASTA sequence (allele_0001.fq.gz, etc.)
    prefix = str(work / "sim")
    cmd = [
        PBSIM,
        "--strategy", "wgs",
        "--genome", ref,
        "--method", "errhmm",
        "--errhmm", PBMODELS / model["errhmm"],
        "--accuracy-mean", str(model["accuracy_mean"]),
        "--difference-ratio", model["difference_ratio"],
        "--depth", "1",        # placeholder; we'll use per-allele depth below
        "--length-mean", str(length_mean),
        "--length-sd",   str(length_sd),
        "--seed", str(scenario.seed),
        "--prefix", prefix,
    ]

    # pbsim --depth applies to ALL sequences equally. To get per-allele depth,
    # run pbsim once per allele on a single-sequence FASTA, then merge.
    merged = work / "reads.fq"
    with open(merged, "wb") as out:
        for i, allele in enumerate(scenario.alleles):
            single_ref = work / f"ref_allele{i}.fa"
            with open(single_ref, "w") as fh:
                seq = LEFT_FLANK + allele.seq + RIGHT_FLANK
                fh.write(f">allele_{i}\n{seq}\n")

            ap = str(work / f"sim_a{i}")
            sim_cmd = [
                PBSIM,
                "--strategy", "wgs",
                "--genome", single_ref,
                "--method", "errhmm",
                "--errhmm", PBMODELS / model["errhmm"],
                "--accuracy-mean", str(model["accuracy_mean"]),
                "--difference-ratio", model["difference_ratio"],
                "--depth", str(allele.depth),
                "--length-mean", str(length_mean),
                "--length-sd",   str(length_sd),
                "--seed", str(scenario.seed + i),
                "--prefix", ap,
            ]
            run(sim_cmd)
            fq_gz = work / f"sim_a{i}_0001.fq.gz"
            if fq_gz.exists():
                with gzip.open(fq_gz) as fq_in:
                    out.write(fq_in.read())

    return merged


def align_and_index(reads: Path, ref: Path, work: Path, mm2_preset: str) -> Path:
    """Align reads to reference with minimap2, sort and index BAM."""
    bam = work / "aligned.bam"
    # Pipe minimap2 → samtools sort in one step
    mm2 = subprocess.Popen(
        ["minimap2", "-a", f"-x{mm2_preset}", "--secondary=no", str(ref), str(reads)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    with open(bam, "wb") as bam_fh:
        subprocess.run(
            ["samtools", "sort", "-"],
            stdin=mm2.stdout, stdout=bam_fh, stderr=subprocess.DEVNULL, check=True,
        )
    mm2.wait()
    run(["samtools", "index", bam])
    return bam


def extract_reads(scenario: Scenario, bam: Path, work: Path) -> Path:
    """Use bedpull to extract reads spanning the STR locus from each allele.

    ANCHOR_PAD bp of flanking sequence are included on each side so that
    poa-consensus sees non-repetitive context, eliminating phase rotation.
    """
    bed = work / "targets.bed"
    with open(bed, "w") as fh:
        for i, allele in enumerate(scenario.alleles):
            start = max(0, FLANK_LEN - ANCHOR_PAD)
            end   = FLANK_LEN + len(allele.seq) + ANCHOR_PAD
            fh.write(f"allele_{i}\t{start}\t{end}\n")

    extracted = work / "extracted.fa"
    run(["bedpull", "-b", bam, "-r", bed, "-o", extracted])
    return extracted


def run_consensus(extracted: Path, work: Path, scenario: Scenario) -> Optional[Path]:
    """Run poa-consensus CLI on extracted reads."""
    consensus = work / "consensus.fa"
    # For single-allele scenarios, use a wider minimum band (50) so that
    # long repetitive alleles (e.g. CAG×100) are aligned correctly.
    # Multi-allele mode uses only the adaptive band because a fixed wider
    # band disrupts bubble detection and allele partitioning.
    band_args = ["--adaptive-band"] + ([] if scenario.multi else ["--band-width", "50"])
    cmd = (
        [str(POA_BIN)] + band_args + ["--semi-global"]
        + (["--multi"] if scenario.multi else [])
        + scenario.poa_flags
        + ["--min-reads", "2", str(extracted)]
    )
    try:
        result = subprocess.run(cmd, capture_output=True, check=True)
        consensus.write_bytes(result.stdout)
        return consensus
    except subprocess.CalledProcessError:
        return None


def evaluate(scenario: Scenario, consensus_path: Optional[Path]) -> dict:
    """Compare consensus output to expected alleles; return metrics dict."""
    result = {
        "scenario": scenario.name,
        "model":    scenario.model,
        "alleles":  "+".join(f"{a.unit}×{a.count}" for a in scenario.alleles),
        "total_depth": sum(a.depth for a in scenario.alleles),
    }

    if consensus_path is None or not consensus_path.exists():
        result["status"] = "FAILED (no consensus)"
        return result

    records = parse_fasta(consensus_path)
    if not records:
        result["status"] = "FAILED (empty output)"
        return result

    if not scenario.is_multi:
        # Single allele: extract the allele region by locating the repeat unit.
        # The unit (e.g. 'CAG') does not appear in the flanking sequences, so
        # all occurrences in the consensus come from the allele.
        cons_seq  = records[0][1]
        truth     = scenario.alleles[0]
        extracted = extract_allele_by_unit(cons_seq, truth.unit)
        result["units_truth"] = truth.count
        if extracted is None:
            # Repeat unit not found at all — completely wrong consensus.
            result["anchor_found"] = False
            result["cons_len"]     = len(cons_seq)
            result["truth_len"]    = len(truth.seq)
            result["delta_len"]    = len(cons_seq) - len(truth.seq)
            result["units_found"]  = 0
            result["delta_units"]  = -truth.count
            ok = False
        else:
            result["anchor_found"] = True
            result["cons_len"]     = len(extracted)
            result["truth_len"]    = len(truth.seq)
            result["delta_len"]    = len(extracted) - len(truth.seq)
            result["units_found"]  = count_repeat(extracted, truth.unit)
            result["delta_units"]  = result["units_found"] - truth.count
            result["edit_dist"]    = levenshtein(extracted.upper(), truth.seq.upper().encode())
            # Tolerance scales with repeat count: ±1 for ≤50 units, ±2 for ≤100, ±4 for ≤200, etc.
            unit_tol = max(1, truth.count // 50)
            edit_tol = max(3, len(truth.seq) // 25)
            ok = (abs(result["delta_units"]) <= unit_tol) and result["edit_dist"] <= edit_tol

        result["status"] = "OK" if ok else "FAIL"

    else:
        # Multi-allele: optimally assign each output record to a truth allele.
        truths = scenario.alleles
        if len(records) != len(truths):
            result["status"] = f"FAIL (expected {len(truths)} alleles, got {len(records)})"
            result["n_alleles_found"] = len(records)
            result["n_alleles_truth"] = len(truths)
            return result

        # Count repeat units directly in the full consensus sequence.
        # Since the unit (e.g. 'CAG') does not appear in the perfect flanking
        # sequences, spurious counts from POA-averaged flank context are rare.
        # This is more robust than boundary extraction for noisy repeat regions.
        unit = truths[0].unit
        found_units_list = [count_repeat(cons_seq, unit) for _, cons_seq in records]

        # Find the optimal one-to-one assignment (minimize total |delta|).
        best_assignment = None
        best_cost = float('inf')
        for perm in permutations(range(len(truths))):
            cost = sum(abs(found_units_list[i] - truths[perm[i]].count) for i in range(len(records)))
            if cost < best_cost:
                best_cost = cost
                best_assignment = perm

        matched = []
        for i, truth_idx in enumerate(best_assignment):
            t = truths[truth_idx]
            found = found_units_list[i]
            matched.append({
                "truth_units": t.count,
                "found_units": found,
                "delta_units": found - t.count,
            })

        result["allele_results"] = matched
        all_ok = all(abs(m["delta_units"]) <= max(1, m["truth_units"] // 50) for m in matched)
        result["status"] = "OK" if all_ok else "FAIL"

    return result


# ── Output formatting ─────────────────────────────────────────────────────────

def print_result(r: dict):
    if "allele_results" in r:
        allele_str = "  ".join(
            f"({m['truth_units']}→{m['found_units']} Δ{m['delta_units']:+d})"
            for m in r["allele_results"]
        )
        print(f"  {r['scenario']:<35} {r['model']:<10} d={r['total_depth']:<3}  "
              f"multi: {allele_str}  {r['status']}")
    else:
        unit_info = f"units: {r.get('units_truth','?')}→{r.get('units_found','?')} (Δ{r.get('delta_units',0):+d})"
        edit_info = f"edit={r.get('edit_dist', '?')}" if "edit_dist" in r else f"Δlen={r.get('delta_len',0):+d}"
        anchor    = "⚓" if r.get("anchor_found") else "~"
        print(f"  {r['scenario']:<35} {r['model']:<10} d={r['total_depth']:<3}  "
              f"{unit_info}  {edit_info}  {anchor}  {r['status']}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", nargs="*", help="Scenario names to run (default: all)")
    ap.add_argument("--workdir", default="bench/work", help="Working directory for temp files")
    ap.add_argument("--keep", action="store_true", help="Keep working files after run")
    args = ap.parse_args()

    # Verify required binaries
    for binary, name in [(POA_BIN, "poa-consensus"), ("minimap2", "minimap2"),
                         ("samtools", "samtools"), ("bedpull", "bedpull"), (PBSIM, "pbsim")]:
        path = shutil.which(str(binary)) or (Path(str(binary)).exists() and str(binary))
        if not path:
            sys.exit(f"error: {name} not found ({binary}); build CLI with: cargo build --release --features cli")

    to_run = SCENARIOS
    if args.scenarios:
        to_run = [s for s in SCENARIOS if s.name in args.scenarios]
        unknown = set(args.scenarios) - {s.name for s in SCENARIOS}
        if unknown:
            print(f"warning: unknown scenarios: {', '.join(sorted(unknown))}", file=sys.stderr)

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*80}")
    print(f"  poa-consensus validation  ({len(to_run)} scenarios)")
    print(f"{'='*80}")
    print(f"  {'SCENARIO':<35} {'MODEL':<10} {'DEPTH':<6}  ACCURACY  STATUS")
    print(f"  {'-'*35} {'-'*10} {'-'*5}  {'─'*30}  ──────")

    results, n_ok, n_fail = [], 0, 0

    for scenario in to_run:
        work = workdir / scenario.name
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)

        model = ERRMODELS[scenario.model]

        try:
            ref        = build_reference(scenario, work)
            reads      = simulate_reads(scenario, ref, work, model)
            bam        = align_and_index(reads, ref, work, model["mm2_preset"])
            extracted  = extract_reads(scenario, bam, work)
            consensus  = run_consensus(extracted, work, scenario)
            r          = evaluate(scenario, consensus)
        except Exception as exc:
            r = {"scenario": scenario.name, "model": scenario.model,
                 "total_depth": sum(a.depth for a in scenario.alleles),
                 "alleles": "+".join(f"{a.unit}×{a.count}" for a in scenario.alleles),
                 "status": f"ERROR: {exc}"}

        results.append(r)
        print_result(r)

        if r["status"] == "OK":
            n_ok += 1
        else:
            n_fail += 1

        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)

    print(f"\n  Result: {n_ok} passed, {n_fail} failed out of {len(to_run)} scenarios")
    print()


if __name__ == "__main__":
    main()
