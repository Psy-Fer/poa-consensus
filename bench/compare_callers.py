#!/usr/bin/env python3
"""
Orthogonal consensus-caller comparison harness.

Runs the same simulated-read pipeline as validate.py (pbsim3 -> minimap2 ->
bedpull) but feeds the resulting extracted reads to THREE independent
consensus callers instead of one:

  - poa-consensus  (this crate's CLI)
  - abPOA          (via the pyabpoa bindings, semi-global mode)
  - SPOA           (via the pyspoa bindings, semi-global mode)

All three see byte-identical input, so disagreement is the signal:
  - poa-consensus agrees with both externals  -> strong correctness evidence
  - poa-consensus alone disagrees              -> a bug worth chasing here
  - an external alone disagrees                -> evidence we're doing better
  - all three disagree with truth              -> a ceiling in the data/
                                                    algorithm family, not a
                                                    poa-consensus-specific bug

Scope: single-allele scenarios only. abPOA/SPOA have no equivalent to this
crate's bubble-based multi-allele splitting, so a fair three-way comparison
is restricted to scenarios where all callers solve the same problem: one
read set in, one consensus out. Multi-allele scenarios in validate.py's
catalogue are skipped (reported, not silently dropped).

Setup (one-time):
    pip install pyabpoa pyspoa

Usage:
    python3 bench/compare_callers.py [--scenarios NAME ...] [--workdir PATH] [--keep]

Requires the same native tools as validate.py: pbsim3, minimap2, samtools, bedpull,
and a release build of the poa-consensus CLI (cargo build --release --features cli).
"""

import argparse
import shutil
import sys
from pathlib import Path
from typing import Optional

import validate as v


def fitting_edit_distance(query: bytes, target: bytes) -> int:
    """Minimum edit distance to align all of `query` against some substring of `target`.

    Standard "fitting" alignment: free leading/trailing gaps in `target`, but
    `query` must be fully consumed. This deliberately avoids any extraction
    heuristic (locating flank anchors, clustering literal repeat-unit matches)
    -- both were tried here and both have real failure modes on this benchmark's
    synthetic data: anchor extraction can lock onto the wrong occurrence of a
    k-mer when the flank is itself periodic (as this benchmark's flanks are, by
    design, to keep the flank free of the repeat unit); literal-unit clustering
    discards an entire side of the repeat when a single small indel disrupts the
    reading frame partway through. Fitting alignment sidesteps both: it just asks
    "how close is the best-matching region of the caller's output to the known
    truth allele," with no assumption about where that region starts or ends.
    """
    n, m = len(query), len(target)
    prev = [0] * (m + 1)
    for i in range(1, n + 1):
        cur = [i] + [0] * m
        qi = query[i - 1]
        for j in range(1, m + 1):
            cost = 0 if qi == target[j - 1] else 1
            cur[j] = min(prev[j - 1] + cost, prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return min(prev)


# ── External caller wrappers ─────────────────────────────────────────────────
# Both run in semi-global mode to match this crate's own documented default
# recommendation for extracted STR reads (see CLAUDE.md "Alignment Mode").

def run_abpoa(seqs: list[bytes]) -> Optional[bytes]:
    try:
        import pyabpoa as pa
    except ImportError:
        return None
    aligner = pa.msa_aligner(aln_mode="e")  # extension alignment == semi-global
    res = aligner.msa([s.decode() for s in seqs], out_cons=True, out_msa=False)
    if not res or not res.cons_seq:
        return None
    return res.cons_seq[0].encode()


def run_spoa(seqs: list[bytes]) -> Optional[bytes]:
    try:
        import spoa
    except ImportError:
        return None
    cons, _msa = spoa.poa([s.decode() for s in seqs], algorithm=2)  # 2 = semi-global
    return cons.encode() if cons else None


CALLERS = [
    ("poa-consensus", None),   # handled separately: it's a subprocess, not a Python call
    ("abPOA", run_abpoa),
    ("SPOA", run_spoa),
]


# ── Per-caller evaluation ────────────────────────────────────────────────────

def evaluate_caller(caller: str, cons_bytes: Optional[bytes], truth: "v.Allele") -> dict:
    """Score one caller's consensus against the known truth allele.

    Scores by `fitting_edit_distance()` against the FULL consensus -- no
    extraction step at all. Two extraction-based approaches were tried and
    both broke on this benchmark's synthetic data: literal repeat-unit
    clustering discards an entire side of the repeat when one small indel
    disrupts the reading frame partway through; anchor-based flank extraction
    can lock onto the wrong occurrence of its k-mer because the flank is
    itself periodic. `units_found` (via literal clustering) is kept purely as
    an informational display field -- it can still misreport on a disrupted
    repeat, so don't trust it over `edit_dist` for the pass/fail call.
    """
    if cons_bytes is None:
        return {"caller": caller, "status": "N/A", "units_found": None, "delta_units": None}

    edit = fitting_edit_distance(truth.seq.upper().encode(), cons_bytes.upper())
    edit_tol = max(3, len(truth.seq) // 25)
    ok = edit <= edit_tol

    cluster = v.extract_allele_by_unit(cons_bytes, truth.unit)
    units_found = v.count_repeat(cluster, truth.unit) if cluster else 0
    delta = units_found - truth.count

    return {
        "caller": caller,
        "status": "OK" if ok else "FAIL",
        "units_found": units_found,
        "delta_units": delta,
        "edit_dist": edit,
    }


def cell(r: dict) -> str:
    if r["status"] == "N/A":
        return "N/A (not installed)"
    if r["units_found"] is None:
        return r["status"]
    return f"units={r['units_found']:>3}(Δ{r['delta_units']:+d}) edit={r['edit_dist']:>2}  {r['status']}"


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", nargs="*", help="Scenario names to run (default: all single-allele scenarios)")
    ap.add_argument("--workdir", default="bench/work_compare", help="Working directory for temp files")
    ap.add_argument("--keep", action="store_true", help="Keep working files after run")
    args = ap.parse_args()

    for binary, name in [(v.POA_BIN, "poa-consensus"), ("minimap2", "minimap2"),
                         ("samtools", "samtools"), ("bedpull", "bedpull"), (v.PBSIM, "pbsim")]:
        path = shutil.which(str(binary)) or (Path(str(binary)).exists() and str(binary))
        if not path:
            sys.exit(f"error: {name} not found ({binary}); build CLI with: cargo build --release --features cli")

    missing_py = [pkg for pkg, mod in [("pyabpoa", "pyabpoa"), ("pyspoa", "spoa")]
                  if __import__("importlib").util.find_spec(mod) is None]
    if missing_py:
        sys.exit(f"error: missing Python packages: {', '.join(missing_py)}; run: pip install {' '.join(missing_py)}")

    all_scenarios = [s for s in v.SCENARIOS if not s.is_multi and not s.long_read]
    skipped_multi = [s.name for s in v.SCENARIOS if s.is_multi]
    to_run = all_scenarios
    if args.scenarios:
        requested = set(args.scenarios)
        to_run = [s for s in all_scenarios if s.name in requested]
        unknown = requested - {s.name for s in all_scenarios}
        if unknown:
            print(f"warning: unknown or unsupported (multi-allele) scenarios skipped: {', '.join(sorted(unknown))}",
                  file=sys.stderr)

    print(f"\n{'='*100}")
    print(f"  orthogonal consensus-caller comparison  ({len(to_run)} scenarios; "
          f"{len(skipped_multi)} multi-allele scenarios not applicable, skipped: {', '.join(skipped_multi)})")
    print(f"{'='*100}")
    header = f"  {'SCENARIO':<28} {'MODEL':<9} {'TRUTH':<6} {'poa-consensus':<30} {'abPOA':<30} {'SPOA':<30}"
    print(header)
    print(f"  {'-'*28} {'-'*9} {'-'*5} {'-'*30} {'-'*30} {'-'*30}")

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    tally = {"all_agree_ok": 0, "poa_better": 0, "poa_worse": 0, "all_fail": 0, "error": 0}
    rows = []

    for scenario in to_run:
        work = workdir / scenario.name
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)

        model = v.ERRMODELS[scenario.model]
        truth = scenario.alleles[0]

        try:
            ref = v.build_reference(scenario, work)
            reads = v.simulate_reads(scenario, ref, work, model)
            bam = v.align_and_index(reads, ref, work, model["mm2_preset"])
            extracted_path = v.extract_reads(scenario, bam, work)
            records = v.parse_fasta(extracted_path)
            seqs = [seq for _, seq in records]

            poa_path = v.run_consensus(extracted_path, work, scenario)
            poa_records = v.parse_fasta(poa_path) if poa_path else []
            poa_cons = poa_records[0][1] if poa_records else None

            results = [
                evaluate_caller("poa-consensus", poa_cons, truth),
                evaluate_caller("abPOA", run_abpoa(seqs), truth),
                evaluate_caller("SPOA", run_spoa(seqs), truth),
            ]
        except Exception as exc:
            print(f"  {scenario.name:<28} ERROR: {exc}")
            tally["error"] += 1
            if not args.keep:
                shutil.rmtree(work, ignore_errors=True)
            continue

        rows.append((scenario, truth, results))
        poa_ok = results[0]["status"] == "OK"
        ext_ok = [r["status"] == "OK" for r in results[1:] if r["status"] != "N/A"]
        if poa_ok and all(ext_ok):
            tally["all_agree_ok"] += 1
        elif poa_ok and not all(ext_ok):
            tally["poa_better"] += 1
        elif not poa_ok and any(ext_ok):
            tally["poa_worse"] += 1
        else:
            tally["all_fail"] += 1

        alleles_str = f"{truth.unit}×{truth.count}"
        print(f"  {scenario.name:<28} {scenario.model:<9} {alleles_str:<6} "
              f"{cell(results[0]):<30} {cell(results[1]):<30} {cell(results[2]):<30}")

        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)

    n = len(rows)
    print(f"\n  Result: {n} scenarios compared"
          + (f", {tally['error']} errored" if tally["error"] else ""))
    print(f"    all three agree (OK):                {tally['all_agree_ok']}")
    print(f"    poa-consensus OK, external(s) FAIL:  {tally['poa_better']}   <- evidence we're doing better")
    print(f"    poa-consensus FAIL, external(s) OK:  {tally['poa_worse']}   <- chase these")
    print(f"    all three FAIL:                      {tally['all_fail']}   <- data/algorithm-family ceiling")
    print()


if __name__ == "__main__":
    main()
