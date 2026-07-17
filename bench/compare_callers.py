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

    `truth.desc is not None` marks a GENERAL_SCENARIOS entry (a non-repeat
    "allele" represented as `Allele(unit=<raw sequence>, count=1, ...)` --
    see GENERAL_SCENARIOS below): `units_found`/`delta_units` are meaningless
    for a non-repeat truth (there is no repeat unit to cluster), so they are
    skipped entirely rather than reporting a misleading number.
    """
    if cons_bytes is None:
        return {"caller": caller, "status": "N/A", "units_found": None, "delta_units": None}

    edit = fitting_edit_distance(truth.seq.upper().encode(), cons_bytes.upper())
    edit_tol = max(3, len(truth.seq) // 25)
    ok = edit <= edit_tol

    if truth.desc is not None:
        return {
            "caller": caller,
            "status": "OK" if ok else "FAIL",
            "units_found": None,
            "delta_units": None,
            "edit_dist": edit,
        }

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
        # General (non-repeat) scenario: no unit-count field applies, but
        # edit_dist (if present) is still the actual scoring signal.
        if "edit_dist" in r:
            return f"edit={r['edit_dist']:>4}  {r['status']}"
        return r["status"]
    return f"units={r['units_found']:>3}(Δ{r['delta_units']:+d}) edit={r['edit_dist']:>2}  {r['status']}"


# ── General (non-repeat) scenarios ───────────────────────────────────────────
#
# Every scenario above (and every scenario in bench/validate.py's own
# SCENARIOS catalogue) is a tandem-repeat allele (Allele(unit, count, depth),
# .seq == unit * count). poa-consensus is a general POA consensus library
# (see CLAUDE.md's Competitive Landscape section, positioned against
# abPOA/SPOA/POASTA generally, not just for STR calling), so this crate's own
# three-way orthogonal comparison against those tools has had zero coverage of
# general, non-repeat "hard" MSA/POA scenarios until now. This section adds
# that coverage, reusing the identical pbsim3 -> minimap2 -> bedpull pipeline
# and fitting_edit_distance() scoring above -- only the truth sequence
# construction differs.
#
# A non-repeat "allele" is represented as `Allele(unit=<raw sequence>,
# count=1, depth=N, desc=<short label>)` -- `.seq` (`unit * count`) is then
# exactly the raw sequence with no other change needed anywhere in
# validate.py's pipeline helpers (build_reference/simulate_reads/
# align_and_index/extract_reads all only ever touch `.seq`/`.depth`).
#
# Scenarios were chosen to be structurally distinct failure-mode categories,
# not "more of the same" (CLAUDE.md's own testing philosophy: "Tests should
# be generic rather than repeat-specific... finds exactly where things break
# as well as where they work"), informed by general MSA/POA literature and
# by abPOA's/SPOA's own test data (both tools' own test_data/test fixtures
# are short (~50bp), non-repetitive, general-ACGT sequences with scattered
# substitution/indel noise across many reads -- confirmed by inspection of
# github.com/yangao07/abPOA/test_data and github.com/rvaser/spoa/test, not
# assumed):
#
#   1. random baseline (ONT R10 and HiFi)   -- sanity: does general sequence
#      work AT ALL through this crate's periodic-repeat-tuned machinery
#      (adaptive band, spine minimizer anchoring, slide-lock) when there is
#      no periodicity to exploit.
#   2. GC-content extremes                  -- error-profile interaction; GC
#      bias affects real sequencer error rates and is a known confounder in
#      general MSA benchmarks, distinct from anything the repeat-only
#      catalogue can probe (every existing flank/unit is already a fixed,
#      moderate-GC synthetic sequence).
#   3. embedded homopolymer run              -- classic ONT/PacBio failure
#      mode (single-base run, not a multi-base tandem unit): distinct
#      structurally from every existing CAG/GAA/CTTTT/AAAAG scenario, which
#      are all >=3bp-period repeats.
#   4. local short-range exact duplication   -- a ~20bp segment repeated
#      immediately adjacent once (not a long, high-copy-number STR): a
#      different scale/shape of alignment ambiguity than this crate's usual
#      focus, closer to a large-scale duplication/CNV-adjacent structure.
#   5. large non-periodic structural insert  -- a substantial (250bp)
#      non-repeating insert as part of the single true allele (not a
#      competing minority allele -- external tools have no bubble-splitting
#      equivalent, so this tests single-consensus reconstruction of a large
#      non-periodic feature, not multi-allele phasing).
#   6. short reads relative to a long region -- semi-global partial-read
#      assembly stress test on non-repeat content (already tested for
#      repeats via lr_*/partial-read unit tests, never against external
#      tools, never on non-periodic sequence).
#   7. shallow depth                         -- low-depth stress test on
#      general sequence, mirroring cag20_d05_r10's repeat equivalent.
#   8. short-period microsatellite (2bp)     -- a different period length
#      class (dinucleotide) from this crate's existing 3-5bp focus; common
#      in real genomes, a distinct edge of the same underlying mechanism.

def _rng_seq(rng, n: int, gc: float = 0.5) -> str:
    """`n` bases drawn i.i.d. from a fixed GC-content distribution."""
    at_each = (1.0 - gc) / 2.0
    gc_each = gc / 2.0
    weights = [at_each, gc_each, gc_each, at_each]  # A, C, G, T
    bases = "ACGT"
    return "".join(rng.choices(bases, weights=weights, k=n))


def _build_general_scenarios() -> list["v.Scenario"]:
    import random

    scenarios = []

    # 1. Random baseline, no structure at all -- ONT R10 and HiFi.
    rng = random.Random(20260001)
    seq = _rng_seq(rng, 300)
    scenarios.append(v.Scenario(
        "gen_random_ont_r10",
        [v.Allele(unit=seq, count=1, depth=20, desc="rand300")],
        "ont_r10",
    ))
    scenarios.append(v.Scenario(
        "gen_random_hifi",
        [v.Allele(unit=seq, count=1, depth=20, desc="rand300")],
        "hifi",
    ))

    # 2. GC-content extremes.
    rng = random.Random(20260002)
    gc_high = _rng_seq(rng, 300, gc=0.80)
    scenarios.append(v.Scenario(
        "gen_gc_high_ont_r10",
        [v.Allele(unit=gc_high, count=1, depth=20, desc="gc80")],
        "ont_r10",
    ))
    rng = random.Random(20260003)
    at_rich = _rng_seq(rng, 300, gc=0.20)
    scenarios.append(v.Scenario(
        "gen_at_rich_ont_r9",
        [v.Allele(unit=at_rich, count=1, depth=20, desc="gc20_r9")],
        "ont_r9",
    ))

    # 3. Embedded homopolymer run (single base, not a multi-base repeat unit)
    # in the middle of otherwise unstructured sequence.
    rng = random.Random(20260004)
    left = _rng_seq(rng, 140)
    right = _rng_seq(rng, 140)
    homopolymer_seq = left + "A" * 20 + right
    scenarios.append(v.Scenario(
        "gen_homopolymer_embedded_ont_r10",
        [v.Allele(unit=homopolymer_seq, count=1, depth=20, desc="homopolymer20A")],
        "ont_r10",
    ))

    # 4. Local short-range exact duplication: one 20bp segment immediately
    # repeated once (not a long high-copy-number STR).
    rng = random.Random(20260005)
    pre = _rng_seq(rng, 100)
    dup_unit = _rng_seq(rng, 20)
    post = _rng_seq(rng, 160)
    dup_seq = pre + dup_unit + dup_unit + post
    scenarios.append(v.Scenario(
        "gen_short_tandem_dup_ont_r10",
        [v.Allele(unit=dup_seq, count=1, depth=20, desc="localdup20x2")],
        "ont_r10",
    ))

    # 5. Large non-periodic structural insert: single true allele = 150bp
    # core + a distinct 250bp non-repeating insert, every read carries it.
    rng = random.Random(20260006)
    core_a = _rng_seq(rng, 75)
    insert = _rng_seq(rng, 250)
    core_b = _rng_seq(rng, 75)
    insert_seq = core_a + insert + core_b
    scenarios.append(v.Scenario(
        "gen_large_structural_insert_ont_r10",
        [v.Allele(unit=insert_seq, count=1, depth=20, desc="core+250ins")],
        "ont_r10",
    ))

    # 6. Short reads relative to a long non-repeat region: force length_mean
    # short enough that reads only partially span the 500bp truth. Also
    # shrinks flank_len (default 2048bp would dominate the reference and
    # make a 300-400bp read miss the truth region almost entirely -- the
    # point here is partial coverage of the TARGET region, not of the whole
    # reference), and sets partial=True: bedpull defaults to spanning-reads-
    # only, which would silently extract ZERO reads here since none of them
    # fully span by design (confirmed empirically -- see investigation
    # notes; this is a harness gap, not a poa-consensus bug).
    rng = random.Random(20260007)
    long_seq = _rng_seq(rng, 500)
    scenarios.append(v.Scenario(
        "gen_short_reads_long_region_ont_r10",
        [v.Allele(unit=long_seq, count=1, depth=20, desc="partial500")],
        "ont_r10",
        length_mean=350,
        flank_len=150,
        partial=True,
    ))

    # 7. Shallow depth on general sequence.
    rng = random.Random(20260008)
    shallow_seq = _rng_seq(rng, 300)
    scenarios.append(v.Scenario(
        "gen_low_depth_ont_r10",
        [v.Allele(unit=shallow_seq, count=1, depth=6, desc="rand300_d6")],
        "ont_r10",
    ))

    # 8. Short-period (dinucleotide) microsatellite -- a different period
    # length class from this crate's usual 3-5bp CAG/GAA/CTTTT/AAAAG focus.
    rng = random.Random(20260009)
    ms_left = _rng_seq(rng, 140)
    ms_right = _rng_seq(rng, 146)
    micro_seq = ms_left + "AT" * 7 + ms_right
    scenarios.append(v.Scenario(
        "gen_microsatellite_at7_ont_r10",
        [v.Allele(unit=micro_seq, count=1, depth=20, desc="AT×7micro")],
        "ont_r10",
    ))

    return scenarios


GENERAL_SCENARIOS = _build_general_scenarios()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", nargs="*", help="Scenario names to run (default: all single-allele scenarios)")
    ap.add_argument("--workdir", default="bench/work_compare", help="Working directory for temp files")
    ap.add_argument("--keep", action="store_true", help="Keep working files after run")
    ap.add_argument("--general", action="store_true",
                    help="Run the general (non-repeat) scenario battery instead of the repeat-based catalogue")
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

    if args.general:
        all_scenarios = GENERAL_SCENARIOS
        skipped_multi = []
    else:
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

    label = "GENERAL (non-repeat)" if args.general else "repeat-based"
    print(f"\n{'='*100}")
    print(f"  orthogonal consensus-caller comparison -- {label}  ({len(to_run)} scenarios"
          + (f"; {len(skipped_multi)} multi-allele scenarios not applicable, skipped: {', '.join(skipped_multi)}"
             if skipped_multi else "")
          + ")")
    print(f"{'='*100}")
    header = f"  {'SCENARIO':<32} {'MODEL':<9} {'TRUTH':<14} {'poa-consensus':<30} {'abPOA':<30} {'SPOA':<30}"
    print(header)
    print(f"  {'-'*32} {'-'*9} {'-'*13} {'-'*30} {'-'*30} {'-'*30}")

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
            print(f"  {scenario.name:<32} ERROR: {exc}")
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

        print(f"  {scenario.name:<32} {scenario.model:<9} {truth.display:<14} "
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
