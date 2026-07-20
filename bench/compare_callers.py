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
import subprocess
import sys
import time
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


# ── Comprehensive matrix ─────────────────────────────────────────────────────
#
# A broad, systematic cross-tool matrix (poa-consensus vs abPOA vs SPOA) that
# deliberately spans sequence-context x depth x error-model x contamination,
# rather than the narrow hand-picked default catalogue. The point is breadth:
# where does poa-consensus stand *generally*, including easy baselines every
# tool should nail, and where are the real gaps.
#
# Representation conventions (all reuse validate.py's pbsim3 -> minimap2 ->
# bedpull pipeline unchanged; only truth construction differs):
#   * A pure tandem repeat is `Allele(unit, count, depth)` (no `desc`), so the
#     unit-count display and `extract_allele_by_unit` reporting work.
#   * Any non-repeat / interrupted / structural truth is a raw sequence carried
#     as `Allele(unit=<raw seq>, count=1, depth, desc=<label>)`; `desc` set means
#     evaluate_caller scores it purely by `fitting_edit_distance` (no unit count).
#   * Contamination and genuine-multi-allele scenarios use `multi=False` with the
#     TRUTH allele first (highest depth) and one or more extra contigs after it.
#     build_reference writes one contig per allele, simulate_reads simulates each
#     at its own depth, and extract_reads pools reads from every contig into one
#     read set -- exactly the established `sv_cag20_out60` shape. All three
#     callers run single-consensus and are scored against `alleles[0]` (the
#     dominant/expected allele). For genuine two-allele cases this measures
#     "recover the dominant allele amid a real second allele"; noted as a caveat.

def _hp(rng, left_n: int, base: str, run: int, right_n: int) -> str:
    """Random flanks around a homopolymer run of `base` x `run`."""
    return _rng_seq(rng, left_n) + base * run + _rng_seq(rng, right_n)


def _build_comprehensive_scenarios() -> list["v.Scenario"]:
    import random
    S: list["v.Scenario"] = []
    A = v.Allele
    Sc = v.Scenario

    # ── A. Tandem repeats: motif size x length x model (clean, depth 20) ──────
    # 3-mer CAG across the full length ladder.
    for c in (5, 10, 20, 50, 100, 200):
        S.append(Sc(f"cmp_cag{c}_r10", [A("CAG", c, 20)], "ont_r10"))
    # 3-mer GAA across the ladder.
    for c in (5, 10, 20, 50, 100, 200):
        S.append(Sc(f"cmp_gaa{c}_r10", [A("GAA", c, 20)], "ont_r10"))
    # Model spread on two representative repeat lengths.
    S.append(Sc("cmp_cag20_r9", [A("CAG", 20, 20)], "ont_r9"))
    S.append(Sc("cmp_cag20_hifi", [A("CAG", 20, 20)], "hifi"))
    S.append(Sc("cmp_gaa50_r9", [A("GAA", 50, 20)], "ont_r9"))
    S.append(Sc("cmp_gaa50_hifi", [A("GAA", 50, 20)], "hifi"))
    # 2-mer AT microsatellite (dinucleotide period).
    for c in (10, 20, 40):
        S.append(Sc(f"cmp_at{c}_r10", [A("AT", c, 20)], "ont_r10"))
    # 4-mer motif.
    for c in (10, 25, 50):
        S.append(Sc(f"cmp_gcat{c}_r10", [A("GCAT", c, 20)], "ont_r10"))
    # 5-mer RFC1-like repeats (Known Bug #4: dense diagonal lattice, hard).
    for c in (10, 30, 60):
        S.append(Sc(f"cmp_aaaag{c}_r10", [A("AAAAG", c, 20)], "ont_r10"))
    S.append(Sc("cmp_ctttt30_r10", [A("CTTTT", 30, 20)], "ont_r10"))
    # 6-mer motif (longer period than this crate's usual 3-5bp focus).
    for c in (10, 30):
        S.append(Sc(f"cmp_motif6_{c}_r10", [A("GCTAGA", c, 20)], "ont_r10"))

    # ── B. Depth stress (low / high) on repeats and random ───────────────────
    S.append(Sc("cmp_cag20_d5_r10", [A("CAG", 20, 5)], "ont_r10"))
    S.append(Sc("cmp_cag20_d10_r10", [A("CAG", 20, 10)], "ont_r10"))
    S.append(Sc("cmp_cag20_d30_r10", [A("CAG", 20, 30)], "ont_r10"))
    S.append(Sc("cmp_cag50_d5_r10", [A("CAG", 50, 5)], "ont_r10"))
    S.append(Sc("cmp_cag50_d30_r10", [A("CAG", 50, 30)], "ont_r10"))

    # ── C. Non-repeat sequence contexts ──────────────────────────────────────
    rng = random.Random(20270101)
    S.append(Sc("cmp_rand150_r10", [A(_rng_seq(rng, 150), 1, 20, desc="rand150")], "ont_r10"))
    S.append(Sc("cmp_rand300_r10", [A(_rng_seq(rng, 300), 1, 20, desc="rand300")], "ont_r10"))
    S.append(Sc("cmp_rand600_r10", [A(_rng_seq(rng, 600), 1, 20, desc="rand600")], "ont_r10"))
    rng = random.Random(20270102)
    _r300 = _rng_seq(rng, 300)
    S.append(Sc("cmp_rand300_hifi", [A(_r300, 1, 20, desc="rand300")], "hifi"))
    S.append(Sc("cmp_rand300_r9", [A(_r300, 1, 20, desc="rand300")], "ont_r9"))
    # Depth stress on random.
    rng = random.Random(20270103)
    _rd = _rng_seq(rng, 300)
    S.append(Sc("cmp_rand300_d5_r10", [A(_rd, 1, 5, desc="rand300_d5")], "ont_r10"))
    S.append(Sc("cmp_rand300_d30_r10", [A(_rd, 1, 30, desc="rand300_d30")], "ont_r10"))
    # GC extremes.
    rng = random.Random(20270104)
    S.append(Sc("cmp_gc80_r10", [A(_rng_seq(rng, 300, gc=0.80), 1, 20, desc="gc80")], "ont_r10"))
    S.append(Sc("cmp_gc20_r10", [A(_rng_seq(rng, 300, gc=0.20), 1, 20, desc="gc20")], "ont_r10"))
    S.append(Sc("cmp_gc80_hifi", [A(_rng_seq(rng, 300, gc=0.80), 1, 20, desc="gc80")], "hifi"))
    # Embedded homopolymer runs of varying length / base / model.
    rng = random.Random(20270105)
    S.append(Sc("cmp_homopolymer_a20_r10", [A(_hp(rng, 140, "A", 20, 140), 1, 20, desc="hpA20")], "ont_r10"))
    S.append(Sc("cmp_homopolymer_a30_r10", [A(_hp(rng, 140, "A", 30, 140), 1, 20, desc="hpA30")], "ont_r10"))
    S.append(Sc("cmp_homopolymer_t15_hifi", [A(_hp(rng, 140, "T", 15, 140), 1, 20, desc="hpT15")], "hifi"))
    # Interrupted repeats: a clean CAG run broken once mid-run by a SNP or indel.
    _cag = "CAG"
    S.append(Sc("cmp_interrupted_cag_snp_r10",
               [A(_cag * 10 + "CAA" + _cag * 10, 1, 20, desc="cagSNPint")], "ont_r10"))
    S.append(Sc("cmp_interrupted_cag_indel_r10",
               [A(_cag * 10 + "CA" + _cag * 10, 1, 20, desc="cagINDELint")], "ont_r10"))

    # ── D. Structural (single-consensus reconstruction) ──────────────────────
    rng = random.Random(20270201)
    S.append(Sc("cmp_struct_insert250_r10",
               [A(_rng_seq(rng, 75) + _rng_seq(rng, 250) + _rng_seq(rng, 75), 1, 20, desc="ins250")], "ont_r10"))
    rng = random.Random(20270202)
    S.append(Sc("cmp_struct_insert100_hifi",
               [A(_rng_seq(rng, 100) + _rng_seq(rng, 100) + _rng_seq(rng, 100), 1, 20, desc="ins100")], "hifi"))
    # Short-range tandem duplication (a ~20bp segment repeated once, adjacent).
    rng = random.Random(20270203)
    _pre, _du, _post = _rng_seq(rng, 100), _rng_seq(rng, 20), _rng_seq(rng, 160)
    S.append(Sc("cmp_tandem_dup_r10", [A(_pre + _du + _du + _post, 1, 20, desc="dup20x2")], "ont_r10"))
    rng = random.Random(20270204)
    _pre2, _du2, _post2 = _rng_seq(rng, 100), _rng_seq(rng, 30), _rng_seq(rng, 160)
    S.append(Sc("cmp_tandem_dup_hifi", [A(_pre2 + _du2 + _du2 + _post2, 1, 20, desc="dup30x2")], "hifi"))

    # ── E. Contamination (multi=False; truth=alleles[0], dominant) ───────────
    # Pure unrelated-sequence noise: contaminant contig is an unrelated random seq.
    rng = random.Random(20270301)
    _truth_u = _rng_seq(rng, 300)
    _noise_a = _rng_seq(rng, 300)
    _noise_b = _rng_seq(rng, 300)
    S.append(Sc("cmp_contam_unrelated_10pct_r10",
               [A(_truth_u, 1, 18, desc="rand300"), A(_noise_a, 1, 2, desc="noise")],
               "ont_r10"))
    S.append(Sc("cmp_contam_unrelated_25pct_r10",
               [A(_truth_u, 1, 15, desc="rand300"), A(_noise_b, 1, 5, desc="noise")],
               "ont_r10"))
    # Near-miss same-locus-wrong-length: contaminant is the same repeat, off by N units.
    S.append(Sc("cmp_contam_nearmiss_cag_10pct_r10",
               [A("CAG", 20, 18), A("CAG", 40, 2)], "ont_r10"))
    S.append(Sc("cmp_contam_nearmiss_cag_25pct_r10",
               [A("CAG", 20, 15), A("CAG", 30, 5)], "ont_r10"))
    S.append(Sc("cmp_contam_nearmiss_gaa_25pct_r10",
               [A("GAA", 50, 15), A("GAA", 70, 5)], "ont_r10"))
    # Heterogeneous one-offs: several distinct single-read contaminants.
    rng = random.Random(20270302)
    _truth_h = _rng_seq(rng, 300)
    _hetero = [A(_rng_seq(rng, 300), 1, 1, desc=f"oneoff{i}") for i in range(4)]
    S.append(Sc("cmp_contam_hetero_oneoffs_r10",
               [A(_truth_h, 1, 18, desc="rand300")] + _hetero, "ont_r10"))
    # Large deletion as contamination: minority carries a shortened (deleted) form.
    rng = random.Random(20270303)
    _full = _rng_seq(rng, 120) + _rng_seq(rng, 120) + _rng_seq(rng, 120)
    _deleted = _full[:120] + _full[240:]  # middle 120bp deleted
    S.append(Sc("cmp_contam_deletion_25pct_r10",
               [A(_full, 1, 15, desc="full360"), A(_deleted, 1, 5, desc="del120")], "ont_r10"))

    # ── F. Partial reads (non-spanning) ──────────────────────────────────────
    rng = random.Random(20270401)
    S.append(Sc("cmp_partial_rand500_r10",
               [A(_rng_seq(rng, 500), 1, 20, desc="partial500")],
               "ont_r10", length_mean=350, flank_len=150, partial=True))
    # Partial reads over a long repeat: no single read spans CAG100.
    S.append(Sc("cmp_partial_cag100_r10",
               [A("CAG", 100, 20)], "ont_r10",
               length_mean=250, flank_len=150, partial=True))

    # ── G. Genuine multi-allele (scored vs dominant allele; caveat) ──────────
    # Single-consensus callers can only emit one sequence; scored against the
    # higher-depth allele (alleles[0]). Not a phasing test -- a "recover the
    # dominant allele amid a real second allele" test.
    S.append(Sc("cmp_multiallele_difflen_cag_r10",
               [A("CAG", 20, 16), A("CAG", 40, 8)], "ont_r10"))
    S.append(Sc("cmp_multiallele_difflen_gaa_r10",
               [A("GAA", 30, 16), A("GAA", 60, 8)], "ont_r10"))
    rng = random.Random(20270501)
    _base = _rng_seq(rng, 150)
    _snp = list(_base)
    _snp[75] = "A" if _base[75] != "A" else "C"  # single-base difference, same length
    S.append(Sc("cmp_multiallele_samelen_snp_r10",
               [A(_base, 1, 16, desc="hapRef"), A("".join(_snp), 1, 8, desc="hapSNP")],
               "ont_r10"))

    return S


COMPREHENSIVE_SCENARIOS = _build_comprehensive_scenarios()


# ── Longer-sequence (multi-kb scale) tier ────────────────────────────────────
#
# Probes the multi-kb regime CLAUDE.md targets (up to ~20kb), where
# poa-consensus's adaptive-band machinery is meant to matter vs SPOA (fully
# unbanded: ~9.6 GB per alignment at 20kb per CLAUDE.md's scale table) and
# abPOA (adaptive/static band). The rest of the matrix tops out ~600bp; this
# scales random sequence and tandem repeats up to ~20kb.
#
# A caller failing at large size is an INFORMATIVE result, not a run-breaker,
# so each caller here runs ISOLATED in its own subprocess under a wall-clock
# timeout and a virtual-memory cap (ulimit -v). One tool OOMing/timing out on a
# huge scenario cannot crash the batch or lose the other tools' results, and we
# record per scenario which completed, which failed (timeout/oom/error),
# runtime, and peak RSS (via /usr/bin/time -v). Caps are printed at run time.

LONG_TIMEOUT_S = 600          # per-caller wall-clock cap
LONG_MEM_KB = 10 * 1024 * 1024  # per-caller virtual-memory cap (10 GB; box has 31 GB)

_CHILD = r'''
import sys
fa, caller = sys.argv[1], sys.argv[2]
seqs, buf = [], []
with open(fa) as fh:
    for line in fh:
        if line.startswith(">"):
            if buf:
                seqs.append("".join(buf)); buf = []
        else:
            buf.append(line.strip())
if buf:
    seqs.append("".join(buf))
if not seqs:
    sys.exit(3)
if caller == "abPOA":
    import pyabpoa as pa
    al = pa.msa_aligner(aln_mode="e")
    r = al.msa(seqs, out_cons=True, out_msa=False)
    cons = r.cons_seq[0] if (r and r.cons_seq) else ""
elif caller == "SPOA":
    import spoa
    cons, _ = spoa.poa(seqs, algorithm=2)
else:
    sys.exit(4)
sys.stdout.write(cons)
'''


def _parse_time_v(path: Path):
    """Parse peak RSS (kbytes) and elapsed wall-clock (s) from /usr/bin/time -v."""
    peak_kb = None
    elapsed_s = None
    if not path.exists():
        return peak_kb, elapsed_s
    for line in path.read_text().splitlines():
        if "Maximum resident set size" in line:
            try:
                peak_kb = int(line.rsplit(":", 1)[-1].strip())
            except ValueError:
                pass
        elif "Elapsed (wall clock)" in line:
            val = line.rsplit(")", 1)[-1].lstrip(": ").strip()
            try:
                secs = 0.0
                for part in val.split(":"):
                    secs = secs * 60 + float(part)
                elapsed_s = secs
            except ValueError:
                pass
    return peak_kb, elapsed_s


def _run_isolated(argv: list, work: Path, tag: str):
    """Run `argv` under ulimit -v + timeout + /usr/bin/time -v.

    Returns (stdout_bytes | None, status, elapsed_s | None, peak_kb | None).
    status is one of: "ok", "timeout", "oom", "error(rc=N)".
    """
    import shlex
    timefile = work / f"time_{tag}.txt"
    inner = " ".join(shlex.quote(str(a)) for a in argv)
    bash = (f"ulimit -v {LONG_MEM_KB}; "
            f"exec /usr/bin/time -v -o {shlex.quote(str(timefile))} "
            f"timeout {LONG_TIMEOUT_S} {inner}")
    proc = subprocess.run(["bash", "-c", bash], capture_output=True)
    rc = proc.returncode
    peak_kb, elapsed_s = _parse_time_v(timefile)
    near_cap = peak_kb is not None and peak_kb >= 0.92 * LONG_MEM_KB
    if rc == 0:
        return proc.stdout, "ok", elapsed_s, peak_kb
    if rc == 124:
        return None, "timeout", elapsed_s, peak_kb
    # 137=SIGKILL(OOM killer), 134=SIGABRT(bad_alloc), 139=SIGSEGV(failed malloc)
    if rc in (134, 137, 139) or near_cap:
        return None, "oom", elapsed_s, peak_kb
    return None, f"error(rc={rc})", elapsed_s, peak_kb


def _mem_str(peak_kb):
    if peak_kb is None:
        return "?"
    return f"{peak_kb / 1024:.0f}MB"


def _rt_str(elapsed_s):
    if elapsed_s is None:
        return "?"
    return f"{elapsed_s:.1f}s"


def _build_longtier_scenarios() -> list["v.Scenario"]:
    import random
    S: list["v.Scenario"] = []
    A = v.Allele
    Sc = v.Scenario

    def _long(name, seq, depth, model, desc):
        # Small flanks + a read length that spans the whole region + anchor, so
        # extracted reads cover the long locus (the scale test wants spanning
        # reads, not partials). length_sd stays validate's mean//5.
        region = len(seq)
        return Sc(name, [A(seq, 1, depth, desc=desc)], model,
                  length_mean=region + 450, flank_len=300)

    def _long_repeat(name, unit, count, depth, model):
        seq_len = len(unit) * count
        return Sc(name, [A(unit, count, depth)], model,
                  length_mean=seq_len + 450, flank_len=300)

    # ── Longer non-repeat random: 1k / 2k / 5k / 10k / 20k ────────────────────
    rng = random.Random(20280001)
    S.append(_long("lng_rand1k_r10", _rng_seq(rng, 1000), 12, "ont_r10", "rand1k"))
    S.append(_long("lng_rand2k_r10", _rng_seq(rng, 2000), 12, "ont_r10", "rand2k"))
    S.append(_long("lng_rand5k_r10", _rng_seq(rng, 5000), 10, "ont_r10", "rand5k"))
    # Fair 5kb comparison: the original lng_rand5k_r10 above starved to ~2
    # spanning reads (read-length distribution vs a tight 5.6kb reference).
    # This variant raises depth and gives read-length headroom (bigger flanks =
    # longer reference, higher length_mean) so bedpull's spanning-only
    # extraction yields an adequate pool (target >= ~10-15 reads). Verified
    # empirically before trusting the result (see the run report).
    rng5 = random.Random(20280011)
    S.append(v.Scenario("lng_rand5k_deep_r10",
                        [A(_rng_seq(rng5, 5000), 1, 60, desc="rand5k")],
                        "ont_r10", length_mean=6400, flank_len=800))
    S.append(_long("lng_rand10k_r10", _rng_seq(rng, 10000), 10, "ont_r10", "rand10k"))
    S.append(_long("lng_rand20k_r10", _rng_seq(rng, 20000), 10, "ont_r10", "rand20k"))
    # One-two HiFi at the smaller long sizes.
    rng = random.Random(20280002)
    S.append(_long("lng_rand1k_hifi", _rng_seq(rng, 1000), 12, "hifi", "rand1k"))
    S.append(_long("lng_rand2k_hifi", _rng_seq(rng, 2000), 12, "hifi", "rand2k"))

    # ── Longer tandem repeats beyond the 200-unit ceiling ─────────────────────
    S.append(_long_repeat("lng_cag300_r10", "CAG", 300, 15, "ont_r10"))   # 900bp
    S.append(_long_repeat("lng_cag500_r10", "CAG", 500, 12, "ont_r10"))   # 1500bp
    S.append(_long_repeat("lng_gaa300_r10", "GAA", 300, 15, "ont_r10"))   # 900bp
    S.append(_long_repeat("lng_gaa500_r10", "GAA", 500, 12, "ont_r10"))   # 1500bp
    # Longer motif forming a multi-kb locus.
    S.append(_long_repeat("lng_motif6_200_r10", "GCTAGA", 200, 12, "ont_r10"))  # 1200bp

    # ── A multi-kb structural case (large insertion within a long region) ──────
    rng = random.Random(20280003)
    struct = _rng_seq(rng, 1000) + _rng_seq(rng, 1000) + _rng_seq(rng, 1000)  # 3kb, mid = insert
    S.append(_long("lng_struct_insert_multikb_r10", struct, 12, "ont_r10", "3kb+ins"))

    return S


LONGTIER_SCENARIOS = _build_longtier_scenarios()


def run_longtier(args):
    """Dedicated runner for the multi-kb scale tier: isolated, resource-capped
    caller invocations, plus a scale/runtime/memory table."""
    to_run = LONGTIER_SCENARIOS
    if args.scenarios:
        requested = set(args.scenarios)
        to_run = [s for s in LONGTIER_SCENARIOS if s.name in requested]
        unknown = requested - {s.name for s in LONGTIER_SCENARIOS}
        if unknown:
            print(f"warning: unknown longtier scenarios skipped: {', '.join(sorted(unknown))}",
                  file=sys.stderr)

    # Distinct workdir so a concurrent --comprehensive run never clobbers this.
    workdir = Path(args.workdir if args.workdir != "bench/work_compare" else "bench/work_longtier")
    workdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*110}")
    print(f"  orthogonal consensus-caller comparison -- LONGER-SEQUENCE (multi-kb) tier  "
          f"({len(to_run)} scenarios)")
    print(f"  per-caller caps: timeout={LONG_TIMEOUT_S}s, virtual-memory={LONG_MEM_KB // (1024*1024)}GB "
          f"(box has 31GB); a caller failing at scale is informative, not a run-breaker")
    print(f"{'='*110}")
    header = f"  {'SCENARIO':<32} {'MODEL':<8} {'TRUTH':<10} {'poa-consensus':<28} {'abPOA':<28} {'SPOA':<28}"
    print(header)
    print(f"  {'-'*32} {'-'*8} {'-'*10} {'-'*28} {'-'*28} {'-'*28}")

    tally = {"all_agree_ok": 0, "poa_better": 0, "poa_worse": 0, "all_fail": 0, "error": 0}
    cat = {"poa_better": [], "poa_worse": [], "all_fail": [], "error": []}
    scale_rows = []   # (name, region_bp, n_reads, per-caller (status,rt,peak))
    t_start = time.time()

    for scenario in to_run:
        work = workdir / scenario.name
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)
        model = v.ERRMODELS[scenario.model]
        truth = scenario.alleles[0]
        region_bp = len(truth.seq)

        try:
            ref = v.build_reference(scenario, work)
            reads = v.simulate_reads(scenario, ref, work, model)
            bam = v.align_and_index(reads, ref, work, model["mm2_preset"])
            extracted_path = v.extract_reads(scenario, bam, work)
            records = v.parse_fasta(extracted_path)
            n_reads = len(records)
        except Exception as exc:
            print(f"  {scenario.name:<32} PIPELINE ERROR: {exc}")
            tally["error"] += 1
            cat["error"].append((scenario.name, f"pipeline: {exc}"))
            if not args.keep:
                shutil.rmtree(work, ignore_errors=True)
            continue

        if n_reads < 2:
            print(f"  {scenario.name:<32} SKIP: only {n_reads} spanning read(s) extracted")
            tally["error"] += 1
            cat["error"].append((scenario.name, f"only {n_reads} spanning reads extracted"))
            if not args.keep:
                shutil.rmtree(work, ignore_errors=True)
            continue

        # poa-consensus: isolated CLI invocation (band-width 50 + adaptive default,
        # matching validate.run_consensus for single-allele).
        poa_argv = [str(v.POA_BIN), "--band-width", "50", "--min-reads", "2", str(extracted_path)]
        poa_out, poa_status, poa_rt, poa_peak = _run_isolated(poa_argv, work, "poa")
        poa_cons = None
        if poa_status == "ok" and poa_out:
            poa_records = v.parse_fasta_bytes(poa_out) if hasattr(v, "parse_fasta_bytes") else None
            if poa_records is None:
                # parse FASTA from stdout bytes inline
                seqbuf, cur = [], []
                for line in poa_out.split(b"\n"):
                    if line.startswith(b">"):
                        if cur:
                            seqbuf.append(b"".join(cur)); cur = []
                    else:
                        cur.append(line.strip())
                if cur:
                    seqbuf.append(b"".join(cur))
                poa_cons = seqbuf[0] if seqbuf else None
            else:
                poa_cons = poa_records[0][1] if poa_records else None

        ab_out, ab_status, ab_rt, ab_peak = _run_isolated(
            [sys.executable, "-c", _CHILD, str(extracted_path), "abPOA"], work, "abpoa")
        ab_cons = ab_out.strip() if (ab_status == "ok" and ab_out) else None

        sp_out, sp_status, sp_rt, sp_peak = _run_isolated(
            [sys.executable, "-c", _CHILD, str(extracted_path), "SPOA"], work, "spoa")
        sp_cons = sp_out.strip() if (sp_status == "ok" and sp_out) else None

        results = [
            evaluate_caller("poa-consensus", poa_cons, truth),
            evaluate_caller("abPOA", ab_cons, truth),
            evaluate_caller("SPOA", sp_cons, truth),
        ]
        # Overlay the isolation status when a caller did not complete, so the
        # cell shows *why* (timeout/oom) rather than a bare N/A.
        statuses = [poa_status, ab_status, sp_status]
        for r, st in zip(results, statuses):
            if st != "ok":
                r["status"] = "N/A"
                r["iso"] = st

        def _cell(r):
            if r.get("iso"):
                return f"{r['iso'].upper()}"
            return cell(r)

        scale_rows.append((scenario.name, region_bp, n_reads,
                           (poa_status, poa_rt, poa_peak),
                           (ab_status, ab_rt, ab_peak),
                           (sp_status, sp_rt, sp_peak)))

        # Tally: a caller that timed out / OOMed counts as a FAIL for it.
        def _ok(r, st):
            return st == "ok" and r["status"] == "OK"
        poa_ok = _ok(results[0], poa_status)
        ext_ok = [_ok(results[1], ab_status), _ok(results[2], sp_status)]
        if poa_ok and all(ext_ok):
            tally["all_agree_ok"] += 1
        elif poa_ok and not all(ext_ok):
            tally["poa_better"] += 1
            cat["poa_better"].append((scenario, truth, results))
        elif not poa_ok and any(ext_ok):
            tally["poa_worse"] += 1
            cat["poa_worse"].append((scenario, truth, results))
        else:
            tally["all_fail"] += 1
            cat["all_fail"].append((scenario, truth, results))

        print(f"  {scenario.name:<32} {scenario.model:<8} {truth.display:<10} "
              f"{_cell(results[0]):<28} {_cell(results[1]):<28} {_cell(results[2]):<28}")

        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)

    # ── Summary ──────────────────────────────────────────────────────────────
    n = len(scale_rows)
    print(f"\n  Result: {n} scenarios compared"
          + (f", {tally['error']} errored/skipped" if tally["error"] else ""))
    print(f"    all three agree (OK):                {tally['all_agree_ok']}")
    print(f"    poa-consensus OK, external(s) FAIL:  {tally['poa_better']}   <- banding advantage at scale")
    print(f"    poa-consensus FAIL, external(s) OK:  {tally['poa_worse']}   <- chase these")
    print(f"    all three FAIL:                      {tally['all_fail']}   <- shared scale ceiling")

    # ── Scale / runtime / memory table ───────────────────────────────────────
    print(f"\n  SCALE / RUNTIME / PEAK-MEMORY table (per caller):")
    print(f"  {'SCENARIO':<32} {'REGION':>7} {'READS':>5}  "
          f"{'poa: st/rt/mem':<26} {'abPOA: st/rt/mem':<26} {'SPOA: st/rt/mem':<26}")
    print(f"  {'-'*32} {'-'*7} {'-'*5}  {'-'*26} {'-'*26} {'-'*26}")
    for name, region_bp, n_reads, poa_t, ab_t, sp_t in scale_rows:
        def _fmt(t):
            st, rt, pk = t
            return f"{st}/{_rt_str(rt)}/{_mem_str(pk)}"
        print(f"  {name:<32} {region_bp:>6}b {n_reads:>5}  "
              f"{_fmt(poa_t):<26} {_fmt(ab_t):<26} {_fmt(sp_t):<26}")

    def _line(scenario, truth, results):
        def _c(r):
            return r.get("iso", "").upper() or cell(r)
        return (f"    - {scenario.name:<34} ({scenario.model}, {truth.display}): "
                f"poa[{_c(results[0])}]  abPOA[{_c(results[1])}]  SPOA[{_c(results[2])}]")

    if cat["poa_better"]:
        print(f"\n  WINS (poa-consensus OK where a competitor failed) -- banding-at-scale advantage:")
        for s, t, r in cat["poa_better"]:
            print(_line(s, t, r))
    if cat["poa_worse"]:
        print(f"\n  LOSSES (poa-consensus FAIL while a competitor OK) -- actionable:")
        for s, t, r in cat["poa_worse"]:
            print(_line(s, t, r))
    if cat["all_fail"]:
        print(f"\n  SHARED SCALE CEILINGS (all three FAIL):")
        for s, t, r in cat["all_fail"]:
            print(_line(s, t, r))
    if cat["error"]:
        print(f"\n  ERRORED / SKIPPED (nothing silently dropped):")
        for name, exc in cat["error"]:
            print(f"    - {name}: {exc}")

    elapsed = time.time() - t_start
    print(f"\n  Total runtime: {elapsed:.1f}s ({elapsed/60:.1f} min) over "
          f"{len(to_run)} scenarios ({tally['error']} errored/skipped)")
    print()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", nargs="*", help="Scenario names to run (default: all single-allele scenarios)")
    ap.add_argument("--workdir", default="bench/work_compare", help="Working directory for temp files")
    ap.add_argument("--keep", action="store_true", help="Keep working files after run")
    ap.add_argument("--general", action="store_true",
                    help="Run the general (non-repeat) scenario battery instead of the repeat-based catalogue")
    ap.add_argument("--comprehensive", action="store_true",
                    help="Run the broad comprehensive matrix (60+ scenarios spanning "
                         "sequence-context x depth x error-model x contamination)")
    ap.add_argument("--longtier", action="store_true",
                    help="Run the longer-sequence (multi-kb, up to ~20kb) scale tier with "
                         "isolated, timeout- and memory-capped caller invocations")
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

    if args.longtier:
        run_longtier(args)
        return

    if args.comprehensive:
        all_scenarios = COMPREHENSIVE_SCENARIOS
        skipped_multi = []
    elif args.general:
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

    label = ("COMPREHENSIVE matrix" if args.comprehensive
             else "GENERAL (non-repeat)" if args.general
             else "repeat-based")
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
    cat = {"poa_better": [], "poa_worse": [], "all_fail": [], "error": []}
    rows = []
    t_start = time.time()

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
            cat["error"].append((scenario.name, str(exc)))
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
            cat["poa_better"].append((scenario, truth, results))
        elif not poa_ok and any(ext_ok):
            tally["poa_worse"] += 1
            cat["poa_worse"].append((scenario, truth, results))
        else:
            tally["all_fail"] += 1
            cat["all_fail"].append((scenario, truth, results))

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

    def _line(scenario, truth, results) -> str:
        return (f"    - {scenario.name:<34} ({scenario.model}, {truth.display}): "
                f"poa[{cell(results[0])}]  abPOA[{cell(results[1])}]  SPOA[{cell(results[2])}]")

    if cat["poa_worse"]:
        print(f"\n  LOSSES (poa-consensus FAIL while abPOA and/or SPOA OK) -- actionable follow-ups:")
        for scenario, truth, results in cat["poa_worse"]:
            print(_line(scenario, truth, results))
    if cat["poa_better"]:
        print(f"\n  WINS (poa-consensus alone OK) -- confirmed differentiators:")
        for scenario, truth, results in cat["poa_better"]:
            print(_line(scenario, truth, results))
    if cat["all_fail"]:
        print(f"\n  SHARED CEILINGS (all three FAIL):")
        for scenario, truth, results in cat["all_fail"]:
            print(_line(scenario, truth, results))
    if cat["error"]:
        print(f"\n  ERRORED / SKIPPED (nothing silently dropped):")
        for name, exc in cat["error"]:
            print(f"    - {name}: {exc}")

    elapsed = time.time() - t_start
    print(f"\n  Total runtime: {elapsed:.1f}s ({elapsed/60:.1f} min) over "
          f"{len(to_run)} scenarios ({tally['error']} errored)")
    print()


if __name__ == "__main__":
    main()
