#!/usr/bin/env python3
"""3-way caller comparison on a CONTROLLED synthetic grid (poa-consensus vs
abPOA vs SPOA) — the "know the edges" tool.

The robustness matrix (tests/single_allele_matrix.rs) shows poa-consensus is
self-consistent with ground truth. This asks the harder question: is it
*competitive*? Read sets + truth are emitted by the Rust generator
(tests/emit_compare_grid.rs) so all three callers see byte-identical input.
Each caller's consensus is scored by fitting edit distance to truth, and every
config is classified:

  POA-WORSE   poa's edit exceeds the best external's by a real margin  -> a
              poa-consensus flaw to fix before release
  FLOOR       all callers are similarly far from truth                 -> the
              problem's inherent difficulty, not our bug
  OK/BEST     poa ties or beats the externals

Setup (one-time):  pip install pyabpoa pyspoa
Emit the grid:     COMPARE_GRID_DIR=/tmp/poa_compare_grid \
                     cargo test --release --features cli --test emit_compare_grid -- --ignored
Run:               python3 bench/compare_grid.py [--dir /tmp/poa_compare_grid]
"""

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from compare_callers import fitting_edit_distance, run_abpoa, run_spoa  # noqa: E402

POA_BIN = Path(__file__).parent.parent / "target/release/poa-consensus"


def read_fasta(path: Path) -> list[bytes]:
    seqs, cur = [], []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur).encode())
                cur = []
        elif line.strip():
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur).encode())
    return seqs


def run_poa(reads_fa: Path) -> bytes | None:
    """Shipping single-allele CLI path (default band 50 + adaptive, semi-global)."""
    proc = subprocess.run([str(POA_BIN), str(reads_fa), "--quiet"], capture_output=True)
    if proc.returncode != 0:
        return None
    seqs = []
    cur = []
    for line in proc.stdout.decode(errors="replace").splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        elif line.strip():
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs[0].encode() if seqs else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="/tmp/poa_compare_grid")
    ap.add_argument("--margin", type=float, default=0.03,
                    help="POA-WORSE if poa_edit - best_external >= max(3, margin*truth_len)")
    args = ap.parse_args()

    d = Path(args.dir)
    manifest = d / "manifest.txt"
    if not manifest.exists():
        print(f"no manifest at {manifest}; emit the grid first (see module docstring)", file=sys.stderr)
        return 2
    if not POA_BIN.exists():
        print(f"build the CLI first: cargo build --release --features cli", file=sys.stderr)
        return 2

    names = [n for n in manifest.read_text().splitlines() if n.strip()]
    have_abpoa = run_abpoa([b"ACGT", b"ACGT"]) is not None
    have_spoa = run_spoa([b"ACGT", b"ACGT"]) is not None
    if not (have_abpoa or have_spoa):
        print("neither pyabpoa nor pyspoa importable — pip install pyabpoa pyspoa", file=sys.stderr)
        return 2

    print(f"{'CONFIG':<28} {'truth':>5} {'poa':>5} {'abPOA':>6} {'SPOA':>5}   verdict")
    print("-" * 70)
    n_worse = n_floor = n_ok = 0
    worst = []
    for name in names:
        reads_fa = d / f"{name}.reads.fa"
        truth = (d / f"{name}.truth").read_bytes()
        seqs = read_fasta(reads_fa)
        tl = len(truth)

        poa = run_poa(reads_fa)
        ab = run_abpoa(seqs) if have_abpoa else None
        sp = run_spoa(seqs) if have_spoa else None

        e_poa = fitting_edit_distance(truth, poa) if poa else 10**9
        e_ab = fitting_edit_distance(truth, ab) if ab else None
        e_sp = fitting_edit_distance(truth, sp) if sp else None
        ext = [e for e in (e_ab, e_sp) if e is not None]
        best_ext = min(ext) if ext else None

        thr = max(3, int(args.margin * tl))
        if best_ext is not None and e_poa - best_ext >= thr:
            verdict = "POA-WORSE"
            n_worse += 1
            worst.append((e_poa - best_ext, name, e_poa, e_ab, e_sp, tl))
        elif best_ext is not None and best_ext >= thr and e_poa >= thr:
            verdict = "floor"
            n_floor += 1
        else:
            verdict = "ok"
            n_ok += 1

        fmt = lambda x: "—" if x is None else str(x)
        print(f"{name:<28} {tl:>5} {e_poa:>5} {fmt(e_ab):>6} {fmt(e_sp):>5}   {verdict}")

    print("-" * 70)
    print(f"{len(names)} configs:  ok/best {n_ok}   floor (all hard) {n_floor}   POA-WORSE {n_worse}")
    if worst:
        print("\nWorst POA-WORSE cases (poa deficit vs best external):")
        for deficit, name, ep, ea, es, tl in sorted(worst, reverse=True)[:15]:
            print(f"  +{deficit:<4} {name:<28} poa={ep} abPOA={ea} SPOA={es} (truth {tl})")
    return 1 if n_worse else 0


if __name__ == "__main__":
    raise SystemExit(main())
