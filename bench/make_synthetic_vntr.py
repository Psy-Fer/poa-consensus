#!/usr/bin/env python3
"""Generate a synthetic homogeneous-VNTR read set that reproduces the greedy
flank-fabrication failure mode WITHOUT any private patient data.

The failure needs three ingredients, all present in MUC1-type loci:
  1. a long, homogeneous tandem repeat (so repeat-unit nodes fold onto each
     other and node coverage stops reflecting biological position);
  2. a unique flank whose last bases coincide with the repeat unit's leading
     bases (so the flank->repeat boundary is positionally ambiguous -- this is
     what lets a phantom insert accrue folded coverage at the boundary);
  3. realistic ONT indel noise so deletions scatter across positions instead of
     consolidating.

Deterministic (fixed seed) so the emitted FASTA is a stable, committable
fixture. Uses only synthetic sequence -- safe to commit.
"""

import argparse
import random


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# A made-up GC-rich 60 bp unit in the spirit of a VNTR repeat (NOT a real MUC1
# motif). Its first 3 bases (GGC) are echoed at the end of the left flank so
# the boundary is ambiguous.
UNIT = "GGCTCAGGTGGCAGCGGAGCCTGGACCTGCAGTGGCAGTCCAGGCTCCACGTGGCAGTGG"
LEFT_FLANK = "ATGACCAGTATGGCTATGGCTGGCAAGGGTGGTAGGAGTATCAGAGTGGTGGC"  # ends ...GGTGGC, echoes unit head
RIGHT_FLANK = "AGCTGAGCCTGATGCAGAGCCTGAGGCCGAGGTGACATTGTGGACTGGAGGGG"


def mutate(seq: str, sub: float, ins: float, dele: float, rng: random.Random) -> str:
    out = []
    for b in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + ins:
            out.append(rng.choice("ACGT"))
            out.append(b)
            continue
        if r < dele + ins + sub:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
            continue
        out.append(b)
    return "".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--units", type=int, default=40, help="repeat-unit count (allele length)")
    ap.add_argument("--depth", type=int, default=45, help="number of reads")
    ap.add_argument("--sub", type=float, default=0.015)
    ap.add_argument("--ins", type=float, default=0.015)
    ap.add_argument("--dele", type=float, default=0.020)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--variant-rate", type=float, default=0.0,
                    help="per-unit-copy substitution rate applied ONCE to the true allele "
                    "(creates variant units = interior unique anchors, i.e. a VNTR rather "
                    "than a pure STR). 0 = identical units (pure STR, the CR-limited case).")
    ap.add_argument("--partial-fraction", type=float, default=0.0,
                    help="fraction of reads truncated to a random sub-span (non-spanning "
                    "partial reads) to stress partial-read handling.")
    ap.add_argument("--rc-fraction", type=float, default=0.75,
                    help="fraction of reads emitted reverse-complemented (mixed strand, like real pulled reads)")
    ap.add_argument("--truth", action="store_true", help="also print the true allele as a comment header")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    # Build the true allele. With --variant-rate, each unit COPY is independently
    # mutated once (fixed in the truth), so units differ -> interior unique k-mers
    # -> a VNTR (scaffold-tractable) instead of a pure STR.
    if args.variant_rate > 0:
        units = []
        for _ in range(args.units):
            u = list(UNIT)
            for i in range(len(u)):
                if rng.random() < args.variant_rate:
                    u[i] = rng.choice([b for b in "ACGT" if b != u[i]])
            units.append("".join(u))
        repeat = "".join(units)
    else:
        repeat = UNIT * args.units
    true_allele = LEFT_FLANK + repeat + RIGHT_FLANK
    if args.truth:
        print(f"; true_allele_len={len(true_allele)} units={args.units}")

    for i in range(args.depth):
        # small per-read unit-count jitter (+-1) as real reads have
        jitter = rng.choice([-1, 0, 0, 0, 1])
        u = args.units + jitter
        if args.variant_rate > 0:
            # keep the variant units; add/drop a copy at the boundary for jitter
            body = repeat
            if jitter == 1:
                body = repeat + units[-1]
            elif jitter == -1:
                body = repeat[:-len(units[-1])] if units else repeat
            seq = LEFT_FLANK + body + RIGHT_FLANK
        else:
            seq = LEFT_FLANK + UNIT * u + RIGHT_FLANK
        seq = mutate(seq, args.sub, args.ins, args.dele, rng)
        if args.partial_fraction and rng.random() < args.partial_fraction:
            # truncate to a random contiguous sub-span (non-spanning partial read)
            L = len(seq)
            a = rng.randint(0, L // 2)
            b = rng.randint(L // 2, L)
            seq = seq[a:b]
        if rng.random() < args.rc_fraction:
            seq = revcomp(seq)
        print(f">read{i}_u{u}")
        print(seq)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
