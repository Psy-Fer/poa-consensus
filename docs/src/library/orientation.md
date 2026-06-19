# Read Orientation

POA assumes all reads are on the same strand. Mixed-strand input produces a garbage graph
silently. Orientation must be handled before calling `consensus`.

## Auto-orientation

```rust
use poa_consensus::{auto_orient, consensus, SeedSelection, select_seed, PoaConfig};

// reads: Vec<Vec<u8>> from your BAM/FASTQ parser
let seed_idx = select_seed(
    &reads.iter().map(|r| r.as_slice()).collect::<Vec<_>>(),
    &SeedSelection::Auto,
)?;

let oriented = auto_orient(&reads, seed_idx);
// oriented: Vec<Cow<[u8]>>
// borrowed if already correct strand, owned (reversed) if flipped

let ref_slices: Vec<&[u8]> = oriented.iter().map(|r| r.as_ref()).collect();
let result = consensus(&ref_slices, seed_idx, &PoaConfig::default())?;
```

`auto_orient` returns `Cow<[u8]>` for each read: a borrowed slice if the read is already on
the correct strand, or an owned `Vec<u8>` (the reverse complement) if it was flipped. This
avoids copying reads that don't need it.

## Per-read strand detection

```rust
use poa_consensus::{orient_to_seed, Strand};

let strand = orient_to_seed(read, seed, 15);  // k = 15
match strand {
    Strand::Forward => { /* use read as-is */ }
    Strand::Reverse => { /* use reverse_complement(read) */ }
}
```

`orient_to_seed` builds a 15-mer frequency set from the seed and counts hits in both
orientations of the read. The orientation with more matches wins. Falls back to alignment
score comparison for ambiguous cases (low k-mer match counts, highly repetitive sequence).

## Reverse complement

```rust
use poa_consensus::reverse_complement;

let rc = reverse_complement(b"ACGTACGT");
assert_eq!(rc, b"ACGTACGT");
```

## When to apply orientation

Always orient reads before POA for reads extracted from a BAM that may contain both strands.
If your upstream pipeline already guarantees strand (e.g. haplotagged BAM with HP tag), you
can skip this step.

For amplicon reads from a symmetric PCR product (where reads from both strands represent the
same sequence), always orient. For strand-specific protocols, check whether your reads are
guaranteed to be on the sense strand.
