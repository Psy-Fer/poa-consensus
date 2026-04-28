# poa-consensus

A pure-Rust banded Partial Order Alignment (POA) library for building consensus sequences from a set of reads.

POA aligns reads into a directed acyclic graph (DAG) using affine-gap dynamic programming (separate scoring for gap open and gap extension) and extracts consensus by following the heaviest (most-supported) path. Insertions and deletions create branches resolved by read support.

**Target use cases:** short tandem repeat (STR) loci, amplicon consensus, structural variants, per-locus nanopore or HiFi read sets (50 bp to ~20 kb with banded DP). Can be paired with [`bedpull`](https://github.com/Psy-Fer/bedpull) for building a consensus of extracted sequences.

> The nearest alternative on crates.io is [`poasta`](https://crates.io/crates/poasta) (Broad Institute, pure Rust, gap-affine A\* alignment), which excels at larger graphs such as bacterial genes and HLA loci. For STR reads where graphs stay small and throughput across many loci matters, banded DP is faster. Both crates are pure Rust with no C dependencies.

## Usage

```toml
[dependencies]
poa-consensus = "0.1"
```

### Functional API

```rust
use poa_consensus::{consensus, consensus_multi, PoaConfig};

let reads: Vec<&[u8]> = vec![
    b"CATCATCAT",
    b"CATCATCAT",
    b"CATCGTCAT",
    b"CATCATCAT",
];

// Single-allele consensus; seed_idx seeds the graph with reads[seed_idx].
let result = consensus(&reads, 0, &PoaConfig::default())?;
println!("{}", String::from_utf8_lossy(&result.sequence));

// Multi-allele: returns one Consensus per detected allele.
let alleles = consensus_multi(&reads, 0, &PoaConfig::default())?;
```

### Stateful API

```rust
use poa_consensus::{PoaGraph, PoaConfig};

let reads: &[&[u8]] = &[b"CATCATCAT", b"CATCATCAT", b"CATCGTCAT"];

let mut graph = PoaGraph::new(reads[0], PoaConfig::default())?;
for read in &reads[1..] {
    graph.add_read(read)?;
}
let consensus = graph.consensus()?;
let stats     = graph.stats();
println!("bubbles: {}", stats.bubble_count);
```

The stateful API lets you inspect graph state between reads and reuse pre-allocated buffers across calls which is important for high-throughput per-locus pipelines.

### Two-pass adaptive mode

```rust
use poa_consensus::{consensus_adaptive, PoaConfig};

// Pass 1 builds the graph and computes GraphStats.
// Pass 2 is chosen automatically: multi-allele split, noise tightening,
// or semi-global switch, depending on what the stats reveal.
let alleles = consensus_adaptive(&reads, 0, &PoaConfig::default())?;
```

### Orientation utilities

```rust
use poa_consensus::{auto_orient, reverse_complement};

// Orient all reads to match the strand of reads[seed_idx] before POA.
// Uses kmer matching
let oriented = auto_orient(&reads, seed_idx);
```

## CLI

```
cargo install poa-consensus --features cli
```

```
poa-consensus reads.fa          # FASTA or FASTQ, auto-detected
poa-consensus reads.fa --multi  # multi-allele mode
```

## Feature flags

| Flag   | Adds                                       |
|--------|--------------------------------------------|
| `cli`  | Binary target; pulls in `clap` + `noodles` |
| `plot` | SVG visualisation helpers via `kuva`       |

Default build: library only, zero external dependencies.

## LLM disclosure

An LLM was used in the creation of this library. I was building bladerunner, and STR detector/genotyper and found the limitations of other POA crates, so after prototyping one inside bladerunner, I spun it out into this crate. Claude was used to help speed that along.

## License

MIT - James Ferguson 2026
