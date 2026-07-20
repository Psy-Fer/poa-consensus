# poa-consensus

[![CI](https://github.com/Psy-Fer/poa-consensus/actions/workflows/ci.yml/badge.svg)](https://github.com/Psy-Fer/poa-consensus/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/poa-consensus.svg)](https://crates.io/crates/poa-consensus)
[![docs.rs](https://docs.rs/poa-consensus/badge.svg)](https://docs.rs/poa-consensus)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.85%2B-orange.svg)](https://www.rust-lang.org)

A pure-Rust banded Partial Order Alignment (POA) library for building consensus sequences from a set of reads.

POA aligns reads into a directed acyclic graph (DAG) using affine-gap dynamic programming and extracts consensus by following the heaviest (most-supported) path. Insertions and deletions create branches resolved by read support.

**Target use cases:** short tandem repeat (STR) loci, amplicon consensus, structural variants, per-locus nanopore or HiFi read sets (50 bp to ~20 kb with banded DP). Can be paired with [`bedpull`](https://github.com/Psy-Fer/bedpull) for building a consensus of extracted sequences.

> The nearest alternative on crates.io is [`poasta`](https://crates.io/crates/poasta) (Broad Institute, pure Rust, gap-affine A\* alignment), which excels at larger graphs such as bacterial genes and HLA loci. For STR reads where graphs stay small and throughput across many loci matters, banded DP is faster. Both crates are pure Rust with no C dependencies.

## Usage

```toml
[dependencies]
poa-consensus = "0.3"
```

### Functional API

```rust
use poa_consensus::{consensus, consensus_multi, consensus_adaptive, PoaConfig};

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

// Adaptive two-pass: inspects graph statistics after pass 1 and automatically
// chooses multi-allele split, noise tightening, or semi-global switch for pass 2.
// Returns AdaptiveResult { consensuses, action } — action records which branch fired.
let result  = consensus_adaptive(&reads, 0, &PoaConfig::default())?;
let alleles = result.consensuses;  // Vec<Consensus>; one or two elements
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

The stateful API lets you inspect graph state between reads and reuse pre-allocated buffers across calls, which is important for high-throughput per-locus pipelines.

### Seed selection

```rust
use poa_consensus::{select_seed, SeedSelection};

// Automatically find the shortest fully-spanning read using terminal k-mer anchors.
// Falls back to longest read when no cluster structure is detected.
// Returns Err(NoSpanningReads) when reads split into two non-overlapping groups.
let seed_idx = select_seed(&reads, &SeedSelection::Auto)?;
let result = consensus(&reads, seed_idx, &PoaConfig::default())?;
```

### Orientation utilities

```rust
use poa_consensus::{auto_orient, reverse_complement};

// Orient all reads to match the strand of reads[seed_idx] before POA.
let oriented = auto_orient(&reads, seed_idx);
```

### Diagnostics

```rust
use poa_consensus::{diagnose, DiagnoseConfig};

let result = consensus(&reads, 0, &PoaConfig::default())?;
let warnings = diagnose(&result, &DiagnoseConfig::default());

if !warnings.is_clean() {
    for (is_warning, msg) in warnings.messages("consensus") {
        let level = if is_warning { "warning" } else { "note" };
        eprintln!("{level}: {msg}");
    }
}
```

`diagnose` checks four independent signals: read depth, coverage gaps, near-zero interior support, and structural competing alleles. All signals are also available as structured fields on `ConsensusWarnings` for programmatic handling.

## CLI

```
cargo install poa-consensus --features cli
```

```
poa-consensus reads.fa                    # FASTA or FASTQ, auto-detected
poa-consensus reads.fa --multi            # multi-allele mode
poa-consensus reads.fa --no-adaptive-band # disable adaptive band (on by default)
poa-consensus reads.fa --quiet            # suppress warnings; errors always printed
```

## Feature flags

| Flag   | Adds                                       |
|--------|--------------------------------------------|
| `cli`  | Binary target; pulls in `clap` + `noodles` |
| `plot` | SVG visualisation helpers via `kuva`       |

Default build: library only, zero external dependencies.

## AI disclosure

This library was developed with AI assistance (Claude). Architecture decisions, testing, validation, and algorithm designs, are the author's own. AI tooling served as an accelerator over existing skill. The library originated from prototyping a POA implementation inside [`bladerunner`](https://github.com/Psy-Fer/bladerunner), a nanopore STR detector and genotyper, and was spun out as a standalone crate when existing POA crates did not meet the throughput and API requirements of that use case.

## License

MIT - James Ferguson 2026
