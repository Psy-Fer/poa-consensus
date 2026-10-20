# Introduction

`poa-consensus` is a pure-Rust library for building a consensus sequence from a set of reads
using Partial Order Alignment (POA). It is designed for short tandem repeat (STR) and VNTR
loci and works on reads from ~50 bp to ~20 kb. It also suits amplicon sequencing and any per-locus use case where reads cover a shared region with length
variation, substitution errors, or true allelic differences.

## What POA does

POA aligns reads into a directed acyclic graph (DAG). Each node in the graph represents a
base; edges represent adjacency. When reads agree they share nodes; when they disagree
(substitution, insertion, or deletion) the graph branches. A deletion in some reads creates
a shortcut edge; an insertion creates an extra node.

The default consensus is best-fit: the engine builds both a heaviest-path consensus (the
path with the highest total read support) and a majority-frequency (MSA-column) consensus,
and keeps whichever the reads support better. Majority-frequency recovers homopolymer and
length-variable runs where the heaviest path over-calls.

## When to use this crate

| Scenario | Recommended |
|---|---|
| STR / VNTR loci (~50 bp to ~20 kb), high throughput across many loci | `poa-consensus` (banded DP) |
| Repeat between known unique flanks; some reads partial | `poa-consensus` `consensus_flanked` |
| Large non-repetitive graphs, bacterial genes, HLA loci | [`poasta`](https://crates.io/crates/poasta) (A\* alignment) |
| Diploid STR with two distinct allele lengths | `poa-consensus` with `consensus_multi` |
| No spanning reads, expansion too long for any single read | `poa-consensus` `bridged_consensus` |

The nearest alternative on crates.io is
[`poasta`](https://crates.io/crates/poasta) (Broad Institute), which uses exact gap-affine
A\* alignment and excels at large graphs. For STR reads where graphs stay small and
throughput across hundreds of loci per sample matters, banded DP is faster and simpler.
Both crates are pure Rust with no C dependencies.

## Feature flags

| Flag | What it adds |
|---|---|
| *(none)* | Library only; zero external dependencies |
| `cli` | Binary target; pulls in `clap` + `noodles` |

## Quick example

```rust
use poa_consensus::{consensus, PoaConfig};

let reads: Vec<&[u8]> = vec![
    b"CATCATCATCAT",
    b"CATCATCATCAT",
    b"CATCATCATCAT",
    b"CATCATCATCATCAT",  // one read with an extra repeat unit
];

let result = consensus(&reads, 0, &PoaConfig::default())?;
println!("{}", String::from_utf8_lossy(&result.sequence));
// => CATCATCATCAT  (majority length wins)
```
