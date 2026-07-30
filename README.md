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
poa-consensus = "0.6"
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

// Single-allele consensus. The engine internally seeds on the median-length
// read and returns the better-fitting of a heaviest-path and a
// majority-frequency (MSA-column) consensus, so `seed_idx` is validated for
// bounds but does not bias the result on periodic repeats.
let result = consensus(&reads, 0, &PoaConfig::default())?;
println!("{}", String::from_utf8_lossy(&result.sequence));

// The result carries per-position coverage, per-base path weights, detected
// bubble sites, and interior coverage gaps for downstream analysis.
let stats = &result.graph_stats;
println!("bubbles: {}", stats.bubble_count);

// Multi-allele: returns one Consensus per detected allele.
let alleles = consensus_multi(&reads, 0, &PoaConfig::default())?;
```

The API is purely functional. (Earlier releases exposed a stateful `PoaGraph`
builder and a two-pass `consensus_adaptive`; both were removed in 0.5.0 when the
engine was rebuilt — see the CHANGELOG. `consensus` now covers the single-allele
case directly, including the seed-robustness the adaptive path used to provide.)

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

The full `PoaConfig` surface and the diagnostic thresholds are exposed as
flags, grouped in `--help` under **Band**, **Scoring**, **Coverage / consensus**,
**Alignment**, **Multi-allele**, and **Diagnostics**. Defaults match the library
defaults, so plain `poa-consensus reads.fa` is unchanged. Examples:

```
# Scoring and band tuning
poa-consensus reads.fa --match 2 --mismatch -2 --gap-open -3 --gap-extend -1
poa-consensus reads.fa --band-width 100 --adaptive-band-b 12 --adaptive-band-f 0.02

# Coverage / consensus
poa-consensus reads.fa --min-reads 5 --min-coverage-fraction 0.6
poa-consensus reads.fa --consensus-mode majority        # default (best-fit) | force MSA-column majority

# Multi-allele (raise min-allele-freq on noisy ONT data)
poa-consensus reads.fa --multi --min-allele-freq 0.4 --phasing-bubble-min-span 8

# Diagnostics thresholds
poa-consensus reads.fa --depth-warn-threshold 15 --truncation-ratio-threshold 0.5

poa-consensus --help                                     # full grouped list
```

## Feature flags

| Flag   | Adds                                       |
|--------|--------------------------------------------|
| `cli`  | Binary target; pulls in `clap` + `noodles` |

Default build: library only, zero external dependencies.

## AI disclosure

This library was developed with AI assistance (Claude). Architecture decisions, testing, validation, and algorithm designs, are the author's own. AI tooling served as an accelerator over existing skill. The library originated from prototyping a POA implementation inside [`bladerunner`](https://github.com/Psy-Fer/bladerunner), a nanopore STR detector and genotyper, and was spun out as a standalone crate when existing POA crates did not meet the throughput and API requirements of that use case.

## License

MIT - James Ferguson 2026
