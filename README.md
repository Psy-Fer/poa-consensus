# poa-consensus

[![CI](https://github.com/Psy-Fer/poa-consensus/actions/workflows/ci.yml/badge.svg)](https://github.com/Psy-Fer/poa-consensus/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/poa-consensus.svg)](https://crates.io/crates/poa-consensus)
[![docs.rs](https://docs.rs/poa-consensus/badge.svg)](https://docs.rs/poa-consensus)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.85%2B-orange.svg)](https://www.rust-lang.org)

A pure-Rust banded Partial Order Alignment (POA) library for building consensus sequences from a set of reads.

POA aligns reads into a directed acyclic graph (DAG) using affine-gap dynamic programming. Insertions and deletions create branches resolved by read support. By default the consensus is best-fit: the engine builds both a heaviest-path and a majority-frequency (MSA-column) consensus and keeps whichever the reads support better.

**Target use cases:** short tandem repeat (STR) and VNTR loci, amplicon consensus, structural variants, per-locus nanopore or HiFi read sets (50 bp to ~20 kb with banded DP). Typical cost with the defaults: HG002 ONT STR loci ~10 ms per call (median); ~5 kb MUC1 VNTR haplotypes ~0.55 s; 20 kb reads ~3.3 s CPU / ~270 MB peak (about 2.4x less CPU with `--no-homopolymer`). Can be paired with [`bedpull`](https://github.com/Psy-Fer/bedpull) for building a consensus of extracted sequences.

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

The API is purely functional (the stateful `PoaGraph` was removed in 0.5.0; see the CHANGELOG).

### Flank-anchored and bridged consensus

```rust
use poa_consensus::{bridged_consensus, consensus_flanked, consensus_multi_flanked, PoaConfig};

let cfg = PoaConfig::default();
let left_flank: &[u8] = b"GATTACAGATTACAGATTAC"; // >= ~20 bp of unique sequence
let right_flank: &[u8] = b"CCGGTTAACCGGTTAACCGG";

// Consensus of the region between the flanks. Partial reads that do not span both
// flanks are auto-excluded when they are a meaningful fraction of the pool.
let repeat = consensus_flanked(&reads, left_flank, right_flank, &cfg)?;
let repeat_alleles = consensus_multi_flanked(&reads, left_flank, right_flank, &cfg)?;

// No read spans the region (e.g. a long expansion): join per-side consensuses,
// with a GapKind::Unknown gap recorded at the junction.
let bridged = bridged_consensus(&left_reads, 0, &right_reads, 0, &cfg)?;
```

### Homopolymer handling (on by default)

Single-allele consensus aligns homopolymer runs in compressed form (runs of 2 or
more identical bases become one symbol, and run lengths are re-estimated afterwards).
It also uses cheaper gap costs for indels inside a homopolymer. This stops noisy
reads from scattering run-length errors as mismatch pairs at homopolymer-rich
repeats such as (AAAAT)n. Each result is spliced into the plain consensus only inside
homopolymer runs, and kept only if it fits the reads strictly better. Defaults: `homopolymer_min_run: 2`,
`homopolymer_gap: Some((-2, -1))` (CLI `--homopolymer-min-run 2 --homopolymer-gap=-2,-1`).
It costs about 2.4x the CPU of the uncompressed engine. Multi-allele does not use it.
To turn it off:

```rust
let cfg = PoaConfig {
    homopolymer_min_run: 0,
    homopolymer_gap: None,
    ..PoaConfig::default()
};
```

CLI: `--no-homopolymer`.

### Seed selection

```rust
use poa_consensus::{select_seed, SeedSelection};

// consensus() seeds itself on the median-length read, so the index mainly serves as the
// orientation reference. select_seed is useful for detecting read sets that split into
// two non-overlapping groups: it returns Err(NoSpanningReads), and bridged_consensus is
// the right call. Otherwise it returns the shortest fully-spanning read (terminal k-mer
// anchors), falling back to the longest read.
let seed_idx = select_seed(&reads, &SeedSelection::Auto)?;
let result = consensus(&reads, seed_idx, &PoaConfig::default())?;
```

### Orientation utilities

```rust
use poa_consensus::auto_orient;

// Orient all reads to match the strand of reads[seed_idx] before POA.
let owned: Vec<Vec<u8>> = reads.iter().map(|r| r.to_vec()).collect();
let oriented = auto_orient(&owned, seed_idx);
let oriented: Vec<&[u8]> = oriented.iter().map(|r| r.as_ref()).collect();
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

`diagnose` checks five independent signals: read depth, coverage gaps, near-zero interior support, structural competing alleles, and suspected truncation (consensus much shorter than the median read). All signals are also available as structured fields on `ConsensusWarnings` for programmatic handling.

Every `Consensus` also carries `support_quality`, a per-base Phred score (capped at 60) for how decisively the reads back each base over its strongest alternative. Below about Q10 the reads disagree there, typically on a repeat length; `analysis::low_support_regions(&result, 10)` returns those stretches, and the CLI writes the scores as the FASTQ quality string with `--fastq`.

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
**Alignment**, **Multi-allele**, **Flanking**, and **Diagnostics**. Defaults match the library
defaults, so plain `poa-consensus reads.fa` is unchanged. Examples:

```
# Scoring and band tuning
poa-consensus reads.fa --match 2 --mismatch -2 --gap-open -3 --gap-extend -1
poa-consensus reads.fa --band-width 100 --adaptive-band-b 12 --adaptive-band-f 0.02

# Coverage / consensus
poa-consensus reads.fa --min-reads 5 --min-coverage-fraction 0.6
poa-consensus reads.fa --consensus-mode majority        # default (best-fit) | force MSA-column majority

# Flank-anchored: consensus of the repeat between two unique flanks
poa-consensus reads.fa --left-flank GATTACAGATTACAGATTAC --right-flank CCGGTTAACCGGTTAACCGG

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
