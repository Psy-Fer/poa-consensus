# Single-Allele Consensus

## Basic usage

```rust
use poa_consensus::{consensus, select_seed, SeedSelection, PoaConfig};

let reads: Vec<&[u8]> = /* your reads */;

// Let the library pick the seed automatically.
let seed_idx = select_seed(&reads, &SeedSelection::Auto)?;
let result   = consensus(&reads, seed_idx, &PoaConfig::default())?;

println!("{}", String::from_utf8_lossy(&result.sequence));
println!("length: {} bp, depth: {}", result.sequence.len(), result.n_reads);
```

## Alignment mode

`AlignmentMode::SemiGlobal` is the default and should be used for extracted STR reads.

**Why not Global?** Global alignment has a subtle failure at homopolymer-flanked loci. If
the seed has more leading identical bases than most reads (e.g. a 17-A seed when the majority
have 15--16 A's), shorter reads must place their deletions somewhere inside the homopolymer
run. The aligner puts them at the end of the run rather than the start, because both
positions are equally penalised. This leaves the extra 1--2 A-nodes at the boundary of the
run on the only spine path, with no shortcut alternative for those reads. The boundary trim
approaches from the left and stops as soon as the first node has adequate coverage -- the
extra nodes are unreachable by the trim. The result is a consensus 1--2 bases longer than
the true majority length.

Semi-global eliminates this: shorter reads take free terminal gaps and never traverse the
extra nodes, giving them genuine low coverage that the boundary trim removes.

## Coverage output

`result.coverage` is the per-position read depth along the consensus path:

```rust
let min_cov = result.coverage.iter().copied().min().unwrap_or(0);
let mean_cov = result.coverage.iter().sum::<u32>() as f64 / result.coverage.len() as f64;
println!("min coverage: {min_cov}, mean coverage: {mean_cov:.1}");
```

Coverage can drop below the mean at coverage gaps (see
[Bridged Consensus](bridged-consensus.md)) and at boundaries where partial reads don't reach.

## Length variation

The heaviest-path consensus picks the allele length supported by the majority of reads. With
20 reads from a heterozygous locus (10 × CAG×20, 10 × CAG×30), the single-allele consensus
will be whichever length happened to have the seed read -- typically the seed length wins the
boundary contest. Use `consensus_multi` for genuinely diploid loci.

With length variation from **errors** (reads that are ±1 repeat unit due to error), the
heaviest path naturally picks the plurality length because the `(weight - 1)` normalisation
penalises single-read detours. A minority of reads with an extra or missing unit do not
shift the consensus.

## Checking result quality

```rust
use poa_consensus::{diagnose, DiagnoseConfig};

let warnings = diagnose(&result, &DiagnoseConfig::default());
if !warnings.is_clean() {
    for (is_warning, msg) in warnings.messages("my_locus") {
        eprintln!("{}: {msg}", if is_warning { "WARN" } else { "NOTE" });
    }
}

// Structured access
if let Some(trunc) = &warnings.truncation_suspected {
    eprintln!("truncation: consensus {}bp vs median read {}bp (ratio {:.2})",
        trunc.consensus_len, trunc.median_read_len, trunc.ratio);
}
```

See [Diagnostics and Warnings](diagnostics.md) for all available checks.
