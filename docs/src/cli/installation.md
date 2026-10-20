# Installation

## From crates.io

```
cargo install poa-consensus --features cli
```

This installs the `poa-consensus` binary to `~/.cargo/bin/`.

## From source

```
git clone https://github.com/Psy-Fer/poa-consensus
cd poa-consensus
cargo install --path . --features cli
```

## Checking the installation

```
poa-consensus --version
poa-consensus --help
```

## Requirements

- Rust 1.88 or later for the CLI (the `noodles` FASTA/FASTQ dependency needs it). The
  library alone builds on Rust 1.85.
- No system libraries required (pure Rust).

Input files must be uncompressed FASTA or FASTQ; gzip is not supported (decompress with
`zcat reads.fa.gz | poa-consensus`). The format is auto-detected from the first byte.
