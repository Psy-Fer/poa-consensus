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

- Rust 1.85 or later
- No system libraries required (pure Rust)

Input files must be FASTA or FASTQ (uncompressed or gzip-compressed). The format is
auto-detected from the first byte of the file.
