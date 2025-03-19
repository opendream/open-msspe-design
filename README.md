# Open-MSSPE-Design

![Build status](https://github.com/opendream/open-msspe-design/actions/workflows/open-msspe-design.yml/badge.svg)

Open-MSSPE-Design is a Rust-based pipeline for designing primers for Metagenomic Sequencing with Spiked Primer Enrichment (MSSPE). This approach supports viral diagnostics and genomic surveillance by enriching viral sequences during sequencing, as described in [Deng et al. (2020)](https://doi.org/10.1038/s41564-019-0637-9), and openly implemented in [nf-msspe](https://github.com/MaestSi/nf-msspe) by Simon Maestri. This implementation introduces significant revisions to optimize and automate the primer design process.

---

## Overview

Key features:
- Fully automated primer design workflow.
- Fully customizable kmer selection process. 
- Uses nearest-neighbor thermodynamic models from the standalone Primer3 package [ntthal](https://manpages.debian.org/testing/primer3/ntthal.1.en.html)
- Enhanced filtering for:
  - nucleotide repeats & homopolymers
  - specific minimum and maximum temperature of melting (tm) values
  - strict tm value ranges (within 2 standard deviations of mean)
  - hairpins (with same primer)
  - cross-dimers (across multiple primers/entire pool)
  - check strength of secondary structures via DeltaG calculations
- Rust implementation for performance and reliability.

---

## Getting Started

### Prerequisites
- [Rust](https://www.rust-lang.org/tools/install) (latest stable version)
- [Cargo](https://crates.io/) (Rust package manager)
- [mafft](https://mafft.cbrc.jp/alignment/server/index.html) for multiple sequence alignment
- [Primer3](https://primer3.org) for designing PCR primers.

### Installation
Clone the repository and navigate to the project directory:
```bash
git clone https://github.com/opendream/open-msspe-design.git
cd open-msspe-design/od-msspe
```

Install dependecies:
#### macOS
```bash
brew install rust
brew install mafft
brew install primer3
```

### Input Requirements
Prepare a FASTA file containing viral genome sequences. This will serve as the input for primer design.

---

## Usage

Build the pipeline:
```bash
cargo build --release
```

To run the pipeline:
```bash
./target/release/od-msspe
```
 OR build & run:
```bash
cargo run
```

### Required Config
- `--input`: Path to the input FASTA file containing viral genome sequences.
- `--output`: Directory where the designed primers will be saved.

### Optional Config
The following arguments control various aspects of the primer design process:

#### Primer Design Parameters
- `--kmer-size`: Size of k-mers used in primer design (default: 13).
- `--window-size`: Window size for genome scanning (default: 500).
- `--overlap-size`: Overlap size between adjacent windows (default: 250).
- `--max-mismatch-segments`: Maximum number of mismatched segments allowed (default: 1).
- `--max-iterations`: Maximum number of iterations for primer optimization (default: 1000).
- `--search-windows-size`: Size of search windows for primer candidates (default: 50).

#### Thermodynamic Parameters
- `--mv-conc`: Monovalent cation concentration in mM (default: 50.0).
- `--dv-conc`: Divalent cation concentration in mM (default: 3.0).
- `--dntp-conc`: dNTP concentration in mM (default: 0.0).
- `--dna-conc`: Primer concentration in nM (default: 250.0).
- `--annealing-temp`: Annealing temperature in °C (default: 25.0).

#### Temperature Thresholds
- `--min-tm`: Minimum melting temperature allowed (default: 30.0).
- `--max-tm`: Maximum melting temperature allowed (default: 60.0).
- `--tm-stddev`: Set the number of standard deviations away from the mean of the tm values (default: 2).
- `--max-self-dimer-any-tm`: Maximum Tm for self-dimer at any position (default: 10°C below max-tm).
- `--max-self-dimer-end-tm`: Maximum Tm for self-dimer at 3' end (default: 10°C below max-tm).
- `--max-hairpin-tm`: Maximum Tm for hairpin structures (default: 10°C below min-tm).
- `--max-delta-g`: Maximum delta G value for secondary structures (default: -9).

#### Boolean Flags
- `--keep-all`: Ignore all filtering criteria and keep all primers.
- `--check-cross-dimers`: Enable cross-dimer checking between all primer pairs.
- `--check-self-dimer`: Enable self-dimer checking for individual primers.
- `--check-hairpin`: Enable hairpin structure checking for individual primers.
- `--disable-tm-stddev`: Turns off tm-stddev config. Use if you do not want strictly similar tm values across all primers.
- `--do-align`: Perform MAFFT multiple sequence alignment if true. Set to false if sequence already aligned.

### Example
```bash
cargo run -- --input data/viral_genomes.fasta --output results/msspe_primers.csv --kmer-size=15
```

Debugging
```bash
RUST_LOG=info cargo run -- --input data/viral_genomes.fasta --output results/msspe_primers.csv
```

---

## Contributions
Contributions, issues, and feature requests are welcome! Please open an issue or submit a pull request to improve the project.

---

## License
This project is licensed under the [MIT License](LICENSE).

Maintainers of this project are grateful for the support of the [Skoll Foundation](https://skoll.org/).
