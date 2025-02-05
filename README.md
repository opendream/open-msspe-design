# Open-MSSPE-Design

Open-MSSPE-Design is a Rust-based pipeline for designing primers for Metagenomic Sequencing with Spiked Primer Enrichment (MSSPE). This approach supports viral diagnostics and genomic surveillance by enriching viral sequences during sequencing, as described in [Deng et al. (2020)](https://doi.org/10.1038/s41564-019-0637-9), and openly implemented in [nf-msspe](https://github.com/MaestSi/nf-msspe) by Simon Maestri. This implementation introduces significant revisions to optimize and automate the primer design process.

---

## Overview

Key features:
- Fully automated primer design workflow.
- Rust implementation for performance and reliability.
- Uses nearest-neighbor thermodynamic models from the standalone Primer3 package [ntthal](https://manpages.debian.org/testing/primer3/ntthal.1.en.html)
- Enhanced filtering for:
  - nucleotide repeats
  - homopolymers
  - extremely high/low Tm values (<30 or >60)
  - hairpins (with same primer)
  - cross-dimers (across multiple primers/entire pool)

---

## Getting Started

### Prerequisites
- [Rust](https://www.rust-lang.org/tools/install) (latest stable version)
- [Cargo](https://crates.io/) (Rust package manager)
- `mafft` for multiple sequence alignment

### Installation
Clone the repository and navigate to the project directory:
```bash
git clone https://github.com/opendream/open-msspe-design.git
cd open-msspe-design
```

Install dependecies:
```bash
brew install rust
brew install mafft
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

### Options
- `--input`: Path to the input FASTA file containing viral genome sequences.
- `--output`: Directory where the designed primers will be saved.
- `--config`: (Optional) Path to a configuration file for advanced settings.

### Example
```bash
cd od-msspe
./target/release/od-msspe --input data/viral_genomes.fasta --output results/msspe_primers.csv
```

Debugging
```bash
cd od-msspe
RUST_LOG=info ./target/release/od-msspe --input data/viral_genomes.fasta --output results/msspe_primers.csv
```

---

## Contributions
Contributions, issues, and feature requests are welcome! Please open an issue or submit a pull request to improve the project.

---

## License
This project is licensed under the [MIT License](LICENSE).

Maintainers of this project are grateful for the support of the [Skoll Foundation](https://skoll.org/).
