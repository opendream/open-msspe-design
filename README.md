# Open-MSSPE-Design

Open-MSSPE-Design is a Rust-based pipeline for designing primers for Metagenomic Sequencing with Spiked Primer Enrichment (MSSPE). This approach supports viral diagnostics and genomic surveillance by enriching viral sequences during sequencing, as described in [Deng et al. (2020)](https://doi.org/10.1038/s41564-019-0637-9), and openly implemented in [nf-msspe](https://github.com/MaestSi/nf-msspe) by Simon Maestri. This implementation introduces significant revisions to optimize and automate the primer design process.

---

## Overview

MSSPE combines metagenomic sequencing with targeted primer enrichment to enhance the detection of specific viral genomes. This repository provides an automated pipeline to generate primers tailored for specific viral targets, ensuring reproducibility and efficiency in primer design.

Key features:
- Fully automated primer design workflow.
- Rust implementation for performance and reliability.
- Enhanced filtering for di-nucleotide repeats, homopolymers, and potential secondary structures via deltaG calculations.

---

## Getting Started

### Prerequisites
- [Rust](https://www.rust-lang.org/tools/install) (latest stable version)
- [Cargo](https://crates.io/) (Rust package manager)
- `bio` crate for handling biological sequence data
- `serde` and `serde_json` crates for configuration parsing
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
cargo install bio serde serde_json
```

Build the pipeline:
```bash
cargo build --release
```

### Input Requirements
Prepare a FASTA file containing viral genome sequences. This will serve as the input for primer design.

---

## Usage

To run the pipeline, use the following command:
```bash
cargo run
```

### Options
- `--input`: Path to the input FASTA file containing viral genome sequences.
- `--output`: Directory where the designed primers will be saved.
- `--config`: (Optional) Path to a configuration file for advanced settings.

### Example
```bash
./target/release/open-msspe-design --input data/viral_genomes.fasta --output results/
```

This generates primers in the specified output directory, optimized for MSSPE.

---

## Contributions
Contributions, issues, and feature requests are welcome! Please open an issue or submit a pull request to improve the project.

---

## License
This project is licensed under the [MIT License](LICENSE).
