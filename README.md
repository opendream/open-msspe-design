# Open-MSSPE-Design

Open-MSSPE-Design is a Rust-based pipeline for designing primers for Metagenomic Sequencing with Spiked Primer Enrichment (MSSPE). This approach supports viral diagnostics and genomic surveillance by enriching viral sequences during sequencing, as described in Deng et al. (2020), and openly implemented in [nf-msspe](https://github.com/MaestSi/nf-msspe). This implementation introduces significant revisions to optimize and automate the primer design process.

---

## Overview

MSSPE combines metagenomic sequencing with targeted primer enrichment to enhance the detection of specific viral genomes. This repository provides an automated pipeline to generate primers tailored for specific viral targets, ensuring reproducibility and efficiency in primer design.

Key features:
- Fully automated primer design workflow.
- Rust implementation for performance and reliability.
- Enhanced filtering for di-nucleotide repeats, homopolymers, and potential secondary structures via deltaG calculations.

Reference:  
Deng, X., Achari, A., Federman, S. et al. *Metagenomic sequencing with spiked primer enrichment for viral diagnostics and genomic surveillance*. *Nat Microbiol* **5**, 443–454 (2020).  
[Read the paper](https://doi.org/10.1038/s41564-019-0637-9)

---

## Getting Started

### Prerequisites
- Rust (latest stable version)
- Cargo (Rust package manager)

### Installation
Clone the repository and navigate to the project directory:
```bash
git clone https://github.com/opendream/open-msspe-design.git
cd open-msspe-design
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
./target/release/open-msspe-design --input <path_to_fasta_file> --output <output_directory>
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
