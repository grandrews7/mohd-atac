# MOHD ATAC-seq Snakemake Pipeline

This repository contains a modular Snakemake pipeline for processing ATAC-seq data from raw FASTQ files through quality control, alignment, peak calling, and downstream analyses (e.g., fragment-length distribution, TSS enrichment, FRiP). The workflow is fully configurable via a `config.yaml` file.

---

## 📋 Table of Contents

1. [Features](#features)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Configuration (`config.yaml`)](#configuration-configyaml)
5. [Usage](#usage)
6. [Pipeline Steps](#pipeline-steps)
7. [Directory Structure](#directory-structure)
8. [Contributing](#contributing)
9. [License](#license)

---

## Features

* Automatic detection of samples and paired-end FASTQ files
* Adapter trimming (Cutadapt) and dual FastQC reports (pre- and post-trimming)
* Bowtie2 alignment (paired-end) with customizable parameters
* BAM filtering, duplicate marking (Picard), and indexing
* Peak calling with MACS3 (multiple q-value thresholds)
* Generation of BigWig coverage tracks
* Quality control plots: fragment-length distribution, TSS enrichment, FRiP score
* Support for ChromBPNet preprocessing and bias modeling
* Pseudoreplicate generation and IDR analysis

---

## Prerequisites

* **Conda** (Miniconda or Anaconda)
* **Snakemake** ≥ 6.0
* **Bowtie2**, **samtools**, **Picard**, **MACS3**, **bedtools**, **bedGraphToBigWig**, **Cutadapt**, **FastQC**
* **Python** 3.7+ with required libraries (listed in environment YAMLs under `envs/`)

---

## Installation

1. **Clone this repository**

   ```bash
   git clone https://github.com/grandrews7/mohd-atac.git
   cd mohd-atac
   ```

2. **Create and activate a Snakemake environment**

   ```bash
   conda create -n snakemake-env snakemake=7.0 python=3.8
   conda activate snakemake-env
   ```

3. **Install additional tools** (example using mamba)

   ```bash
   mamba install -n snakemake-env -c bioconda \
     bowtie2 samtools picard macs3 bedtools \
     bedGraphToBigWig cutadapt fastqc
   ```

4. **Verify installation**

   ```bash
   snakemake --version
   ```

---

## Configuration (`config.yaml`)

All pipeline parameters and file locations are controlled in `config.yaml`. A minimal example:

```yaml
# config.yaml

data_dir: "data/merged_data"
fastq_pattern: "{sample}_{read}.fastq.gz"
samples: []
reads:
  - "R1"
  - "R2"
genomes:
  - "GRCh38"
qvals:
  - 0.05
resources_dir: "resources"
```

Modify paths or add additional genomes, q-values, or samples as needed.

---

## Usage

Run the full pipeline with four cores (adjust `--cores` as needed):

```bash
snakemake \
  --cores 4 \
  --use-conda \
  --configfile config.yaml \
  --printshellcmds
```

Or submit via a cluster profile:

```bash
snakemake \
  --profile slurm \
  --configfile config.yaml
```

---

## Pipeline Steps

1. **FastQC** (raw)
2. **Cutadapt** trimming
3. **FastQC** (trimmed)
4. **Bowtie2** alignment (paired-end)
5. **BAM** filtering, mate fixing, deduplication (Picard)
6. **BEDPE** & tagAlign conversion
7. **MACS3** peak calling (fixed & variable q-values)
8. **bedGraph** → **BigWig** conversion
9. **QC metrics:** fragment-length distribution, TSS enrichment, FRiP
10. **ChromBPNet** preprocessing & bias modeling
11. **Pseudoreplicate** splitting & **IDR** analysis

---

## Directory Structure

```
├── Snakefile
├── config.yaml
├── envs/                # conda environment YAMLs
├── scripts/             # custom Python scripts
├── resources/           # genome FASTA, sizes, blacklist, etc.
├── results/             # pipeline outputs
├── logs/                # per-rule log files
└── README.md            # this file
```

---

## Contributing

Feel free to open issues or pull requests. For substantial changes, please discuss first via GitHub Issues.

---

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.
