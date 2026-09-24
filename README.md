# BRB-seq Snakemake pipeline

Snakemake-based pipeline for Bulk RNA Barcoding and sequencing (BRB-seq). This pipeline handles raw read quality control, STARsolo mapping, UMI deduplication, sample demultiplexing, count matrix generation, BAM file splitting (per-sample BAMs), and comprehensive QC reporting.

For every new BRB-seq run, you only need to update the `config/config.yaml` to point to your new big FASTQ lane files and update your run-specific sample demultiplexing TSV file. 

The final outputs include a dense gene-by-sample count matrix, per-sample BAM files, a custom QC PDF, and a MultiQC report.

---

## Repository contents

- `Snakefile` — main workflow. 
- `run_pipeline.sh` — wrapper for dry-run / unlock / touch / SLURM submission. Automatically handles Apptainer cache mapping.
- `slurm_status.sh` — helps Snakemake track job status so it can deem timed-out / OOM jobs as failures and mark them for auto-restart with more time / memory.
- `scripts/` — helper scripts executed by pipeline rules:
  - `prepare_starsolo_inputs.py` — matches sample barcodes against the BRB library to generate a STAR whitelist and mapping file.
  - `reformat_starsolo.py` — converts STARsolo's sparse matrix into a final, dense `all_samples_counts.tsv`, strips empty droplets, and renames columns to your Sample IDs.
  - `plot_qc.R` — R script for generating custom BRB-seq QC plots (UMI depth, Genes Detected, Saturation, PCA).
- `config/`
  - `config.yaml` — main configuration file containing reference paths, tool parameters, and pointers to the current run's FASTQs.
  - `demultiplex_tsv` — run-specific TSV mapping `Sample_name` to `Barcode` (e.g., `Elke1` -> `G16`).
  - `alitheia_brb_barcodes_v5D_384_brb.txt` — static reference dictionary of BRB-seq barcode names to nucleotide sequences.
- `rules/` — modular Snakemake rule files (`prepare_inputs.smk`, `qc.smk`, `starsolo.smk`, `reformat_counts.smk`).

---

## Workflow overview (rules)

| Stage | Rule(s) | What it does | Main outputs |
|------:|---------|--------------|--------------|
| 1 | `prepare_starsolo_inputs` | Extracts valid BRB barcode sequences for the current run to create a whitelist and mapping file. | `whitelist.txt`, `barcode_to_sample.tsv` |
| 2 | `subset_fastq` | Instantly extracts the first 1 million reads (4M lines) to speed up FastQC. | `qc/subset/*_subset.fq.gz` |
| 3 | `fastqc_lane` | Runs FastQC on the subsetted reads (finishes in ~60 seconds). | `qc/fastqc/*_subset_fastqc.html` |
| 4 | `starsolo` | Aligns reads to the genome, demultiplexes via CB/UMI tags, performs UMI deduplication, and generates a unified BAM file and sparse count matrix. | `starsolo/Solo.out/Gene/raw/matrix.mtx`, `bams/Aligned.sortedByCoord.out.bam` |
| 5 | `samtools_index` | Indexes the massive unified BAM file. | `bams/Aligned.sortedByCoord.out.bam.bai` |
| 6 | `split_bams_by_sample` | Splits the massive STAR BAM by Cell Barcode (`CB`) and renames the resulting files to your actual Sample IDs. | `bams/split_by_sample/<Sample_name>.bam` (+ `.bai`) |
| 7 | `reformat_counts` | Converts STAR's sparse matrix into a clean, dense TSV. Renames columns to Sample IDs and filters out unassigned droplets. | `counts/all_samples_counts.tsv` |
| 8 | `plot_brb_qc` | Generates a custom PDF with UMI distribution, genes detected, saturation curves, and PCA plots. | `qc/BRB_seq_QC_Plots.pdf` |
| 9 | `multiqc` | Aggregates FastQC and STARsolo metrics into a final HTML report. | `qc/multiqc_report.html` |

---

## Requirements

- Conda (Miniconda/Anaconda)
- A working Snakemake installation with **Apptainer** (formerly Singularity).
- Cluster submission is handled by `run_pipeline.sh` via `sbatch` (SLURM).
- All bioinformatics tools (STAR, Samtools, R, Python, FastQC, MultiQC) are provided automatically via containers defined in the config.

---

## 1) First-time installation (do once)

### 1.1 Environment configuration

Create a lightweight environment that runs Snakemake, Pandas (for the duplicate check), and Apptainer. 
```bash
conda create -n snakemake-c -c bioconda -c conda-forge snakemake pandas apptainer scipy
```

> **Note on Containers:** The workflow uses `--use-singularity` to run per-rule containers. The wrappers are pre-configured to cache images in a shared HPC directory (`/hpc/umc_kaaij/...`) so they only need to be downloaded once.

### 1.2 Clone the github repo

On your cluster or workstation, clone the repository and enter the folder:
```bash
git clone https://github.com/YourUsername/brb_seq_snakemake.git
cd brb_seq_snakemake
```

---

## 2) Configure a new BRB-seq run

For every new sequencing run, you only need to update two files.

### 2.1 Edit the Demultiplexing TSV
Open or replace your demultiplexing TSV (e.g., `config/demultiplex_tsv`). It must contain at least two columns: `Sample_name` and `Barcode`.
```tsv
Sample_name	Barcode
Elke1	G16
Elke2	H16
Elke3	I16
```
*(The pipeline will automatically cross-reference these Barcodes with the static Alitheia BRB barcode library).*

### 2.2 Edit the Config File
Open `config/config.yaml` and update the run-specific paths:
1. **`lane_r1`**: Path to the huge R1 FASTQ.
2. **`lane_r2`**: Path to the huge R2 FASTQ.
3. **`demultiplex_tsv`**: Path to the TSV you just edited.
4. *(Optional)* **`generate_bam`**: Set to `True` to generate per-sample BAMs, or `False` to skip BAM generation and save time/disk space.
5. *(Optional)* **`result_dir`**: The output directory name for this specific run.

---

## 3) Run the pipeline

### 3.1 Activate Snakemake environment & Start a screen session (recommended)
Before running the pipeline, it is strongly recommended to start a `screen` session. This keeps the Snakemake process alive if your SSH connection drops.
```bash
conda activate snakemake-c
screen -R snakemake
```

### 3.2 Dry-run (recommended before launching)
```bash
./run_pipeline.sh -n
```
*During the dry-run, check that the pipeline parses your demultiplex TSV correctly. If you accidentally included duplicate sample names, the pipeline will instantly crash and tell you which samples to fix.*

### 3.3 Run the pipeline
```bash
./run_pipeline.sh
```
The script will automatically submit the jobs to SLURM. You can safely detach from your screen (`Ctrl+A` then `D`).

### 3.4 Unlock (if Snakemake crashed previously)
```bash
./run_pipeline.sh -u
```

---

## 4) Outputs

All results are written to the `result_dir` specified in your `config.yaml`. 

The most important outputs are:
- `counts/all_samples_counts.tsv` — The final dense Gene-by-Sample count matrix ready for DESeq2 / Seurat.
- `qc/BRB_seq_QC_Plots.pdf` — Custom R-generated PDF showing UMI depth, Genes Detected, sequencing saturation, and a PCA plot (color-coded by `Experiment` if provided in metadata).
- `qc/multiqc_report.html` — Interactive QC report aggregating STAR mapping rates and FastQC metrics.
- `bams/split_by_sample/` — (If enabled) Individual BAM files cleanly renamed to your Sample IDs, complete with `.bai` indexes.

---
