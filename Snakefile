configfile: "config/config.yaml"

import pandas as pd
from snakemake.exceptions import WorkflowError

# Read the file
_demux = pd.read_csv(config["demultiplex_tsv"], sep="\t")

# ── 1. IDIOT-PROOF VALIDATION CHECK ──────────────────────────────────────────
# Find all rows where 'Sample_name' is duplicated
_duplicates = _demux[_demux.duplicated(subset=["Sample_name"], keep=False)]

if not _duplicates.empty:
    # Get a list of the offending sample names
    bad_samples = _duplicates["Sample_name"].unique().tolist()
    
    # Halt the pipeline instantly with a massive error message
    raise WorkflowError(
        f"\n\n"
        f"🛑 PIPELINE HALTED: Duplicate Sample Names Detected!\n"
        f"The demultiplex TSV contains duplicate sample names. This is not allowed.\n"
        f"The following samples appear more than once: {bad_samples}\n"
        f"Please fix the file before running the pipeline:\n"
        f"File: {config['demultiplex_tsv']}\n"
        f"-------------------------------------------------------------------\n"
    )
# ─────────────────────────────────────────────────────────────────────────────

# If we pass the check, proceed as normal
_demux["sample_id"] = (
    _demux["Sample_name"]
    .str.strip()
    .str.replace(r"[\s/\\]", "_", regex=True)
)
SAMPLES = _demux["sample_id"].tolist()

RD = config["result_dir"]

# Order of include is important! In some cases the output of one rule is an input to another, so we need to define the rules in the correct order.
include: "rules/prepare_inputs.smk"
include: "rules/starsolo.smk"
include: "rules/qc.smk"
include: "rules/reformat_counts.smk"

rule all:
    input:
        RD + "/counts/all_samples_counts.tsv",
        RD + "/qc/multiqc_report.html",
        RD + "/qc/BRB_seq_QC_Plots.pdf",
        RD + "/bams/split_by_sample" if config.get("generate_bam", True) else []