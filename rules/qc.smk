import os

RD = config["result_dir"]

# Dynamically extract the base filename (removes .fq.gz or .fastq.gz)
R1_BASE = os.path.basename(config["lane_r1"]).replace(".fastq.gz", "").replace(".fq.gz", "")
R2_BASE = os.path.basename(config["lane_r2"]).replace(".fastq.gz", "").replace(".fq.gz", "")

rule subset_fastq:
    """Instantly subset the first 4M lines (1M reads) to make FastQC finish in 60 seconds."""
    input:
        r1 = config["lane_r1"],
        r2 = config["lane_r2"]
    output:
        r1_sub = RD + f"/qc/subset/{R1_BASE}_subset.fq.gz",
        r2_sub = RD + f"/qc/subset/{R2_BASE}_subset.fq.gz"
    threads: 2
    # Note: On macOS, zcat sometimes acts weird with .gz files, use `gzip -cd` if zcat fails
    shell:
        """
        # Disable pipefail because `head` closes the pipe early, 
        # which intentionally causes a harmless SIGPIPE in gzip.
        set +o pipefail
        
        gzip -cd {input.r1} | head -n 4000000 | gzip > {output.r1_sub}
        gzip -cd {input.r2} | head -n 4000000 | gzip > {output.r2_sub}
        """

rule fastqc_lane:
    """Run FastQC only on the subsetted reads."""
    input:
        r1 = rules.subset_fastq.output.r1_sub,
        r2 = rules.subset_fastq.output.r2_sub
    output:
        # FastQC automatically appends '_fastqc' to the input filename
        html_r1 = RD + f"/qc/fastqc/{R1_BASE}_subset_fastqc.html",
        zip_r1  = RD + f"/qc/fastqc/{R1_BASE}_subset_fastqc.zip",
        html_r2 = RD + f"/qc/fastqc/{R2_BASE}_subset_fastqc.html",
        zip_r2  = RD + f"/qc/fastqc/{R2_BASE}_subset_fastqc.zip"
    params:
        outdir = RD + "/qc/fastqc"
    resources:
        runtime = 15 # Will only take ~1 minute now!
    log:
        RD + "/logs/fastqc/fastqc_lane.log"
    benchmark:
        RD + "/benchmarks/fastqc/fastqc_lane.tsv"
    threads: config["threads"]["fastqc"]
    container: config["containers"]["fastqc"]
    shell:
        """
        fastqc {input.r1} {input.r2} -o {params.outdir} --threads {threads} > {log} 2>&1
        """

rule plot_brb_qc:
    """Generate custom summary PDFs for the researcher."""
    input:
        counts = RD + "/counts/all_samples_counts.tsv",
        metadata = config["demultiplex_tsv"]
    output:
        report = RD + "/qc/BRB_seq_QC_Plots.pdf"
    log:
        RD + "/logs/qc/plot_brb_qc.log"
    container: 
        config["containers"]["r_env"]
    script:
        "../scripts/plot_qc.R"

rule multiqc:
    """Aggregate FastQC and STARsolo metrics."""
    input:
        fastqc_zip_r1 = rules.fastqc_lane.output.zip_r1,
        fastqc_zip_r2 = rules.fastqc_lane.output.zip_r2,
        star_log      = rules.starsolo.output.star_log,
        star_summary  = rules.starsolo.output.star_summary, # Ensure MultiQC waits for STARsolo!
        custom_plots  = rules.plot_brb_qc.output.report # Ensure custom plots run before multiqc finishes
    output:
        RD + "/qc/multiqc_report.html"
    params:
        outdir = RD + "/qc",
        dirs = RD + "/qc/fastqc " + RD + "/starsolo"
    resources:
        runtime = 30
    log:
        RD + "/logs/multiqc/multiqc.log"
    benchmark:
        RD + "/benchmarks/multiqc/multiqc.tsv"
    container: config["containers"]["multiqc"]
    shell:
        """
        multiqc {params.dirs} -o {params.outdir} --force > {log} 2>&1
        """