RD = config["result_dir"]

# Check the config file for the BAM setting
GENERATE_BAM = config.get("generate_bam", True)
# Check the config file for UMI dedup
UMIS_DEDUP = config.get("UMIs_dedup", True)

rule starsolo:
    input:
        r1        = config["lane_r1"],
        r2        = config["lane_r2"],
        whitelist = RD + "/starsolo/whitelist.txt",
        genome    = config["star_genome_dir"],
        gtf       = config["gtf"]
    output:
        matrix   = RD + "/starsolo/Solo.out/Gene/raw/matrix.mtx",
        barcodes = RD + "/starsolo/Solo.out/Gene/raw/barcodes.tsv",
        features = RD + "/starsolo/Solo.out/Gene/raw/features.tsv",
        star_log = RD + "/starsolo/Log.final.out",
        star_summary = RD + "/starsolo/Solo.out/Gene/Summary.csv",
        # Conditionally require the BAM file as an output
        bam      = [RD + "/bams/Aligned.sortedByCoord.out.bam"] if GENERATE_BAM else []
    params:
        prefix   = RD + "/starsolo/",
        bam_dir  = RD + "/bams/",
        cb_len   = config["barcode_len"],
        umi_len  = config["umi_len"],
        
        umi_dedup_arg = "1MM_Directional" if UMIS_DEDUP else "NoDedup",
        make_bam = "True" if GENERATE_BAM else "False",
        
        # BULLETPROOF: directly read the thread count from the config, not from Snakemake's dynamic allocation
        bam_args = (
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outSAMattributes NH HI AS NM CR UR CB UB GX GN sS sQ sM "
            f"--limitBAMsortRAM 27112659435 "
            f"--outBAMsortingThreadN {config['threads']['star']}"
        ) if GENERATE_BAM else "--outSAMtype None"
        
    log:
        RD + "/logs/starsolo/starsolo.log"
    benchmark:
        RD + "/benchmarks/starsolo/starsolo.tsv"
    retries: 3
    resources:
        mem_mb = lambda wildcards, attempt: 80000 + ((attempt - 1) * 30000),
        runtime = lambda wildcards, attempt: 720 + ((attempt - 1) * 1440)
    threads:
        config["threads"]["star"]
    container:
        config["containers"]["star"]
    shell:
        """
        umi_start=$(( {params.cb_len} + 1 ))

        # BULLETPROOF: --runThreadN explicitly uses the config file value!
        STAR \
            --runMode alignReads \
            --runThreadN {config[threads][star]} \
            --genomeDir {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --readFilesIn {input.r2} {input.r1} \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {input.whitelist} \
            --soloCBstart 1 --soloCBlen {params.cb_len} \
            --soloUMIstart $umi_start --soloUMIlen {params.umi_len} \
            --clipAdapterType CellRanger4 \
            --soloFeatures Gene \
            --soloUMIdedup {params.umi_dedup_arg} \
            --soloBarcodeReadLength 0 \
            --soloCellFilter EmptyDrops_CR \
            --outFileNamePrefix {params.prefix} \
            {params.bam_args} \
        > {log} 2>&1

        # Only move the BAM file if we actually told STAR to create it
        if [ "{params.make_bam}" == "True" ]; then
            mv {params.prefix}Aligned.sortedByCoord.out.bam {params.bam_dir}
        fi
        """

rule samtools_index:
    input:
        RD + "/bams/Aligned.sortedByCoord.out.bam"
    output:
        RD + "/bams/Aligned.sortedByCoord.out.bam.bai"
    log:
        RD + "/logs/samtools_index/samtools_index.log"
    benchmark:
        RD + "/benchmarks/samtools_index/samtools_index.tsv"
    threads: 4                           
    resources:                           
        mem_mb = 16000,                   
        runtime = 120                    
    container:
        config["containers"]["samtools"]
    shell:
        "samtools index -@ {threads} {input} {output} > {log} 2>&1"
        
        
rule split_bams_by_sample:
    input:
        bam = RD + "/bams/Aligned.sortedByCoord.out.bam",
        bai = RD + "/bams/Aligned.sortedByCoord.out.bam.bai", # Ensure it's indexed first
        barcode_map = RD + "/starsolo/barcode_to_sample.tsv"  # <--- Added this to get the names!
    output:
        split_dir = directory(RD + "/bams/split_by_sample")
    log:
        RD + "/logs/samtools_split/samtools_split.log"
    threads: 
        8
    resources:
        mem_mb = 32000,
        runtime = 240
    container:
        config["containers"]["samtools"]
    shell:
        """
        mkdir -p {output.split_dir}
        
        echo "Splitting massive BAM file by Cell Barcode (CB)..." > {log}
        samtools split -@ {threads} -d CB \
            -f "{output.split_dir}/%!.bam" \
            {input.bam} >> {log} 2>&1
            
        echo "Renaming BAMs to sample names..." >> {log}
        
        # Read the TSV file (Format: barcode_sequence \t Sample_name)
        while IFS=$'\\t' read -r barcode sample; do
            # Skip the header row if there is one
            if [[ "$barcode" == "barcode"* || -z "$barcode" ]]; then continue; fi
            
            src="{output.split_dir}/$barcode.bam"
            dest="{output.split_dir}/$sample.bam"
            
            # If the sequence BAM exists, rename it to the Sample Name
            if [[ -f "$src" ]]; then
                mv "$src" "$dest"
                # Create an index (.bai) for the new BAM so it's ready for IGV!
                samtools index -@ 2 "$dest"
            fi
        done < {input.barcode_map}
        
        # Rename the "unassigned" file so it makes sense
        if [[ -f "{output.split_dir}/-.bam" ]]; then
            mv "{output.split_dir}/-.bam" "{output.split_dir}/Unassigned_Reads.bam"
            samtools index -@ 2 "{output.split_dir}/Unassigned_Reads.bam"
        fi
        
        echo "Done!" >> {log}
        """