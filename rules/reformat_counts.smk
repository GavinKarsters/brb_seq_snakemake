RD = config["result_dir"]

rule reformat_counts:
    input:
        matrix   = RD + "/starsolo/Solo.out/Gene/raw/matrix.mtx",
        barcodes = RD + "/starsolo/Solo.out/Gene/raw/barcodes.tsv",
        features = RD + "/starsolo/Solo.out/Gene/raw/features.tsv",
        mapping  = RD + "/starsolo/barcode_to_sample.tsv"
    output:
        RD + "/counts/all_samples_counts.tsv"
    params:
        solo_dir = RD + "/starsolo/Solo.out/Gene/raw/"
    log:
        RD + "/logs/reformat_counts/reformat_counts.log"
    benchmark:
        RD + "/benchmarks/reformat_counts/reformat_counts.tsv"
    container:
        config["containers"]["anaconda"]
    shell:
        """
        python scripts/reformat_starsolo.py \
            {params.solo_dir} \
            {input.mapping} \
            {output} \
        > {log} 2>&1
        """