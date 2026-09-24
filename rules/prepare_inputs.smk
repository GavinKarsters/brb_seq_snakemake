RD = config["result_dir"]

rule prepare_starsolo_inputs:
    input:
        brb   = config["brb_barcodes"],
        demux = config["demultiplex_tsv"]
    output:
        whitelist = RD + "/starsolo/whitelist.txt",
        mapping   = RD + "/starsolo/barcode_to_sample.tsv"
    log:
        RD + "/logs/prepare_inputs/prepare_inputs.log"
    benchmark:
        RD + "/benchmarks/prepare_inputs/prepare_inputs.tsv"
    container:
        config["containers"]["python"]
    shell:
        """
        python scripts/prepare_starsolo_inputs.py \
            {input.brb} \
            {input.demux} \
            {output.whitelist} \
            {output.mapping} \
        > {log} 2>&1
        """