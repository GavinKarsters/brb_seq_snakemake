#!/usr/bin/env python3
"""
Reads the BRB barcode library file and the sample-to-barcode mapping.
Produces:
  1. whitelist.txt          — one barcode sequence per line (STARsolo --soloCBwhitelist)
  2. barcode_to_sample.tsv  — sequence <TAB> sample_id  (used for count matrix renaming)
"""
import sys
import csv
import re


def sanitize(name):
    return re.sub(r"[\s/\\]", "_", name.strip())
    #return name.strip().replace(" ", "_").replace("/", "_").replace("\\", "_")


def main():
    brb_file      = sys.argv[1]
    demux_file    = sys.argv[2]
    whitelist_out = sys.argv[3]
    mapping_out   = sys.argv[4]

    # Parse BRB barcode reference: Name → sequence
    brb_seqs = {}
    with open(brb_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            brb_seqs[row["Name"].strip()] = row["B1"].strip()

    with open(whitelist_out, "w") as wl, open(mapping_out, "w") as mp:
        mp.write("barcode_sequence\tsample_id\n")
        with open(demux_file) as fin:
            reader = csv.DictReader(fin, delimiter="\t")
            for row in reader:
                sample_id    = sanitize(row["Sample_name"])
                barcode_name = row["Barcode"].strip()
                if barcode_name not in brb_seqs:
                    raise ValueError(
                        f"Barcode '{barcode_name}' for sample '{sample_id}' "
                        f"not found in BRB barcode reference file."
                    )
                seq = brb_seqs[barcode_name]
                wl.write(seq + "\n")
                mp.write(f"{seq}\t{sample_id}\n")

    print(f"[prepare_starsolo_inputs] Done.", flush=True)


if __name__ == "__main__":
    main()