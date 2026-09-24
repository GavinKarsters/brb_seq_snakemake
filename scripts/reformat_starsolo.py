import sys
import os
import pandas as pd
import scipy.io

def main():
    if len(sys.argv) != 4:
        print("Usage: python reformat_starsolo.py <solo_dir> <mapping_file> <output_file>")
        sys.exit(1)

    solo_dir = sys.argv[1]
    mapping_file = sys.argv[2]
    out_file = sys.argv[3]

    # 1. Load the sparse matrix
    matrix_file = os.path.join(solo_dir, "matrix.mtx")
    print(f"Reading {matrix_file}...")
    # STARsolo mtx is: rows=genes, cols=cells (barcodes)
    mat = scipy.io.mmread(matrix_file).tocsc()
    
    # Convert to dense pandas DataFrame
    df = pd.DataFrame.sparse.from_spmatrix(mat).sparse.to_dense()

    # 2. Load Barcodes and Features (Genes)
    barcodes_file = os.path.join(solo_dir, "barcodes.tsv")
    features_file = os.path.join(solo_dir, "features.tsv")

    print(f"Reading {barcodes_file} and {features_file}...")
    # STARsolo doesn't have headers for these files
    barcodes = pd.read_csv(barcodes_file, sep="\t", header=None, names=["barcode"])
    # Features usually has 3 columns: Ensembl_ID, Gene_Symbol, Type
    features = pd.read_csv(features_file, sep="\t", header=None)

    # Assign row names (Genes) and column names (Barcodes)
    # Using column 1 (Gene Symbols) if available, otherwise column 0 (Ensembl IDs)
    if features.shape[1] > 1:
        df.index = features[1].values
    else:
        df.index = features[0].values
        
    df.columns = barcodes["barcode"].values

    # 3. Load the mapping file
    print(f"Reading {mapping_file}...")
    # We explicitly tell pandas to read the header
    mapping = pd.read_csv(mapping_file, sep="\t")
    
    # Strip any whitespace from column names just in case
    mapping.columns = mapping.columns.str.strip()

    # Create a dictionary to map barcode -> sample_id
    # Ensure we use the exact column names from your file: 'barcode_sequence' and 'sample_id'
    if 'barcode_sequence' not in mapping.columns or 'sample_id' not in mapping.columns:
         print(f"ERROR: Mapping file must contain 'barcode_sequence' and 'sample_id' columns. Found: {list(mapping.columns)}")
         sys.exit(1)

    barcode_to_sample = dict(zip(mapping["barcode_sequence"], mapping["sample_id"]))

    # 4. Rename the columns (barcodes) to Sample IDs
    # If a barcode in the matrix isn't in our mapping file, we keep the barcode sequence
    new_columns = [barcode_to_sample.get(bc, bc) for bc in df.columns]
    df.columns = new_columns

    # 5. Filter the DataFrame to only include columns that were in our mapping file
    # This removes empty droplets / barcodes that aren't our actual samples
    valid_samples = mapping["sample_id"].tolist()
    
    # Find columns that are in our valid samples list
    cols_to_keep = [col for col in df.columns if col in valid_samples]
    
    df_filtered = df[cols_to_keep]

    # 6. Save the final counts table
    print(f"Saving formatted counts to {out_file}...")
    # Write to TSV, keeping the index (Gene names)
    df_filtered.to_csv(out_file, sep="\t", index=True, index_label="Gene")
    print("Done!")

if __name__ == "__main__":
    main()