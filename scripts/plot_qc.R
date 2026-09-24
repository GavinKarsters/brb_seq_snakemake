# scripts/plot_qc.R
suppressPackageStartupMessages({
    library(ggplot2)
})

# ── Logging ──────────────────────────────────────────────────────────────────
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ── I/O ──────────────────────────────────────────────────────────────────────
counts_file <- snakemake@input[["counts"]]
output_pdf  <- snakemake@output[["report"]]

# Safely check if metadata was provided in the Snakefile
meta_file <- NULL
if ("metadata" %in% names(snakemake@input)) {
    meta_file <- snakemake@input[["metadata"]]
}

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# ── 1. Read & Prep Counts ────────────────────────────────────────────────────
counts_raw <- read.table(counts_file, header=TRUE, sep="\t", check.names=FALSE)
gene_col <- colnames(counts_raw)[1]

# Remove duplicate genes, keeping only the first occurrence
counts_unique <- counts_raw[!duplicated(counts_raw[[gene_col]]), ]
row.names(counts_unique) <- counts_unique[[gene_col]]
counts_mat <- counts_unique[ , -1]

# Calculate basic QC metrics
qc_df <- data.frame(
    Sample = colnames(counts_mat),
    Total_UMIs = colSums(counts_mat),
    Genes_Detected = colSums(counts_mat > 0),
    stringsAsFactors = FALSE
)

# ── 2. Read & Prep Metadata ──────────────────────────────────────────────────
has_experiment <- FALSE
if (!is.null(meta_file) && file.exists(meta_file)) {
    # Read the TSV
    meta <- read.table(meta_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE, quote="")
    
    # Check if we have the necessary columns
    if ("Experiment" %in% colnames(meta) && "Sample_name" %in% colnames(meta)) {
        has_experiment <- TRUE
        # Merge QC data with Metadata
        qc_df <- merge(qc_df, meta[, c("Sample_name", "Experiment")], by.x="Sample", by.y="Sample_name", all.x=TRUE)
        # If any samples were missing from metadata, label them "Unknown"
        qc_df$Experiment[is.na(qc_df$Experiment)] <- "Unknown"
    }
}

# If no experiment column exists, create a dummy column so plotting code is identical
if (!has_experiment) {
    qc_df$Experiment <- "All Samples"
}

# ── 3. Calculate PCA ─────────────────────────────────────────────────────────
do_pca <- ncol(counts_mat) >= 3 # Need at least 3 samples for PCA
if (do_pca) {
    # Normalize to Counts per Million (CPM), then log2 transform
    cpm <- sweep(counts_mat, 2, qc_df$Total_UMIs, FUN="/") * 1e6
    logcpm <- log2(cpm + 1)
    
    # Filter out genes with zero variance
    var_genes <- apply(logcpm, 1, var)
    logcpm_sub <- logcpm[var_genes > 0, , drop=FALSE]
    
    if (nrow(logcpm_sub) > 10) {
        pca <- prcomp(t(logcpm_sub), scale. = TRUE)
        
        # Calculate % variance explained
        var_exp <- pca$sdev^2 / sum(pca$sdev^2)
        pc1_label <- paste0("PC1 (", round(var_exp[1] * 100, 1), "%)")
        pc2_label <- paste0("PC2 (", round(var_exp[2] * 100, 1), "%)")
        
        # Create PCA dataframe
        pca_df <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2], stringsAsFactors=FALSE)
        pca_df <- merge(pca_df, qc_df[, c("Sample", "Experiment")], by="Sample", all.x=TRUE)
    } else {
        do_pca <- FALSE
    }
}

# ── 4. Plotting Function ─────────────────────────────────────────────────────
# A helper function to add commas to large numbers on the axes
comma_fmt <- function(x) format(x, big.mark = ",", scientific = FALSE)

# Function to generate a suite of beautiful plots for any given dataframe
generate_plots <- function(q_df, p_df, title_suffix="") {
    
    # Order by Total UMIs (Highest to Lowest)
    q_df$Sample <- factor(q_df$Sample, levels = q_df$Sample[order(q_df$Total_UMIs, decreasing = TRUE)])
    
    # Base aesthetic theme
    my_theme <- theme_light(base_size = 14) +
                theme(plot.title = element_text(face = "bold"),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
                      panel.grid.minor = element_blank())

    plots <- list()
    
    # 1. Histogram
    plots[[1]] <- ggplot(q_df, aes(x = Total_UMIs)) +
        geom_histogram(bins = 30, fill = "#3A7CA5", color = "black", alpha = 0.8) +
        scale_x_continuous(labels = comma_fmt) +
        my_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=10)) +
        labs(title = paste0("1. Distribution of Total Counts", title_suffix),
             x = "Total UMIs per Sample", y = "Number of Samples")

    # 2. Bar Plot (Total UMIs)
    plots[[2]] <- ggplot(q_df, aes(x = Sample, y = Total_UMIs, fill = Experiment)) +
        geom_col(color = "black", linewidth = 0.1) +
        scale_y_continuous(labels = comma_fmt) +
        my_theme + 
        labs(title = paste0("2. Total UMIs per Sample", title_suffix),
             x = "Sample ID", y = "Total UMIs")
             
    # 3. Bar Plot (Genes Detected)
    plots[[3]] <- ggplot(q_df, aes(x = Sample, y = Genes_Detected, fill = Experiment)) +
        geom_col(color = "black", linewidth = 0.1) +
        scale_y_continuous(labels = comma_fmt) +
        my_theme + 
        labs(title = paste0("3. Unique Genes Detected", title_suffix),
             x = "Sample ID", y = "Genes Detected")

    # 4. Saturation
    plots[[4]] <- ggplot(q_df, aes(x = Total_UMIs, y = Genes_Detected, color = Experiment)) +
        geom_point(size = 3, alpha = 0.8) +
        scale_x_continuous(labels = comma_fmt) +
        scale_y_continuous(labels = comma_fmt) +
        my_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=10)) +
        labs(title = paste0("4. Depth vs Genes Detected", title_suffix),
             x = "Total UMIs", y = "Genes Detected")
             
    # 5. PCA
    if (!is.null(p_df)) {
        plots[[5]] <- ggplot(p_df, aes(x = PC1, y = PC2, color = Experiment)) +
            geom_point(size = 4, alpha = 0.8) +
            my_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=10)) +
            labs(title = paste0("5. Principal Component Analysis", title_suffix),
                 x = pc1_label, y = pc2_label)
    }
    
    return(plots)
}

# ── 5. Generate PDF ──────────────────────────────────────────────────────────
pdf(output_pdf, width = 12, height = 7)

# Print Global Overview (All Samples)
global_plots <- generate_plots(qc_df, if(do_pca) pca_df else NULL, " (Overview)")
for(p in global_plots) print(p)

# If experiments are defined, loop through and print a set for each experiment
if (has_experiment) {
    unique_exps <- unique(qc_df$Experiment)
    
    # Only create per-experiment plots if there is more than 1 experiment
    if (length(unique_exps) > 1) {
        for (exp in unique_exps) {
            # Subset data
            q_sub <- qc_df[qc_df$Experiment == exp, ]
            p_sub <- NULL
            if (do_pca) p_sub <- pca_df[pca_df$Experiment == exp, ]
            
            # Generate & Print
            sub_plots <- generate_plots(q_sub, p_sub, paste0(" (", exp, ")"))
            for(p in sub_plots) print(p)
        }
    }
}

dev.off()

# ── Close log ────────────────────────────────────────────────────────────────
sink(type = "message")
sink()
close(log)