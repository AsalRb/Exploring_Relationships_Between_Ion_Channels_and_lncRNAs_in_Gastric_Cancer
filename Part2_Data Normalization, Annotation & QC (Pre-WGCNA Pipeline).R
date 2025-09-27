
# Part 2:
# ================= Data Normalization, Annotation & QC (Pre-WGCNA Pipeline) =================

# ----- Normalization -----
library(DESeq2)

# 1. Create full DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = data.frame(condition = sample_labels),
  design = ~ condition
)

# 2. Estimate size factors for normalization
dds <- estimateSizeFactors(dds)

# 3. Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# 4. Optional: check size factors
sizeFactors(dds)

# 5. Checking result
# First few genes across first 3 samples
head(normalized_counts[, 1:3])
# Check for NA values
any(is.na(normalized_counts))  # Should be FALSE
# Check for non-finite values (Inf, -Inf, NaN)
any(!is.finite(as.matrix(normalized_counts)))  # Should also be FALSE
# Check for rows with all zero counts
sum(rowSums(normalized_counts) == 0)  # Should be a small number or 0
# Summary statistics
summary(as.vector(as.matrix(normalized_counts)))

# 6. Save for later
saveRDS(dds, "dds_full.rds")
saveRDS(normalized_counts, "normalized_counts.rds")

# Load normalized DESeqDataSet
dds <- readRDS("dds_full.rds")
# Load normalized count matrix
normalized_counts <- readRDS("normalized_counts.rds")


# ----- Renaming columns from UUID-filenames to short barcodes -----
# 1. Checking: Duplicate short_barcodes in sample_info
sum(duplicated(sample_info$short_barcode))

# 2. Remove the ".tsv" from file_name to match column names
sample_info$file_name_short <- gsub(".tsv$", "", sample_info$file_name)

# 3. Now create a named vector: names = old column names, values = short_barcode
name_map <- setNames(sample_info$short_barcode, sample_info$file_name_short)

# 4. Check match
all(colnames(raw_counts) %in% names(name_map))  # Should be TRUE

# 5. Checking: One-to-one match between file_name_short and expression matrix columns
length(unique(sample_info$file_name_short)) == ncol(raw_counts) # Should be TRUE

# 6. Rename columns
colnames(raw_counts) <- name_map[colnames(raw_counts)]
colnames(normalized_counts) <- name_map[colnames(normalized_counts)]

# 7. Check results
head(colnames(raw_counts))                   # Should now be TCGA-style names
head(colnames(normalized_counts))            # Same
any(duplicated(colnames(raw_counts)))        # Should return FALSE
any(duplicated(colnames(normalized_counts))) # Should return FALSE

# 8. Save for later
saveRDS(raw_counts, "raw_counts_named.rds")
saveRDS(normalized_counts, "normalized_counts_named.rds")

# Load raw counts with renamed sample columns
raw_counts <- readRDS("raw_counts_named.rds")
# Load normalized counts with renamed sample columns
normalized_counts <- readRDS("normalized_counts_named.rds")


# ----- Gene Annotation -----
# 1. Checking: Count how many duplicated Ensembl IDs
sum(duplicated(rownames(raw_counts)))

# 2. Finding suitable GENCODE version
library(rtracklayer)
library(dplyr)

# Versions to check (adjust as needed)
versions <- 22:43
matches <- data.frame(version = integer(), matched_genes = integer(), missing_genes = integer())

options(timeout = 600)  # increase to 10 minutes

for (ver in versions) {
  message("Checking GENCODE v", ver, "...")
  
  # Download GTF
  gtf_url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_%d/gencode.v%d.annotation.gtf.gz", ver, ver)
  tmpfile <- tempfile(fileext = ".gtf.gz")
  tryCatch({
    download.file(gtf_url, tmpfile, quiet = TRUE)
    gtf <- import(tmpfile)
    
    # Extract gene annotation
    gene_annotation <- as.data.frame(gtf) %>%
      filter(type == "gene") %>%
      select(ensembl_gene_id = gene_id) %>%
      mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id))
    
    # Compare to raw_counts
    matched <- sum(rownames(raw_counts) %in% gene_annotation$ensembl_gene_id)
    missing <- sum(!rownames(raw_counts) %in% gene_annotation$ensembl_gene_id)
    
    matches <- rbind(matches, data.frame(version = ver, matched_genes = matched, missing_genes = missing))
  }, error = function(e) {
    message("Failed for v", ver)
  })
}

# Show best match
matches <- matches[order(-matches$matched_genes), ]
print(matches)

# 3. Creating annotation file
library(rtracklayer)
library(dplyr)

# Download GENCODE v36 GTF
gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz"
download.file(gtf_url, destfile = "gencode.v36.annotation.gtf.gz")

# Import the GTF
gtf <- import("gencode.v36.annotation.gtf.gz")

# Keep only gene entries
gene_annotation <- as.data.frame(gtf) %>%
  filter(type == "gene") %>%
  select(ensembl_gene_id = gene_id,
         gene_name = gene_name,
         gene_biotype = gene_type) %>%
  mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id)) %>%
  filter(ensembl_gene_id %in% rownames(raw_counts))  # keep only genes in counts

# 4. Check coverage
matched_genes <- sum(rownames(raw_counts) %in% gene_annotation$ensembl_gene_id)
missing_genes <- setdiff(rownames(raw_counts), gene_annotation$ensembl_gene_id)

cat("Matched genes:", matched_genes, "\n")
cat("Missing genes:", length(missing_genes), "\n")

# 5. Save to file
write.table(gene_annotation,
            file = "TCGA_STAD_gene_annotation_v36.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Annotation file saved as: TCGA_STAD_gene_annotation_v36.tsv\n")

# 6. Load Gene Annotation file
gene_annotation <- read.delim("TCGA_STAD_gene_annotation_v36.tsv", stringsAsFactors = FALSE)


# ----- Quality Control (QC) -----
library(DESeq2)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(WGCNA)

# QC Data Preparation (Labeling, Ordering, and VST Transformation):
# 1. Create display names without changing sample_info
display_names <- sample_labels %>%
  as.character() %>%
  { paste0(., "_", ave(seq_along(.), ., FUN = seq_along)) }

# 2. Order indices so Normals first, Tumors second (for display)
order_idx <- order(sample_labels)
sample_labels_ordered <- sample_labels[order_idx]
display_names_ordered <- display_names[order_idx]

# 3. Create DESeq2 object from *raw* counts (ordered for display)
dds_qc <- DESeqDataSetFromMatrix(
  countData = raw_counts[, order_idx],
  colData = data.frame(condition = sample_labels_ordered),
  design = ~ condition
)

# 4. Variance stabilizing transformation (all genes)
vsd_qc <- vst(dds_qc, blind = TRUE)


# PCA BEFORE outlier removal:
# 1. Pretty display names (doesn't change downstream objects)
display_names <- sample_labels |>
  as.character() |>
  (\(.) paste0(., "_", ave(seq_along(.), ., FUN = seq_along)))()

# 2. Order for display (Normals first)
order_idx <- order(sample_labels)
sample_labels_ordered <- sample_labels[order_idx]
display_names_ordered <- display_names[order_idx]

# 3. Build DESeq2 object from RAW counts (ordered for display only)
dds_qc <- DESeqDataSetFromMatrix(
  countData = raw_counts[, order_idx],
  colData   = data.frame(condition = sample_labels_ordered),
  design    = ~ condition
)

# 4. VST (QC-style, all genes)
vsd_qc <- vst(dds_qc, blind = TRUE)

# 5. PCA data
pca_before_df <- plotPCA(vsd_qc, intgroup = "condition", returnData = TRUE)
pv_before     <- round(100 * attr(pca_before_df, "percentVar"))

# 6. Write x/y limits to reuse in the AFTER plot for identical axes
xlim_all <- range(pca_before_df$PC1)
ylim_all <- range(pca_before_df$PC2)
saveRDS(list(xlim = xlim_all, ylim = ylim_all), file = "PCA_limits_before.rds")

# 7. Consistent colors across both figures
pca_cols <- c(Normal = "#F17BB0", Tumor = "#17BECF")  # pink / teal

# 8. Plot
PCA_before_outlier_removal <- ggplot(pca_before_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = pca_cols) +
  # coord_cartesian(xlim = xlim_all, ylim = ylim_all)  
  xlab(paste0("PC1 (", pv_before[1], "%)")) +
  ylab(paste0("PC2 (", pv_before[2], "%)")) +
  ggtitle("PCA (All genes) — BEFORE outlier removal") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

PCA_before_outlier_removal

# 9. Save
ggsave("QC_PCA_before_outlier_removal.pdf", PCA_before_outlier_removal, width = 7, height = 5)


# Heatmap BEFORE outlier removal:
# 1. Select top 50 most variable genes on the full cohort (from VST matrix)
mat_before <- assay(vsd_qc)
top50_ids  <- head(order(rowVars(mat_before), decreasing = TRUE), 50)

# 2. Z-score per gene (center/scale rows)
z_before <- t(scale(t(mat_before[top50_ids, ])))

# 3. Fix color scale limit from BEFORE (99th percentile of |z|)
zlim <- as.numeric(quantile(abs(z_before), 0.99, na.rm = TRUE))
zlim <- max(2.5, zlim)  # keep a sensible minimum range
breaks <- seq(-zlim, zlim, length.out = 101)

# 4. Map Ensembl IDs → symbols for nicer row names
ens_ids <- rownames(z_before)
symb <- gene_annotation$gene_name[match(ens_ids, gene_annotation$ensembl_gene_id)]
rownames(z_before) <- ifelse(is.na(symb) | symb == "", ens_ids, symb)

# 5. Column annotation (Group) and pretty labels
annotation_col <- data.frame(Group = sample_labels_ordered)
rownames(annotation_col) <- colnames(z_before)

annotation_colors <- list(
  Group = c(Normal = "#F17BB0", Tumor = "#17BECF")  # pink / teal
)

# 6. Save items needed for AFTER plot (gene set + color range)
saveRDS(list(top50_ids = rownames(mat_before)[top50_ids], zlim = zlim),
        file = "heatmap_comparability.rds")

# 7. Plot & save PDF
pdf("QC_Heatmap_Top50_before_outlier_removal.pdf", width = 8, height = 10)
pheatmap(
  z_before,
  color            = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(100),
  breaks           = breaks,
  annotation_col   = annotation_col,
  annotation_colors= annotation_colors,
  labels_col       = display_names_ordered,
  show_rownames    = TRUE,
  fontsize_col     = 6,
  main             = "Top 50 variable genes — BEFORE outlier removal (z-score)"
)
dev.off()


# Sample–Sample Correlation (before outlier removal):
library(pheatmap)
library(matrixStats)
library(ggplot2)

# 1. Build correlation matrix on vst expression 
mat_qc <- assay(vsd_qc)                # genes x samples
# Choose correlation type: "pearson" (default) or "spearman"
cor_type <- "pearson"
C <- cor(mat_qc, method = cor_type, use = "pairwise.complete.obs")  # samples x samples

# 2. Order the matrix by hierarchical clustering on correlation distance ---
D <- as.dist(1 - C)                    # distance = 1 - correlation
hc <- hclust(D, method = "average")
ord <- hc$order

C_ord <- C[ord, ord]
labs_barcode <- colnames(C)[ord]       # original barcodes in order
labs_pretty  <- display_names_ordered[ord]  # pretty labels aligned to order

# Safety check: lengths must match
stopifnot(length(labs_pretty) == ncol(C_ord))

# 3. Column annotations (Normal/Tumor) for the heatmap ---
anno <- data.frame(Group = sample_labels_ordered)
rownames(anno) <- colnames(C)          # annotation rows must match ORIGINAL colnames
anno <- anno[labs_barcode, , drop = FALSE]  # reorder to dendrogram order

# 4. Fixed colors for clarity
annotation_colors <- list(
  Group = c(Normal = "#2b6cb0", Tumor = "#d53f8c")  # blue, pink
)

# 5. Pretty heatmap (no row labels; show pretty sample labels on columns) ---
pdf("QC_Sample–Sample Correlation_before_outlier_removal.pdf", width = 8, height = 10)
pheatmap(
  C_ord,
  # keep the same ordering as hc object
  clustering_distance_rows = D,
  clustering_distance_cols = D,
  clustering_method = "average",
  
  annotation_col    = anno,
  annotation_colors = annotation_colors,
  
  labels_col     = labs_pretty,   # pretty column labels
  show_rownames  = FALSE,         # <- hide row labels (instead of labels_row)
  treeheight_row = 0,             # <- remove row dendrogram 
  
  main = sprintf("Sample–Sample %s Correlation (vst, pre-outlier)",
                 tools::toTitleCase(cor_type)),
  display_numbers = FALSE,
  legend          = TRUE,
  # setting to FALSE if it's needed to hide the annotation legend box
  annotation_legend = TRUE,
  fontsize_col   = 6,
  border_color   = NA
)
dev.off()

# 6. Quick automatic flags for low-correlation samples
# Flag samples whose median correlation to others is low.
med_cor <- apply(C, 2, function(x) median(x[!is.na(x) & x < 0.999]))  # exclude self-corr
flag_thresh <- quantile(med_cor, 0.05)   # bottom 5% as a heuristic
flags <- med_cor <= flag_thresh

message("Low-correlation threshold (5th pct): ", round(flag_thresh, 3))
message("Flagged samples (potential QC concerns):")
print(names(med_cor)[flags])

# 7. Save flagged samples only (1 per line, plain text)
writeLines(names(med_cor)[flags],
           "QC_flagged_samples_correlation.txt")

# 8. Load flagged samples
readLines("QC_flagged_samples_correlation.txt")

# 9. Scatter of median correlation
df_med <- data.frame(
  sample = seq_along(med_cor),
  median_cor = med_cor,
  label = display_names_ordered,
  group = sample_labels_ordered
)
p_scatter <- ggplot(df_med, aes(sample, median_cor, color = group)) +
  geom_point() +
  geom_hline(yintercept = flag_thresh, linetype = 2, color = "grey40") +
  scale_x_continuous(breaks = NULL) +
  labs(title = "Per-sample median correlation (vst)",
       y = sprintf("Median %s correlation", cor_type), x = "Samples") +
  theme_bw()

p_scatter

ggsave("QC_MedianCorrelationScatter_preOutlier.pdf",
       p_scatter, width = 8, height = 5)


# Find and remove outlier samples by Hierarchical clustering (before DEG):
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 1. Load inputs
dds <- readRDS("dds_full.rds")  # built earlier with ~ condition

# Nice, human-friendly sample labels 
grp <- colData(dds)$condition
display_names <- paste0(as.character(grp), "_",
                        ave(seq_along(grp), grp, FUN = seq_along))

# 2. Build expression matrix for clustering 
# Use vst with blind=TRUE for QC-style clustering
vsd <- vst(dds, blind = TRUE)
expr_all <- assay(vsd)                # genes x samples (numeric)

# 3. Use all genes or top variable genes
top_n <- 5000                         
rv <- rowVars(expr_all)
keep_idx <- order(rv, decreasing = TRUE)[seq_len(min(top_n, nrow(expr_all)))]
expr_matrix <- expr_all[keep_idx, , drop = FALSE]

# Apply friendly labels only to this QC matrix
colnames(expr_matrix) <- display_names

# 4. Build sample dendrogram 
sampleTree <- hclust(dist(t(expr_matrix)), method = "average")

# 5. Sweep several cut heights and summarize
cand_cutH <- sort(unique(c(
  quantile(sampleTree$height, c(0.90, 0.92, 0.95, 0.97, 0.99), na.rm = TRUE),
  seq(from = floor(min(sampleTree$height)/10)*10,
      to   = ceiling(max(sampleTree$height)/10)*10,
      by   = 20)
)))

res_list <- lapply(cand_cutH, function(h) {
  cl <- cutreeStatic(sampleTree, cutHeight = as.numeric(h), minSize = 10)
  outliers <- colnames(expr_matrix)[cl == 0]
  data.frame(
    cutH = as.numeric(h),
    n_outliers = length(outliers),
    outlier_examples = paste(head(outliers, 5), collapse = ", ")
  )
})

cutH_summary <- do.call(rbind, res_list)
cutH_summary <- cutH_summary[order(cutH_summary$cutH), ]
print(cutH_summary, row.names = FALSE)

# 6. Save the sweep summary for the record
write.table(cutH_summary, "wgcna_outlier_cutH_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 7. Load cutH_summary
cutH_summary <- read.table("wgcna_outlier_cutH_summary.tsv",
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# Auto-pick cut height + Plot + Save + Finalize:
stopifnot(inherits(sampleTree, "hclust"))
stopifnot(is.matrix(expr_matrix))
stopifnot(all(sampleTree$labels %in% colnames(expr_matrix)))

# 1. Auto-pick a reasonable cut height (≤ 5% outliers)
target_frac <- 0.05
n_samples   <- ncol(expr_matrix)

best <- cutH_summary[cutH_summary$n_outliers <= target_frac * n_samples, ]
if (nrow(best)) {
  cutH_auto <- min(best$cutH)
} else {
  cutH_auto <- max(cutH_summary$cutH)  # fallback if none meet the target
}
cat("Auto-selected cut height:", cutH_auto, "\n")
writeLines(sprintf("cutH_auto\t%s", cutH_auto), "wgcna_cutH_auto.txt")

# 2. Make outlier calls at the auto-selected height
clust_auto <- cutreeStatic(sampleTree, cutHeight = cutH_auto, minSize = 10)

# 3. Define colors
col_keep   <- "#1f77b4"  # blue
col_out    <- "#d62728"  # red
col_normal <- "#2ca02c"  # green
col_tumor  <- "#17becf"  # teal

# 4. Build pretty-label -> condition mapping from dds
group <- colData(dds)$condition
display_names <- paste0(as.character(group), "_",
                        ave(seq_along(group), group, FUN = seq_along))
label_map <- data.frame(
  label     = display_names,        # pretty label (matches expr_matrix colnames)
  barcode   = colnames(dds),
  condition = as.character(group),
  stringsAsFactors = FALSE
)

stopifnot(all(colnames(expr_matrix) %in% label_map$label))

# condition colors in dendrogram order
labs <- sampleTree$labels
grp_lookup <- setNames(label_map$condition, label_map$label) #label -> condition
grp_for_expr <- grp_lookup[colnames(expr_matrix)]
names(grp_for_expr) <- colnames(expr_matrix)
grp_for_expr <- grp_for_expr[labs]
group_col_vec <- ifelse(grp_for_expr == "Tumor", col_tumor, col_normal)

# outlier bar in dendrogram order
outlier_vec <- ifelse(clust_auto == 0, col_out, col_keep)
names(outlier_vec) <- colnames(expr_matrix)
outlier_vec <- outlier_vec[labs]

# 5. Save ONE publication figure (tree + bars + separate legend row)
pdf("SampleClustering_withOutliers.pdf", width = 9, height = 7)  # vector output
layout(matrix(c(1,2), nrow=2), heights=c(8,1))  # 8:1 ratio

# (top) dendrogram + bars
par(mar = c(1.2, 4, 3.2, 1))
plotDendroAndColors(
  sampleTree,
  colors          = cbind(Outlier = outlier_vec, Group = group_col_vec),
  groupLabels     = NULL,           # no bar labels inside plot, cleaner
  main            = "Sample clustering (auto-selected cut height)",
  dendroLabels    = FALSE,          # hide dense sample names
  cex.colorLabels = 0.9,
  addGuide        = TRUE,
  guideHang       = 0.05
)

# (bottom) clean legend
par(mar = c(0,0,0,0))
plot.new()
legend("center",
       legend = c("Outlier samples", "Kept samples", "Normal samples", "Tumor samples",
                  sprintf("Cut height = %s", round(cutH_auto, 1))),
       fill   = c(col_out, col_keep, col_normal, col_tumor, NA),
       border = NA, bty = "n", horiz = TRUE, cex = 0.9,
       text.col = c("black","black","black","black","darkgreen"))
dev.off()


# Finalize the cut height, save outlier calls:
cutH <- cutH_auto   # or e.g., cutH <- 200
writeLines(sprintf("cutH_used\t%s", cutH), "wgcna_cutH_used.txt")

clust <- cutreeStatic(sampleTree, cutHeight = cutH, minSize = 10)
outlier_samples <- colnames(expr_matrix)[clust == 0]

cat("Cut height:", cutH, "\n")
cat("Number of outliers:", length(outlier_samples), "\n")
if (length(outlier_samples)) print(outlier_samples)

# Save label-level outlier calls (pretty labels)
write.table(
  data.frame(sample_label = colnames(expr_matrix),
             cluster = clust,
             is_outlier = clust == 0),
  file = sprintf("wgcna_outliers_labels_cutH_%s.tsv", cutH),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Save matched barcodes (handy later)
outlier_barcodes <- label_map$barcode[match(outlier_samples, label_map$label)]
write.table(outlier_barcodes,
            "QC_flagged_samples_HC.tsv",
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# Load outliers_labels and outlier_barcodes
outlier_calls <- read.table(
  sprintf("wgcna_outliers_labels_cutH_%s.tsv", cutH),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

outlier_barcodes <- read.table(
  "QC_flagged_samples_HC.tsv",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)[,1]


# Reconcile QC outliers (Correlation vs Hierarchical Clustering):
cor_outliers_file <- "QC_flagged_samples_correlation.txt"
hc_outliers_file  <- "QC_flagged_samples_HC.tsv"

# 1. Normalize IDs (trim/uppercase) so set ops are reliable
.norm <- function(x) unique(sort(toupper(trimws(x))))

# 2. Load files (empty if missing)
cor_outliers <- if (file.exists(cor_outliers_file)) {
  .norm(scan(cor_outliers_file, what = character(), quiet = TRUE))
} else character(0)

hc_outliers <- if (file.exists(hc_outliers_file)) {
  .norm(scan(hc_outliers_file, what = character(), quiet = TRUE))
} else character(0)

cat("Loaded:", length(cor_outliers), "correlation-flagged samples\n")
cat("Loaded:", length(hc_outliers),  "hierarchical clustering outliers\n")

# 3. Sets
consensus_outliers <- intersect(cor_outliers, hc_outliers)   # BOTH
union_outliers     <- union(cor_outliers, hc_outliers)       # EITHER
only_cor           <- setdiff(cor_outliers, hc_outliers)     # correlation-only
only_hc            <- setdiff(hc_outliers,  cor_outliers)    # HC-only

# 4. Choose removal policy: "consensus" (safer), "union" (strictest), or "hc_only"
REMOVE_MODE <- "consensus"

final_outliers <- switch(
  REMOVE_MODE,
  "consensus" = consensus_outliers,
  "union"     = union_outliers,
  "hc_only"   = hc_outliers,
  { warning("Unknown REMOVE_MODE; defaulting to 'consensus'"); consensus_outliers }
)

cat(sprintf("\nRemoval policy: %s\n", REMOVE_MODE))
cat("n(consensus) =", length(consensus_outliers),
    "| n(union) =", length(union_outliers),
    "| n(HC-only) =", length(only_hc),
    "| n(COR-only) =", length(only_cor), "\n")
cat("n(final_outliers) =", length(final_outliers), "\n\n")

# 5. Map to pretty labels if dds exists; otherwise skip gracefully
if (exists("dds")) {
  grp <- colData(dds)$condition
  display_names <- paste0(as.character(grp), "_",
                          ave(seq_along(grp), grp, FUN = seq_along))
  label_map <- data.frame(
    barcode   = .norm(colnames(dds)),
    label     = display_names,
    condition = as.character(grp),
    stringsAsFactors = FALSE
  )
} else {
  label_map <- NULL
}

# 6. Writer that includes pretty labels when available
.write_set <- function(vec, file, label_map = NULL) {
  vec <- .norm(vec)
  if (length(vec) == 0) {
    write.table(data.frame(barcode = character(0), pretty_label = character(0)),
                file, sep = "\t", quote = FALSE, row.names = FALSE)
    return(invisible(NULL))
  }
  if (!is.null(label_map)) {
    # normalize to match label_map$barcode
    pretty <- label_map$label[match(vec, label_map$barcode)]
    out <- data.frame(barcode = vec, pretty_label = ifelse(is.na(pretty), vec, pretty))
  } else {
    out <- data.frame(barcode = vec, pretty_label = vec)
  }
  write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# 7. Save all reports in the current working directory
.write_set(consensus_outliers, "consensus_outliers.tsv",                label_map)
.write_set(union_outliers,     "union_outliers.tsv",                    label_map)
.write_set(only_hc,            "hc_only_outliers.tsv",                  label_map)
.write_set(only_cor,           "correlation_only_outliers.tsv",         label_map)

# Final list for actual removal (barcodes only, one per line)
writeLines(.norm(final_outliers), "final_outliers_to_remove.txt")

cat("QC reports saved under:", getwd(), "\n",
    "- consensus_outliers.tsv\n",
    "- union_outliers.tsv\n",
    "- hc_only_outliers.tsv\n",
    "- correlation_only_outliers.tsv\n",
    "- final_outliers_to_remove.txt\n")

# 8. load files
consensus_outliers_df <- read.table("consensus_outliers.tsv",
                                    sep = "\t", header = TRUE, stringsAsFactors = FALSE)
union_outliers_df  <- read.table("union_outliers.tsv", sep = "\t", header = TRUE)
hc_only_df         <- read.table("hc_only_outliers.tsv", sep = "\t", header = TRUE)
cor_only_df        <- read.table("correlation_only_outliers.tsv", sep = "\t", header = TRUE)
final_outliers <- readLines("final_outliers_to_remove.txt")


# Create clean objects using reconciled outliers:
# 1. Load inputs
dds       <- readRDS("dds_full.rds")
raw_counts <- readRDS("raw_counts_named.rds")
normalized_counts <- readRDS("normalized_counts_named.rds")

# Load the final list (In fresh session condition):
if (!exists("final_outliers")) {
  final_outliers <- if (file.exists("final_outliers_to_remove.txt")) {
    scan("final_outliers_to_remove.txt", what = character(), quiet = TRUE)
  } else character(0)
}

# 2. Remove outlier samples and build clean objects
dds_clean <- dds[, !colnames(dds) %in% final_outliers]
colData(dds_clean)$condition <- droplevels(colData(dds_clean)$condition)

raw_counts_clean        <- raw_counts[, colnames(raw_counts) %in% colnames(dds_clean), drop = FALSE]
normalized_counts_clean <- normalized_counts[, colnames(normalized_counts) %in% colnames(dds_clean), drop = FALSE]

cat("Original samples:", ncol(dds), "\n")
cat("Removed samples :", length(final_outliers), "\n")
cat("Kept samples    :", ncol(dds_clean), "\n")
print(table(colData(dds_clean)$condition))

# 3. save clean objects
saveRDS(dds_clean, "dds_clean.rds")
saveRDS(raw_counts_clean, "raw_counts_clean.rds")
saveRDS(normalized_counts_clean, "normalized_counts_clean.rds")

# 4. Load clean objects
dds_clean <- readRDS("dds_clean.rds")
raw_counts_clean <- readRDS("raw_counts_clean.rds")
normalized_counts_clean <- readRDS("normalized_counts_clean.rds")


# AFTER outlier removal: PCA & Heatmap (comparable to BEFORE)
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap); library(matrixStats)
})

# 1) Load cleaned data (or skip if already in memory)
if (!exists("dds_clean")) dds_clean <- readRDS("dds_clean.rds")

# 2) VST on cleaned set
vsd_clean <- vst(dds_clean, blind = TRUE)


# PCA (AFTER outlier removal):
pca_after_df <- plotPCA(vsd_clean, intgroup = "condition", returnData = TRUE)
pv_after     <- round(100 * attr(pca_after_df, "percentVar"))

# Reuse BEFORE limits if available (ensures identical axes)
pca_limits <- if (file.exists("PCA_limits_before.rds")) readRDS("PCA_limits_before.rds") else NULL
xlim_all <- if (!is.null(pca_limits)) pca_limits$xlim else range(pca_after_df$PC1)
ylim_all <- if (!is.null(pca_limits)) pca_limits$ylim else range(pca_after_df$PC2)

pca_cols <- c(Normal = "#F17BB0", Tumor = "#17BECF")  # same colors as BEFORE

PCA_after_outlier_removal <- ggplot(pca_after_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = pca_cols) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  xlab(paste0("PC1 (", pv_after[1], "%)")) +
  ylab(paste0("PC2 (", pv_after[2], "%)")) +
  ggtitle("PCA (All genes) — AFTER outlier removal") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

print(PCA_after_outlier_removal)
ggsave("QC_PCA_after_outlier_removal.pdf", PCA_after_outlier_removal, width = 7, height = 5)


# Heatmap (AFTER outlier removal, same genes & scale):
mat_after <- assay(vsd_clean)

# 1. Load the BEFORE comparability info (same gene set + same z-limit)
cmp <- readRDS("heatmap_comparability.rds")
top50_ids_before <- cmp$top50_ids
zlim             <- cmp$zlim
breaks           <- seq(-zlim, zlim, length.out = 101)

# 2. Keep the same genes; drop any missing
keep    <- intersect(top50_ids_before, rownames(mat_after))
z_after <- t(scale(t(mat_after[keep, ])))

# 3. Map Ensembl -> symbols (optional, same as BEFORE)
gene_annotation <- read.delim("TCGA_STAD_gene_annotation_v36.tsv", stringsAsFactors = FALSE)
ens_ids <- rownames(z_after)
symb    <- gene_annotation$gene_name[match(ens_ids, gene_annotation$ensembl_gene_id)]
rownames(z_after) <- ifelse(is.na(symb) | symb == "", ens_ids, symb)

# 4. Column annotation & pretty labels for cleaned set
grp_clean   <- droplevels(colData(dds_clean)$condition)
plot_labels <- paste0(as.character(grp_clean), "_",
                      ave(seq_along(grp_clean), grp_clean, FUN = seq_along))
anno <- data.frame(Group = grp_clean)
rownames(anno) <- colnames(z_after)

annotation_colors <- list(Group = c(Normal = "#F17BB0", Tumor = "#17BECF"))

pdf("QC_Heatmap_Top50_after_outlier_removal.pdf", width = 8, height = 10)
pheatmap(
  z_after,
  color              = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(100),
  breaks             = breaks,                  # SAME scale as BEFORE
  annotation_col     = anno,
  annotation_colors  = annotation_colors,
  labels_col         = plot_labels,
  show_rownames      = TRUE,
  fontsize_col       = 6,
  main               = "Top 50 variable genes — AFTER outlier removal (same genes & scale)"
)
dev.off()