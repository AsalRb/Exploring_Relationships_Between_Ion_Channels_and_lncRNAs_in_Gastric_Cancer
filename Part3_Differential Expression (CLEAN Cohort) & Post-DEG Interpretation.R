
# Part 3:
# ================= Differential Expression (CLEAN Cohort) & Post-DEG Interpretation =================

# ----- Differential Expression Analysis (CLEAN cohort) -----
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})

# 1) Load cleaned DESeqDataSet
dds <- readRDS("dds_clean.rds")

# Safety: ensure two groups and "Normal" is the reference
stopifnot("condition" %in% colnames(colData(dds)))
dds$condition <- droplevels(dds$condition)
if (nlevels(dds$condition) < 2) stop("Only one group left; cannot run DEG.")
dds$condition <- relevel(dds$condition, ref = "Normal")

# 2) Run DESeq on the CLEAN data
dds <- DESeq(dds)

# 3) Results: Tumor vs Normal (positive log2FC = higher in Tumor)
alpha <- 0.05
res <- results(dds,
               contrast = c("condition", "Tumor", "Normal"),
               alpha = alpha)
res <- res[order(res$padj), ]

# 4) LFC shrinkage for volcano & ranking
if (requireNamespace("apeglm", quietly = TRUE)) {
  res <- lfcShrink(dds, contrast = c("condition","Tumor","Normal"),
                   type = "apeglm")
}

# 5) Add gene symbols
gene_annotation <- read.delim("TCGA_STAD_gene_annotation_v36.tsv", stringsAsFactors = FALSE)
symb <- gene_annotation$gene_name[match(rownames(res), gene_annotation$ensembl_gene_id)]
res$gene_symbol <- ifelse(is.na(symb) | symb == "", rownames(res), symb)

# 6) Save results (both RDS and CSV)
saveRDS(res, "DEG_Tumor_vs_Normal_clean.rds")
write.csv(as.data.frame(res), "DEG_Tumor_vs_Normal_clean.csv")

# 7) Also save a significant set for downstream plots
sig <- as.data.frame(res) %>%
  filter(!is.na(padj), padj < alpha, abs(log2FoldChange) > 1)
write.table(sig, "significant_genes_clean.tsv",
            sep = "\t", quote = FALSE, row.names = TRUE)

cat("DEG done on CLEAN data.\n",
    "Significant genes (padj<", alpha, " & |LFC|>1): ", nrow(sig), "\n", sep="")

# 8) Load objects
res <- readRDS("DEG_Tumor_vs_Normal_clean.rds")
res_csv <- read.csv("DEG_Tumor_vs_Normal_clean.csv", row.names = 1)
sig <- read.delim("significant_genes_clean.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 9) Quick summary to console
summary(res)

# 10) Create and save DESeq2 MA-plot (log2 fold change vs mean expression)
pdf("MAplot_DESeq2.pdf", width = 7, height = 5)  # open PDF device
plotMA(res, ylim = c(-5, 5))
dev.off()  

# 11) Volcano plot (CLEAN cohort)
volc_df <- as.data.frame(res)
volc_df$sig <- with(volc_df, !is.na(padj) & padj < alpha & abs(log2FoldChange) > 1)

gg_volcano <- ggplot(volc_df,
                     aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.8, na.rm = TRUE) +
  scale_color_manual(values = c(`TRUE` = "#d62728", `FALSE` = "grey70"),
                     labels = c("Not sig.", "padj<0.05 & |LFC|>1"),
                     name = "") +
  xlab("log2 Fold Change (Tumor / Normal)") +
  ylab(expression(-log[10]~padj)) +
  ggtitle("Volcano: Tumor vs Normal (CLEAN cohort)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(gg_volcano)
ggsave("Volcano_Tumor_vs_Normal_CLEAN.pdf", gg_volcano, width = 7, height = 5)


# ----- Count gene biotypes in significant genes -----
library(dplyr)
library(ggplot2)

# 1) Inputs (fresh session):
sig <- read.delim("significant_genes_clean.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Make sure rownames carry the Ensembl IDs for the join:
if (is.null(rownames(sig))) stop("Row names of 'sig' must contain Ensembl IDs.")
gene_annotation <- read.delim("TCGA_STAD_gene_annotation_v36.tsv", stringsAsFactors = FALSE)

# 2) Be robust to different column names in the annotation:
if (!"gene_biotype" %in% names(gene_annotation) && "gene_type" %in% names(gene_annotation)) {
  gene_annotation$gene_biotype <- gene_annotation$gene_type
}
stopifnot(all(c("ensembl_gene_id","gene_biotype") %in% names(gene_annotation)))

# 3) Strip Ensembl version (e.g., ENSG000001.5 -> ENSG000001)
sig$ensembl_gene_id <- sub("\\..*$", "", rownames(sig))

# 4) Join and count
sig_genes_annot <- sig %>%
  left_join(gene_annotation[, c("ensembl_gene_id","gene_biotype")], by = "ensembl_gene_id")

biotype_counts <- sig_genes_annot %>%
  mutate(gene_biotype = ifelse(is.na(gene_biotype), "unknown", gene_biotype)) %>%
  count(gene_biotype, name = "count") %>%
  arrange(desc(count))

print(biotype_counts)

# 5) Save counts table
write.table(biotype_counts, "significant_genes_biotype_counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 6) Plot + save
p_biotype <- ggplot(biotype_counts, aes(x = reorder(gene_biotype, -count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal(base_size = 14) +
  xlab("Gene Biotype") +
  ylab("Number of Significant Genes") +
  ggtitle("Distribution of Gene Biotypes in Significant Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("DEG_BiotypeCounts.pdf", p_biotype, width = 9, height = 5)


# ----- Interpretation Stage (After DEG) -----
# PCA AFTER DEG
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

# Inputs
dds <- readRDS("dds_clean.rds")
sig <- read.delim("significant_genes_clean.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Map sig genes (strip Ensembl versions) to rows in dds
sig_ids_core <- sub("\\..*$", "", rownames(sig))
dds_rows     <- rownames(dds)
dds_core     <- sub("\\..*$", "", dds_rows)
match_idx    <- match(sig_ids_core, dds_core)
sig_dds_rows <- dds_rows[na.omit(match_idx)]
stopifnot(length(sig_dds_rows) > 1)

# VST, subset to sig genes
vsd <- vst(dds, blind = FALSE)
vsd_sig <- assay(vsd)[sig_dds_rows, , drop = FALSE]

# PCA
pca <- prcomp(t(vsd_sig), scale. = TRUE)
grp <- droplevels(colData(dds)$condition)
pv  <- round(summary(pca)$importance[2, 1:2] * 100, 1)

# Reuse the AFTER-QC PCA limits that it's saved earlier
qc_limits_path <- "PCA_limits_before.rds"   
pca_limits <- if (file.exists(qc_limits_path)) readRDS(qc_limits_path) else NULL
xlim_all   <- if (!is.null(pca_limits)) pca_limits$xlim else range(pca$x[,1])
ylim_all   <- if (!is.null(pca_limits)) pca_limits$ylim else range(pca$x[,2])

# Same colors used in QC plots
pca_cols <- c(Normal = "#F17BB0", Tumor = "#17BECF")

p_pca_interp <- ggplot(
  data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], condition = grp),
  aes(PC1, PC2, color = condition)
) +
  geom_point(size = 3) +
  scale_color_manual(values = pca_cols) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +    
  xlab(paste0("PC1 (", pv[1], "%)")) +
  ylab(paste0("PC2 (", pv[2], "%)")) +
  ggtitle("PCA (Significant DEGs)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

print(p_pca_interp)
ggsave("INTERP_PCA_SigDEGs.pdf", p_pca_interp, width = 7, height = 5)


# Heatmap AFTER DEG
suppressPackageStartupMessages({ library(pheatmap); library(dplyr); library(matrixStats) })

dds <- readRDS("dds_clean.rds")
sig <- read.delim("significant_genes_clean.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Map sig genes to dds rownames
sig_ids_core <- sub("\\..*$", "", rownames(sig))
dds_rows     <- rownames(dds)
dds_core     <- sub("\\..*$", "", dds_rows)
sig_dds_rows <- dds_rows[na.omit(match(sig_ids_core, dds_core))]
stopifnot(length(sig_dds_rows) > 0)

# Normalized counts and subset
nc <- counts(dds, normalized = TRUE)[sig_dds_rows, , drop = FALSE]

# Pick top-N most variable sig genes 
N <- min(50, nrow(nc))
top_idx <- head(order(rowVars(as.matrix(nc)), decreasing = TRUE), N)
nc_top  <- nc[top_idx, , drop = FALSE]

# Z-score per gene
z <- t(scale(t(nc_top)))
z[is.na(z)] <- 0

# Reuse the SAME z-limit used in QC heatmap
cmp <- if (file.exists("heatmap_comparability.rds")) readRDS("heatmap_comparability.rds") else NULL
if (!is.null(cmp) && !is.null(cmp$zlim)) {
  zlim <- cmp$zlim
} else {
  # fallback: compute a robust limit here
  zlim <- as.numeric(quantile(abs(z), 0.99, na.rm = TRUE))
  zlim <- max(2.5, zlim)
}
breaks <- seq(-zlim, zlim, length.out = 101)
z[z >  zlim] <-  zlim
z[z < -zlim] <- -zlim

# Gene symbol mapping
if (file.exists("TCGA_STAD_gene_annotation_v36.tsv")) {
  ga <- read.delim("TCGA_STAD_gene_annotation_v36.tsv", stringsAsFactors = FALSE)
  if (!"gene_name" %in% names(ga) && "gene_symbol" %in% names(ga)) ga$gene_name <- ga$gene_symbol
  sym <- ga$gene_name[match(sub("\\..*$","", rownames(z)), ga$ensembl_gene_id)]
  rownames(z) <- ifelse(is.na(sym) | sym == "", rownames(z), sym)
}

# Annotation & labels, same colors as QC
grp <- droplevels(colData(dds)$condition)
lab <- paste0(as.character(grp), "_", ave(seq_along(grp), grp, FUN = seq_along))
anno <- data.frame(Group = grp); rownames(anno) <- colnames(z)
ann_colors <- list(Group = c(Normal = "#F17BB0", Tumor = "#17BECF"))

p <- pheatmap(
  z,
  color              = colorRampPalette(c("#3b4cc0", "white", "#b40426"))(100),
  breaks             = breaks,                   # <- identical color scale to QC
  annotation_col     = anno,
  annotation_colors  = ann_colors,
  labels_col         = lab,
  show_rownames      = TRUE,
  fontsize_col       = 6,
  main               = sprintf("Top %d Significant DEGs", nrow(z))
)
print(p)
ggsave("INTERP_Heatmap_TopSigDEGs.pdf", p, width = 8, height = 10)
