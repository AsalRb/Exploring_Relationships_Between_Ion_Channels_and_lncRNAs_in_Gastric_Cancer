
# Part 1:
# ================= Data Acquisition & Assembly: GDC Manifest, Counts Merge, and Tumor/Normal Labeling =================

#----- Preparing manifest file for downloading data -----
# Set the working directory to the project folder (TCGA-STAD gastric cancer analysis)
setwd("C:/Users/Asal/Desktop/Gastric cancer/tcga_stad_project")
getwd()

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

query <- GDCquery(project = "TCGA-STAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

manifest_df <- getManifest(query)
write.table(manifest_df, file = "TCGA_STAD_manifest.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#----- Load Datas into R -----
# Loading required package
library(dplyr)

# List all .tsv files
# 1. Find only real RNA-seq expression files (not parcel logs)
files <- list.files(
  path = ".",
  pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",  # Matches only final expression tables
  recursive = TRUE,
  full.names = TRUE
)

# 2. Checking result
# How many were found
length(files)  # Should match the number of samples (around 300+ for TCGA-STAD)


#----- Read all files and merge them into a matrix -----
# Loading required package
library(data.table)

# 1. Preallocate a list
expression_list <- vector("list", length(files))

# 2. Read files one by one (no metadata, only gene_id and unstranded)
for (i in seq_along(files)) {
  dt <- fread(files[i], header = TRUE, sep = "\t", skip = "#", select = c("gene_id", "unstranded"))
  dt <- dt[!grepl("^N_", gene_id)]
  colname <- tools::file_path_sans_ext(basename(files[i]))
  setnames(dt, "unstranded", colname)
  expression_list[[i]] <- dt
}

# 3. Merge all with Reduce
merged_counts <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = FALSE), expression_list)

# 4. Clean up
merged_counts[, gene_id := sub("\\..*", "", gene_id)]
merged_counts <- merged_counts[!duplicated(gene_id)]
setkey(merged_counts, gene_id)
raw_counts <- as.data.frame(merged_counts)
rownames(raw_counts) <- raw_counts$gene_id
raw_counts <- raw_counts[, -1]

# 5. Checking result
dim(raw_counts)
nrow(raw_counts)  # Should be around ~60,000 genes
ncol(raw_counts)  # Should be 448 
head(raw_counts[, 1:3])  # First few genes across first 3 samples
summary(raw_counts)
any(is.na(raw_counts))          # Should be FALSE
any(!is.finite(as.matrix(raw_counts)))  # Should be FALSE
head(rownames(raw_counts))  # Should return things like "ENSG00000000003"
raw_counts["ENSG00000141510", ]
summary(as.vector(as.matrix(raw_counts)))
sum(rowSums(raw_counts) == 0)  # 2381!
sum(rowSums(raw_counts) < 10)  # 3136!

# 6. Filter out genes with low overall expression (e.g., row sum ≤ 10)
raw_counts_filtered <- raw_counts[rowSums(raw_counts) > 10, ]
# Check dimensions
dim(raw_counts_filtered)

# 7. Overwriting raw_counts with filtered version (genes with rowSums > 0)
raw_counts <- raw_counts_filtered

# 8. Saving raw_counts to file
saveRDS(raw_counts, "raw_counts.rds")

# Load raw counts
raw_counts <- readRDS("raw_counts.rds")


# ----- Tumor vs. Normal Labeling -----
# Loading required packages
library(TCGAbiolinks)
library(dplyr)
library(httr)

# Set SSL + timeout settings (if connection is slow)
set_config(config(ssl_verifypeer = 0L))
set_config(timeout(300))  # increase timeout globally

# 1. Extract sample barcodes from full dataset
sample_barcodes <- colnames(raw_counts)

# 2. Extract UUIDs from barcodes
sample_uuids <- gsub("\\.rna_seq\\.augmented_star_gene_counts.*", "", sample_barcodes)

# 3. Extract results from stored GDC query object
query_results <- query$results[[1]]

# 4. Match UUIDs to file names in GDC query result
expected_file_names <- paste0(sample_uuids, ".rna_seq.augmented_star_gene_counts.tsv")

# Build metadata: UUID, file name, case ID, short barcode
sample_info <- query_results %>%
  filter(file_name %in% expected_file_names) %>%
  select(id, file_name, cases) %>%
  mutate(short_barcode = substr(cases, 1, 16))  # Extract TCGA barcodes

# 5. Retrieve biospecimen metadata (contains sample types)
sample_types <- GDCquery_clinic(project = "TCGA-STAD", type = "biospecimen")

# 6. Merge metadata with sample types
sample_info <- merge(
  sample_info,
  sample_types[, c("submitter_id", "sample_type")],
  by.x = "short_barcode", by.y = "submitter_id",
  all.x = TRUE
)

# 7. Create tumor/normal labels
sample_labels <- ifelse(sample_info$sample_type == "Solid Tissue Normal", "Normal", "Tumor")
sample_labels <- factor(sample_labels, levels = c("Normal", "Tumor"))

# 8. Checking result
# First few genes across first 6 samples
head(sample_info[, 1:5])
# How many samples are tumor/normal
print(table(sample_labels))  # Should show something like: Normal XX, Tumor XXX
# Check all Sample Types
table(sample_info$sample_type)   # Primary Tumor = 412, Solid Tissue Normal = 36
# Check for Duplicates
sum(duplicated(sample_info$short_barcode))  # Should be 0

# 9. Save for later use
saveRDS(sample_info, "sample_metadata.rds")
saveRDS(sample_labels, "sample_labels.rds")

# Load metadata 
sample_info <- readRDS("sample_metadata.rds")
# Load sample labels (Tumor/Normal)
sample_labels <- readRDS("sample_labels.rds")
