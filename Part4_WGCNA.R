
# Part 4:
# ================= WGCNA on All genes =================

# ----- Input preparation (VST, most variable genes) -----
suppressPackageStartupMessages({
  library(DESeq2); library(WGCNA); library(matrixStats); library(dplyr)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Load inputs
dds_clean <- readRDS("dds_clean.rds")

# VST that respects design
vsd <- vst(dds_clean, blind = FALSE)
expr <- assay(vsd)                        # genes x samples

# Keep top-N most variable genes (adjust N)
N_top <- min(10000, nrow(expr))          # default: up to 10k
rv <- rowVars(expr)
keep <- order(rv, decreasing = TRUE)[seq_len(N_top)]
expr_top <- expr[keep, , drop = FALSE]

# Transpose for WGCNA: samples x genes
datExpr <- t(expr_top)

# Save for reproducibility
saveRDS(datExpr, file = "WGCNA_datExpr_vst_topVar.rds")

# Load 
datExpr <- readRDS("WGCNA_datExpr_vst_topVar.rds")

# ----- Sample clustering sanity check -----
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

# datExpr should be samples × genes (from your VST, top-variable set)
stopifnot(is.matrix(datExpr) || is.data.frame(datExpr))

gsg <- goodSamplesGenes(datExpr, verbose = 3)

cat("\n[goodSamplesGenes] allOK =", gsg$allOK, "\n")
cat("Bad genes   :", sum(!gsg$goodGenes), "\n")
cat("Bad samples :", sum(!gsg$goodSamples), "\n\n")

# Save reports (empty files if nothing is flagged)
bad_genes   <- colnames(datExpr)[!gsg$goodGenes]
bad_samples <- rownames(datExpr)[!gsg$goodSamples]

write.table(data.frame(bad_gene = bad_genes),
            "WGCNA_goodSamplesGenes_bad_genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(bad_sample = bad_samples),
            "WGCNA_goodSamplesGenes_bad_samples.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# NOTE: We are NOT removing anything here because QC already handled outliers.

# Hierarchical clustering of samples for visualization
sampleTree <- hclust(dist(datExpr), method = "average")

pdf("WGCNA_SampleClustering_postQC.pdf", width = 10, height = 5)
plot(sampleTree,
     main = "WGCNA: Sample clustering (post-QC; for visualization only)",
     xlab = "", sub = "", labels = FALSE)
dev.off()

# ----- Choose soft-thresholding power -----
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

# Plot scale-free fit
pdf("WGCNA_SoftThreshold.pdf", width = 10, height = 5)
par(mfrow = c(1,2))

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.9, col="blue")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=powers, col="red")

dev.off()

# ----- Build network & detect modules -----
suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
set.seed(1)

# Inputs
stopifnot(exists("datExpr"))
softPower     <- 10          # chosen from SFT
minModuleSize <- 30          # increase (e.g., 40–50) to avoid tiny modules
mergeCut      <- 0.25        # height threshold to merge similar modules (0.2–0.35 typical)

# 1) Compute TOM directly from expression (signed, bicor)
#   - robust to outliers; uses biweight midcorrelation
TOM <- TOMsimilarityFromExpr(
  datExpr,
  power        = softPower,
  networkType  = "signed",
  corType      = "bicor",
  maxPOutliers = 0.1
)
dissTOM  <- 1 - TOM

# 2) Gene tree (hierarchical clustering on TOM distance)
geneTree <- hclust(as.dist(dissTOM), method = "average")

pdf("WGCNA_GeneDendrogram_PreMerge.pdf", width = 12, height = 6)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering dendrogram (pre-merge)", labels = FALSE)
dev.off()

# 3) Dynamic tree cut to define initial modules
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM  = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize   = minModuleSize
)

moduleColors_before <- labels2colors(dynamicMods)

pdf("WGCNA_Dendrogram_DynamicTreeCut.pdf", width = 12, height = 6)
plotDendroAndColors(
  geneTree,
  moduleColors_before,
  "Dynamic Tree Cut",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
)
dev.off()

# 4) Merge very similar modules (recommended)
# Compute module eigengenes, then merge modules with highly correlated eigengenes
MEs0  <- moduleEigengenes(datExpr, moduleColors_before)$eigengenes
MEs0  <- orderMEs(MEs0)

merge <- mergeCloseModules(
  datExpr,
  moduleColors_before,
  cutHeight = mergeCut,
  verbose = 3
)

moduleColors_after <- merge$colors
MEs_after          <- merge$newMEs

# 5) Document before vs after merge
pdf("WGCNA_Dendrogram_Before_AfterMerge.pdf", width = 12, height = 6)
plotDendroAndColors(
  geneTree,
  cbind(`Before merge` = moduleColors_before, `After merge` = moduleColors_after),
  groupLabels = c("Before merge", "After merge"),
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
)
dev.off()

# Eigengene clustering & correlation heatmap (after merge)
pdf("WGCNA_ModuleEigengeneTree_AfterMerge.pdf", width = 8, height = 5)
MEtree <- hclust(as.dist(1 - cor(MEs_after)), method = "average")
plot(MEtree, main = "Clustering of module eigengenes (after merge)", xlab = "", sub = "")
dev.off()

pdf("WGCNA_MEcorr_Heatmap_AfterMerge.pdf", width = 7, height = 8)
labeledHeatmap(
  Matrix    = cor(MEs_after),
  xLabels   = names(MEs_after),
  yLabels   = names(MEs_after),
  ySymbols  = names(MEs_after),
  colorLabels = FALSE,
  colors    = blueWhiteRed(50),
  zlim      = c(-1, 1),
  main      = "Module eigengene correlation (after merge)"
)
dev.off()

# 6) Save artifacts for downstream steps
tab_before <- sort(table(moduleColors_before), decreasing = TRUE)
tab_after  <- sort(table(moduleColors_after),  decreasing = TRUE)

write.table(data.frame(module = names(tab_before), size = as.integer(tab_before)),
            "WGCNA_ModuleSizes_PreMerge.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(module = names(tab_after),  size = as.integer(tab_after)),
            "WGCNA_ModuleSizes_AfterMerge.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Final objects to use later
moduleColors <- moduleColors_after
MEs          <- MEs_after

saveRDS(list(
  moduleColors_before = moduleColors_before,
  moduleColors_after  = moduleColors_after,
  MEs_after           = MEs_after,
  geneTree            = geneTree,
  softPower           = softPower,
  minModuleSize       = minModuleSize,
  mergeCutHeight      = mergeCut
), "WGCNA_modules.rds")

# Also export per-gene color table
write.table(data.frame(Gene = colnames(datExpr), ModuleColor = moduleColors),
            "WGCNA_Gene2Module.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ----- Module–trait correlations (Tumor vs Normal) -----
suppressPackageStartupMessages(library(WGCNA))
library(dplyr)

# Use the merged eigengenes from Step 4
# MEs should be a data.frame with rows = samples, cols = ME<color>
stopifnot(exists("MEs"), exists("datExpr"))

# Trait: 0 = Normal, 1 = Tumor
grp <- droplevels(colData(dds_clean)$condition)
trait_df <- data.frame(Condition = as.numeric(grp == "Tumor"),
                       row.names = colnames(dds_clean))

# Align rows of traits and eigengenes to datExpr (samples)
common_samples <- intersect(rownames(trait_df), rownames(datExpr))
MEs_aligned    <- MEs[common_samples, , drop = FALSE]
traits_aligned <- trait_df[common_samples, , drop = FALSE]

# ensure ME columns are ordered canonically
MEs_aligned <- orderMEs(MEs_aligned)

# Correlations (robust bicor) + p-values
moduleTraitCor <- bicor(MEs_aligned, traits_aligned,
                        use = "pairwise.complete.obs", maxPOutliers = 0.1)
moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples = nrow(MEs_aligned))

# Summary table (with FDR)
MT_df <- data.frame(
  Module        = rownames(moduleTraitCor),
  Cor_Condition = moduleTraitCor[, "Condition"],
  P_Condition   = moduleTraitP[, "Condition"],
  FDR_Condition = p.adjust(moduleTraitP[, "Condition"], method = "BH"),
  stringsAsFactors = FALSE
) %>% arrange(desc(abs(Cor_Condition)))

write.table(MT_df, "WGCNA_ModuleTrait_Tumor.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap (values + p-values)
txt <- paste0(signif(moduleTraitCor, 2), "\n(",
              signif(moduleTraitP, 1), ")")

pdf("WGCNA_ModuleTrait_Heatmap.pdf", width = 7, height = 8)
labeledHeatmap(
  Matrix       = moduleTraitCor,
  xLabels      = colnames(traits_aligned),
  yLabels      = rownames(moduleTraitCor),
  ySymbols     = rownames(moduleTraitCor),
  colors       = blueWhiteRed(50),
  colorLabels  = FALSE,
  textMatrix   = txt,
  cex.text     = 0.7,
  zlim         = c(-1, 1),
  main         = "Module–Trait Relationships (bicor)"
)
dev.off()

#  ----- Relating WGCNA modules to clinical traits and survival -----
# Download & tidy TCGA clinical traits:
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
  library(stringr)
})

# 1) Pull TCGA-STAD clinical (patient-level)
clin_raw <- GDCquery_clinic(project = "TCGA-STAD", type = "clinical")

# 2) Select traits of interest (edit as you like)
clin <- clin_raw %>%
  transmute(
    patient_id          = submitter_id,                 # e.g., TCGA-XX-XXXX
    age_at_diagnosis    = as.numeric(age_at_diagnosis),
    gender              = ifelse(tolower(gender) %in% c("male","female"), tolower(gender), NA),
    ajcc_pathologic_stage = tolower(ajcc_pathologic_stage),  # e.g., stage iiia
    vital_status        = tolower(vital_status),        # alive / dead
    days_to_death       = as.numeric(days_to_death),
    days_to_last_follow_up = as.numeric(days_to_last_follow_up),
    # add other fields if present for your cohort:
    # tumor_grade = tolower(tumor_grade),
    # residual_disease = tolower(residual_disease),
    # ...
  )

# 3) Harmonize stage to I/II/III/IV where possible
stage_collapse <- function(x) {
  x <- tolower(x)
  x <- gsub("^stage[[:space:]]*", "", x)
  x <- gsub("[abcd]+$", "", x)      # drop sub-stages (iia → ii)
  x <- trimws(x)
  x[!x %in% c("i","ii","iii","iv")] <- NA
  x
}
clin$stage_simple <- stage_collapse(clin$ajcc_pathologic_stage)

# 4) Build a sample-level table keyed to your expression samples
# Your columns are TCGA *sample* barcodes (length 16); map → patient (length 12)
sample_barcodes <- colnames(dds_clean)
patient_ids <- substr(sample_barcodes, 1, 12)

sample_map <- data.frame(
  sample_barcode = sample_barcodes,
  patient_id     = patient_ids,
  stringsAsFactors = FALSE
)

# 5) Merge patient clinical → sample level
traits_sample <- sample_map %>%
  left_join(clin, by = "patient_id")
rownames(traits_sample) <- traits_sample$sample_barcode

# Build a clean traits matrix aligned to MEs / datExpr:
# We will create numeric columns suitable for correlation:
# - numeric traits kept as is (e.g., age)
# - binary factors → 0/1
# - multi-level factors (e.g., stage I/II/III/IV) → one-hot dummies

# 5.1 Start with age
trait_df <- data.frame(row.names = rownames(traits_sample))
if ("age_at_diagnosis" %in% names(traits_sample))
  trait_df$Age <- traits_sample$age_at_diagnosis

# 5.2 Gender (binary 0/1 if both present)
if ("gender" %in% names(traits_sample)) {
  g <- traits_sample$gender
  if (any(!is.na(g))) {
    # Make 'female'=0, 'male'=1 by convention (change if you prefer)
    trait_df$Male <- ifelse(g == "male", 1L, ifelse(g == "female", 0L, NA))
  }
}

# 5.3 SAFE one-hot encoding for stage (I/II/III/IV)
if ("stage_simple" %in% names(traits_sample)) {
  # Normalize and lock levels
  stage_levels <- c("i","ii","iii","iv")
  st <- factor(traits_sample$stage_simple, levels = stage_levels)
  
  # Preallocate (nSamples x 4)
  mm <- matrix(NA_integer_, nrow = nrow(traits_sample), ncol = length(stage_levels))
  colnames(mm) <- c("Stage_I","Stage_II","Stage_III","Stage_IV")
  rownames(mm) <- rownames(traits_sample)
  
  # Fill one-hot: 1 for the matching stage, 0 for other stages, NA if stage is NA
  for (k in seq_along(stage_levels)) {
    col <- as.integer(st == stage_levels[k])
    # st==... yields NA where st is NA; change FALSE(0)/TRUE(1) but keep NAs
    mm[, k] <- ifelse(is.na(col), NA_integer_, col)
  }
  
  # Bind to your trait_df (which already has the same rownames)
  trait_df <- cbind(trait_df, mm[rownames(trait_df), , drop = FALSE])
}

stopifnot(nrow(trait_df) == nrow(traits_sample))
stopifnot(identical(rownames(trait_df), rownames(traits_sample)))
colSums(is.na(trait_df))  # just to see missingness

# 5.4 Condition (Tumor/Normal) from your clean DESeq2 object (already used before)
cond <- droplevels(colData(dds_clean)$condition)
trait_df$Condition <- as.numeric(cond == "Tumor")  # 0=Normal, 1=Tumor
rownames(trait_df) <- colnames(dds_clean)

# 5.5 Align to datExpr / MEs (rows = samples)
common_samples <- Reduce(intersect, list(rownames(trait_df), rownames(datExpr), rownames(MEs)))
trait_df <- trait_df[common_samples, , drop = FALSE]
MEs_aligned <- MEs[common_samples, , drop = FALSE]

# 6) Module–trait correlations (auto: Pearson for binary/zero-MAD; bicor otherwise)
suppressPackageStartupMessages(library(WGCNA))

me_names <- colnames(MEs_aligned)
tr_names <- colnames(trait_df)

moduleTraitCor <- matrix(NA_real_, nrow = ncol(MEs_aligned), ncol = ncol(trait_df),
                         dimnames = list(me_names, tr_names))
moduleTraitP   <- moduleTraitCor

# 6.1 Correlate each trait with module eigengenes
for (j in seq_len(ncol(trait_df))) {
  y <- trait_df[, j]
  # Use bicor only if trait has non-zero MAD and >2 unique numeric values
  use_bicor <- is.numeric(y) && is.finite(mad(y, na.rm = TRUE)) &&
    mad(y, na.rm = TRUE) > 0 && length(unique(na.omit(y))) > 2
  if (use_bicor) {
    cmat <- bicor(MEs_aligned, y, use = "pairwise.complete.obs", maxPOutliers = 0.1)
  } else {
    cmat <- cor(MEs_aligned, y, use = "pairwise.complete.obs", method = "pearson")
  }
  moduleTraitCor[, j] <- cmat[, 1]
  moduleTraitP[, j]   <- corPvalueStudent(moduleTraitCor[, j], nSamples = nrow(MEs_aligned))
}

# 6.2 FDR per trait
moduleTraitFDR <- apply(moduleTraitP, 2, function(p) p.adjust(p, method = "BH"))

# 6.3 Record N per trait (pairwise complete obs)
trait_N <- sapply(trait_df, function(y) sum(complete.cases(MEs_aligned[,1,drop=FALSE], y)))
write.table(data.frame(Trait = names(trait_N), N = as.integer(trait_N)),
            "WGCNA_ModuleTrait_TraitNs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 6.4 Tidy long table (with correct per-trait FDR)
corr_tbl <- do.call(rbind, lapply(seq_len(nrow(moduleTraitCor)), function(i) {
  data.frame(
    Module = rownames(moduleTraitCor)[i],
    Trait  = colnames(moduleTraitCor),
    Cor    = as.numeric(moduleTraitCor[i, ]),
    P      = as.numeric(moduleTraitP[i, ]),
    stringsAsFactors = FALSE
  )
})) |>
  dplyr::group_by(Trait) |>
  dplyr::mutate(FDR = p.adjust(P, method = "BH")) |>
  dplyr::arrange(dplyr::desc(abs(Cor)), .by_group = TRUE) |>
  dplyr::ungroup()

write.table(corr_tbl, "WGCNA_ModuleTrait_AllTraits.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 6.5 Top-10 modules per trait by |cor| (with FDR)
top_by_trait <- corr_tbl |>
  dplyr::group_by(Trait) |>
  dplyr::slice_max(order_by = abs(Cor), n = 10, with_ties = FALSE) |>
  dplyr::ungroup()
write.table(top_by_trait, "WGCNA_ModuleTrait_Top10PerTrait.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 6.6 Heatmap (correlation with P in parentheses)
txt <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitP, 1), ")")
pdf("WGCNA_ModuleTrait_Heatmap_AllTraits.pdf", width = 8, height = 0.45*nrow(moduleTraitCor) + 4)
labeledHeatmap(
  Matrix      = moduleTraitCor,
  xLabels     = colnames(trait_df),
  yLabels     = rownames(moduleTraitCor),
  ySymbols    = rownames(moduleTraitCor),
  colors      = blueWhiteRed(50),
  colorLabels = FALSE,
  textMatrix  = txt,
  cex.text    = 0.7,
  zlim        = c(-1, 1),
  main        = "Module–Trait Relationships"
)
dev.off()


# 7) Overall survival (OS) association per module (univariate Cox)
suppressPackageStartupMessages({
  library(survival); library(survminer)
  library(dplyr)
})

# 7.1 Build OS time/event at *patient* level
os_df <- clin %>%
  transmute(
    patient_id,
    time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    event = ifelse(!is.na(days_to_death) & days_to_death > 0, 1L, 0L)
  )

# 7.2 Map OS to samples (if multiple samples per patient, keep Tumor first; otherwise first)
sample_os <- sample_map %>%
  left_join(os_df, by = "patient_id") %>%
  mutate(
    condition = as.character(colData(dds_clean)$condition[match(sample_barcode, colnames(dds_clean))])
  ) %>%
  arrange(patient_id, desc(condition == "Tumor")) %>%
  group_by(patient_id) %>%
  slice(1) %>%                   # pick one sample per patient
  ungroup()

sample_os <- as.data.frame(sample_os)
rownames(sample_os) <- sample_os$sample_barcode

# 7.3 Align with MEs_aligned
keep <- intersect(rownames(MEs_aligned), rownames(sample_os))
MEs_surv <- MEs_aligned[keep, , drop = FALSE]
surv_time  <- sample_os[keep, "time",  drop = TRUE]
surv_event <- sample_os[keep, "event", drop = TRUE]

# 7.4 Remove NAs
ok <- is.finite(surv_time) & !is.na(surv_event)
MEs_surv   <- MEs_surv[ok, , drop = FALSE]
surv_time  <- surv_time[ok]
surv_event <- surv_event[ok]

# 7.5 Cox per module
cox_results <- lapply(colnames(MEs_surv), function(me) {
  df <- data.frame(time = surv_time, event = surv_event, ME = MEs_surv[, me])
  fit <- try(coxph(Surv(time, event) ~ ME, data = df), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  sm <- summary(fit)
  data.frame(
    Module = me,
    HR     = exp(coef(fit)),
    HR_low = exp(confint(fit))[1],
    HR_up  = exp(confint(fit))[2],
    Wald_p = sm$wald["pvalue"],
    stringsAsFactors = FALSE
  )
})

# 7.6 Save Cox results (+FDR)
cox_tbl <- do.call(rbind, cox_results)
if (!is.null(cox_tbl)) {
  cox_tbl$FDR <- p.adjust(cox_tbl$Wald_p, method = "BH")
  write.table(cox_tbl, "WGCNA_Module_OS_Cox.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
}


# ----- Focus modules → lncRNAs & ion channels -----

suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(readr)
})

options(stringsAsFactors = FALSE)

# 1) Choose modules of interest 
# pick by hand
modules_of_interest <- c()  

# 2) Prepare per-gene table (annotation + module + hub metrics) 
# Quick pre-checks
stopifnot(nrow(MEs) == nrow(datExpr))
stopifnot(all(rownames(MEs) == rownames(datExpr)))
if (is.null(names(moduleColors))) names(moduleColors) <- colnames(datExpr)
stopifnot(length(moduleColors) == ncol(datExpr))

# Align object dimensions
stopifnot(ncol(datExpr) == length(moduleColors))
stopifnot(all(rownames(MEs) == rownames(datExpr)))  # samples aligned

# kME (module membership): correlation of gene with module eigengenes
kME_all <- signedKME(datExpr, MEs, outputColumnName = "kME")

# Gene significance (GS) for Condition (Tumor=1 vs Normal=0)
cond_vec <- trait_df[rownames(datExpr), "Condition"]
GS_cond  <- as.numeric(cor(datExpr, cond_vec, use = "pairwise.complete.obs", method = "pearson"))

# Build a master gene table
genes_df <- data.frame(
  EnsemblID = colnames(datExpr),
  Module    = moduleColors,
  GS_Condition = GS_cond,
  stringsAsFactors = FALSE
)

# Add per-module kME (columns kMEturquoise, kMEblue, ...)
genes_df <- cbind(genes_df, kME_all[genes_df$EnsemblID, , drop = FALSE])

# Attach symbols/biotypes from your annotation
# (annotation Ensembl IDs should be version-stripped as you already did)
ann <- gene_annotation[, c("ensembl_gene_id","gene_name","gene_biotype")]
colnames(ann) <- c("EnsemblID","Symbol","Biotype")
genes_df <- genes_df %>%
  left_join(ann, by = "EnsemblID") %>%
  relocate(Symbol, Biotype, .after = EnsemblID)

# Save full per-gene hub table
write.table(genes_df, "WGCNA_AllGenes_with_kME_GS.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 3) Enrichment filter for ion channels / lncRNAs 
suppressPackageStartupMessages({ library(dplyr); library(readr) })

lnc_like <- c("lncRNA","lincRNA","antisense","processed_transcript",
              "sense_intronic","sense_overlapping","3prime_overlapping_ncRNA",
              "bidirectional_promoter_lncRNA","non_coding")

genes_df$is_lncRNA <- ifelse(tolower(genes_df$Biotype) %in% tolower(lnc_like), 1L, 0L)

ion_file <- "ion_channel_genes.txt"
if (file.exists(ion_file)) {
  ion_syms <- unique(trimws(read_lines(ion_file)))
  genes_df$is_ionChannel <- ifelse(!is.na(genes_df$Symbol) & genes_df$Symbol %in% ion_syms, 1L, 0L)
} else {
  patt <- "^(KCN|KCNA|KCNB|KCNC|KCND|KCNE|KCNF|KCNG|KCNH|KCNJ|KCNK|KCNQ|KCNS|SCN|CACN|CACNA|CACNB|CACNG|TRP|TRPC|TRPM|TRPV|TRPA|TRPN|CLCN|HCN|GRIN|GRIA|GRIK|GRM)"
  genes_df$is_ionChannel <- ifelse(!is.na(genes_df$Symbol) & grepl(patt, toupper(genes_df$Symbol)), 1L, 0L)
}

total_genes <- nrow(genes_df)
total_ion   <- sum(genes_df$is_ionChannel)

mod_stats <- genes_df %>%
  group_by(Module) %>%
  summarise(
    n_genes     = n(),
    ion_hits    = sum(is_ionChannel),
    lnc_hits    = sum(is_lncRNA),
    lnc_frac    = lnc_hits / n_genes,
    .groups = "drop"
  ) %>%
  mutate(
    ion_enrich_p = phyper(q = ion_hits - 1, m = total_ion, n = total_genes - total_ion,
                          k = n_genes, lower.tail = FALSE),
    ion_enrich_fdr = p.adjust(ion_enrich_p, method = "BH")
  ) %>%
  arrange(ion_enrich_fdr, desc(lnc_frac))

write.table(mod_stats, "WGCNA_Module_IonChannel_LncRNA_Summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Pick modules of interest: combine clinical & content-based relevance
clin_hits <- c("turquoise", "cyan")   # from module–trait heatmap
if (file.exists("WGCNA_Module_OS_Cox.tsv")) {
  cox_tbl <- read.delim("WGCNA_Module_OS_Cox.tsv")
  os_sig  <- subset(cox_tbl, FDR < 0.10)$Module
  clin_hits <- unique(c(clin_hits, sub("^ME","", os_sig)))
}

selected <- mod_stats %>%
  filter(Module %in% clin_hits) %>%
  filter(ion_enrich_fdr < 0.25 | lnc_frac >= 0.20) %>%
  arrange(ion_enrich_fdr, desc(lnc_frac)) %>%
  pull(Module)

if (length(selected) == 0) selected <- c("turquoise", "cyan")

modules_of_interest <- unique(selected)
cat("Modules_of_interest:", paste(modules_of_interest, collapse = ", "), "\n")

# 4) Extract tables for selected modules
# Quick pre-checks
stopifnot(exists("genes_df"), exists("kME_all"), exists("modules_of_interest"))
stopifnot(all(c("EnsemblID","Module") %in% names(genes_df)))
stopifnot(all(paste0("kME", unique(genes_df$Module)) %in% colnames(kME_all)))

sel_genes <- genes_df %>% filter(Module %in% modules_of_interest)

# For each selected module, compute per-gene module membership (kME for THAT module)
# Match the kME column name pattern: kME<moduleColor>
kme_colnames <- colnames(kME_all)
mod_kme_col <- function(color) {
  # usually "kMEturquoise", "kMEblue", etc.
  ix <- grep(paste0("^kME", color, "$"), kme_colnames, ignore.case = FALSE)
  if (length(ix) == 1) kme_colnames[ix] else NA_character_
}

sel_genes <- sel_genes %>%
  rowwise() %>%
  mutate(kME_inModule = {
    col <- mod_kme_col(Module)
    if (is.na(col)) NA_real_ else kME_all[EnsemblID, col]
  }) %>%
  ungroup()

# Save per-module gene lists with hub metrics
for (mod in modules_of_interest) {
  out <- sel_genes %>% filter(Module == mod) %>%
    arrange(desc(abs(kME_inModule)))
  fn <- paste0("WGCNA_ModuleGenes_", mod, ".tsv")
  write.table(out, fn, sep = "\t", quote = FALSE, row.names = FALSE)
}

# 5) Build lncRNA ↔ ion-channel edge list (within-module correlations)
# Quick pre-check
ann <- read.delim("TCGA_STAD_gene_annotation_v36.tsv", stringsAsFactors = FALSE)
colnames(ann) <- c("EnsemblID","Symbol","Biotype")

# Parameters for edge extraction
edge_cor_method <- "pearson"   # or "bicor" if you prefer (use WGCNA::bicor)
edge_min_absCor <- 0.40        # threshold for absolute correlation
edge_p_adj      <- 0.05        # BH-FDR cutoff (optional)

all_edges <- list()

for (mod in modules_of_interest) {
  mod_ids <- sel_genes$EnsemblID[sel_genes$Module == mod]
  expr_mod <- datExpr[, mod_ids, drop = FALSE]   # samples × genes in module
  
  # Partition lncRNAs and ion channels within this module
  is_lnc   <- sel_genes$is_lncRNA[sel_genes$Module == mod] == 1
  is_ion   <- sel_genes$is_ionChannel[sel_genes$Module == mod] == 1
  
  lnc_ids  <- mod_ids[is_lnc]
  ion_ids  <- mod_ids[is_ion]
  
  if (length(lnc_ids) == 0 || length(ion_ids) == 0) next
  
  # Compute pairwise correlations (lncRNAs × ion channels)
  E_lnc <- expr_mod[, lnc_ids, drop = FALSE]
  E_ion <- expr_mod[, ion_ids, drop = FALSE]
  
  if (edge_cor_method == "pearson") {
    C <- cor(E_lnc, E_ion, use = "pairwise.complete.obs", method = "pearson")
    # p-values (Student) per pair
    nS <- nrow(datExpr)
    P <- corPvalueStudent(C, nS)
  } else {
    C <- bicor(E_lnc, E_ion, use = "pairwise.complete.obs", maxPOutliers = 0.1)
    nS <- nrow(datExpr)
    P <- corPvalueStudent(C, nS)
  }
  
  # Multiple-testing (BH) across all pairs in this module
  Padj <- matrix(p.adjust(as.vector(P), method = "BH"), nrow = nrow(P), ncol = ncol(P),
                 dimnames = dimnames(P))
  
  # Keep edges passing thresholds
  keep <- which(abs(C) >= edge_min_absCor & Padj <= edge_p_adj, arr.ind = TRUE)
  if (nrow(keep) > 0) {
    edges_mod <- data.frame(
      Module = mod,
      lnc_Ensembl = rownames(C)[keep[,1]],
      ion_Ensembl = colnames(C)[keep[,2]],
      Cor        = as.numeric(C[keep]),
      P          = as.numeric(P[keep]),
      FDR        = as.numeric(Padj[keep]),
      stringsAsFactors = FALSE
    )
    # Add symbols for readability
    sym_map <- ann[, c("EnsemblID","Symbol")]
    edges_mod <- edges_mod %>%
      left_join(sym_map, by = c("lnc_Ensembl" = "EnsemblID")) %>%
      rename(lnc_Symbol = Symbol) %>%
      left_join(sym_map, by = c("ion_Ensembl" = "EnsemblID")) %>%
      rename(ion_Symbol = Symbol)
    all_edges[[mod]] <- edges_mod
  }
}

edge_tbl <- if (length(all_edges)) bind_rows(all_edges) else data.frame()
if (nrow(edge_tbl)) {
  write.table(edge_tbl, "WGCNA_lncRNA_IonChannel_Edges.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved lncRNA↔ion-channel edges:", nrow(edge_tbl), "\n")
} else {
  cat("No lncRNA↔ion-channel edges passed the thresholds in selected modules.\n")
}

# 6) Preparing Clean Cytoscape Input Tables from WGCNA Results 
# (Deduplication, Annotation, and Hub Identification)
suppressPackageStartupMessages(library(dplyr))

# 6.a Harmonize Ensembl IDs (drop version) everywhere 
genes_df$EnsemblID <- sub("\\..*$", "", genes_df$EnsemblID)

# If your edge table still has versions, strip them too
edge_tbl$lnc_Ensembl <- sub("\\..*$", "", edge_tbl$lnc_Ensembl)
edge_tbl$ion_Ensembl <- sub("\\..*$", "", edge_tbl$ion_Ensembl)

# Prepare a clean 1:1 annotation map: EnsemblID -> Symbol_fallback
ann$EnsemblID <- sub("\\..*$", "", ann$ensembl_gene_id)
ann$Symbol    <- ann$gene_name

ann2 <- ann %>%
  filter(!is.na(EnsemblID), EnsemblID != "") %>%
  group_by(EnsemblID) %>%
  summarize(
    Symbol_fallback = {
      x <- Symbol[!is.na(Symbol) & Symbol != ""]
      if (length(x) == 0) NA_character_
      else names(sort(table(x), decreasing = TRUE))[1]
    },
    .groups = "drop"
  )
stopifnot(!anyDuplicated(ann2$EnsemblID))

# 6.b Deduplicate genes_df to 1 row per EnsemblID 
# Strategy:
#  - compute a per-row "own-module kME" if available (kME<color>)
#  - if that column is missing, use max |kME*| across all kME columns
#  - within each EnsemblID, keep the row with the largest |own kME|

kme_cols <- grep("^kME", colnames(genes_df), value = TRUE)

# helper to grab own-module kME for a row
own_kME <- function(module, row_idx) {
  col_name <- paste0("kME", as.character(module[row_idx]))
  if (col_name %in% kme_cols) {
    return(genes_df[row_idx, col_name, drop = TRUE])
  } else {
    # fallback: max |kME*| across all modules
    vals <- as.numeric(genes_df[row_idx, kme_cols, drop = TRUE])
    return(if (length(vals)) vals[which.max(abs(vals))] else NA_real_)
  }
}

genes_df$kME_inOwn_tmp <- vapply(seq_len(nrow(genes_df)),
                                 function(i) own_kME(genes_df$Module, i),
                                 numeric(1))

# Collapse duplicates: keep the row with max |kME_inOwn_tmp|
genes_df_dedup <- genes_df %>%
  group_by(EnsemblID) %>%
  arrange(desc(abs(kME_inOwn_tmp)), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Fill missing symbols from annotation fallback (1:1 map)
genes_df_dedup <- genes_df_dedup %>%
  left_join(ann2, by = "EnsemblID") %>%
  mutate(Symbol = dplyr::coalesce(Symbol, Symbol_fallback)) %>%
  select(-Symbol_fallback, -kME_inOwn_tmp)

# Sanity checks
stopifnot(!anyDuplicated(genes_df_dedup$EnsemblID))
cat("genes_df rows (before -> after):", nrow(genes_df), "->", nrow(genes_df_dedup), "\n")
cat("Remaining NA symbols:", sum(is.na(genes_df_dedup$Symbol) | genes_df_dedup$Symbol == ""), "\n")

# 6.c Build the Cytoscape tables using the de-duplicated genes_df 

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

# 0. Preconditions & ID Harmonization for Cytoscape Build (genes_df_dedup + edge_tbl)
# - genes_df_dedup: per-gene table (EnsemblID, Symbol, Biotype, Module, GS_Condition, kME* columns)
# - edge_tbl: lncRNA↔ion-channel edges (lnc_Ensembl, ion_Ensembl, Cor, P, FDR, Module, lnc_Symbol, ion_Symbol)

stopifnot(exists("genes_df_dedup"), exists("edge_tbl"))

# Harmonize Ensembl IDs (strip versions) in both objects
genes_df <- genes_df_dedup %>%
  mutate(EnsemblID = sub("\\..*$", "", EnsemblID))

edge_tbl <- edge_tbl %>%
  mutate(
    lnc_Ensembl = sub("\\..*$", "", lnc_Ensembl),
    ion_Ensembl = sub("\\..*$", "", ion_Ensembl)
  )

# Define lncRNA biotypes in LOWERCASE to match tolower(Biotype)
lnc_set <- tolower(c(
  "lncRNA","lincRNA","antisense","processed_transcript",
  "sense_intronic","sense_overlapping","3prime_overlapping_ncRNA",
  "bidirectional_promoter_lncRNA","non_coding"
))

# 1. Node table: only genes that appear in edges
nodes_in_edges <- unique(c(edge_tbl$lnc_Ensembl, edge_tbl$ion_Ensembl))

# helper: fetch kME for the row's own module; fallback to max-|kME*|
kme_cols <- grep("^kME", colnames(genes_df), value = TRUE)

node_df <- genes_df %>%
  filter(EnsemblID %in% nodes_in_edges) %>%
  rowwise() %>%
  transmute(
    id       = EnsemblID,            # Cytoscape node key
    symbol   = Symbol,
    biotype  = Biotype,
    module   = Module,
    is_lncRNA = as.integer(tolower(Biotype) %in% lnc_set),
    is_ion    = as.integer(!is.na(Symbol) & grepl(
      "^(KCN|KCNA|KCNB|KCNC|KCND|KCNE|KCNF|KCNG|KCNH|KCNJ|KCNK|KCNQ|KCNS|SCN|CACN|CACNA|CACNB|CACNG|TRP|TRPC|TRPM|TRPV|TRPA|TRPN|CLCN|HCN|GRIN|GRIA|GRIK|GRM)",
      toupper(Symbol)
    )),
    GS_Condition,
    kME_inModule = {
      col <- paste0("kME", Module)
      if (col %in% kme_cols) {
        genes_df[[col]][match(EnsemblID, genes_df$EnsemblID)]
      } else if (length(kme_cols)) {
        vals <- as.numeric(genes_df[match(EnsemblID, genes_df$EnsemblID), kme_cols, drop = TRUE])
        vals[which.max(abs(vals))]
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()

# 2. Edge table: Cytoscape-friendly columns
edge_out <- edge_tbl %>%
  transmute(
    source        = lnc_Ensembl,
    target        = ion_Ensembl,
    weight        = Cor,
    pval          = P,
    fdr           = FDR,
    module        = Module,
    source_symbol = lnc_Symbol,
    target_symbol = ion_Symbol
  )

# 3. Write files
write.table(node_df, "Cytoscape_nodes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(edge_out, "Cytoscape_edges.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 4. Quick hubs (degree per node)
deg_tbl <- edge_out %>%
  select(source, target) %>%
  pivot_longer(cols = everything(), values_to = "id") %>%
  count(id, name = "degree")

hub_tbl <- node_df %>%
  left_join(deg_tbl, by = "id") %>%
  mutate(degree = ifelse(is.na(degree), 0L, degree)) %>%
  arrange(desc(degree))

top_lnc <- hub_tbl %>% filter(is_lncRNA == 1) %>% slice_head(n = 30)
top_ion <- hub_tbl %>% filter(is_ion == 1)    %>% slice_head(n = 30)

write.table(hub_tbl, "WGCNA_lncIon_HubSummary_All.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_lnc, "WGCNA_lncIon_TopLncRNAHubs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_ion, "WGCNA_lncIon_TopIonHubs.tsv",     sep = "\t", quote = FALSE, row.names = FALSE)

# 5. Final console check
cat("n(nodes) =", nrow(node_df), " | n(edges) =", nrow(edge_out), "\n")
cat("lncRNA nodes =", sum(node_df$is_lncRNA), " | ion-channel nodes =", sum(node_df$is_ion), "\n")

# Quick checks
# Every edge endpoint appears in node table?
all(edge_out$source %in% node_df$id) && all(edge_out$target %in% node_df$id)

# Any zero-degree nodes? (should be none because nodes were taken from edges)
sum(!(node_df$id %in% unique(c(edge_out$source, edge_out$target))))

# lncRNA/ion breakdown by module (handy before Cytoscape)
table(node_df$module, node_df$is_lncRNA)
table(node_df$module, node_df$is_ion)


# 7) Creating cytoscape for all genes of turquoise module
# Install once: BiocManager::install("RCy3")
suppressPackageStartupMessages(library(RCy3))

# Connect to Cytoscape (make sure Cytoscape app is running)
cytoscapePing()

# 7.a Load your exported tables 
nodes <- read.delim("Cytoscape_nodes.tsv", sep="\t", header=TRUE, check.names=FALSE)
edges <- read.delim("Cytoscape_edges.tsv", sep="\t", header=TRUE, check.names=FALSE)

# 7.b Focus on one module (e.g., turquoise)
mod <- "turquoise"
nodes_mod <- subset(nodes, module == mod)
edges_mod <- subset(edges, module == mod)

# 7.c Create network from data frames
# just turquoise
createNetworkFromDataFrames(nodes_mod, edges_mod,
                            title=paste("lnc↔ion (", mod, ")", sep=""),
                            collection="WGCNA")

# 7.d Style the turquoise network you just created
suppressPackageStartupMessages({
  library(RCy3)
  library(RColorBrewer)
})

# 0. Make sure we're styling the network you just created
net.suid <- getNetworkSuid()
setCurrentNetwork(net.suid)

# 1. (Optional) add a single "type" column so shapes don't collide
#    lncRNA=diamond, ion=triangle, other=ellipse
nodes_mod$type <- ifelse(nodes_mod$is_lncRNA==1 & nodes_mod$is_ion==0, "lncRNA",
                         ifelse(nodes_mod$is_ion==1 & nodes_mod$is_lncRNA==0, "ion_channel", "other"))
# Push the new column to the Cytoscape node table (key column in Cytoscape is "name")
loadTableData(nodes_mod[, c("id","type","kME_inModule","GS_Condition","symbol")],
              data.key.column="id", table.key.column="name")

# 2. Create & apply a clean visual style
style.name <- "WGCNA_Turquoise"
if (!(style.name %in% getVisualStyleNames())) {
  createVisualStyle(style.name,
                    defaults = list(
                      "NODE_FILL_COLOR"        = "#bdbdbd",
                      "NODE_SHAPE"             = "ELLIPSE",
                      "NODE_SIZE"              = 30,
                      "NODE_BORDER_WIDTH"      = 0.5,
                      "NODE_TRANSPARENCY"      = 220,
                      "NODE_LABEL_COLOR"       = "#222222",
                      "NODE_LABEL_FONT_SIZE"   = 12,
                      "EDGE_TRANSPARENCY"      = 160,
                      "EDGE_STROKE_UNSELECTED_PAINT" = "#999999",
                      "EDGE_WIDTH"             = 2
                    ))
}
setVisualStyle(style.name)

# 3. Node label = gene symbol (passthrough)
mapVisualProperty("NODE_LABEL", "symbol", "p")

# 4. Make sure we're styling your custom style
style.name <- "WGCNA_Turquoise"
setVisualStyle(style.name)

# 5. Node color by type
setNodeColorMapping(
  table.column        = "type",
  table.column.values = c("lncRNA","ion_channel","other"),
  colors              = c("#8da0cb", "#fc8d62", "#bdbdbd"),
  mapping.type        = "d",
  style.name          = style.name
)

# 6. Quick check
getCurrentStyle()

tbl <- getTableColumns(table = "node")  # Cytoscape node table
cat("Distinct 'type' values in Cytoscape:\n")
print(sort(unique(tbl$type)))

# 7. Node shape by type
setNodeShapeMapping(
  table.column        = "type",
  table.column.values = c("lncRNA","ion_channel","other"),
  shapes              = c("DIAMOND","TRIANGLE","ELLIPSE"),
  style.name          = style.name
)

# 8. Node size by kME
setNodeSizeMapping(
  table.column        = "kME_inModule",
  table.column.values = c(0, 1),
  sizes               = c(22, 80),
  style.name          = style.name
)

# 9. Edge width by correlation strength
wmin <- min(edges_mod$weight, na.rm = TRUE)
wmax <- max(edges_mod$weight, na.rm = TRUE)

setEdgeLineWidthMapping(
  table.column        = "weight",
  table.column.values = c(wmin, wmax),
  widths              = c(1, 8),
  style.name          = style.name
)

# 10. Edge color by correlation sign
setEdgeColorMapping(
  table.column        = "weight",
  table.column.values = c(wmin, 0, wmax),
  colors              = c("#2c7bb6", "#dddddd", "#d7191c"),
  style.name          = style.name
)

# 11. See available layouts
getLayoutNames()

# 12. Apply a nice layout
layoutNetwork("force-directed")
fitContent()

# 13. Export
exportImage(filename="Cytoscape_turquoise_network", type="PNG", resolution=300)
saveSession("WGCNA_turquoise.cys")
