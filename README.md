
# Exploring Relationships Between Ion Channels and lncRNAs in Gastric Cancer

---

##  Overview

This repository explores the role of ion channels and long non-coding RNAs (lncRNAs) in gastric cancer (GC).

Using **RNA-seq data from TCGA-STAD**, the project applies **bioinformatics and systems biology approaches** such as differential expression analysis, weighted gene co-expression network analysis (WGCNA), and clinical trait correlations. The ultimate goal is to uncover potential **lncRNA–ion channel interactions** relevant to gastric cancer biology.

 **Note:** This is a **work in progress**. The repository currently provides the pipeline and scripts. Outputs and finalized results will be added later.

---

## Importance of the Project

* **Gastric cancer** is one of the leading causes of cancer-related deaths worldwide, with poor survival rates in advanced stages.
* **Ion channels** are key regulators of cellular signaling, proliferation, and apoptosis, and their dysregulation has been linked to cancer progression.
* **lncRNAs** are emerging as crucial regulators of gene expression and tumor biology, often acting through complex gene networks.
* Despite their individual importance, **the interplay between ion channels and lncRNAs in gastric cancer remains poorly understood**.
* By identifying **lncRNA–ion channel co-expression modules** and linking them to **clinical traits and survival**, this project may contribute to:

  * Improved understanding of gastric cancer mechanisms
  * Discovery of potential **biomarkers** for diagnosis or prognosis
  * Highlighting new **therapeutic targets** for future studies

---

##  Objectives

* Download and preprocess **TCGA-STAD RNA-seq data**
* Perform **differential expression analysis (DEG)**
* Construct co-expression networks using **WGCNA**
* Detect modules enriched in **lncRNAs and ion channels**
* Correlate modules with **clinical traits and survival**
* Export networks for visualization in **Cytoscape**

---

##  Repository Structure

```
.
├── Part1_DataAcquisition/                                                              # TCGA query, manifest, counts merging
├── Part2_Normalization_QC/                                                             # DESeq2 normalization, QC scripts
├── Part3_DEG_Analysis/                                                                 # Differential expression analysis
├── Part4_WGCNA/                                                                        # WGCNA scripts & module detection
├── Exploring Relationships Between Ion Channels and lncRNAs in Gastric Cancer.pdf      # Project report draft (methods + notes, in progress)
└── README.md
```

---

##  Workflow

### Requirements

* **R ≥ 4.2**
* Key packages:
  `TCGAbiolinks`, `SummarizedExperiment`, `DESeq2`,
  `WGCNA`, `matrixStats`, `dplyr`, `survival`, `survminer`, `RCy3`

### Pipeline

1. **Data Acquisition**

   * Query and download **TCGA-STAD RNA-seq raw counts** using `TCGAbiolinks`.
   * Save the manifest and prepare a clean dataset for downstream analysis.

2. **Normalization & Quality Control**

   * Normalize raw counts using **DESeq2** (variance-stabilizing transformation).
   * Perform sample quality control: PCA, clustering, detection of outliers.
   * Annotate genes with GENCODE v36.

3. **Differential Expression Analysis (DEG)**

   * Identify **tumor vs normal differentially expressed genes**.
   * Separate lncRNAs and protein-coding genes.
   * Generate summary DEG tables and volcano plots.

4. **Weighted Gene Co-expression Network Analysis (WGCNA)**

   * Construct a signed network using variance-stabilized counts.
   * Detect gene modules and merge highly similar ones.
   * Correlate module eigengenes with clinical traits (tumor/normal, stage, age, survival).
   * Export **module–trait heatmaps and dendrograms**.

5. **Focus on Ion Channels & lncRNAs**

   * Enrich modules for **ion channel genes** and **lncRNAs**.
   * Identify lncRNA–ion channel co-expression pairs.
   * Build **Cytoscape-ready networks** (nodes & edges tables).

6. **Future Steps (Planned)**

   * Perform **functional enrichment** of selected modules.
   * Validate findings with external datasets and literature.
   * Integrate survival associations with hub genes.

---

##  Project Status

*  Pipeline and scripts implemented
*  Analysis in progress (outputs will be added after refinement)
*  Next steps: validation, enrichment, and biological interpretation
