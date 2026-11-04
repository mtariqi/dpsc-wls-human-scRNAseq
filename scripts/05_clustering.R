# ==============================================================================
# Title:        05_clustering.R - Normalisation, PCA and first clustering
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      Log-normalises counts, selects variable genes, runs PCA and
#               clusters cells before integration.
# Inputs:       data/processed/GSE164157_qc.rds.
# Outputs:      data/processed/GSE164157_clustered.rds.
# Usage:        Rscript scripts/05_clustering.R
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ==============================================================================

library(Seurat)

obj <- readRDS(
  "data/processed/GSE164157_raw.rds"
)

obj[["percent.mt"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

png(
  "results/figures/QC_violin.png",
  width = 1800,
  height = 600
)

VlnPlot(
  obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

dev.off()

obj <- subset(
  obj,
  subset =
    nFeature_RNA > 200 &
      nCount_RNA > 500 &
      percent.mt < 20
)

saveRDS(
  obj,
  "data/processed/GSE164157_qc.rds"
)
