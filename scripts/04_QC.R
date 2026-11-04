# ==============================================================================
# Title:        04_QC.R - Quality control of GSE164157
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      Filters cells with > 200 genes, > 500 UMIs and < 20%
                mitochondrial reads.
# Inputs:       data/processed/GSE164157_raw.rds.
# Outputs:      data/processed/GSE164157_qc.rds; QC tables and figures.
# Usage:        Rscript scripts/04_QC.R
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ==============================================================================

library(Seurat)

samples <- c(
  "GSM4998457_Pulp1",
  "GSM4998458_Pulp2",
  "GSM4998459_Pulp3",
  "GSM4998460_Pulp4",
  "GSM4998461_Pulp5"
)

objs <- list()

for(i in seq_along(samples)){

  prefix <- samples[i]

  counts <- ReadMtx(
    mtx = paste0("data/GSE164157/", prefix, "_matrix.mtx.gz"),
    features = paste0("data/GSE164157/", prefix, "_genes.tsv.gz"),
    cells = paste0("data/GSE164157/", prefix, "_barcodes.tsv.gz")
  )

  objs[[i]] <- CreateSeuratObject(
    counts = counts,
    project = prefix
  )

}

combined <- merge(
  objs[[1]],
  y = objs[-1]
)

saveRDS(
  combined,
  "data/processed/GSE164157_raw.rds"
)

print(combined)
