# ==============================================================================
# Title:        02_download_GEO.R - Download GSE164157 from GEO
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      Downloads the GSE164157 count matrices for five donors with
#               GEOquery.
# Inputs:       None.
# Outputs:      data/raw/GSE164157/.
# Usage:        Rscript scripts/02_download_GEO.R
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ==============================================================================

library(GEOquery)

datasets <- c(
  "GSE164157",
  "GSE202476",
  "GSE185222"
)

for (ds in datasets) {
  cat("\nDownloading:", ds, "\n")

  tryCatch(
    {
      getGEOSuppFiles(
        GEO = ds,
        makeDirectory = TRUE,
        baseDir = "data/raw"
      )
    },
    error = function(e) {
      cat("Failed:", ds, "-", conditionMessage(e), "\n")
    }
  )
}

cat("\nDone. Downloaded files are in data/raw/\n")
