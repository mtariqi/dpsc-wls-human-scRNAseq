# ==============================================================================
# Title:        install_stage06.R - Analysis script
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      See Details below.
# Inputs:       See Details below.
# Outputs:      See Details below.
# Usage:        Rscript
                /home/mtariq/dpsc-wls-human-scRNAseq/scripts/install_stage06.R
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ------------------------------------------------------------------------------
# Details:
# Packages needed for stage 06 (run once inside the project so renv records them)
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("harmony", "patchwork", "dplyr", "tidyr", "ggplot2"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in c("scDblFinder", "SingleCellExperiment", "BiocParallel"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
if (!requireNamespace("presto", quietly = TRUE))   # makes FindAllMarkers much faster (optional)
  try(remotes::install_github("immunogenomics/presto"), silent = TRUE)
if (requireNamespace("renv", quietly = TRUE)) renv::snapshot(prompt = FALSE)
message("Stage 06 packages ready.")
