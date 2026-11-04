# ==============================================================================
# Title:        00_install_packages.R - Install R packages and record the environment
# Project:      Gpr177/Wntless (WLS) and Wnt/Hedgehog genes in human DPSCs
# Author:       Md Tariqul Islam
# Affiliation:  Laboratory of Professor Wei Hsu, The ADA Forsyth Institute,
#               Somerville, MA, USA
# Last updated: 2026-11-04
# ------------------------------------------------------------------------------
# Purpose:      Installs the CRAN and Bioconductor packages used by the
                pipeline and records exact versions in renv.lock.
# Inputs:       None.
# Outputs:      renv.lock; project package library (renv/).
# Usage:        Rscript scripts/00_install_packages.R
# Dependencies: R 4.6; Seurat 5, harmony, scDblFinder, edgeR, GEOquery, ggplot2,
#               patchwork, dplyr (exact versions in renv.lock)
# Development:  Written with AI assistance (Claude, Anthropic); reviewed, run and
#               validated by the author.
# ==============================================================================

cran_packages <- c(
  "Seurat",
  "tidyverse",
  "patchwork",
  "cowplot",
  "data.table",
  "remotes"
)

for(pkg in cran_packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
  }
}

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}

BiocManager::install("GEOquery")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
