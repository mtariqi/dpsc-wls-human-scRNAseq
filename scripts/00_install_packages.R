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
