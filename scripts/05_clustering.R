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
