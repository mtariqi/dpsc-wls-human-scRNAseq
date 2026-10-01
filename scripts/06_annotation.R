library(Seurat)

obj <- readRDS(
  "data/processed/GSE164157_qc.rds"
)

obj <- NormalizeData(obj)

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 3000
)

obj <- ScaleData(obj)

obj <- RunPCA(
  obj,
  npcs = 50
)

png(
  "results/figures/PCA_Elbow.png",
  width = 1500,
  height = 1200
)

ElbowPlot(obj)

dev.off()

obj <- RunUMAP(
  obj,
  dims = 1:30
)

obj <- FindNeighbors(
  obj,
  dims = 1:30
)

obj <- FindClusters(
  obj,
  resolution = 0.5
)

png(
  "results/figures/UMAP_clusters.png",
  width = 1800,
  height = 1500
)

DimPlot(
  obj,
  reduction = "umap",
  label = TRUE
)

dev.off()

saveRDS(
  obj,
  "data/processed/GSE164157_clustered.rds"
)

print(obj)
