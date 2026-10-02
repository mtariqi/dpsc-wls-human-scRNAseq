# =============================================================================
# 06a_batch_doublets.R
# Stage 6a: batch inspection, doublet removal and Harmony integration
# Project: WLS (GPR177) expression in human dental pulp (GSE164157)
#
# Input : clustered Seurat object from 05_clustering.R
# Output: data/processed/GSE164157_integrated.rds
#         results/06a_batch/  (figures + tables)
#
# Run from the repo root:  Rscript scripts/06a_batch_doublets.R
# =============================================================================

## ---- 0. Settings ------------------------------------------------------------
IN_RDS   <- "data/processed/GSE164157_clustered.rds"
OUT_RDS  <- "data/processed/GSE164157_integrated.rds"
OUT_DIR  <- "results/06a_batch"
N_PCS    <- 20      # elbow plot flattened around PC15-20
RES      <- 0.5     # clustering resolution after integration
set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
  library(harmony); library(scDblFinder); library(SingleCellExperiment)
})
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_RDS), recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, w = 9, h = 6) ggsave(file.path(OUT_DIR, name), p, width = w, height = h, dpi = 200)

## ---- 1. Load ----------------------------------------------------------------
if (!file.exists(IN_RDS)) {
  found <- list.files(".", pattern = "\\.rds$", recursive = TRUE, ignore.case = TRUE)
  stop("Could not find ", IN_RDS, ".\nRDS files in this repo:\n  ",
       paste(found, collapse = "\n  "),
       "\nSet IN_RDS at the top of this script to the clustered object.")
}
obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)          # Seurat v5: merged samples are stored as separate layers
cat("Loaded:", ncol(obj), "cells x", nrow(obj), "features\n")

## ---- 2. Make sure every cell has a sample label ------------------------------
if (!"sample" %in% colnames(obj[[]])) {
  if (length(unique(obj$orig.ident)) > 1) {
    obj$sample <- obj$orig.ident
  } else {
    # merge() without add.cell.ids appends _1, _2, ... to duplicated barcodes;
    # this recovers the sample index only if every sample got a suffix.
    stop("No 'sample' column and orig.ident has one value. ",
         "Re-run 04_QC.R with CreateSeuratObject(project = 'Pulp1') etc. ",
         "or merge(..., add.cell.ids = c('Pulp1', ...)).")
  }
}
obj$sample <- factor(obj$sample)
print(table(obj$sample))

mt_col <- intersect(c("percent.mt", "percent_mt", "pct_mt"), colnames(obj[[]]))[1]
if (is.na(mt_col)) obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-") else
  if (mt_col != "percent.mt") obj$percent.mt <- obj[[mt_col]][, 1]

## ---- 3. Batch inspection on the existing (unintegrated) clusters ------------
if (!"umap" %in% Reductions(obj)) obj <- RunUMAP(obj, dims = 1:N_PCS)
obj$clusters_unintegrated <- if ("seurat_clusters" %in% colnames(obj[[]])) obj$seurat_clusters else Idents(obj)

p <- DimPlot(obj, group.by = "sample", shuffle = TRUE, raster = FALSE) + ggtitle("Before integration: by sample") |
     DimPlot(obj, group.by = "clusters_unintegrated", label = TRUE, raster = FALSE) + NoLegend() + ggtitle("Before integration: clusters")
save_plot(p, "01_UMAP_by_sample_before.png", w = 14)
save_plot(DimPlot(obj, split.by = "sample", group.by = "clusters_unintegrated", ncol = 3, raster = FALSE) + NoLegend(),
          "02_UMAP_split_by_sample_before.png", w = 14, h = 9)

comp <- as.data.frame.matrix(table(obj$clusters_unintegrated, obj$sample))
comp_frac <- sweep(comp, 1, rowSums(comp), "/")
batch_tab <- data.frame(cluster = rownames(comp), n_cells = rowSums(comp),
                        dominant_sample = colnames(comp_frac)[apply(comp_frac, 1, which.max)],
                        dominant_fraction = round(apply(comp_frac, 1, max), 3),
                        n_samples_over_5pct = rowSums(comp_frac > 0.05))
batch_tab$flag_sample_specific <- batch_tab$dominant_fraction > 0.8
write.csv(cbind(batch_tab, comp), file.path(OUT_DIR, "cluster_by_sample_before.csv"), row.names = FALSE)
cat("\nClusters with >80% of cells from one sample (possible batch clusters):\n")
print(subset(batch_tab, flag_sample_specific))

long <- data.frame(cluster = obj$clusters_unintegrated, sample = obj$sample)
save_plot(ggplot(long, aes(cluster, fill = sample)) + geom_bar(position = "fill") +
            labs(y = "Fraction of cells", title = "Sample composition of each cluster (before integration)") +
            theme_minimal(), "03_cluster_composition_before.png", w = 12, h = 5)

save_plot(VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "clusters_unintegrated",
                  pt.size = 0, ncol = 1), "04_QC_by_cluster.png", w = 14, h = 10)

## ---- 4. Doublet detection (scDblFinder, run per sample) ----------------------
sce <- SingleCellExperiment(list(counts = LayerData(obj, assay = "RNA", layer = "counts")),
                            colData = obj[[]])
sce <- scDblFinder(sce, samples = "sample", BPPARAM = BiocParallel::SerialParam())
obj$doublet_class <- sce$scDblFinder.class
obj$doublet_score <- sce$scDblFinder.score
dbl <- as.data.frame.matrix(table(obj$sample, obj$doublet_class))
dbl$doublet_pct <- round(100 * dbl$doublet / rowSums(dbl), 1)
write.csv(dbl, file.path(OUT_DIR, "doublets_by_sample.csv"))
cat("\nDoublets by sample:\n"); print(dbl)

dbl_cl <- obj[[]] %>% group_by(clusters_unintegrated) %>%
  summarise(n = n(), doublet_pct = round(100 * mean(doublet_class == "doublet"), 1)) %>% arrange(desc(doublet_pct))
write.csv(dbl_cl, file.path(OUT_DIR, "doublets_by_cluster.csv"), row.names = FALSE)
save_plot(DimPlot(obj, group.by = "doublet_class", cols = c(singlet = "grey80", doublet = "red"),
                  order = "doublet", raster = FALSE) + ggtitle("scDblFinder calls"), "05_UMAP_doublets.png", w = 8, h = 6)

n_before <- ncol(obj)
obj <- subset(obj, subset = doublet_class == "singlet")
cat("\nRemoved", n_before - ncol(obj), "doublets;", ncol(obj), "singlets remain\n")

## ---- 5. Re-process singlets and integrate with Harmony ----------------------
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
save_plot(ElbowPlot(obj, ndims = 50), "06_elbow_50PCs.png", w = 7, h = 5)

obj <- RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:N_PCS, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:N_PCS, reduction.name = "umap.harmony", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
obj <- FindClusters(obj, resolution = RES, cluster.name = "harmony_clusters", verbose = FALSE)
Idents(obj) <- "harmony_clusters"

p <- DimPlot(obj, reduction = "umap.harmony", group.by = "sample", shuffle = TRUE, raster = FALSE) + ggtitle("After Harmony: by sample") |
     DimPlot(obj, reduction = "umap.harmony", label = TRUE, raster = FALSE) + NoLegend() + ggtitle(paste0("After Harmony: ", nlevels(Idents(obj)), " clusters"))
save_plot(p, "07_UMAP_after_harmony.png", w = 14)

comp2 <- as.data.frame.matrix(table(obj$harmony_clusters, obj$sample))
frac2 <- sweep(comp2, 1, rowSums(comp2), "/")
after <- data.frame(cluster = rownames(comp2), n_cells = rowSums(comp2),
                    dominant_sample = colnames(frac2)[apply(frac2, 1, which.max)],
                    dominant_fraction = round(apply(frac2, 1, max), 3))
write.csv(cbind(after, comp2), file.path(OUT_DIR, "cluster_by_sample_after.csv"), row.names = FALSE)
cat("\nAfter integration, clusters still >80% from one sample:\n")
print(subset(after, dominant_fraction > 0.8))

## ---- 6. Gene name check before annotation -----------------------------------
for (g in c("WLS", "GPR177", "THY1", "MCAM", "COL1A1")) {
  hits <- grep(paste0("^", g, "($|[.-])"), rownames(obj), value = TRUE)
  cat(sprintf("%-7s -> %s\n", g, if (length(hits)) paste(hits, collapse = ", ") else "NOT FOUND"))
}

saveRDS(obj, OUT_RDS)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nSaved", OUT_RDS, "\nNext: Rscript scripts/06_annotation.R\n")
