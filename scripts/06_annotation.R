# =============================================================================
# 06_annotation.R
# Stage 6b: cell type annotation and candidate DPSC identification
# Project: WLS (GPR177) expression in human dental pulp (GSE164157)
#
# Input : data/processed/GSE164157_integrated.rds   (from 06a_batch_doublets.R)
# Output: data/processed/GSE164157_annotated.rds
#         results/06_annotation/  (figures + tables)
#
# Annotation is two-pass:
#   Pass 1 (automatic): marker module scores suggest a label per cluster and
#          write results/06_annotation/cluster_labels_suggested.csv
#   Pass 2 (manual):    copy that file to data/annotation/cluster_labels.csv,
#          check/edit the 'label' column against the DotPlot and marker table,
#          then rerun this script. Your labels then take priority.
# =============================================================================

IN_RDS     <- "data/processed/GSE164157_integrated.rds"
OUT_RDS    <- "data/processed/GSE164157_annotated.rds"
OUT_DIR    <- "results/06_annotation"
MANUAL_CSV <- "data/annotation/cluster_labels.csv"
RUN_FINDALLMARKERS <- FALSE   # already saved in the first run; set TRUE to recompute

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, w = 10, h = 7) ggsave(file.path(OUT_DIR, name), p, width = w, height = h, dpi = 200)

obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "harmony_clusters"
red <- if ("umap.harmony" %in% Reductions(obj)) "umap.harmony" else "umap"
cat("Cells:", ncol(obj), " Clusters:", nlevels(Idents(obj)), "\n")

## ---- 1. Marker panels for human dental pulp --------------------------------
panels <- list(
  Pulp_fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "TWIST2"),
  MSC_progenitor    = c("THY1", "NT5E", "ENG", "PRRX1", "LEPR", "CXCL12", "NGFR", "FRZB"),
  Perivascular_MCAM = c("MCAM", "RGS5", "PDGFRB", "KCNJ8", "NOTCH3", "ACTA2", "MYH11"),
  Odontoblast       = c("DSPP", "DMP1", "PHEX", "SP7", "MEPE"),
  Endothelial       = c("PECAM1", "CDH5", "VWF", "PLVAP", "CLDN5"),
  Lymphatic         = c("PROX1", "LYVE1", "CCL21", "TFF3"),
  Schwann_glia      = c("SOX10", "PLP1", "MPZ", "S100B", "CDH19"),
  T_NK              = c("CD3E", "CD3D", "IL7R", "NKG7", "GNLY"),
  B_plasma          = c("MS4A1", "CD79A", "IGHG1", "JCHAIN"),
  Myeloid           = c("CD14", "LYZ", "CD68", "C1QA", "CSF1R"),
  Mast              = c("TPSAB1", "CPA3", "KIT"),
  Epithelial        = c("KRT14", "KRT5", "EPCAM", "KRT19"),
  Erythrocyte       = c("HBB", "HBA1"),
  Proliferating     = c("MKI67", "TOP2A", "CENPF")
)
present <- function(g) intersect(g, rownames(obj))
panels  <- lapply(panels, present)
missing <- setdiff(c("COL1A1", "THY1", "MCAM", "DSPP", "PECAM1", "SOX10", "PTPRC", "WLS"), rownames(obj))
if (length(missing)) cat("WARNING: markers not found:", paste(missing, collapse = ", "), "\n")
panels <- panels[lengths(panels) >= 2]

dot_genes <- unique(c(unlist(panels), "PTPRC", "WLS"))
dot_genes <- present(dot_genes)
save_plot(DotPlot(obj, features = dot_genes, cluster.idents = TRUE) + RotatedAxis() +
            theme(axis.text.x = element_text(size = 7)) + ggtitle("Lineage markers by cluster"),
          "01_DotPlot_lineage_markers.png", w = 22, h = 9)

## ---- 2. Module scores and suggested labels ----------------------------------
for (nm in names(panels)) {
  obj <- AddModuleScore(obj, features = list(panels[[nm]]), name = paste0("score_", nm), seed = 2026)
  colnames(obj@meta.data)[colnames(obj@meta.data) == paste0("score_", nm, "1")] <- paste0("score_", nm)
}
score_cols <- paste0("score_", names(panels))
cl_scores <- obj[[]] %>% group_by(cluster = harmony_clusters) %>%
  summarise(n_cells = n(), across(all_of(score_cols), mean), .groups = "drop")
m <- as.matrix(cl_scores[, score_cols]); colnames(m) <- names(panels)
top2 <- t(apply(m, 1, function(r) { o <- order(r, decreasing = TRUE)[1:2]; c(names(r)[o], r[o]) }))
suggest <- data.frame(cluster = cl_scores$cluster, n_cells = cl_scores$n_cells,
                      suggested = top2[, 1], score = round(as.numeric(top2[, 3]), 3),
                      runner_up = top2[, 2], margin = round(as.numeric(top2[, 3]) - as.numeric(top2[, 4]), 3))
suggest$confidence <- ifelse(suggest$margin > 0.3, "high", ifelse(suggest$margin > 0.1, "medium", "low: check manually"))
suggest$label <- suggest$suggested
write.csv(suggest, file.path(OUT_DIR, "cluster_labels_suggested.csv"), row.names = FALSE)
print(suggest)

hm <- as.data.frame(m); hm$cluster <- cl_scores$cluster
hm <- pivot_longer(hm, -cluster, names_to = "panel", values_to = "score")
save_plot(ggplot(hm, aes(panel, cluster, fill = score)) + geom_tile() +
            scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
            theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = "Mean module score per cluster", x = NULL, y = "Cluster"),
          "02_module_score_heatmap.png", w = 10, h = 9)

## ---- 3. Exploratory cluster markers (cell-level Wilcoxon) -------------------
if (RUN_FINDALLMARKERS) {
  mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE)
  write.csv(mk, file.path(OUT_DIR, "cluster_markers_all.csv"), row.names = FALSE)
  top10 <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
  write.csv(top10, file.path(OUT_DIR, "cluster_markers_top10.csv"), row.names = FALSE)
}

## ---- 4. Apply labels (manual file wins if present) --------------------------
labs <- if (file.exists(MANUAL_CSV)) { cat("Using manual labels from", MANUAL_CSV, "\n"); read.csv(MANUAL_CSV) } else {
  cat("No manual labels yet: using suggested labels (provisional)\n"); suggest }
obj$cell_type <- labs$label[match(as.character(obj$harmony_clusters), as.character(labs$cluster))]
obj$cell_type[is.na(obj$cell_type)] <- "Unassigned"

save_plot(DimPlot(obj, reduction = red, group.by = "cell_type", label = TRUE, repel = TRUE, raster = FALSE) +
            ggtitle(if (file.exists(MANUAL_CSV)) "Cell types" else "Cell types (provisional, automatic)"),
          "03_UMAP_cell_types.png", w = 11, h = 8)
ct <- as.data.frame.matrix(table(obj$cell_type, obj$sample))
ct$total <- rowSums(ct); ct$pct <- round(100 * ct$total / sum(ct$total), 1)
write.csv(ct[order(-ct$total), ], file.path(OUT_DIR, "cell_type_by_sample.csv"))

## ---- 5. Candidate DPSC / progenitor cells within the mesenchyme -------------
mes_types <- c("Pulp_fibroblast", "MSC_progenitor", "Perivascular_MCAM", "Odontoblast")
if (!any(obj$cell_type %in% mes_types))
  stop("No cells labelled as ", paste(mes_types, collapse = "/"),
       ". If you renamed mesenchymal labels in ", MANUAL_CSV, ", edit mes_types to match.")
mes <- subset(obj, subset = cell_type %in% mes_types)
cat("\nMesenchymal cells for sub-clustering:", ncol(mes), "\n")
mes <- FindVariableFeatures(mes, nfeatures = 2000, verbose = FALSE)
mes <- ScaleData(mes, verbose = FALSE)
mes <- RunPCA(mes, npcs = 30, verbose = FALSE)
mes <- harmony::RunHarmony(mes, group.by.vars = "sample", dims.use = 1:15, verbose = FALSE)
mes <- RunUMAP(mes, reduction = "harmony", dims = 1:15, reduction.name = "umap.mes", verbose = FALSE)
mes <- FindNeighbors(mes, reduction = "harmony", dims = 1:15, verbose = FALSE)
mes <- FindClusters(mes, resolution = 0.4, cluster.name = "mes_clusters", verbose = FALSE)

# Two separate scores so the DPSC call is not built from perivascular markers:
#   stem_core    = mesenchymal stem/progenitor markers
#   perivascular = mural/pericyte markers (reported alongside, not used for the call)
stem_core <- present(c("THY1", "NT5E", "ENG", "NGFR", "FRZB", "LEPR", "PRRX1", "CXCL12"))
periv     <- present(c("MCAM", "PDGFRB", "RGS5", "KCNJ8", "NOTCH3", "ACTA2"))
mes <- AddModuleScore(mes, features = list(stem_core), name = "stem_score", seed = 2026)
mes <- AddModuleScore(mes, features = list(periv),     name = "periv_score", seed = 2026)
mes$stem_score  <- mes$stem_score1;  mes$stem_score1  <- NULL
mes$periv_score <- mes$periv_score1; mes$periv_score1 <- NULL

pct_pos <- function(g) tapply(FetchData(mes, g)[, 1] > 0, mes$mes_clusters, mean)
marker_pct <- sapply(stem_core, pct_pos)
mes_sum <- mes[[]] %>% group_by(mes_clusters) %>%
  summarise(n_cells = n(), stem_score = round(mean(stem_score), 3),
            periv_score = round(mean(periv_score), 3),
            main_type = names(which.max(table(cell_type))), .groups = "drop")
mes_sum <- cbind(mes_sum, round(100 * marker_pct[as.character(mes_sum$mes_clusters), , drop = FALSE], 1))
mes_sum <- mes_sum[order(-mes_sum$stem_score), ]
mes_sum$candidate_DPSC <- mes_sum$stem_score > 0     # enriched above background genes
write.csv(mes_sum, file.path(OUT_DIR, "mesenchymal_subclusters_DPSC_score.csv"), row.names = FALSE)
print(mes_sum[, c("mes_clusters", "n_cells", "stem_score", "periv_score", "main_type", "candidate_DPSC")])

p <- DimPlot(mes, reduction = "umap.mes", group.by = "mes_clusters", label = TRUE) + ggtitle("Mesenchymal sub-clusters") |
     FeaturePlot(mes, "stem_score", reduction = "umap.mes", order = TRUE) + ggtitle("Stem/progenitor score") |
     FeaturePlot(mes, "periv_score", reduction = "umap.mes", order = TRUE) + ggtitle("Perivascular score")
save_plot(p, "04_mesenchyme_DPSC_score.png", w = 20, h = 6)
save_plot(DotPlot(mes, features = unique(c(stem_core, periv, present(c("COL1A1", "DSPP", "MKI67", "WLS"))))) +
            RotatedAxis() + ggtitle("Stem, perivascular and other markers in mesenchymal sub-clusters"),
          "05_DotPlot_mesenchyme_markers.png", w = 13, h = 6)

top_cl <- as.character(mes_sum$mes_clusters[mes_sum$candidate_DPSC])
mes$candidate_DPSC <- as.character(mes$mes_clusters) %in% top_cl
obj$mes_cluster <- NA; obj$mes_cluster[colnames(mes)] <- as.character(mes$mes_clusters)
obj$candidate_DPSC <- FALSE; obj$candidate_DPSC[colnames(mes)] <- mes$candidate_DPSC
obj$stem_score <- NA;  obj$stem_score[colnames(mes)]  <- mes$stem_score
cat("Candidate DPSC sub-clusters (stem score > 0):", paste(top_cl, collapse = ", "),
    " (", sum(obj$candidate_DPSC), "cells )\n")

## ---- 6. First look at WLS (full analysis in 07_WLS_expression.R) ------------
wls <- grep("^WLS($|[.-])", rownames(obj), value = TRUE)[1]
if (!is.na(wls)) {
  save_plot(FeaturePlot(obj, wls, reduction = red, order = TRUE, raster = FALSE) | 
              VlnPlot(obj, wls, group.by = "cell_type", pt.size = 0) + NoLegend(),
            "06_WLS_preview.png", w = 16, h = 6)
  expr <- FetchData(obj, c(wls, "cell_type", "sample", "candidate_DPSC"))
  names(expr)[1] <- "WLS"
  expr$group <- ifelse(expr$candidate_DPSC, "Candidate_DPSC", expr$cell_type)
  prev <- expr %>% group_by(group) %>%
    summarise(n_cells = n(), pct_WLS_pos = round(100 * mean(WLS > 0), 1),
              mean_logexpr = round(mean(WLS), 3), .groups = "drop") %>% arrange(desc(pct_WLS_pos))
  write.csv(prev, file.path(OUT_DIR, "WLS_preview_by_cell_type.csv"), row.names = FALSE)
  print(prev)
} else cat("WLS not found in feature names: check gene naming before stage 07\n")

saveRDS(obj, OUT_RDS)
saveRDS(mes, "data/processed/GSE164157_mesenchyme.rds")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nSaved", OUT_RDS, "and data/processed/GSE164157_mesenchyme.rds\n")
if (!file.exists(MANUAL_CSV)) cat("Next: review", file.path(OUT_DIR, "cluster_labels_suggested.csv"),
                                  "then copy it to", MANUAL_CSV, "with corrected labels and rerun.\n")
