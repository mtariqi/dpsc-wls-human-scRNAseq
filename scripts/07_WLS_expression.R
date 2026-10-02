# =============================================================================
# 07_WLS_expression.R
# Stage 7: WLS (GPR177) expression across cell types and in candidate DPSCs,
#          with donor-level statistics
# Project: WLS expression in human dental pulp (GSE164157)
#
# Input : data/processed/GSE164157_annotated.rds  (from 06_annotation.R)
# Output: results/07_WLS/  (figures + tables)
#
# Statistics: donors (n = 5), not cells, are the unit of replication.
#   (1) per-donor % WLS+ and mean expression for every cell group
#   (2) pseudobulk edgeR QL tests with donor as a blocking factor
#   (3) mixed-effects logistic model of WLS detection with donor random effect
# =============================================================================

IN_RDS  <- "data/processed/GSE164157_annotated.rds"
OUT_DIR <- "results/07_WLS"
MIN_CELLS <- 20          # minimum cells per donor x group for a pseudobulk sample
EXCLUDE   <- c("Erythrocyte", "Epithelial")   # too few cells / not pulp-resident
set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(edgeR); library(Matrix)
})
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, w = 10, h = 6) ggsave(file.path(OUT_DIR, name), p, width = w, height = h, dpi = 200)

obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"
obj <- obj[, !(obj$cell_type %in% EXCLUDE)]
obj$group <- ifelse(obj$candidate_DPSC, "Candidate_DPSC", obj$cell_type)
obj$group[obj$group == "Perivascular_MCAM"] <- "Perivascular_other"
red <- if ("umap.harmony" %in% Reductions(obj)) "umap.harmony" else "umap"
WLS <- grep("^WLS($|[.-])", rownames(obj), value = TRUE)[1]
stopifnot(!is.na(WLS))
cat("Cells:", ncol(obj), "| WLS feature:", WLS, "\n"); print(table(obj$group))

## ---- 1. Visualisation -------------------------------------------------------
ord <- obj[[]] %>% mutate(w = FetchData(obj, WLS)[, 1] > 0) %>% group_by(group) %>%
  summarise(p = mean(w)) %>% arrange(desc(p)) %>% pull(group)
obj$group <- factor(obj$group, levels = ord)

save_plot(FeaturePlot(obj, WLS, reduction = red, order = TRUE, raster = FALSE) + ggtitle("WLS expression") |
          DimPlot(obj, reduction = red, group.by = "group", label = TRUE, repel = TRUE, raster = FALSE) + NoLegend(),
          "01_UMAP_WLS_and_groups.png", w = 16, h = 7)
save_plot(VlnPlot(obj, WLS, group.by = "group", pt.size = 0) + NoLegend() + ggtitle("WLS by cell group"),
          "02_Violin_WLS.png", w = 10, h = 5)

wnts  <- intersect(c(paste0("WNT", c("1","2","2B","3","3A","4","5A","5B","6","7A","7B","8A","8B","9A","9B","10A","10B","11","16"))), rownames(obj))
recep <- intersect(c(paste0("FZD", 1:10), "LRP5", "LRP6", "PORCN"), rownames(obj))
save_plot(DotPlot(obj, features = c(WLS, "PORCN", wnts), group.by = "group") + RotatedAxis() +
            ggtitle("WLS, PORCN and WNT ligands"), "03_DotPlot_WLS_WNT_ligands.png", w = 13, h = 5)
save_plot(DotPlot(obj, features = recep, group.by = "group") + RotatedAxis() +
            ggtitle("Wnt receptors and co-receptors"), "04_DotPlot_Wnt_receptors.png", w = 12, h = 5)

## ---- 2. Per-donor summaries ----------------------------------------------------
df <- data.frame(group = obj$group, sample = obj$sample, WLS = FetchData(obj, WLS)[, 1])
per_donor <- df %>% group_by(group, sample) %>%
  summarise(n_cells = n(), pct_pos = 100 * mean(WLS > 0), mean_expr = mean(WLS), .groups = "drop")
write.csv(per_donor, file.path(OUT_DIR, "WLS_per_donor_per_group.csv"), row.names = FALSE)

overall <- 100 * mean(df$WLS > 0)
summ <- per_donor %>% filter(n_cells >= MIN_CELLS) %>% group_by(group) %>%
  summarise(n_donors = n(), cells = sum(n_cells),
            pct_pos_median = round(median(pct_pos), 1), pct_pos_min = round(min(pct_pos), 1),
            pct_pos_max = round(max(pct_pos), 1), mean_expr = round(mean(mean_expr), 3),
            fold_enrichment_vs_all = round(median(pct_pos) / overall, 2), .groups = "drop") %>%
  arrange(desc(pct_pos_median))
write.csv(summ, file.path(OUT_DIR, "WLS_summary_by_group.csv"), row.names = FALSE)
cat("\nWLS by group (donor-level, groups with >=", MIN_CELLS, "cells per donor):\n"); print(summ)

save_plot(ggplot(filter(per_donor, n_cells >= MIN_CELLS), aes(group, pct_pos)) +
            geom_boxplot(outlier.shape = NA, fill = "grey92") +
            geom_point(aes(colour = sample), size = 2.5, position = position_jitter(width = .15, seed = 1)) +
            theme_minimal() + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
            labs(x = NULL, y = "% WLS+ cells", title = "WLS detection per donor"),
          "05_WLS_pct_per_donor.png", w = 11, h = 6)

## ---- 3. Pseudobulk (sum counts per donor x group) ------------------------------
counts <- LayerData(obj, assay = "RNA", layer = "counts")
key <- factor(paste(obj$sample, obj$group, sep = "|"))
pb  <- counts %*% sparse.model.matrix(~ 0 + key)
colnames(pb) <- levels(key)
ncell <- table(key)[colnames(pb)]
meta <- data.frame(id = colnames(pb), sample = sub("\\|.*", "", colnames(pb)),
                   group = sub(".*\\|", "", colnames(pb)), n_cells = as.integer(ncell))
keep <- meta$n_cells >= MIN_CELLS
pb <- pb[, keep]; meta <- meta[keep, ]

y_all <- DGEList(as.matrix(pb), samples = meta)
y_all <- normLibSizes(y_all[filterByExpr(y_all, group = meta$group), , keep.lib.sizes = FALSE])
lc <- cpm(y_all, log = TRUE)
if (WLS %in% rownames(lc)) {
  pb_wls <- cbind(meta, logCPM_WLS = round(lc[WLS, ], 3))
  write.csv(pb_wls, file.path(OUT_DIR, "WLS_pseudobulk_logCPM.csv"), row.names = FALSE)
  save_plot(ggplot(pb_wls, aes(reorder(group, -logCPM_WLS, median), logCPM_WLS)) +
              geom_boxplot(outlier.shape = NA, fill = "grey92") + geom_point(aes(colour = sample), size = 2.5) +
              theme_minimal() + theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
              labs(x = NULL, y = "WLS log2 CPM (pseudobulk)", title = "WLS pseudobulk expression per donor"),
            "06_WLS_pseudobulk_logCPM.png", w = 11, h = 6)
}

# Paired tests: candidate DPSC vs each comparator, donor as blocking factor
compare <- function(ref) {
  m <- meta[meta$group %in% c("Candidate_DPSC", ref), ]
  donors <- names(which(table(m$sample) == 2))           # donors with both groups
  m <- m[m$sample %in% donors, ]
  if (length(donors) < 3) return(data.frame(comparison = ref, n_donors = length(donors), note = "too few donors"))
  y <- DGEList(as.matrix(pb[, m$id]), samples = m)
  m$group <- relevel(factor(m$group), ref = ref)
  design <- model.matrix(~ sample + group, data = m)
  y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
  if (!WLS %in% rownames(y)) return(data.frame(comparison = ref, n_donors = length(donors), note = "WLS filtered out (low counts)"))
  y <- normLibSizes(y); y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- glmQLFTest(fit, coef = "groupCandidate_DPSC")
  tt <- topTags(res, n = Inf)$table
  data.frame(comparison = paste("Candidate_DPSC vs", ref), n_donors = length(donors),
             WLS_log2FC = round(tt[WLS, "logFC"], 3), WLS_PValue = signif(tt[WLS, "PValue"], 3),
             WLS_FDR_genomewide = signif(tt[WLS, "FDR"], 3),
             WLS_rank = which(rownames(tt) == WLS), genes_tested = nrow(tt), note = "")
}
# "Rest" = all non-DPSC cells pooled
meta_rest <- meta; pb_rest <- pb
rest_ids <- meta$id[meta$group != "Candidate_DPSC"]
for (s in unique(meta$sample)) {
  ids <- intersect(rest_ids, meta$id[meta$sample == s]); if (!length(ids)) next
  pb_rest <- cbind(pb_rest, Matrix::rowSums(pb[, ids, drop = FALSE]))
  colnames(pb_rest)[ncol(pb_rest)] <- paste0(s, "|All_other_cells")
  meta_rest <- rbind(meta_rest, data.frame(id = paste0(s, "|All_other_cells"), sample = s,
                                           group = "All_other_cells", n_cells = sum(meta$n_cells[meta$id %in% ids])))
}
pb <- pb_rest; meta <- meta_rest
comparators <- intersect(c("Pulp_fibroblast", "Endothelial", "Schwann_glia", "Myeloid", "All_other_cells"), meta$group)
tests <- bind_rows(lapply(comparators, compare))
write.csv(tests, file.path(OUT_DIR, "WLS_pseudobulk_edgeR_tests.csv"), row.names = FALSE)
cat("\nPseudobulk edgeR (positive log2FC = higher in candidate DPSCs):\n"); print(tests)

## ---- 4. Mixed-effects model of WLS detection ----------------------------------
if (requireNamespace("lme4", quietly = TRUE)) {
  mes_df <- df[df$group %in% c("Candidate_DPSC", "Pulp_fibroblast", "Perivascular_other"), ]
  mes_df$is_DPSC <- mes_df$group == "Candidate_DPSC"; mes_df$pos <- mes_df$WLS > 0
  fitm <- lme4::glmer(pos ~ is_DPSC + (1 | sample), data = mes_df, family = binomial)
  co <- summary(fitm)$coefficients["is_DPSCTRUE", ]
  ci <- exp(co["Estimate"] + c(-1.96, 1.96) * co["Std. Error"])
  mm <- data.frame(model = "WLS+ ~ candidate DPSC + (1|donor), mesenchymal cells",
                   odds_ratio = round(exp(co["Estimate"]), 3), CI95_low = round(ci[1], 3),
                   CI95_high = round(ci[2], 3), p_value = signif(co["Pr(>|z|)"], 3))
  write.csv(mm, file.path(OUT_DIR, "WLS_mixed_model.csv"), row.names = FALSE)
  cat("\nMixed-effects logistic model (within mesenchyme):\n"); print(mm)
} else cat("\nlme4 not installed: skipping mixed model (install.packages('lme4'))\n")

## ---- 5. Co-expression with WNT ligands ----------------------------------------
if (length(wnts)) {
  W <- FetchData(obj, c(WLS, wnts)) > 0
  any_wnt <- rowSums(W[, -1, drop = FALSE]) > 0
  co_tab <- data.frame(group = obj$group, wls = W[, 1], wnt = any_wnt) %>% group_by(group) %>%
    summarise(n = n(), pct_WLS = round(100 * mean(wls), 1), pct_anyWNT = round(100 * mean(wnt), 1),
              pct_both = round(100 * mean(wls & wnt), 2),
              pct_both_expected = round(100 * mean(wls) * mean(wnt), 2), .groups = "drop") %>%
    mutate(obs_over_expected = round(pct_both / pct_both_expected, 2))
  write.csv(co_tab, file.path(OUT_DIR, "WLS_WNT_coexpression.csv"), row.names = FALSE)
  cat("\nWLS / WNT ligand co-expression:\n"); print(co_tab)
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nDone. Results in", OUT_DIR, "\n")
