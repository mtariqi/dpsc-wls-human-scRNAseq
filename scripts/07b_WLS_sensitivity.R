# =============================================================================
# 07b_WLS_sensitivity.R
# Stage 7b: does the WLS result depend on how candidate DPSCs are defined?
#           + sequencing-depth-adjusted WLS / WNT co-expression
#           + identification of the low-WLS donor
#
# Input : data/processed/GSE164157_annotated.rds
# Output: results/07b_sensitivity/
# =============================================================================

IN_RDS  <- "data/processed/GSE164157_annotated.rds"
OUT_DIR <- "results/07b_sensitivity"
MIN_CELLS <- 20
set.seed(2026)

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(edgeR); library(Matrix)
})
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, w = 10, h = 6) ggsave(file.path(OUT_DIR, name), p, width = w, height = h, dpi = 200)

obj <- readRDS(IN_RDS)
DefaultAssay(obj) <- "RNA"
WLS <- grep("^WLS($|[.-])", rownames(obj), value = TRUE)[1]
mes_types <- c("Pulp_fibroblast", "Perivascular_MCAM", "MSC_progenitor", "Odontoblast")
mes <- obj[, obj$cell_type %in% mes_types]
cat("Mesenchymal cells:", ncol(mes), "\n")

## ---- 1. Alternative DPSC definitions (all within the mesenchyme) -------------
stem_core <- intersect(c("THY1", "NT5E", "ENG", "NGFR", "FRZB", "LEPR", "PRRX1", "CXCL12"), rownames(mes))
mes <- AddModuleScore(mes, features = list(stem_core), name = "stemcell", seed = 2026)
mes$stemcell <- mes$stemcell1; mes$stemcell1 <- NULL
X <- FetchData(mes, c(stem_core, "MCAM", WLS))
n_markers <- rowSums(X[, stem_core] > 0)

defs <- list(
  A_perivascular_subclusters = mes$candidate_DPSC,                                     # stage 06 definition
  B_top20pct_stem_score      = mes$stemcell >= quantile(mes$stemcell, 0.80),          # cell-level, no perivascular genes
  C_3plus_stem_markers       = n_markers >= 3,                                         # marker count, cluster-free
  D_MCAM_positive            = X$MCAM > 0                                              # classic CD146+ DPSC
)
overlap <- sapply(defs, function(a) sapply(defs, function(b) round(sum(a & b) / sum(a), 2)))
write.csv(overlap, file.path(OUT_DIR, "definition_overlap.csv"))
cat("\nFraction of row-definition cells also in column-definition:\n"); print(overlap)

pseudobulk_test <- function(is_dpsc, label) {
  grp <- ifelse(is_dpsc, "DPSC", "Other_mesenchyme")
  key <- factor(paste(mes$sample, grp, sep = "|"))
  pb  <- LayerData(mes, layer = "counts") %*% sparse.model.matrix(~ 0 + key)
  colnames(pb) <- levels(key)
  meta <- data.frame(id = colnames(pb), sample = sub("\\|.*", "", colnames(pb)),
                     group = sub(".*\\|", "", colnames(pb)), n = as.integer(table(key)[colnames(pb)]))
  meta <- meta[meta$n >= MIN_CELLS, ]
  donors <- names(which(table(meta$sample) == 2)); meta <- meta[meta$sample %in% donors, ]
  pct <- data.frame(sample = mes$sample, grp = grp, pos = X[, WLS] > 0) %>%
    group_by(sample, grp) %>% summarise(p = 100 * mean(pos), .groups = "drop")
  pct_d <- median(pct$p[pct$grp == "DPSC"]); pct_o <- median(pct$p[pct$grp == "Other_mesenchyme"])
  out <- data.frame(definition = label, cells = sum(is_dpsc), pct_of_mesenchyme = round(100 * mean(is_dpsc), 1),
                    n_donors = length(donors), pctWLS_DPSC_median = round(pct_d, 1),
                    pctWLS_other_median = round(pct_o, 1))
  if (length(donors) < 3) return(cbind(out, log2FC = NA, PValue = NA, note = "too few donors"))
  meta$group <- relevel(factor(meta$group), ref = "Other_mesenchyme")
  y <- DGEList(as.matrix(pb[, meta$id]), samples = meta)
  design <- model.matrix(~ sample + group, data = meta)
  y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
  if (!WLS %in% rownames(y)) return(cbind(out, log2FC = NA, PValue = NA, note = "WLS filtered"))
  y <- normLibSizes(y); y <- estimateDisp(y, design)
  tt <- topTags(glmQLFTest(glmQLFit(y, design), coef = "groupDPSC"), n = Inf)$table
  cbind(out, log2FC = round(tt[WLS, "logFC"], 3), PValue = signif(tt[WLS, "PValue"], 3), note = "")
}
sens <- bind_rows(Map(pseudobulk_test, defs, names(defs)))
write.csv(sens, file.path(OUT_DIR, "WLS_sensitivity_by_definition.csv"), row.names = FALSE)
cat("\nWLS in candidate DPSCs vs other mesenchymal cells, by definition\n(positive log2FC = higher in DPSCs):\n")
print(sens)

save_plot(ggplot(sens, aes(reorder(definition, log2FC), log2FC)) + geom_col(fill = "#3b6ea5") +
            geom_hline(yintercept = 0) + coord_flip() + theme_minimal() +
            geom_text(aes(label = paste0("p=", signif(PValue, 2))), hjust = ifelse(sens$log2FC > 0, -0.1, 1.1), size = 3.5) +
            labs(x = NULL, y = "WLS log2FC, candidate DPSC vs other mesenchyme (pseudobulk)",
                 title = "Sensitivity of the WLS result to the DPSC definition"),
          "01_sensitivity_log2FC.png", w = 10, h = 4.5)

## ---- 2. Depth-adjusted WLS / WNT co-expression -------------------------------
wnts <- intersect(paste0("WNT", c("1","2","2B","3","3A","4","5A","5B","6","7A","7B","8A","8B","9A","9B","10A","10B","11","16")), rownames(obj))
obj$group <- ifelse(obj$candidate_DPSC, "Candidate_DPSC", obj$cell_type)
keep_groups <- c("Candidate_DPSC", "Pulp_fibroblast", "Endothelial", "Schwann_glia", "Myeloid")
D <- FetchData(obj, c(WLS, wnts, "nFeature_RNA", "group", "sample"))
D <- D[D$group %in% keep_groups, ]
D$wls <- D[[WLS]] > 0; D$wnt <- rowSums(D[, wnts] > 0) > 0
coex <- bind_rows(lapply(keep_groups, function(g) {
  d <- D[D$group == g, ]
  f0 <- glm(wls ~ wnt + sample, data = d, family = binomial)
  f1 <- glm(wls ~ wnt + log(nFeature_RNA) + sample, data = d, family = binomial)
  c0 <- coef(summary(f0))["wntTRUE", ]; c1 <- coef(summary(f1))["wntTRUE", ]
  data.frame(group = g, cells = nrow(d), OR_unadjusted = round(exp(c0[1]), 2),
             OR_depth_adjusted = round(exp(c1[1]), 2),
             CI95 = sprintf("%.2f-%.2f", exp(c1[1] - 1.96 * c1[2]), exp(c1[1] + 1.96 * c1[2])),
             p_adjusted = signif(c1[4], 3))
}))
write.csv(coex, file.path(OUT_DIR, "WLS_WNT_coexpression_depth_adjusted.csv"), row.names = FALSE)
cat("\nWLS-WNT co-expression odds ratio, before and after adjusting for genes detected per cell:\n"); print(coex)

lig <- D %>% filter(group %in% c("Candidate_DPSC", "Pulp_fibroblast")) %>% group_by(group) %>%
  summarise(across(all_of(wnts), ~ round(100 * mean(.x > 0), 2)), .groups = "drop")
lig <- lig[, c(TRUE, colSums(lig[, -1]) > 0.5)]
write.csv(lig, file.path(OUT_DIR, "WNT_ligand_pct_by_group.csv"), row.names = FALSE)
cat("\n% cells expressing each WNT ligand (ligands >0.5% in either group):\n"); print(as.data.frame(lig))

## ---- 3. Which donor has low WLS in candidate DPSCs? -------------------------
donor <- obj[[]] %>% mutate(wls = FetchData(obj, WLS)[, 1] > 0) %>%
  filter(candidate_DPSC) %>% group_by(sample) %>%
  summarise(cells = n(), pct_WLS = round(100 * mean(wls), 1),
            median_genes = median(nFeature_RNA), median_UMI = median(nCount_RNA),
            median_pct_mt = round(median(percent.mt), 1), .groups = "drop")
write.csv(donor, file.path(OUT_DIR, "candidate_DPSC_by_donor_QC.csv"), row.names = FALSE)
cat("\nCandidate DPSCs per donor, with QC:\n"); print(donor)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nDone. Results in", OUT_DIR, "\n")
