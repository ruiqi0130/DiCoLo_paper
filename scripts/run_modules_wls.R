# =============================================================================
# run_modules_wls.R   (Fig S5, Fig S6)
#
# Wls-KO vs wild-type control (E14.5). Wnt secretion is lost in the mutant, so
# the co-localized programs are the ones that exist in CTL and disappear in
# Wls; the operator is therefore computed for CTL against Wls.
#
#   Fig S6 A-D   derivation of the CTL-specific modules
#   Fig S5 A     merged UMAP, coloured by condition
#   Fig S5 B     per-condition UMAP, by cell type and by module activity
#   Fig S5 C     MSigDB hallmark pathway enrichment of the modules
#
# Also writes CTL_dicolo_modules.rds for this dataset.
#
# Workflow follows the DiCoLo package vignette:
#   https://klugerlab.github.io/DiCoLo/articles/DiCoLo_demo.html
#
# Usage (from the repository root):
#   Rscript scripts/run_modules_wls.R
# =============================================================================

suppressPackageStartupMessages({
  library(DiCoLo)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()

data.path <- file.path(DATA_DIR, "wlsko")
emd.path  <- file.path(data.path, "parameters")

conditions <- c(MUT = "data_S_wlsko_dermal_E14.5_MUT.rds",
                CTL = "data_S_wlsko_dermal_E14.5_CTL.rds")
COMP <- "CTL"   # direction of the contrast: CTL-specific co-localization

# =============================================================================
# 1. Load data
# =============================================================================
data_S_ls <- lapply(conditions, function(f) readRDS(file.path(data.path, f)))
names(data_S_ls) <- names(conditions)

# =============================================================================
# 2. Fig S5A - merged embedding, coloured by condition
# =============================================================================
data_S_merge <- merge(data_S_ls[[1]], data_S_ls[-1])
data_S_merge <- data_S_merge %>% NormalizeData() %>%
  FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(dims = 1:15, seed.use = GLOBAL_SEED)

p <- DimPlot(data_S_merge, group.by = "orig.ident", shuffle = TRUE) +
  NoAxes() + ggtitle(NULL)
ggsave(file.path(FIGURE_DIR, "figS5A_condition_umap.png"), p,
       width = 5, height = 3, dpi = 300)

# =============================================================================
# 3. Step 1 - gene-gene OT distances per condition
# =============================================================================
common_genes <- SelectCommonGenes(data_S_ls[["MUT"]], data_S_ls[["CTL"]],
                                  ngenes = 500)

emd_paths <- setNames(
  file.path(emd.path, paste0("GeneTrajectory_", names(conditions))),
  names(conditions))

for (cond in names(emd_paths)) {
  if (!file.exists(file.path(emd_paths[[cond]], "emd.csv"))) {
    message("Computing gene-gene OT distances for ", cond)
    ComputeGeneEMD(data_S_ls[[cond]], common_genes,
                   dir.path = emd_paths[[cond]], script_dir = OT_SCRIPT_DIR)
  }
}

emd_list  <- LoadAlignedGeneEMD(emd_paths)
tsne_list <- lapply(emd_list, function(x) ObtainGeneTSNE(as.dist(x)))

# =============================================================================
# 4. Fig S6 - steps 2-3 and module assembly (CTL direction)
# =============================================================================
res <- DeriveDiCoLoModules(emd_list, tsne_list, comp = COMP, 
                           min_gene = 10, lfdr_thresh = 0.01)

png(file.path(FIGURE_DIR, "figS6A_eigenvalue_spectrum.png"),
    width = 300 * 5, height = 300 * 5, res = 300)
res$panels$A()
dev.off()

ggsave(file.path(FIGURE_DIR, "figS6B_eigenvector_loadings.png"),
       res$panels$B, width = 4 * res$n_eigvec, height = 6, dpi = 300,
       limitsize = FALSE)
ggsave(file.path(FIGURE_DIR, "figS6C_candidate_genes.png"),
       res$panels$C, width = 4 * res$n_eigvec, height = 3.5, dpi = 300,
       limitsize = FALSE)
ggsave(file.path(FIGURE_DIR, "figS6D_gene_modules.png"),
       res$panels$D, width = 10, height = 4, dpi = 300)

gene_modules <- res$modules
saveRDS(gene_modules,
        file.path(data.path, sprintf("%s_dicolo_modules.rds", COMP)))

gene_modules_ls <- split(names(gene_modules), gene_modules)

# =============================================================================
# 5. Fig S5B - cell types and module activity on the per-condition manifolds
# =============================================================================
color_ct <- setNames(c("darkred", "darkblue", "#AC9362"), c("DC", "LD", "UD"))
module_prefix <- "Module_"

data_S_ls <- lapply(data_S_ls, function(data_S) {
  DefaultAssay(data_S) <- "RNA"
  data_S <- AddModuleScore(object = data_S, features = gene_modules_ls,
                           ctrl = 100, name = module_prefix)
  module_score <- data_S@meta.data[, grepl(module_prefix,
                                           colnames(data_S@meta.data))]
  module_score <- apply(module_score, 2, function(x) scale(x)[, 1])
  module_score <- apply(module_score, 2,
                        function(x) (x - min(x)) / (max(x) - min(x)))
  AddMetaData(data_S, module_score)
})

p_celltype <- wrap_plots(lapply(data_S_ls, function(data_S) {
  DimPlot(data_S, group.by = "celltype", cols = color_ct) + NoAxes() +
    ggtitle(NULL)
}), nrow = 2) + plot_layout(guides = "collect")

p_modules <- wrap_plots(lapply(data_S_ls, function(data_S) {
  FeaturePlot(data_S,
              features = paste0(module_prefix, seq_along(gene_modules_ls)),
              ncol = length(gene_modules_ls)) &
    scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) &
    NoAxes()
}), nrow = 2) + plot_layout(guides = "collect")

p <- (p_celltype | p_modules) +
  plot_layout(widths = c(1, length(gene_modules_ls)))
ggsave(file.path(FIGURE_DIR, "figS5B_celltype_and_modules.png"), p,
       width = 4 * (1 + length(gene_modules_ls)), height = 8, dpi = 300,
       limitsize = FALSE)

# =============================================================================
# 6. Fig S5C - MSigDB hallmark pathway enrichment
# =============================================================================
suppressPackageStartupMessages({
  library(msigdbr)
  library(clusterProfiler)
})

m_df <- msigdbr(species = "Mus musculus", category = "H")
m_df$gs_name <- gsub("HALLMARK_", "", m_df$gs_name)

msig_result <- do.call(rbind, lapply(names(gene_modules_ls), function(i) {
  df <- enricher(gene = gene_modules_ls[[i]],
                 TERM2GENE = m_df[, c("gs_name", "gene_symbol")],
                 pAdjustMethod = "none", pvalueCutoff = 0.1)@result
  df$modules <- i
  df <- df %>% filter(p.adjust < 0.2)
  df <- df[1:5, ]
  df[complete.cases(df), ]
}))

msig_result$GeneRatio <- unlist(lapply(msig_result$GeneRatio,
                                       function(x) eval(parse(text = x))))
msig_result$logP <- -log(msig_result$pvalue, base = 10)
msig_result$Description <- factor(
  msig_result$Description,
  levels = rev(unique(msig_result$Description[order(msig_result$pvalue)])))
msig_result$modules <- factor(msig_result$modules,
                              levels = names(gene_modules_ls))
levels(msig_result$modules) <- gsub("CTL", "CTL_M",
                                    levels(msig_result$modules))

p <- ggplot(msig_result, aes(x = modules, y = Description)) +
  geom_point(aes(size = GeneRatio, color = logP), alpha = 0.7) +
  scale_color_gradientn(colors = colorRampPalette(c("blue", "red"))(99),
                        name = "-log10_pvalue") +
  scale_size_continuous(name = "Gene Ratio", limits = c(0, 1)) +
  theme_minimal() +
  labs(size = "Gene Ratio", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10))

ggsave(file.path(FIGURE_DIR, "figS5C_pathway_enrichment.png"), p,
       width = 8, height = 4, dpi = 300)

message("Done. CTL modules: ", paste(levels(gene_modules), collapse = ", "))
