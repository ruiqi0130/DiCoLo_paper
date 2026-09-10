# =============================================================================
# run_modules_smom2.R   (Fig S2, Fig S3, Supplemental Table 2)
#
# Derivation of the condition-specific co-localized gene modules for the SmoM2
# contrast (SmoM2 mutant vs wild-type control, E13.5), in both directions:
#
#   Fig S2 A-D   SmoM2-specific modules (operator computed for MUT vs CTL)
#   Fig S3 A-D   CTL-specific modules   (operator computed for CTL vs MUT)
#
# Also writes MUT_dicolo_modules.rds and CTL_dicolo_modules.rds, which are the
# upstream dependency of Fig 4B, 4C, 4D, Fig S4 and run_concordance_analysis.R,
# and the gene lists behind Supplemental Table 2.
#
# Workflow follows the DiCoLo package vignette:
#   https://klugerlab.github.io/DiCoLo/articles/DiCoLo_demo.html
#
# Usage (from the repository root):
#   Rscript scripts/run_modules_smom2.R
# =============================================================================

suppressPackageStartupMessages({
  library(DiCoLo)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()

data.path <- file.path(DATA_DIR, "smom2")
emd.path  <- file.path(data.path, "parameters")

conditions <- c(MUT = "data_S_smom2_dermal_E13.5_MUT.rds",
                CTL = "data_S_smom2_dermal_E13.5_CTL.rds")

# =============================================================================
# 1. Load data
# =============================================================================
data_S_ls <- lapply(conditions, function(f) readRDS(file.path(data.path, f)))
names(data_S_ls) <- names(conditions)

# =============================================================================
# 2. Step 1 - gene-gene OT distances per condition
#
# Skipped when emd.csv already exists, which is the case for the copies
# distributed in DiCoLo_data. Computing from scratch takes ~5-10 min per
# condition and needs the Python POT library (see README - Dependencies).
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
# 3. Steps 2-3 and module assembly, once per direction
# =============================================================================
directions <- list(MUT = "figS2", CTL = "figS3")
results <- list()

for (comp in names(directions)) {
  prefix <- directions[[comp]]
  message("Deriving ", comp, "-specific modules -> ", prefix)

  res <- DeriveDiCoLoModules(emd_list, tsne_list, comp = comp)
  results[[comp]] <- res

  # Panel A - eigenvalue spectrum with the knee-point selection
  png(file.path(FIGURE_DIR, sprintf("%sA_eigenvalue_spectrum.png", prefix)),
      width = 300 * 5, height = 300 * 5, res = 300)
  res$panels$A()
  dev.off()

  # Panel B - eigenvector loadings on both gene embeddings
  ggsave(file.path(FIGURE_DIR, sprintf("%sB_eigenvector_loadings.png", prefix)),
         res$panels$B, width = 4 * res$n_eigvec, height = 6, dpi = 300,
         limitsize = FALSE)

  # Panel C - locFDR candidate genes
  ggsave(file.path(FIGURE_DIR, sprintf("%sC_candidate_genes.png", prefix)),
         res$panels$C, width = 4 * res$n_eigvec, height = 3.5, dpi = 300,
         limitsize = FALSE)

  # Panel D - modules on both gene embeddings
  ggsave(file.path(FIGURE_DIR, sprintf("%sD_gene_modules.png", prefix)),
         res$panels$D, width = 10, height = 4, dpi = 300)

  # Module assignments, consumed by Fig 4B/4C/4D, Fig S4 and the
  # concordance analysis
  saveRDS(res$modules,
          file.path(data.path, sprintf("%s_dicolo_modules.rds", comp)))
}

# =============================================================================
# 4. Supplemental Table 2 - module gene lists
# =============================================================================
table.path <- file.path(FIGURE_DIR, "tables")
if (!dir.exists(table.path)) dir.create(table.path, recursive = TRUE)

module_table <- do.call(rbind, lapply(names(results), function(comp) {
  m <- results[[comp]]$modules
  data.frame(direction = comp,
             module    = as.character(m),
             gene      = names(m),
             row.names = NULL)
}))
module_table <- module_table[order(module_table$module, module_table$gene), ]
write.csv(module_table,
          file.path(table.path, "supplemental_table2_smom2_modules.csv"),
          row.names = FALSE)

message("Done. Modules: ",
        paste(vapply(results, function(r) paste(levels(r$modules), collapse = ", "),
                     character(1)), collapse = " | "))
