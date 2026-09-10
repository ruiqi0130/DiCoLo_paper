# =============================================================================
# run_covariance_diagnosis.R   (Fig S8F, S8G)
#
# Design diagnostic for the K-factor covariance model (R3-Q1). It re-uses the
# params.rds saved by run_benchmark_on_simulations.R for the cov_strength test,
# re-injects each one (deterministic, NO EMD), and measures the empirical
# gene-gene correlation among the injected genes across the injected cells.
#
# Two outputs:
#   (1) Violin: pooled off-diagonal correlations per cov_strength (all 10 reps).
#       -> shows the distribution spreads and dips NEGATIVE as strength rises,
#          flat ~0 at cov_strength = 0. This is the main claim.
#   (2) Heatmap: one rep per cov_strength, shared gene order.
#       -> shows the covariance is STRUCTURED (pair-varying, both signs), not
#          just estimation noise. Set DRAW_HEATMAP <- FALSE to skip.
#
# The injected matrices themselves are NOT saved by the driver, but injection
# is fully seeded, so reloading params.rds + re-injecting reproduces exactly the
# data the benchmark used. Requires the NEW (K-factor) simulation_functions.R.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
})

# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()
dir.path    <- DATA_DIR
figure.path <- FIGURE_DIR
source("./code/simulation_functions.R")     # <-- must be the K-factor version

# ---- settings: match your S7-D run ---------------------------------------
data.path   <- file.path(dir.path, "pbmc10k")   # dataset S7-D was run on
data.file   <- "pbmc10k.rds"
sample_ls   <- c("monocyte1", "monocyte2")
para_test   <- "cov_strength"
test_id     <- "1vs0"                            # "1vs0" => sample 1 injected
sample_diff <- 1                                 # the injected sample
cov_range   <- seq(0, 0.8, 0.2)                  # exactly as in the driver
n_reps      <- 10
cov_K       <- 3
heatmap_rep <- 1                                 # which rep to show as heatmaps
DRAW_HEATMAP <- TRUE

res.path <- file.path(data.path, "parameters", para_test, test_id)
out.path <- file.path(data.path, "parameters", para_test, "diagnosis")
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

# ---- reload + preprocess EXACTLY as the driver (so re-injection matches) ---
data_S_merge <- readRDS(file.path(data.path, data.file))
data_S_ls <- SplitObject(data_S_merge, split.by = "batch")
names(data_S_ls) <- sample_ls
data_S_ls <- lapply(data_S_ls, function(s) {
  s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>%
    RunPCA(npcs = 50, verbose = FALSE) %>% RunUMAP(dims = 1:10, seed = 42)
})
srat_base <- data_S_ls[[sample_diff]]

# ---- reproduce an injected object from a saved params.rds -----------------
reinject <- function(params) {
  inj <- with(as.list(params),
    Inject_local_signals_v2(srat_base, assay = "RNA", layer = "counts",
      reduc = "pca", npc = 10, min_knn = max(ncell_per_neighbors_ls),
      ncell_per_neighbors_ls = ncell_per_neighbors_ls,
      zinb_prob = zinb_prob, diff.pct = diff.pct,
      gene_params = gene_params, seed = seed,
      cov_strength = cov_strength, cov_K = cov_K,
      diffuse = diffuse, overlap_frac = overlap_frac))
  # Defensive: restore signal_genes on the assay in case SetAssayData dropped
  # meta.features (signal_cells is cell-level meta and always persists).
  gg <- params$gene_params$selected_gene_groups
  sg <- unlist(lapply(names(gg), function(nm) setNames(rep(nm, length(gg[[nm]])), gg[[nm]])))
  inj[["RNA"]]@meta.features$signal_genes <- NA
  inj[["RNA"]]@meta.features[names(sg), "signal_genes"] <- sg
  inj
}

# ---- loop over all reps x cov_strength ------------------------------------
cor_long  <- list()   # for the violin
heat_mats <- list()   # for the heatmap (heatmap_rep only)

for (cv in cov_range) {
  for (r in seq_len(n_reps)) {
    pf <- file.path(res.path, paste0("rep", r),
                    sprintf("%s_params%d.rds", cv, sample_diff))
    if (!file.exists(pf)) { warning("missing params: ", pf); next }
    params <- readRDS(pf)
    diag   <- Summarize_injected_covariance(reinject(params))

    od <- diag$pooled$offdiag
    if (length(od))
      cor_long[[paste(cv, r, sep = "_")]] <-
        data.frame(cov_strength = cv, rep = r, correlation = od)

    if (r == heatmap_rep && length(diag$cor_matrices))
      heat_mats[[as.character(cv)]] <- diag$cor_matrices[[1]]   # single module
  }
}

cor_df <- do.call(rbind, cor_long)
cor_df$cov_strength <- factor(cor_df$cov_strength, levels = as.character(cov_range))
write.csv(cor_df, file.path(out.path, "cov_offdiag_correlations.csv"), row.names = FALSE)

# quick numeric summary for the point-to-point reply
summ <- cor_df %>% group_by(cov_strength) %>%
  summarise(frac_negative = mean(correlation < 0),
            min = min(correlation), mean = mean(correlation), max = max(correlation),
            .groups = "drop")
write.csv(summ, file.path(out.path, "cov_correlation_summary.csv"), row.names = FALSE)
print(summ)

# ---- Plot 1: violin of off-diagonal correlations per cov_strength ---------
p_violin <- ggplot(cor_df, aes(cov_strength, correlation)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(aes(fill = cov_strength), scale = "width", alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.12, outlier.size = 0.25, alpha = 0.85) +
  scale_fill_viridis_d(guide = "none") +
  labs(x = "cov strength",
       y = "Designed pairwise cov\n(log-scale)") +
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"))
ggsave(file.path(figure.path, "figS8F_covariance_violin.png"), p_violin,
       width = 6, height = 4.5, dpi = 300)

# ---- Plot 2: heatmaps, one rep per cov_strength, shared gene order ---------
if (DRAW_HEATMAP && length(heat_mats) >= 1) {
  common_genes <- Reduce(intersect, lapply(heat_mats, rownames))
  if (length(common_genes) >= 2) {
    strongest <- heat_mats[[as.character(max(cov_range))]]
    strongest <- strongest[common_genes, common_genes]
    ord <- hclust(as.dist(1 - strongest))$order
    gene_order <- common_genes[ord]
    
    # global color range
    all_vals <- unlist(lapply(heat_mats, function(m) m[upper.tri(m)]))
    clim <- max(abs(all_vals))
    
    plot_list <- lapply(names(heat_mats), function(cv) {
      m <- heat_mats[[cv]][gene_order, gene_order]
      d <- reshape2::melt(m)
      colnames(d) <- c("gene1", "gene2", "value")
      d$value[d$gene1 == d$gene2] <- NA
      d$gene1 <- factor(d$gene1, levels = gene_order)
      d$gene2 <- factor(d$gene2, levels = rev(gene_order))
      
      ggplot(d, aes(gene1, gene2, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, limits = c(-clim, clim)) +
        coord_fixed() +
        labs(title = paste0("cov strength = ", cv), x = NULL, y = NULL, fill = "cov") +
        theme(
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank())    
      })
    
    p_heat <- wrap_plots(plot_list, nrow = 2) + plot_layout(guides = "collect")
    ggsave(file.path(figure.path, "figS8G_covariance_heatmap.png"), p_heat,
           width = 3 * 2, height = 3 * 2 + 1, dpi = 300)
  }
}

message("Done. Wrote to ", out.path,
        "\n  figS8F_covariance_violin.png\n  figS8G_covariance_heatmap.png\n",
        "  cov_offdiag_correlations.csv\n  cov_correlation_summary.csv")
