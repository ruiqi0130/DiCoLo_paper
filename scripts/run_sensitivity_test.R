# =============================================================================
# run_sensitivity_test.R
#
# Fig S7 - Robustness of DiCoLo to cell-graph construction parameters.
#
# For the SmoM2 dataset, we test the differential operator of SmoM2 AGAINST CTL
# (comp = SmoM2, proj = CTL) and measure how stable the set of significant
# differentially co-localized genes is when the CELL graph is built differently:
#   (Left)   varying number of diffusion-map components (npc), K = 10, dm
#   (Middle) varying cell-graph neighborhood size (K),        npc = 10, dm
#   (Right)  top-10 principal components (pca) vs the default top-10 dm
# Stability is the Jaccard index of the selected genes vs the default
# (npc = 10, K = 10, reduction = "dm"), coloured by the number of leading
# eigenvectors of the differential operator used for gene selection.
#
# The graph knobs live in ComputeGeneEMD(); the gene-graph operator knn is kept
# at its default (10) throughout, since S6 probes the cell graph, not the gene
# graph.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DiCoLo)          # ComputeGeneEMD, LoadGeneEMD, SelectCommonGenes,
                           # ComputeGraphOperator, ComputeDifferentialOperator,
                           # RunSVD, SelectSignificantGenes
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()
dir.path    <- DATA_DIR
figure.path <- FIGURE_DIR
# If DiCoLo is not installed as a package, point ComputeGeneEMD at the bundled
# OT backend by setting script_dir below (dir containing
# gene_distance_cal_parallel.py). Left NULL => use the installed package copy.
script_dir <- OT_SCRIPT_DIR

# ----------------------------------------------------------------------------
# 0. User settings -- ADJUST THESE to match your SmoM2 object
# ----------------------------------------------------------------------------
data.path     <- file.path(DATA_DIR, "smom2")
res.path      <- file.path(data.path, "sensitivity")   # outputs go here
condition_col <- "condition"                # metadata column holding condition
cond_levels   <- c("SmoM2", "CTL")          # comp first, proj second

lfdr_thresh   <- 0.2                         # locFDR threshold (package default)
n_ev_range    <- 1:5                         # leading eigenvectors (colour var)
npc_range     <- c(5, 10, 15, 20, 30)        # Left panel; 10 = default anchor
K_range       <- c(5, 10, 20, 30, 50)        # Middle panel; 10 = default anchor
wait_timeout_min <- 240                      # give up waiting on OT after this

dir.create(res.path, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 1. Load SmoM2 object, split into the two conditions, preprocess
#    (PCA is required so reduction = "pca" is available for the Right panel)
# ----------------------------------------------------------------------------
data_S <- readRDS(file.path(data.path, "smom2.rds"))
Idents(data_S) <- data_S[[condition_col]][, 1]
data_S_ls <- SplitObject(data_S, split.by = condition_col)
stopifnot(all(cond_levels %in% names(data_S_ls)))
data_S_ls <- data_S_ls[cond_levels]          # enforce comp = SmoM2, proj = CTL

data_S_ls <- lapply(data_S_ls, function(s) {
  s %>% NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
})

# Common gene universe: depends only on expression, not on the cell graph,
# so it is fixed across all configurations (makes the Jaccard well defined).
common_genes <- SelectCommonGenes(data_S_ls[[1]], data_S_ls[[2]], ngenes = 500)
saveRDS(common_genes, file.path(res.path, "common_genes.rds"))

# ----------------------------------------------------------------------------
# 2. Enumerate configurations. The default (npc=10, K=10, dm) is computed once
#    and re-used as the anchor of all three panels.
# ----------------------------------------------------------------------------
configs <- c(
  list(list(id = "default", reduction = "dm", npc = 10, K = 10)),
  lapply(setdiff(npc_range, 10), function(p)
    list(id = sprintf("npc%d", p), reduction = "dm", npc = p,  K = 10)),
  lapply(setdiff(K_range, 10), function(k)
    list(id = sprintf("K%d", k),   reduction = "dm", npc = 10, K = k)),
  list(list(id = "pca", reduction = "pca", npc = 10, K = 10))
)

emd_dir <- function(cfg, cond) file.path(res.path, "emd", cfg$id, cond)

# ----------------------------------------------------------------------------
# 3. Launch gene-EMD (OT) for every config x condition. ComputeGeneEMD writes
#    inputs then fires the Python OT backend in the background (nohup ... &),
#    so we launch all jobs, then wait for every emd.csv. Existing results are
#    skipped, so the script is resumable.
#    NOTE: these OT jobs are the expensive part -- ~2 x length(configs) runs on
#    the real SmoM2 data. Throttle by launching in batches if the box is small.
# ----------------------------------------------------------------------------
for (cfg in configs) {
  for (ci in seq_along(cond_levels)) {
    outdir <- emd_dir(cfg, cond_levels[ci])
    if (!file.exists(file.path(outdir, "emd.csv"))) {
      message(sprintf("Launching OT: config=%s cond=%s (npc=%d, K=%d, reduction=%s)",
                      cfg$id, cond_levels[ci], cfg$npc, cfg$K, cfg$reduction))
      ComputeGeneEMD(data_S_ls[[ci]], common_genes, dir.path = outdir,
                     script_dir = script_dir,
                     npc = cfg$npc, K = cfg$K, reduction = cfg$reduction)
    }
  }
}

expected <- unlist(lapply(configs, function(cfg)
  vapply(cond_levels, function(cond) file.path(emd_dir(cfg, cond), "emd.csv"),
         character(1))))

t0 <- Sys.time()
repeat {
  missing <- expected[!file.exists(expected)]
  if (length(missing) == 0) { message("All OT jobs finished."); break }
  if (as.numeric(difftime(Sys.time(), t0, units = "mins")) > wait_timeout_min) {
    warning(sprintf("Timeout: %d OT output(s) never appeared; affected configs will be skipped:\n%s",
                    length(missing), paste(missing, collapse = "\n")))
    break
  }
  Sys.sleep(10)
}

# ----------------------------------------------------------------------------
# 4. Per config: build gene graphs -> differential operator (SmoM2 vs CTL) ->
#    SVD -> significant genes for each number of leading eigenvectors.
# ----------------------------------------------------------------------------
select_genes_for_config <- function(cfg) {
  emd_ls <- lapply(cond_levels, function(cond)
    LoadGeneEMD(file.path(emd_dir(cfg, cond), "")))
  if (any(vapply(emd_ls, is.null, logical(1)))) {
    warning(sprintf("Config %s: missing EMD, skipping.", cfg$id)); return(NULL)
  }
  # align to shared genes across the two conditions
  g <- Reduce(intersect, lapply(emd_ls, rownames))
  emd_ls <- lapply(emd_ls, function(x) x[g, g])

  op_ls <- lapply(emd_ls, ComputeGraphOperator)   # gene-graph knn = 10 (default)
  Pdiff <- ComputeDifferentialOperator(op_ls[[1]], op_ls[[2]])  # comp=SmoM2, proj=CTL

  # Pdiff is PSD, so leading (largest) eigenvectors = SmoM2-against-CTL direction
  keep  <- min(max(n_ev_range), nrow(Pdiff) - 1)
  E     <- RunSVD(Pdiff, eig_keep = keep)
  V     <- E$vectors                              # genes x keep, gene-named rows

  setNames(lapply(n_ev_range, function(n) {
    Vn  <- V[, seq_len(n), drop = FALSE]
    ind <- SelectSignificantGenes(Vn, lfdr_thresh = lfdr_thresh)
    rownames(ind)[rowSums(ind) > 0]
  }), as.character(n_ev_range))
}

gene_sets <- setNames(lapply(configs, select_genes_for_config),
                      vapply(configs, `[[`, character(1), "id"))
saveRDS(gene_sets, file.path(res.path, "gene_sets.rds"))
stopifnot(!is.null(gene_sets[["default"]]))     # anchor must exist

# ----------------------------------------------------------------------------
# 5. Jaccard index of each config vs the default, per number of eigenvectors.
# ----------------------------------------------------------------------------
jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

df <- do.call(rbind, lapply(configs, function(cfg) {
  gs <- gene_sets[[cfg$id]]
  if (is.null(gs)) return(NULL)
  do.call(rbind, lapply(n_ev_range, function(n) {
    data.frame(config    = cfg$id,
               reduction = cfg$reduction,
               npc       = cfg$npc,
               K         = cfg$K,
               n_ev      = n,
               jaccard   = jaccard(gs[[as.character(n)]],
                                   gene_sets[["default"]][[as.character(n)]]),
               stringsAsFactors = FALSE)
  }))
}))
write.csv(df, file.path(res.path, "sensitivity_jaccard.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 6. Plot the three panels. Colour = number of leading eigenvectors.
# ----------------------------------------------------------------------------
df$n_ev <- factor(df$n_ev)

base_theme <- theme(
  legend.title  = element_text(size = 14),
  legend.text   = element_text(size = 12),
  axis.title    = element_text(size = 15),
  axis.text     = element_text(size = 12),
  panel.grid    = element_blank(),
  panel.background = element_blank(),
  axis.line     = element_line(colour = "black"))

# Left: diffusion-map components (dm, K = 10). npc = 10 is the self-anchor (=1).
p_npc <- df %>% filter(reduction == "dm", K == 10) %>%
  ggplot(aes(x = npc, y = jaccard, colour = n_ev, group = n_ev)) +
  geom_line() + geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "# diffusion map components", y = "Jaccard index", colour = "# eigenvectors") +
  base_theme

# Middle: cell-graph neighborhood size (dm, npc = 10). K = 10 is the anchor (=1).
p_K <- df %>% filter(reduction == "dm", npc == 10) %>%
  ggplot(aes(x = K, y = jaccard, colour = n_ev, group = n_ev)) +
  geom_line() + geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "cell-graph neighborhood size (K)", y = "Jaccard index", colour = "# eigenvectors") +
  base_theme

# Right: dm (default) vs pca, both top-10 components, K = 10.
p_pca <- df %>% filter(config %in% c("default", "pca")) %>%
  mutate(reduction = factor(reduction, levels = c("dm", "pca"))) %>%
  ggplot(aes(x = reduction, y = jaccard, colour = n_ev, group = n_ev)) +
  geom_line() + geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "cell embedding", y = "Jaccard index", colour = "# eigenvectors") +
  base_theme

p <- (p_npc | p_K | p_pca) + plot_layout(guides = "collect")
ggsave(file.path(res.path, "figS7_sensitivity.png"), p, width = 15, height = 5, dpi = 300)

message("Done. Wrote:\n  ", file.path(res.path, "sensitivity_jaccard.csv"),
        "\n  ", file.path(res.path, "figS7_sensitivity.png"))
