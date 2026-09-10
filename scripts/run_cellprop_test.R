# =============================================================================
# run_cellprop_test.R
# -----------------------------------------------------------------------------
# Composition-only control for DiCoLo (addresses Reviewer 1, Weakness 6):
# "Can DiCoLo distinguish TRUE differential co-localization from apparent
#  signals driven by shifts in cell-type / cell-state PROPORTIONS?"
#
# Design (your scheme):
#   * Use a SINGLE real condition (e.g. SmoM2-CTL, which has only LD/UD -> no DC
#     to be accidentally sampled away). Split its cells 50/50, STRATIFIED by
#     cell type, into two pseudo-conditions A and B. Because both halves are the
#     same real cells, the within-type transcriptional program is IDENTICAL by
#     construction -- this is what makes the "program unchanged" claim a
#     demonstration rather than a synthetic assumption.
#   * Sweep the proportion of a TARGET cell type (default UD) in B from 100% -> 0%
#     by stratified downsampling (only the target type is thinned; the LOCKED
#     type, default LD, is untouched in B).
#   * To remove the total-cell-number nuisance, A is downsampled to match B's
#     current TOTAL size while PRESERVING A's baseline type ratio (so A's UD:LD
#     stays fixed; only its absolute N shrinks). e.g. start 100/100 at 1:1;
#     B -> 50 (UD:LD = 0:1); A -> 50 keeping 1:1 (= 25:25).
#
# Read-out ("signal", y-axis):
#   * number of TARGET-type marker genes with FDR < cutoff (primary)
#   * number of ALL genes with FDR < cutoff (secondary)
#   * leading differential eigenvalue (secondary)
#
# Expected curve: signal ~ 0 for mild proportion changes (100% -> ~50%), and only
# rises as the target proportion approaches 0 (true presence/absence = a genuine
# co-localization change that SHOULD be flagged). This turns R1's assumed
# dichotomy into an empirical boundary.
#
# NOTE ON SIGNIFICANCE: the FDR procedure below (dicolo_gene_significance) is a
# self-contained default. For a figure that is numerically consistent with your
# R2 false-positive panel, REPLACE its body with the exact null/FDR function you
# used there (see the clearly-marked hook).
# =============================================================================

# ---- Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tibble)
  library(dplyr)
  library(igraph)
  library(pdist)
  library(reshape2)
  require(DiCoLo)
  require(ggplot2)
  require(Seurat)
  require(parallel)
})

# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()
dir.path    <- DATA_DIR
figure.path <- FIGURE_DIR
source("./code/simulation_functions.R")
source("./code/benchmarking_functions.R")

# ======================= CONFIG ==============================================
cfg <- list(
  data_rds     = file.path(DATA_DIR, "smom2", "data_S_smom2_dermal_E13.5_CTL.rds"),  # single real condition (CTL only)
  celltype_col = "celltype",                  # metadata column with LD/UD/(DC)
  target_type  = "UD",                         # the type whose proportion we sweep
  lock_type    = "LD",                         # locked (untouched in B)
  exclude_types = c("DC"),                      # drop rare types from THIS test
  fractions    = c(1.0, 0.75, 0.5, 0.25, 0.10, 0.05, 0.03, 0.0),  # target retained in B
  n_reps       = 10,
  base_seed    = 2026,
  ngenes_common = 500,       # HVGs for the DiCoLo gene set (marker genes added in)
  n_marker = 50,
  eig_agg      = "leading",  # "leading" (1st eigvec) or "topk" (knee, weighted)
  npc          = 10,
  out_dir      = file.path(DATA_DIR, "smom2", "cellprop_test")
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ======================= HELPERS =============================================
preprocess_for_emd <- function(srat) {
  srat %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
}

# Stratified 50/50 split of a single real condition into pseudo-conditions A/B.
# Returns cell-name vectors per type per side, so downstream sampling is exact.
make_stratified_halves <- function(cell_types) {
  by_type <- split(names(cell_types), cell_types)
  A <- list(); B <- list()
  for (ct in names(by_type)) {
    cells <- sample(by_type[[ct]])          # shuffle
    half  <- floor(length(cells) / 2)
    A[[ct]] <- cells[seq_len(half)]
    B[[ct]] <- cells[(half + 1):length(cells)]
  }
  list(A = A, B = B)
}

# Build the A/B cell sets for one proportion point.
#   mode = "sweep":    B keeps all lock_type, retains frac of target_type;
#                      A downsampled to |B| total, preserving A's baseline ratio.
#   mode = "nuisance": B keeps target:lock ratio fixed, total downsampled to the
#                      SAME |B| as the matched sweep point (isolates total-N effect).
build_AB_cells <- function(halves, target_type, lock_type, frac) {
  A0 <- halves$A; B0 <- halves$B
  
  n_tgt_B0 <- length(B0[[target_type]])
  n_lck_B0 <- length(B0[[lock_type]])
  
  n_tgt_B <- round(frac * n_tgt_B0)
  B_cells <- c(B0[[lock_type]],
               sample(B0[[target_type]], n_tgt_B))
  N_B <- length(B_cells)

  # A: downsample to N_B, preserving A's baseline type ratio (all types).
  a_counts0 <- sapply(A0, length)
  r_A <- a_counts0 / sum(a_counts0)
  n_A <- round(r_A * N_B)
  # fix rounding drift
  drift <- N_B - sum(n_A); if (drift != 0) n_A[which.max(n_A)] <- n_A[which.max(n_A)] + drift
  A_cells <- unlist(lapply(names(A0), function(ct)
    sample(A0[[ct]], min(n_A[[ct]], a_counts0[[ct]]))))

  list(A = A_cells, B = B_cells,
       n_target_B = length(intersect(B_cells, B0[[target_type]])),
       n_B = N_B,
       prop_target_B = length(intersect(B_cells, B0[[target_type]])) / N_B)
}

# ======================= LOAD DATA & DEFINE MARKERS ==========================
srat_ctl <- readRDS(cfg$data_rds)
stopifnot(cfg$celltype_col %in% colnames(srat_ctl@meta.data))

cell_types <- setNames(as.character(srat_ctl@meta.data[[cfg$celltype_col]]),
                       colnames(srat_ctl))
message("Cell-type counts used:"); print(table(cell_types))
stopifnot(cfg$target_type %in% cell_types, cfg$lock_type %in% cell_types)

# Fixed target-type marker set (computed once on the full condition -> stable
# across the sweep). These are the genes we monitor as "should NOT be flagged".
Idents(srat_ctl) <- cfg$celltype_col
srat_ctl <- NormalizeData(srat_ctl, verbose = FALSE)
tgt_cells <- colnames(srat_ctl)[cell_types == cfg$target_type]
expr_pct  <- Matrix::rowSums(
  GetAssayData(srat_ctl, assay = "RNA", layer = "data")[, tgt_cells] > 0
) / length(tgt_cells)
eligible  <- names(expr_pct)[expr_pct > 0.005 & expr_pct < 0.5]
mk <- FindMarkers(srat_ctl, ident.1 = cfg$target_type,
                  only.pos = TRUE, verbose = FALSE)
mk <- mk[rownames(mk) %in% eligible, ]
target_markers <- rownames(mk[order(mk$p_val_adj, -mk$avg_log2FC), ])[1:cfg$n_marker]
target_markers <- rownames(mk[order(mk$p_val_adj, -mk$avg_log2FC), ])[1:500]

message(sprintf("Monitoring %d %s markers.", length(target_markers), cfg$target_type))

# ======================= MAIN SWEEP ==========================================
results <- NULL
sample_diff = 1
for(frac in cfg$fractions){
  for(rep_id in seq_len(cfg$n_reps)){
    set.seed(cfg$base_seed + rep_id * 1000 + round(frac * 100))
    halves <- make_stratified_halves(cell_types)
    ab <- build_AB_cells(halves, cfg$target_type, cfg$lock_type, frac)
    if (length(ab$A) < 30 || length(ab$B) < 30) return(NULL)  # too few cells
    
    srat_A <- subset(srat_ctl, cells = ab$A)
    srat_B <- subset(srat_ctl, cells = ab$B)
    
    # Run EMD
    tmp_root <- file.path(cfg$out_dir, "emd_tmp",
                          sprintf("f%.2f_r%d", frac, rep_id))
    emd_paths <- file.path(tmp_root, paste0("GeneTrajectory", 1:2))
    
    if (all(file.exists(file.path(emd_paths, "emd.csv")))) {
      gene_emd_ls <- lapply(emd_paths, function(p) LoadGeneEMD(file.path(p, "")))
    } else{
      dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
      common_genes <- union(SelectCommonGenes(srat_A, srat_B, ngenes = cfg$ngenes_common),
                            target_markers[1:cfg$n_marker])
      data_S_ls <- list(preprocess_for_emd(srat_A), preprocess_for_emd(srat_B))
      emd_paths <- lapply(1:2, function(i) {
        tmp_path <- file.path(tmp_root, paste0("GeneTrajectory", i))
        if (!file.exists(file.path(tmp_path, "emd.csv"))) {
          ComputeGeneEMD(data_S_ls[[i]], common_genes, dir.path = tmp_path)
        }
        tmp_path
      })
      emd_paths <- unlist(emd_paths)
      
      repeat {
        still_running <- !sapply(emd_paths, function(p) file.exists(file.path(p, "emd.csv")))
        if (!any(still_running)) break
        Sys.sleep(10)
      }
      gene_emd_ls <- lapply(emd_paths, function(p) LoadGeneEMD(file.path(p, "")))
    }
    
    if (any(sapply(gene_emd_ls, is.null))) next
    
    # ---- align genes -> graph -> diff operator -> SVD ----
    g <- Reduce(intersect, lapply(gene_emd_ls, rownames))
    gene_emd_ls <- lapply(gene_emd_ls, function(x) x[g, g])
    gene_graph_ls <- lapply(gene_emd_ls, ComputeGraphOperator)
    # sample_diff = 2 -> look for structure present in B (condition 2) not A
    diff.op <- ComputeDifferentialOperator(gene_graph_ls[[sample_diff]],
                                           gene_graph_ls[[3 - sample_diff]])
    E.list <- RunSVD(diff.op, eig_keep = nrow(diff.op))
    rownames(E.list$vectors) <- g
   
    # ---- significance --------------------------------
    # robust z-scores on v1 (for reporting / loading-vs-null scatter)
    v   <- E.list$vectors[, 1]
    med <- median(v, na.rm = TRUE)
    s   <- mad(v, center = med, constant = 1.4826, na.rm = TRUE)
    if (s < 1e-12) s <- sd(v)
    z_scores <- (v - med) / s
    
    sig_mat <- tryCatch(
      SelectSignificantGenes(E.list$vectors[, 1, drop = FALSE],max_genes = nrow(E.list$vectors)),
      error = function(e) {
        cat("  locfdr failed (null-like distribution):", conditionMessage(e), "\n")
        matrix(0, nrow = nrow(E.list$vectors), ncol = 1)
      }
    )
    signf_genes <- rownames(sig_mat)[sig_mat[, 1] > 0]
    markers_in <- intersect(target_markers, g)
    
    results <- c(results,list(list(
      frac            = frac,
      prop_target_B   = ab$prop_target_B,     # UD as fraction of B (reviewer-facing x)
      n_B             = ab$n_B,
      rep             = rep_id,
      signf_genes = signf_genes,
      marker_tested = markers_in,
      E.list = E.list,
      z_scores = z_scores
    )))
  }
}

saveRDS(results, file = file.path(cfg$out_dir,"results.rds"))
# ==============STAT Table ==============================
null_bounds <- function(n, level = 0.975) {
  a_n <- sqrt(2 * log(n))
  b_n <- a_n - (log(log(n)) + log(4 * pi)) / (2 * a_n)
  maxz_env <- b_n - log(-log(level)) / a_n          # upper envelope for max|z|
  kurt_se  <- sqrt(24 / n)
  kurt_ub  <- qnorm(level) * kurt_se                # upper bound for excess kurtosis
  list(maxz_env = maxz_env, kurt_ub = kurt_ub)
}

excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x); s <- mean((x - m)^2)
  mean((x - m)^4) / s^2 - 3
}


top_ns <- c(10:50,100,300,500)
stat_table <- do.call(rbind, lapply(results, function(r) {
  z <- r$z_scores
  z <- z[is.finite(z)]
  n <- length(z)
  if (n < 3) return(NULL)
  
  nb        <- null_bounds(n, level = 0.975)
  max_absz  <- max(abs(z))
  exc_kurt  <- excess_kurtosis(z)
  
  # --- AUPRC for each top-N marker set ---
  V     <- r$E.list$vectors
  de <- data.frame(gene = rownames(V), score = abs(V[, 1]),
                   stringsAsFactors = FALSE)
  de <- de %>% arrange(desc(score)) %>% distinct(gene, .keep_all = TRUE)
  de$rank <- seq_len(nrow(de))
  real_score <- setNames(-de$rank, de$gene)     # higher = better, matches benchmark
  ng <- nrow(de)
  rank_pct <- setNames(de$rank / ng, de$gene)
  
  # per top-N: median rank percentile of markers present in this gene set
  rankpct_vals <- sapply(top_ns, function(N) {
    gt <- head(target_markers, N); gt <- gt[gt %in% de$gene]
    if (length(gt) < 1) return(NA_real_)
    median(rank_pct[gt])
  })
  names(rankpct_vals) <- paste0("rankpct_top", top_ns)
  
  # AUPRC / AUROC
  auc_vals <- sapply(top_ns, function(N) {
    gt <- head(target_markers, N)
    gt <- gt[gt %in% de$gene]
    if (length(gt) < 1) return(NA_real_)
    get_auc(real_score = real_score, gt_gene_ls = gt, metric = "auprc", plot = FALSE)
    # get_auc(real_score = real_score, gt_gene_ls = gt, metric = "auroc", plot = FALSE)
  })
  names(auc_vals) <- paste0("auc_top", top_ns)

  data.frame(
    frac          = r$frac,
    prop_target_B = r$prop_target_B,
    rep           = r$rep,
    n_B           = r$n_B,
    n_genes       = n,
    n_sig         = length(r$signf_genes),
    n_sig_marker  = sum(r$signf_genes %in% r$marker_tested),
    max_absz      = max_absz,
    excess_kurt   = exc_kurt,
    n_exceed_env  = sum(abs(z) > nb$maxz_env),
    maxz_env95    = nb$maxz_env,
    kurt_ub95     = nb$kurt_ub,
    first_eigen_val = r$E.list$values[1],
    as.list(rankpct_vals),
    as.list(auc_vals),
    stringsAsFactors = FALSE
  )
}))

write.csv(stat_table, file = file.path(cfg$out_dir,"stat_table.csv"),row.names = FALSE)
# ======================= PLOT ================================================
stat_table$prop_bin <- factor(round(stat_table$prop_target_B, 3))  
p <- ggplot(stat_table, aes(x = prop_bin, y = rankpct_top100)) +
  geom_boxplot(outlier.size = 0.8, width = 0.6) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.4) +
  scale_y_reverse(limits = c(1, 0)) +
  labs(x = sprintf("%s proportion in condition B", cfg$target_type),
       y = "Median rank percentile of top50 UD markers\n(top = ranked first, bottom = ranked last)") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))
ggsave(file.path(FIGURE_DIR, "figS11_celltype_abundance.png"),
       p, width = 10, height = 8)

long_df <- stat_table %>%
  select(prop_target_B, rep, starts_with("rankpct_top")) %>%
  pivot_longer(starts_with("rankpct_top"), names_to = "topN", values_to = "rankpct") %>%
  mutate(topN = factor(as.integer(sub("rankpct_top", "", topN)), levels = top_ns))

summ <- long_df %>% group_by(prop_target_B, topN) %>%
  summarise(med = median(rankpct, na.rm = TRUE), .groups = "drop")

p2 <- ggplot(summ, aes(prop_target_B, med, color = topN, group = topN)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
  scale_x_reverse() + scale_y_reverse(limits = c(1, 0)) +
  scale_color_viridis_d(end = 0.9) +
  labs(x = "UD proportion in condition B",
       y = "Median rank percentile of UD markers", color = "top-N") +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 15), axis.text = element_text(size = 12))


