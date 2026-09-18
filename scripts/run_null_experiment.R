# =============================================================================
# Null experiment: DiCoLo on pbmc10k with no injected signal
# Goal: show that under null (same condition vs itself), locfdr returns ~0 genes
# =============================================================================

suppressPackageStartupMessages(
  {
    library(DiCoLo)
    library(SingleCellExperiment)
    library(tibble)
    library(dplyr)
    library(igraph)
    library(pdist)
    library(reshape2)
    require(GeneTrajectory)
    require(plotly)
    require(ggplot2)
    library(ggplotify)
    require(patchwork)
    require(viridis)
    require(scales)
    require(Seurat)
    require(parallel)
    require(SeuratWrappers)
    require(RColorBrewer)
  }
)

# --- EDIT THESE PATHS TO MATCH YOUR ENVIRONMENT ---
# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()
dir.path    <- DATA_DIR
figure.path <- FIGURE_DIR
source("./code/simulation_functions.R")
# =============================================================================
# NULL TYPE 1: Random split of ONE sample into two halves------
# This is the cleanest null: same cells, same batch, randomly partitioned.
# EMD must be computed per split since cells differ.
# =============================================================================

# --- Load pbmc10k ---
data.path <- file.path(dir.path, "pbmc10k")
res.path <- file.path(data.path, "null_experiment")
if (!dir.exists(res.path)) dir.create(res.path, recursive = TRUE)

data_S_merge <- readRDS(file.path(data.path, "pbmc10k.rds"))

n_splits <- 10
null_results_split <- list()

for (split_id in seq_len(n_splits)) {
  cat("\n=== Random split", split_id, "===\n")
  set.seed(split_id)
  # Use monocyte1 (or whichever is larger)
  srat <- data_S_merge
  cells <- colnames(srat)
  half <- sample(seq_along(cells), floor(length(cells) / 2))
  
  srat_A <- subset(srat, cells = cells[half])
  srat_B <- subset(srat, cells = cells[-half])
  
  # Re-preprocess each half
  srat_A <- srat_A %>% NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
  srat_B <- srat_B %>% NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
  common_genes = SelectCommonGenes(srat_A, srat_B)
  # Compute gene EMD for each half
  emd_path_A <- file.path(res.path, paste0("split", split_id, "_A"))
  emd_path_B <- file.path(res.path, paste0("split", split_id, "_B"))
  
  if (!file.exists(file.path(emd_path_A, "emd.csv"))) {
    ComputeGeneEMD(srat_A, common_genes, dir.path = emd_path_A)
  }
  if (!file.exists(file.path(emd_path_B, "emd.csv"))) {
    ComputeGeneEMD(srat_B, common_genes, dir.path = emd_path_B)
  }
  message("Run on backend\n")
  
  # Wait for all jobs
  repeat {
    still_running <- !sapply(c(emd_path_A,emd_path_B), function(p) file.exists(file.path(p, "emd.csv")))
    if (!any(still_running)) break
    Sys.sleep(10)
  }
  message("All jobs finished. Proceeding...")
  
  # Load EMD
  gene_emd_A <- LoadGeneEMD(file.path(emd_path_A, ""))
  gene_emd_B <- LoadGeneEMD(file.path(emd_path_B, ""))
  
  # Align genes
  g <- intersect(rownames(gene_emd_A), rownames(gene_emd_B))
  gene_emd_A <- gene_emd_A[g, g]
  gene_emd_B <- gene_emd_B[g, g]
  
  # Graph operators
  graph_A <- ComputeGraphOperator(gene_emd_A)
  graph_B <- ComputeGraphOperator(gene_emd_B)
  
  # Differential operator: A vs B (should be ~null)
  diff.op <- ComputeDifferentialOperator(graph_A, graph_B)
  E.list <- RunSVD(diff.op, eig_keep = nrow(diff.op))
  sig_mat <- SelectSignificantGenes(E.list$vectors[, 1, drop = FALSE])
  signf_genes = rownames(sig_mat)[sig_mat[,1] > 0]
  
  v = E.list$vectors[, 1]
  med <- median(v, na.rm = TRUE)
  s <- mad(v, center = med, constant = 1.4826, na.rm = TRUE)
  if (s < 1e-12) s <- sd(v)
  z_scores <- (v - med)/s
  
  null_results_split[[split_id]] <- list(
    dataset = split_id,
    top_eigenvalues = E.list$values[1:min(10, length(E.list$values))],
    z_scores_v1 = z_scores,
    signf_genes_v1 = signf_genes
  )
  
  # cat("  Top 5 eigenvalues:", round(E.list$values[1:5], 4), "\n")
  # cat("  Knee selected components:", n_eigvec, "\n")
  # cat("  Significant genes (v1 locfdr):", n_sig_genes, "\n")
}


# =============================================================================
# Batch-only Null (pbmcsca, pancreas) --------
# =============================================================================

datasets <- list(
  pbmcsca = list(
    path = file.path(dir.path, "pbmcsca"),
    file = "pbmcsca.rds",
    samples = c("Drop-seq", "inDrops")
  ),
  panc8 = list(
    path = file.path(dir.path, "panc8"),
    file = "panc8.rds",
    samples = c("human2", "human3")
  )
)

batch_null_results <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  cat("\n==========", ds_name, "==========\n")
  
  res.path <- file.path(ds$path, "null_experiment")
  if (!dir.exists(res.path)) dir.create(res.path, recursive = TRUE)
  
  # Load and preprocess
  data_S_merge <- readRDS(file.path(ds$path, ds$file))
  data_S_ls <- SplitObject(data_S_merge, split.by = "batch")
  names(data_S_ls) <- ds$samples
  
  data_S_ls <- lapply(data_S_ls, function(data_S) {
    data_S %>% NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(npcs = 50, verbose = FALSE)
  })
  
  # Common genes (computed fresh for this pair)
  common_genes <- SelectCommonGenes(data_S_ls[[1]], data_S_ls[[2]], ngenes = 500)
  cat("Common genes:", length(common_genes), "\n")
  
  # Compute EMD for each sample
  emd_paths <- lapply(seq_along(data_S_ls), function(i) {
    emd_path <- file.path(res.path, paste0("GeneTrajectory_", ds$samples[i]))
    if (!file.exists(file.path(emd_path, "emd.csv"))) {
      cat("  Computing EMD for", ds$samples[i], "...\n")
      ComputeGeneEMD(data_S_ls[[i]], common_genes, dir.path = emd_path)
    }
    emd_path
  })
  
  # Check if EMD is ready
  repeat {
    still_running <- !sapply(emd_paths, function(p) file.exists(file.path(p, "emd.csv")))
    if (!any(still_running)) break
    Sys.sleep(10)
  }
  message("All jobs finished. Proceeding...")
  
  # Load EMD
  gene_emd_ls <- lapply(emd_paths, function(p) LoadGeneEMD(file.path(p, "")))
  
  if (any(sapply(gene_emd_ls, is.null))) {
    cat("  Failed to load EMD for", ds_name, "\n")
    next
  }
  
  # Align genes
  g <- Reduce(intersect, lapply(gene_emd_ls, rownames))
  gene_emd_ls <- lapply(gene_emd_ls, function(x) x[g, g])
  cat("  Genes after alignment:", length(g), "\n")
  
  # Graph operators
  gene_graph_ls <- lapply(gene_emd_ls, ComputeGraphOperator)
  
  # Differential operator in both directions
  for (direction in c("1vs2", "2vs1")) {
    if (direction == "1vs2") {
      diff.op <- ComputeDifferentialOperator(gene_graph_ls[[1]], gene_graph_ls[[2]])
      label <- paste0(ds$samples[1], " vs ", ds$samples[2])
    } else {
      diff.op <- ComputeDifferentialOperator(gene_graph_ls[[2]], gene_graph_ls[[1]])
      label <- paste0(ds$samples[2], " vs ", ds$samples[1])
    }
    
    E.list <- RunSVD(diff.op, eig_keep = nrow(diff.op))
    sig_mat <- SelectSignificantGenes(E.list$vectors[, 1, drop = FALSE])
    signf_genes = rownames(sig_mat)[sig_mat[,1] > 0]
    
    v = E.list$vectors[, 1]
    med <- median(v, na.rm = TRUE)
    s <- mad(v, center = med, constant = 1.4826, na.rm = TRUE)
    if (s < 1e-12) s <- sd(v)
    z_scores <- (v - med)/s
    
    batch_null_results[[paste0(ds_name, "_", direction)]] <- list(
      dataset = ds_name,
      direction = label,
      top_eigenvalues = E.list$values[1:min(10, length(E.list$values))],
      z_scores_v1 = z_scores,
      signf_genes_v1 = signf_genes
    )
  }
}

# ---------------------------
# Inject signal --------------
# ---------------------------
# default parameters (ngene = 15); one representative replicate per dataset
rep_id = 5
inject_config <- list(
  list(dataset = "pbmc10k",  dir_id = "1vs0", rep = rep_id),
  list(dataset = "pbmcsca",  dir_id = "1vs0", rep = rep_id),
  list(dataset = "pbmcsca",  dir_id = "0vs1", rep = rep_id),
  list(dataset = "panc8",    dir_id = "1vs0", rep = rep_id),
  list(dataset = "panc8",    dir_id = "0vs1", rep = rep_id)
)

inject_results <- list()
ngene = 50
for (cfg in inject_config) {
  cat("\n===", cfg$dataset, cfg$dir_id, "rep", cfg$rep, "===\n")
  
  rep_path <- file.path(dir.path, cfg$dataset, "parameters", "ngene",
                        cfg$dir_id, paste0("rep", cfg$rep))
  
  # --- load injected (ground truth) genes ---
  inj_idx <- if (cfg$dir_id == "1vs0") 1 else 2
  params_file <- file.path(rep_path, sprintf("%s_params%d.rds", ngene,inj_idx))
  
  # injected sample carries the "50_" prefix
  if (cfg$dir_id == "1vs0") {
    emd_path_inj <- file.path(rep_path, paste0(ngene,"_GeneTrajectory1"))
    emd_path_ref <- file.path(rep_path, "GeneTrajectory2")
  } else {
    emd_path_inj <- file.path(rep_path, paste0(ngene,"_GeneTrajectory2"))
    emd_path_ref <- file.path(rep_path, "GeneTrajectory1")
  }
  
  gene_emd_inj <- LoadGeneEMD(file.path(emd_path_inj, ""))
  gene_emd_ref <- LoadGeneEMD(file.path(emd_path_ref, ""))
  
  if (is.null(gene_emd_inj) || is.null(gene_emd_ref)) {
    cat("  EMD not found - check paths:\n  ", emd_path_inj, "\n  ", emd_path_ref, "\n")
    next
  }
  
  g <- intersect(rownames(gene_emd_inj), rownames(gene_emd_ref))
  graph_inj <- ComputeGraphOperator(gene_emd_inj[g, g])
  graph_ref <- ComputeGraphOperator(gene_emd_ref[g, g])
  
  # injected condition as condition 1, so Q^(1) highlights the injected structure
  diff.op <- ComputeDifferentialOperator(graph_inj, graph_ref)
  E.list  <- RunSVD(diff.op, eig_keep = nrow(diff.op))
  
  sig_mat <- tryCatch(
    SelectSignificantGenes(E.list$vectors[, 1, drop = FALSE], max_genes = 100),
    error = function(e) {
      cat("  locfdr failed (null-like distribution):", conditionMessage(e), "\n")
      matrix(0, nrow = nrow(E.list$vectors), ncol = 1)
    }
  )
  signf_genes <- rownames(sig_mat)[sig_mat[, 1] > 0]
  
  v   <- E.list$vectors[, 1]
  med <- median(v, na.rm = TRUE)
  s   <- mad(v, center = med, constant = 1.4826, na.rm = TRUE)
  if (s < 1e-12) s <- sd(v)
  z_scores <- (v - med) / s
  
  injected_genes <- NULL
  if (file.exists(params_file)) {
    params <- readRDS(params_file)
    injected_genes <- unlist(params$gene_params$selected_gene_groups)
    cat("  injected genes:", length(injected_genes),
        "| present in loadings:", sum(injected_genes %in% names(z_scores)), "\n")
  } else {
    cat("  params file not found:", params_file, "\n")
  }
  
  key <- paste0(cfg$dataset, "_", cfg$dir_id)
  inject_results[[key]] <- list(
    dataset        = cfg$dataset,
    direction      = cfg$dir_id,
    rep            = cfg$rep,
    n_genes        = length(g),
    z_scores_v1    = z_scores,
    signf_genes_v1 = signf_genes,
    injected_genes = injected_genes
  )
  
}

# =============================================================================
# Real signal: SmoM2 and WlsKO -----
# =============================================================================
datasets <- list(
  smom2 = list(
    path = file.path(dir.path, "smom2"),
    file = c("data_S_smom2_dermal_E13.5_MUT.rds","data_S_smom2_dermal_E13.5_CTL.rds"),
    samples = c("MUT", "CTL")
  ),
  wlsko = list(
    path = file.path(dir.path, "wlsko"),
    file = c("data_S_wlsko_dermal_E14.5_MUT.rds","data_S_wlsko_dermal_E14.5_CTL.rds"),
    samples = c("MUT", "CTL")
  )
)

real_signal_results <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  cat("\n==========", ds_name, "==========\n")
  
  res.path <- file.path(ds$path, "parameters")
  if (!dir.exists(res.path)) dir.create(res.path, recursive = TRUE)
  
  # Load and preprocess
  data_S_ls <- lapply(ds$file,function(name){
    readRDS(file.path(ds$path, name))
  })
  names(data_S_ls) <- ds$samples
  
  # Common genes (computed fresh for this pair)
  common_genes <- SelectCommonGenes(data_S_ls[[1]], data_S_ls[[2]], ngenes = 500)
  cat("Common genes:", length(common_genes), "\n")
  
  # Compute EMD for each sample
  emd_paths <- lapply(seq_along(data_S_ls), function(i) {
    emd_path <- file.path(res.path, paste0("GeneTrajectory_", ds$samples[i]))
    if (!file.exists(file.path(emd_path, "emd.csv"))) {
      cat("  Computing EMD for", ds$samples[i], "...\n")
      ComputeGeneEMD(data_S_ls[[i]], common_genes, dir.path = emd_path)
    }
    emd_path
  })
  
  # Check if EMD is ready
  repeat {
    still_running <- !sapply(emd_paths, function(p) file.exists(file.path(p, "emd.csv")))
    if (!any(still_running)) break
    Sys.sleep(10)
  }
  message("All jobs finished. Proceeding...")
  
  # Load EMD
  gene_emd_ls <- lapply(emd_paths, function(p) LoadGeneEMD(file.path(p, "")))
  
  if (any(sapply(gene_emd_ls, is.null))) {
    cat("  Failed to load EMD for", ds_name, "\n")
    next
  }
  gene_emd_mut = gene_emd_ls[[1]]
  gene_emd_ctl = gene_emd_ls[[2]]

  # Align genes
  g <- intersect(rownames(gene_emd_mut), rownames(gene_emd_ctl))
  gene_emd_mut <- gene_emd_mut[g, g]
  gene_emd_ctl <- gene_emd_ctl[g, g]
  cat("  Genes:", length(g), "\n")
  
  # Graph operators
  graph_mut <- ComputeGraphOperator(gene_emd_mut)
  graph_ctl <- ComputeGraphOperator(gene_emd_ctl)
  
  # Both directions
  for (direction in c("MUT_vs_CTL", "CTL_vs_MUT")) {
    if (direction == "MUT_vs_CTL") {
      diff.op <- ComputeDifferentialOperator(graph_mut, graph_ctl)
    } else {
      diff.op <- ComputeDifferentialOperator(graph_ctl, graph_mut)
    }
    
    E.list <- RunSVD(diff.op, eig_keep = nrow(diff.op))
    pi1 <- compute_pi1(E.list, n_ev = 1)
    sig_mat <- tryCatch(
      SelectSignificantGenes(E.list$vectors[, 1, drop = FALSE]),
      error = function(e) {
        cat("  locfdr failed (null-like distribution):", conditionMessage(e), "\n")
        matrix(0, nrow = nrow(E.list$vectors), ncol = 1)
      }
    )
    signf_genes = rownames(sig_mat)[sig_mat[,1] > 0]
    
    v = E.list$vectors[, 1]
    med <- median(v, na.rm = TRUE)
    s <- mad(v, center = med, constant = 1.4826, na.rm = TRUE)
    if (s < 1e-12) s <- sd(v)
    z_scores <- (v - med)/s
    
    real_signal_results[[paste0(ds_name, "_", direction)]] <- list(
      dataset = ds_name,
      direction = direction,
      top_eigenvalues = E.list$values[1:min(10, length(E.list$values))],
      z_scores_v1 = z_scores,
      pi1_v1 = pi1,
      signf_genes_v1 = signf_genes
    )
  }
}


# =============================================================================
# Supplementary Figure: z-score panels -----
# =============================================================================
plot_zscore_panel <- function(z_scores, title = NULL, 
                              injected_genes = NULL, selected_genes = NULL, 
                              conf = 0.95) {
  n <- length(z_scores)
  
  # reference for max|Z| under iid N(0,1); max|Z| over n vars ~ max Z over 2n vars
  m <- 2 * n
  a_m <- sqrt(2 * log(m))
  b_m <- a_m - (log(log(m)) + log(4 * pi)) / (2 * a_m)
  max_z_hi <- b_m - log(-log(conf)) / a_m
  
  max_z <- max(abs(z_scores), na.rm = TRUE)
  n_outside <- sum(abs(z_scores) > max_z_hi)
  
  subtitle <- sprintf("max|z|: %.2f\nselected: %d / %d", max_z,length(selected_genes),length(z_scores))
  
  inj_df <- NULL
  if (!is.null(injected_genes)) {
    z_inj <- z_scores[names(z_scores) %in% injected_genes]
    if (length(z_inj) > 0) inj_df <- data.frame(z = z_inj)
    subtitle <- sprintf("%s\nof which %d / %d injected",
                        subtitle,sum(injected_genes %in% selected_genes),length(injected_genes))
  }
  
  df <- data.frame(z = z_scores)
  z_range <- seq(min(z_scores) - 0.5, max(z_scores) + 0.5, length.out = 300)
  null_df <- data.frame(z = z_range, density = dnorm(z_range))
  
  p <- ggplot(df, aes(x = z)) +
    # annotate("rect", xmin = -max_z_hi, xmax = max_z_hi, ymin = -Inf, ymax = Inf,
    #          fill = "grey90", alpha = 0.5) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.3,
                   fill = "grey60", color = "white", linewidth = 0.2) +
    geom_line(data = null_df, aes(x = z, y = density), color = "#2ca02c", linewidth = 1) 
  #+  geom_vline(xintercept = c(-max_z_hi, max_z_hi), color = "#d62728",
    #            linetype = "dashed", linewidth = 0.4)
  
  if (!is.null(inj_df) && nrow(inj_df) > 0) {
    p <- p + geom_rug(data = inj_df, aes(x = z), inherit.aes = FALSE,
                      color = "#d62728", linewidth = 0.8, length = unit(0.10, "npc"))
  }
  
  p +
    labs(title = title, subtitle = subtitle, x = "z-score", y = "density") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 8),
          plot.subtitle = element_text(size = 6.5, color = "grey40"))
}

# Row 1: Null splits (10 panels)
null_panels <- lapply(seq_along(null_results_split), function(i) {
  res <- null_results_split[[i]]
  plot_zscore_panel(res$z_scores_v1, 
                    sprintf("pbmc10k (random split %s)", i),
                    selected_genes = res$signf_genes_v1)
})

# Row 2: Batch only
batch_panels <- lapply(names(batch_null_results), function(nm) {
  res <- batch_null_results[[nm]]
  ds_label <- ifelse(res$dataset == "panc8", "pancreas", "pbmcsca")
  plot_zscore_panel(res$z_scores_v1, 
                    paste0(ds_label, ": ", res$direction),
                    selected_genes = res$signf_genes_v1)
})

# Row 3: Signal (3 panels)
signal_panels <- lapply(names(inject_results), function(nm) {
  res <- inject_results[[nm]]
  ds_label <- ifelse(res$dataset == "panc8", "pancreas", res$dataset)
  if(ds_label == "pancreas"){
    res$'direction' = ifelse(res$'direction' == "1vs0",
                             "human2 vs human3","human3 vs human2")
  }
  if(ds_label == "pbmcsca"){
    res$'direction' = ifelse(res$'direction' == "1vs0",
                             "Drop-seq vs inDrops","inDrops vs Drop-seq")
  }
  title = ifelse(ds_label == "pbmc10k","pbmc10k (random split)",paste0(ds_label, ": ", res$'direction'))
  plot_zscore_panel(res$z_scores_v1, 
                    title,
                    res$injected_genes,
                    selected_genes = res$signf_genes_v1
                    )
})

# Assemble
row1 <- wrap_plots(null_panels, ncol = 5)
row2 <- wrap_plots(batch_panels, ncol = 5)
row3 <- wrap_plots(signal_panels, ncol = 5)

row_label <- function(txt) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = txt, angle = 90,
             size = 3.2, fontface = "bold") +
    theme_void()
}


lab_w <- 0.03  # label column width relative to plots

fig <- (
  (row_label("Random splits")  | row1) + plot_layout(widths = c(lab_w, 1))
) / (
  (row_label("Real batches")           | row2) + plot_layout(widths = c(lab_w, 1))
) / (
  (row_label("Injected signal")      | row3) + plot_layout(widths = c(lab_w, 1))
) + plot_layout(heights = c(2, 1, 1))


# fig <- (row1 / row2 / row3)  +
#   patchwork::plot_layout(heights = c(2, 1, 1))

ggsave(file.path(FIGURE_DIR, "figS10_null_and_negative_control.png"),
       fig, width = 13, height = 8)


