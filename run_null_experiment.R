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
setwd("/data/ruiqi/DiCoLo_paper")
source("./code/simulation_functions.R")
compute_pi1 <- function(E.list, n_ev = 1) {
  sapply(seq_len(n_ev), function(j) {
    v <- E.list$vectors[, j]
    med <- median(v, na.rm = TRUE)
    s <- mad(v, center = med, constant = 1.4826, na.rm = TRUE)
    if (s < 1e-12) s <- sd(v)
    z <- (v - med) / s
    lf <- tryCatch(
      locfdr(z, nulltype = 0, plot = 0),
      error = function(e) NULL
    )
    if (is.null(lf)) return(NA)
    1 - lf$fp0["mlest", "p0"]
  })
}

dir.path <- "DiCoLo_data"

# =============================================================================
# NULL TYPE 1: Random split of ONE sample into two halves (repeat 10 times) ---
# This is the cleanest null: same cells, same batch, randomly partitioned.
# EMD must be computed per split since cells differ.
# =============================================================================

data.path <- file.path(dir.path, "pbmc10k")
res.path <- file.path(data.path, "null_experiment")
if (!dir.exists(res.path)) dir.create(res.path, recursive = TRUE)

# --- Load pbmc10k ---
data_S_merge <- readRDS(file.path(data.path, "pbmc10k.rds"))

n_splits <- 10
null_results_split <- list()

for (split_id in seq_len(n_splits)) {
  cat("\n=== Random split", split_id, "===\n")
  set.seed(42 + split_id)
  
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
  pi1 <- compute_pi1(E.list, n_ev = 1)
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
    pi1_v1 = pi1,
    signf_genes_v1 = signf_genes
  )
  
  # cat("  Top 5 eigenvalues:", round(E.list$values[1:5], 4), "\n")
  # cat("  Knee selected components:", n_eigvec, "\n")
  # cat("  Significant genes (v1 locfdr):", n_sig_genes, "\n")
}


# =============================================================================
# NULL TYPE 2: Batch-only null (pbmcsca, pancreas) --------
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
    pi1 <- compute_pi1(E.list, n_ev = 1)
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
      pi1_v1 = pi1,
      signf_genes_v1 = signf_genes
    )
  }
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

