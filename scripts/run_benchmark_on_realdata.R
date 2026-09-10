# =============================================================================
# Real-data benchmark: Sennett DC signature recovery
# Uses independent biological reference (Sennett et al. 2015 DC signature).
#
# Contrast:
#   SmoM2 vs CTL (E13.5): rank toward SmoM2 (precocious DC formation)
# =============================================================================

# Load libraries ----
suppressPackageStartupMessages({
  library(DiCoLo)
  library(SingleCellExperiment)
  library(tibble)
  library(dplyr)
  library(igraph)
  library(reshape2)
  require(GeneTrajectory)
  require(ggplot2)
  require(patchwork)
  require(viridis)
  require(scales)
  require(Seurat)
  require(parallel)
  require(RColorBrewer)
})

# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()
dir.path    <- DATA_DIR
figure.path <- FIGURE_DIR
source("./R/benchmarking_functions.R")
# =============================================================================
# 0. Independent reference: Sennett et al. 2015 DC signature
# =============================================================================
# Data downloaded from https://ars.els-cdn.com/content/image/1-s2.0-S153458071500430X-mmc2.xlsx
sennett_dc_genes <- readxl::read_excel(file.path(DATA_DIR, "smom2","Sennett_gene_list.xlsx"),sheet = "DC") %>% as.data.frame()
sennett_dc_genes <- sennett_dc_genes[5:nrow(sennett_dc_genes),2] 
cat(sprintf("Sennett DC signature: %d genes loaded\n", length(sennett_dc_genes)))

# =============================================================================
# 1. Define datasets
# =============================================================================
# Each entry: list of two Seurat objects + direction for ranking + label
datasets <- list(
  smom2 = list(
    file1 = "data_S_smom2_dermal_E13.5_MUT.rds",  # condition1 = MUT
    file2 = "data_S_smom2_dermal_E13.5_CTL.rds",   # condition2 = CTL
    direc = "condition1",  # rank toward MUT (precocious DC)
    label = "SmoM2 vs CTL",
    n_eig = 1
  )
  # wlsko = list(
  #   file1 = "data_S_wlsko_dermal_E14.5_MUT.rds", # condition1 = WlsKO
  #   file2 = "data_S_wlsko_dermal_E14.5_CTL.rds",    # condition2 = CTL
  #   direc = "condition2",  # rank toward CTL (DC present in CTL, absent in Wls)
  #   label = "Wls vs CTL",
  #   n_eig = 3
  # )
)

dir.path.data <- DATA_DIR
method_ls <- c("DiCoLo", "milode", "lemur", "DGCA", "memento")

# =============================================================================
# 2. Run all methods on each dataset
# =============================================================================
all_results <- lapply(names(datasets), function(dataset_name) {
  cat("\n========================================\n")
  cat(sprintf("Processing: %s\n", datasets[[dataset_name]]$label))
  cat("========================================\n")
  
  ds <- datasets[[dataset_name]]
  # --- Check if method results already exist ---
  res.path <- file.path(dir.path.data, dataset_name, "parameters")
  res.file <- file.path(res.path, "method_res.rds")
  if (file.exists(res.file)) {
    cat("Loading existing method results...\n")
    method_res <- readRDS(res.file)
    res_DiCoLo <- method_res$res_DiCoLo
    res_milode <- method_res$res_milode
    res_lemur <- method_res$res_lemur
    res_DGCA <- method_res$res_DGCA
    res_memento <- method_res$res_memento
    common_genes_dicolo <- method_res$common_genes
    gt_genes <- method_res$gt_genes
  } else {
    # --- Load data ---
    data_S1 <- readRDS(file.path(dir.path.data, dataset_name, ds$file1))
    data_S2 <- readRDS(file.path(dir.path.data, dataset_name, ds$file2))
    
    # --- Select common genes (same as DiCoLo demo) ---
    common_genes <- SelectCommonGenes(data_S1, data_S2, ngenes = 500)
    cat(sprintf("Common genes: %d\n", length(common_genes)))
    
    # Restrict Sennett signature to common gene universe
    gt_genes <- intersect(sennett_dc_genes, common_genes)
    cat(sprintf("Sennett DC genes in common universe: %d / %d\n", 
                length(gt_genes), length(sennett_dc_genes)))
    
    if (length(gt_genes) < 5) {
      warning(sprintf("Too few reference genes (%d) in %s — skipping", 
                      length(gt_genes), dataset_name))
      return(NULL)
    }
    
    # ===================
    # DiCoLo
    # ===================
    cat("Running DiCoLo...\n")
    emd.path1 <- file.path(dir.path.data, dataset_name, "parameters", "GeneTrajectory_MUT")
    emd.path2 <- file.path(dir.path.data, dataset_name, "parameters", "GeneTrajectory_CTL")
    
    # Compute gene OT distances (runs in background)
    if (!file.exists(file.path(emd.path1, "emd.csv"))) {
      ComputeGeneEMD(data_S1, common_genes, dir.path = emd.path1)
    }
    if (!file.exists(file.path(emd.path2, "emd.csv"))) {
      ComputeGeneEMD(data_S2, common_genes, dir.path = emd.path2)
    }
    
    # Wait for OT computation to finish
    emd_paths <- c(emd.path1, emd.path2)
    repeat {
      still_running <- !sapply(emd_paths, function(p) file.exists(file.path(p, "emd.csv")))
      if (!any(still_running)) break
      cat("  Waiting for OT computation...\n")
      Sys.sleep(10)
    }
    cat("  OT computation complete.\n")
    
    # Load and align gene distances
    gene_emd_ls <- list(
      LoadGeneEMD(file.path(emd.path1, "")),
      LoadGeneEMD(file.path(emd.path2, ""))
    )
    g <- Reduce(intersect, lapply(gene_emd_ls, rownames))
    gene_emd_ls <- lapply(gene_emd_ls, function(x) x[g, g])
    
    # Differential operator + eigendecomposition
    gene_graph_ls <- lapply(gene_emd_ls, ComputeGraphOperator)
    
    # Direction: condition1 - condition2
    #   For SmoM2: diff.op highlights SmoM2-specific co-localization
    #   For WlsKO: diff.op_CTL highlights CTL-specific co-localization
    if (ds$direc == "condition1") {
      diff.op <- ComputeDifferentialOperator(gene_graph_ls[[1]], gene_graph_ls[[2]])
    } else {
      diff.op <- ComputeDifferentialOperator(gene_graph_ls[[2]], gene_graph_ls[[1]])
    }
    E.list <- RunSVD(diff.op, eig_keep = nrow(diff.op))
    res_DiCoLo <- E.list$vectors
    
    # Update common_genes to genes present in DiCoLo output
    common_genes_dicolo <- rownames(res_DiCoLo)
    
    # ===================
    # MiloDE
    # ===================
    cat("Running miloDE...\n")
    sample_diff <- ifelse(ds$direc == "condition1", 1, 2)
    res_milode <- tryCatch({
      RunMiloDE(data_S1, data_S2, input_genes = common_genes_dicolo,
                query_id = sample_diff)
    }, error = function(e) {
      message(sprintf("miloDE failed: %s", conditionMessage(e)))
      NULL
    })
    
    # ===================
    # LEMUR
    # ===================
    cat("Running LEMUR...\n")
    res_lemur <- tryCatch({
      RunLEMUR(data_S1, data_S2, input_genes = common_genes_dicolo)
    }, error = function(e) {
      message(sprintf("LEMUR failed: %s", conditionMessage(e)))
      NULL
    })
    
    # ===================
    # DGCA
    # ===================
    cat("Running DGCA...\n")
    res_DGCA <- tryCatch({
      RunDGCA(data_S1, data_S2, input_genes = common_genes_dicolo)
    }, error = function(e) {
      message(sprintf("DGCA failed: %s", conditionMessage(e)))
      NULL
    })
    
    # memento
    # ===================
    cat("Running memento...\n")
    res_memento <- tryCatch({
      RunMemento(data_S1, data_S2, input_genes = common_genes_dicolo)
    }, error = function(e) {
      message(sprintf("DGCA failed: %s", conditionMessage(e)))
      NULL
    })
    
    # ===================
    # Save method results
    # ===================
    res.path <- file.path(dir.path.data, dataset_name, "parameters")
    if (!dir.exists(res.path)) dir.create(res.path, recursive = TRUE)
    
    saveRDS(list(
      res_DiCoLo = res_DiCoLo,
      res_milode = res_milode,
      res_lemur = res_lemur,
      res_DGCA = res_DGCA,
      res_memento = res_memento,
      common_genes = common_genes_dicolo,
      gt_genes = gt_genes
    ), file = file.path(res.path, "method_res.rds"))
    
  }
  
  # ===================
  # Compute AUPRC
  # ===================
  cat("Computing AUPRC...\n")
  
  # Direction for Generate_rank_table
  rank_direc <- ds$direc
  
  results <- lapply(method_ls, function(method) {
    res_name <- paste0("res_", method)
    de_res <- get(res_name)
    if (is.null(de_res)) return(NA_real_)
    
    tryCatch({
      df <- Generate_rank_table(de_res = de_res,
                                method = method,
                                input_genes = common_genes_dicolo,
                                direc = rank_direc,
                                n_eig = ds$n_eig)
      get_auc(real_score = setNames(-df$rank, df$gene),
              gt_gene_ls = gt_genes,
              metric = "auprc", plot = TRUE)
    }, error = function(e) {
      message(sprintf("  AUPRC failed for %s: %s", method, conditionMessage(e)))
      NA_real_
    })
  })
  names(results) = method_ls
  scores <- sapply(results, function(x) if(is.null(x)) NA else x$auc_score)
  curve_df <- do.call(rbind, lapply(method_ls, function(m) {
    if (is.null(results[[m]])) return(NULL)
    cbind(results[[m]]$curve, method = m)
  }))
  result_df <- data.frame(
    method = method_ls,
    auprc = scores,
    dataset = dataset_name,
    label = ds$label,
    n_gt_genes = length(gt_genes),
    n_common_genes = length(common_genes_dicolo),
    stringsAsFactors = FALSE
  )
  
  cat("\n--- Results ---\n")
  print(result_df)
  return(list(curve_df = curve_df,
         result_df = result_df))
})

# =============================================================================
# 3. Combine and save results
# =============================================================================
df_all <- do.call(rbind, lapply(all_results,function(x) x[["result_df"]]))
write.csv(df_all, file = file.path(dir.path.data,"smom2","benchmark_result.csv"))

# =============================================================================
# 4. Plot: Figure 4E
# =============================================================================
df_all$method <- factor(df_all$method, levels = c("DiCoLo", "DGCA", "lemur", "milode", "memento"))
levels(df_all$method) = c("DiCoLo","DGCA","LEMUR","miloDE","Memento")
p <- ggplot(df_all, aes(x = method, y = auprc, fill = method)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", auprc)), vjust = -0.5, size = 4) +
  # facet_wrap(~ label) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Method",
       y = "Normalized AUPRC") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 13),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )
ggsave(file.path(figure.path, "fig4E_signature_recovery.png"), p, 
       width = 5, height = 5)


baseline <- unique(all_results[[1]]$result_df$n_gt_genes) / 
  unique(all_results[[1]]$result_df$n_common_genes)

ggplot(all_results[[1]]$curve_df, aes(x = fpr, y = tpr, color = method)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = baseline, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.9, y = baseline + 0.01, label = "random baseline") +
  labs(x = "Recall", y = "Precision") +
  theme_classic()
ggplot(all_results[[1]]$cumul_df, aes(x = rank_frac, y = recall, color = method)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey50") +  # random baseline
  labs(x = "Fraction of genes (ranked)", y = "Fraction of DC signature recovered") +
  theme_classic()