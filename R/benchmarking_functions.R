# Benchmarking =======
# remotes::install_github("MarioniLab/miloDE") 
# remotes::install_github("const-ae/lemur")
# remotes::install_github("csglab/GEDI")
# devtools::install_github("andymckenzie/DGCA")
library(reticulate)
# Python interpreter is selected via config.R::use_dicolo_python();
# RunMemento() calls it lazily so that sourcing this file has no side effects.

# supervised embedding is more suitable for sensitive DE detection
add_azimuth_supervised = function(sce , genes , split.by = "sample", ref_samples , query_samples, nPC = 30, reducedDim.name , bpparam){
  require(Seurat)
  require(SeuratObject)
  require(SingleCellExperiment)
  sce_seurat <- CreateSeuratObject(counts = counts(sce[genes,]))
  sce_seurat@assays$RNA$data = logcounts(sce[genes,])
  sce_seurat = AddMetaData(sce_seurat, as.data.frame(colData(sce)), col.name = NULL)
  sce_seurat.list <- SplitObject(sce_seurat, split.by = split.by)
  # normalise and hvgs
  for (i in 1:length(sce_seurat.list)) {
    VariableFeatures(sce_seurat.list[[i]]) = genes
  }
  sce_reference.list = sce_seurat.list[ref_samples]
  
  # integrate reference
  sce_reference.merged <- merge(sce_reference.list[[1]], sce_reference.list[-1])
  sce_reference.merged <- sce_reference.merged %>% NormalizeData() %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = nPC, verbose = FALSE)
  sce_reference.integrated <- IntegrateLayers(
    object = sce_reference.merged,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    dims = 1:nPC,
    features = genes
  )
  
  # map query
  pca_proj_query = bplapply(query_samples , function(current.sample){
    current.sce = sce_seurat.list[[which(names(sce_seurat.list) == current.sample)]]
    query.anchors <- FindTransferAnchors(reference = sce_reference.integrated, query = current.sce, dims = 1:nPC,
                                            reference.reduction = "integrated.cca")
    current.sce = MapQuery(
      anchorset = query.anchors,
      reference = sce_reference.integrated,
      query = current.sce,
      reference.reduction = "integrated.cca"
    )
    out = Embeddings(current.sce[["ref.integrated.cca"]])
    colnames(out) = paste0("PC_" , c(1:nPC))
    return(out)
  } , BPPARAM = bpparam)
  
  # combine
  pca_proj_query = do.call(rbind , pca_proj_query)
  pca_ref = Embeddings(sce_reference.integrated[["integrated.cca"]])
  colnames(pca_ref) = paste0("PC_" , c(1:nPC))
  pca_joint = rbind(pca_proj_query , pca_ref)
  pca_joint = pca_joint[order(rownames(pca_joint)) , ]
  sce = sce[, order(colnames(sce))]
  reducedDim(sce , reducedDim.name) = pca_joint
  return(sce)
}

check_integration_neighbors <- function(labels, emb, k = 20, min_mix = 0.05) {
  library(FNN)
  nn <- get.knn(emb, k = k)$nn.index
  
  mix_scores <- sapply(1:nrow(emb), function(i) {
    mean(labels[nn[i, ]] != labels[i])
  })
  
  avg_mix <- mean(mix_scores)
  cat("Average neighbor mixing =", round(avg_mix, 3), "\n")
  
  if (avg_mix < min_mix) {
    stop("❌ Conditions are completely separated (no mixing).")
  } else {
    message("✅ Conditions show overlap (not completely separated).")
  }
}

library(miloDE)
RunMiloDE <- function(srat1, srat2, seed = 42, input_genes, query_id = 1, ncores = 6){
  require(BiocParallel)
  srat_ls = list(srat1,srat2)
  srat_ls = lapply(1:length(srat_ls),function(i){
    srat = srat_ls[[i]]
    srat$'condition' = paste0("condition",i)
    assays_to_remove <- setdiff(names(srat@assays), "RNA")
    for(assay in assays_to_remove) {
      srat[[assay]] <- NULL
    }
    srat
  })
  query_sample = paste0("condition",query_id)
  # Pseudo-replicate
  srat_ls = lapply(1:length(srat_ls),function(i){
    srat = srat_ls[[i]]
    set.seed(seed+i)
    srat$'replicate' = sample(c(rep(1, floor(ncol(srat)/2)), rep(2, ncol(srat) - floor(ncol(srat)/2))))
    srat$'replicate' = paste0(srat$'condition',"_rep",srat$'replicate')
    srat
  })
  
  srat = merge(srat_ls[[1]],srat_ls[-1])
  srat = srat %>%
    NormalizeData(verbose = FALSE)
  current.genes <- split(row.names(srat@meta.data), srat@meta.data$replicate) %>% lapply(function(cells_use) {
    srat[,cells_use] %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
      VariableFeatures()
  }) %>% unlist %>% unique
  # azimuth supervised
  sce <- as.SingleCellExperiment(srat)
  sce = scuttle::logNormCounts(sce)
  ref_samples = unique(sce$replicate[!sce$condition == query_sample])
  query_samples = unique(sce$replicate[sce$condition == query_sample])
  mcparam = MulticoreParam(workers = ncores)
  register(mcparam)
  sce = add_azimuth_supervised(sce, genes = current.genes, 
                               split.by = "replicate",
                               ref_samples = ref_samples , 
                               query_samples = query_samples , 
                               reducedDim.name = "azimuth_supervised" , 
                               bpparam = mcparam)
  # # harmony
  # srat <- srat %>%
  #   ScaleData(verbose = FALSE) %>%
  #   RunPCA(features = current.genes, npcs = 20, verbose = FALSE)
  # srat <- srat %>% harmony::RunHarmony("replicate")
  # harm = Embeddings(srat[["harmony"]])
  # SingleCellExperiment::reducedDim(sce, "harmony") <- harm
  check_integration_neighbors(labels = sce$condition, emb = reducedDim(sce, "azimuth_supervised"))
  
  sce_milo <- miloDE::assign_neighbourhoods(sce, 
                                            # k = 20, order = 2,
                                            # filtering = TRUE, 
                                            reducedDim_name = "azimuth_supervised",
                                            verbose = FALSE)
  system.time(
    de_milo <- miloDE::de_test_neighbourhoods(sce_milo[input_genes,], sample_id = "replicate", design = ~ condition, covariates = c("condition"))
  )
  # define direction
  tmp_row = de_milo[which.max(de_milo$logFC),]
  tmp_cells = sce_milo@nhoods[,tmp_row[,"Nhood_center"]]
  tmp_cells = names(tmp_cells)[tmp_cells==1]
  tmp_expr = FetchData(subset(srat, cells = tmp_cells),c(tmp_row[,"gene"],"condition"))
  colnames(tmp_expr)[1] = "value"
  tmp_expr = tmp_expr %>%
    group_by(condition) %>%
    summarise(mean_value = mean(value))
  
  pos_sample = ifelse(tmp_expr$mean_value[tmp_expr$condition == "condition1"] > tmp_expr$mean_value[tmp_expr$condition == "condition2"], "condition1","condition2")
  neg_sample = setdiff(tmp_expr$condition,pos_sample)
  de_milo$'direction' = ifelse(de_milo$logFC > 0, pos_sample, ifelse(de_milo$logFC < 0, neg_sample, NA))
  # de_milo$'direction' = ifelse(de_milo$logFC > 0, "condition1", ifelse(de_milo$logFC < 0, "condition2", NA))
  
  return(de_milo)
}

library(lemur)
RunLEMUR <- function(srat1,srat2,seed = 42, input_genes){
  srat_ls = list(srat1,srat2)
  srat_ls = lapply(1:length(srat_ls),function(i){
    srat = srat_ls[[i]]
    srat$'condition' = paste0("condition",i)
    assays_to_remove <- setdiff(names(srat@assays), "RNA")
    for(assay in assays_to_remove) {
      srat[[assay]] <- NULL
    }
    srat
  })
  # Pseudo-replicate
  srat_ls = lapply(1:length(srat_ls),function(i){
    srat = srat_ls[[i]]
    set.seed(seed+i)
    srat$'replicate' = sample(c(rep(1, floor(ncol(srat)/2)), rep(2, ncol(srat) - floor(ncol(srat)/2))))
    srat$'replicate' = paste0(srat$'condition',"_rep",srat$'replicate')
    srat
  })
  srat = merge(srat_ls[[1]],srat_ls[-1])
  sce <- as.SingleCellExperiment(srat)
  fit <- lemur(sce, design = ~ condition, n_embedding = 30)
  fit <- align_harmony(fit)
  fit <- test_de(fit, contrast = cond(condition = "condition1") - cond(condition = "condition2"))
  nei <- find_de_neighborhoods(fit[input_genes,], group_by = vars(replicate,condition))
  nei <- nei %>% mutate(direction = ifelse(lfc > 0, "condition1", ifelse(lfc < 0, "condition2", NA)))
  return(nei)
}


# Differential co-expression
library(DGCA)
RunDGCA <- function(srat1, srat2, input_genes){
  data_S = merge(srat1, srat2)
  data_S$'condition' = ifelse(colnames(data_S) %in% colnames(srat1),"condition1","condition2")
  data_mtx = GetAssayData(data_S, assay = "RNA", layer = "data")[input_genes,]
  design_mat = data_S@meta.data %>% select(condition)
  design_mat <- model.matrix(~ 0 + condition, data = design_mat)
  colnames(design_mat) = c("condition1","condition2")
  
  ddcor_res = ddcorAll(inputMat = data.frame(data_mtx), design = design_mat,
                       compare = c("condition1", "condition2"))
  
  ddcor_res = ddcor_res %>% mutate(direction = ifelse(zScoreDiff < 0, "condition1", "condition2"))
  return(ddcor_res)
}


RunMemento <- function(srat1, srat2, input_genes,
                       capture_rate = 0.25, num_boot = 1000L, num_cpus = 4L){
  # Bind the Python interpreter here rather than at load time (see config.R).
  if (exists("use_dicolo_python")) use_dicolo_python()

  data_S = merge(srat1, srat2)
  data_S$'condition' = ifelse(colnames(data_S) %in% colnames(srat1), "condition1", "condition2")

  counts = as.matrix(GetAssayData(data_S, assay = "RNA", layer = "counts"))
  meta   = data_S@meta.data["condition"]
  meta$'stim' = as.integer(ifelse(meta$condition == "condition1", 0L, 1L))  # cond1->0, cond2->1
  
  anndata   <- import("anndata")
  sparse    <- import("scipy.sparse")
  itertools <- import("itertools")
  builtins  <- import_builtins()
  np        <- import("numpy")
  memento_pkg   <- import("memento")
  
  adata <- anndata$AnnData(
    X   = sparse$csr_matrix(np$array(t(counts))),
    obs = meta,
    var = data.frame(gene = rownames(counts), row.names = rownames(counts))
  )
  adata$obs['capture_rate'] = capture_rate
  
  # expand binary_test_2d
  memento_pkg$setup_memento(adata, q_column = 'capture_rate')
  memento_pkg$create_groups(adata, label_columns = list('stim'))
  memento_pkg$compute_1d_moments(adata, min_perc_group = 0.9)
  
  surviving_genes <- py_to_r(adata$uns[['memento']][['gene_list']])
  test_genes <- intersect(input_genes, surviving_genes)
  gene_pairs <- builtins$list(itertools$combinations(test_genes, 2L))
  memento_pkg$compute_2d_moments(adata, gene_pairs)
  
  sample_meta <- memento_pkg$get_groups(adata)[c('stim')]
  memento_pkg$ht_2d_moments(adata, treatment = sample_meta,
                        num_boot = as.integer(num_boot),
                        num_cpus = as.integer(num_cpus), verbose = 1L)
  
  ht  <- py_to_r(memento_pkg$get_2d_ht_result(adata))
  mom <- py_to_r(memento_pkg$get_2d_moments(adata)[[1]])
  
  c1_col <- "sg^0"; c2_col <- "sg^1"
  stopifnot(all(c(c1_col, c2_col) %in% colnames(mom)))
  
  res <- dplyr::inner_join(ht, mom[, c("gene_1","gene_2",c1_col,c2_col)],
                           by = c("gene_1","gene_2"))
  res$'corr_c1' = res[[c1_col]]
  res$'corr_c2' = res[[c2_col]]
  res$'delta'   = res$corr_c2 - res$corr_c1
  res$'padj'    = p.adjust(res$corr_pval, method = "BH")
  res$'direction' = ifelse(res$delta > 0, "condition2", "condition1")
  res$'winner_pos' = ifelse(res$delta > 0, res$corr_c2 > 0, res$corr_c1 > 0)
  attr(res, "n_genes_tested") = length(union(res$gene_1, res$gene_2))
  attr(res, "n_genes_input")  = length(input_genes)
  return(res)
  
  # res <- memento_pkg$binary_test_2d(
  #   adata         = adata,
  #   gene_pairs    = gene_pairs,
  #   capture_rate  = capture_rate,
  #   treatment_col = "stim",
  #   num_cpus      = as.integer(num_cpus),
  #   num_boot      = as.integer(num_boot)
  # )
  # res <- res$reset_index()
  # res <- py_to_r(res)
  # res$'padj'      = p.adjust(res$corr_pval, method = "BH")
  # res$'direction' = ifelse(res$corr_coef > 0, "condition2", "condition1")
  # return(res)
}


Generate_rank_table <- function(de_res, method, input_genes, direc = "condition1", n_eig = 1){
  if(method == "milode"){
    de_res = de_res %>% mutate(padj = pval_corrected_across_genes) %>% 
      filter(direction == direc) %>% 
      arrange(.,padj,desc(abs(logFC))) %>% distinct(.,gene, .keep_all = TRUE)
    de_res$'rank' = 1:nrow(de_res)
    de_res$'score' = -log10(de_res$padj)
  }
  if(method == "lemur"){
    de_res = de_res %>% mutate(padj = adj_pval) %>% 
      filter(direction == direc) %>%
      arrange(.,padj,desc(lfc)) %>% distinct(.,name, .keep_all = TRUE)
    de_res$'rank' = 1:nrow(de_res)
    de_res$'score' = -log10(de_res$padj)
    de_res$'gene' = de_res$name
  }
  
  if(method == "DGCA"){
    de_res = de_res %>% filter(direction == direc) %>% filter(Classes != "NonSig") %>%
      mutate(condition1 = unlist(lapply(strsplit(as.character(Classes),"/"),function(x) x[[1]])), 
             condition2 = unlist(lapply(strsplit(as.character(Classes),"/"),function(x) x[[2]]))) %>%
      mutate(padj = pValDiff_adj)
    de_res = de_res[de_res[,direc] == "+",]
    de_res = de_res %>%
      tidyr::pivot_longer(cols = c(Gene1, Gene2), values_to = "gene") %>% select(gene,padj,zScoreDiff)
    de_res = de_res %>% arrange(.,padj,desc(abs(zScoreDiff))) %>% distinct(.,gene, .keep_all = TRUE)
    
    de_res$'rank' = 1:nrow(de_res)
    de_res$'score' = -log10(de_res$padj)
  }
  
  if(method == "DiCoLo"){ 
    if(n_eig == 1){
      de_res = data.frame(gene = rownames(de_res), score = abs(de_res[,1]))
    } else {
      de_res = data.frame(gene = rownames(de_res),
                          score = apply(abs(de_res[, 1:n_eig, drop = FALSE]), 1, max))
    }
    de_res = de_res %>% arrange(.,desc(score)) %>% 
      distinct(.,gene, .keep_all = TRUE)
    de_res$'rank' = 1:nrow(de_res)
  }
  
  if(method == "memento"){
    de_res = de_res %>% filter(direction == direc) %>%
      tidyr::pivot_longer(cols = c(gene_1, gene_2), values_to = "gene") %>%
      select(gene, corr_pval, corr_coef)
    if(nrow(de_res) == 0){
      de_res = data.frame(gene = character(0), corr_pval = numeric(0), corr_coef = numeric(0),
                          rank = nrow(de_res), score = numeric(0))
    }else{
      de_res = de_res %>% arrange(., corr_pval, desc(abs(corr_coef))) %>%
        distinct(., gene, .keep_all = TRUE)
      de_res$'rank'  = 1:nrow(de_res)
      de_res$'score' = -log10(de_res$corr_pval)
    }
  }
  
  rank_df = data.frame(gene = input_genes)
  rank_df <- dplyr::left_join(rank_df,de_res,by = "gene") %>% select(gene, rank, score)
  rank_df[is.na(rank_df$rank),"rank"] = nrow(rank_df)
  rank_df[is.na(rank_df$score),"score"] = 0
  
  return(rank_df)
}

get_auc = function(gene_list = NULL, real_score, 
                   gt_gene_ls, metric = "auprc", plot = FALSE){
  if(!is.null(names(real_score))){
    df = data.frame(real_score)
  }else{
    df = data.frame(real_score, row.names = gene_list)
  }
  df$'marker_binary' = 0
  gt_gene_ls = gt_gene_ls[gt_gene_ls %in% rownames(df)]
  df[gt_gene_ls,"marker_binary"] = 1
  if(metric == "auprc"){
    pr <- PRROC::pr.curve(scores.class0 = df$real_score[df$marker_binary == 1], 
                          scores.class1 = df$real_score[df$marker_binary == 0], 
                          curve = plot)
    auprc = pr$auc.integral
    
    # balance auprc
    baseline = sum(df$marker_binary == 1) / nrow(df)
    auprc_balanced = (auprc - baseline) / (1 - baseline)
    if(plot){
      auc_coords <- as.data.frame(pr$'curve')
      colnames(auc_coords) <- c("recall","precision","threshold")
      auc_coords <- auc_coords %>% arrange(recall)
      return(list(auc_score = auprc_balanced, 
                  curve = auc_coords))
    }else{
      return(auprc_balanced)
    }

  }else if(metric == "auroc"){
    cat("Rank in decreasing order\n")
    df[order(df$real_score,decreasing = TRUE),'rank'] = 1:nrow(df)
    roc_obj = pROC::roc(response = df$marker_binary, 
                        predictor = df$rank, direction = ">",
                        levels = c(0, 1),quiet = TRUE)
    auroc = as.numeric(pROC::auc(roc_obj))
    
    if(plot){
      roc_coords <- data.frame(
        fpr = 1 - roc_obj$specificities,
        tpr = roc_obj$sensitivities
      )
      return(list(auc_score = auroc, 
                  curve = roc_coords))
    }else{
      return(auroc)
    }
  }
}


# ---------------------------------------------------------------------------
# plot_module_enrichment
#
# GSEA-style running-enrichment curve for one method against one module.
# Used by scripts/run_concordance_analysis.R (Supplemental Tables 3-4).
# ---------------------------------------------------------------------------
plot_module_enrichment <- function(ranks, module_genes, method_name, 
                                   module_name = "", color = "#2196F3") {
  # Single method, single module — GSEA-style running enrichment
  rank_df <- ranks[ranks$method == method_name, ]
  rank_df <- rank_df[order(rank_df$rank), ]
  n <- nrow(rank_df)
  n_hit <- sum(rank_df$gene %in% module_genes)
  
  # Running sum
  rank_df$hit <- ifelse(rank_df$gene %in% module_genes, 1, 0)
  rank_df$running_sum <- cumsum(rank_df$hit / n_hit - (1 - rank_df$hit) / (n - n_hit))
  
  ggplot(rank_df, aes(x = 1:n, y = running_sum)) +
    geom_line(color = color, linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    # Barcode ticks at bottom
    geom_segment(data = rank_df[rank_df$hit == 1, ],
                 aes(x = which(rank_df$hit == 1), xend = which(rank_df$hit == 1),
                     y = -0.05, yend = -0.15),
                 color = "black", linewidth = 0.3) +
    labs(x = "Gene rank", y = "Enrichment score",
         title = paste0(method_name, " — ", module_name)) +
    theme_minimal()
}
