library(MatrixGenerics)
library(rdist)
make_seed_grid <- function(n_reps, n_sample, base_seed = 233){
  expand.grid(rep = seq_len(n_reps), sample = seq_len(n_sample)) |>
    transform(seed = base_seed + seq_len(n_reps * n_sample))
}

Generate_gene_parameters <- function(srat1, srat2, 
                                assay = "RNA", layer = "counts",
                                n_neighbors,
                                ngene_per_neighbors_ls, seed, mean_quantile = 0.95){
  set.seed(seed)
  # Candidate gene list
  count_mtx = GetAssayData(srat1, assay = assay, layer = layer)
  expr_percent = apply(count_mtx > 0, 1, sum)/ncol(srat1)
  expr_mean = apply(count_mtx,1,mean)
  expr_var = apply(count_mtx,1,var)
  candidate_genes = names(expr_percent)[which(expr_percent > 0.005 & 
                                                      expr_percent < 0.5)]
  candidate_genes = intersect(candidate_genes,names(expr_mean)[expr_mean < expr_var])
  
  expr_compare = apply(GetAssayData(srat2, assay = assay, layer = layer) > 0, 1, sum)
  candidate_genes = candidate_genes[candidate_genes %in% names(expr_compare)]
  candidate_genes = candidate_genes[expr_compare[candidate_genes]>0]
  
  # Pick genes for each neighborhood
  genes_picked <- sample(candidate_genes, sum(ngene_per_neighbors_ls), replace = FALSE)
  selected_gene_groups <- split(genes_picked, rep(seq_along(ngene_per_neighbors_ls), ngene_per_neighbors_ls))
  
  # Signal & noise parameters
  gene_raw_mean = MatrixGenerics::rowMeans(count_mtx)
  gene_raw_disp = (MatrixGenerics::rowVars(count_mtx) - gene_raw_mean)/(gene_raw_mean^2)
  gene_raw_disp = gene_raw_disp[gene_raw_disp>0]; gene_raw_disp = gene_raw_disp[!is.na(gene_raw_disp)]
  
  mean_signal = quantile(gene_raw_mean,mean_quantile)
  dispersion_signal = quantile(gene_raw_disp,0.5)
  mean_noise = quantile(gene_raw_mean,0.5)
  dispersion_noise = quantile(gene_raw_disp,0.5)

  return(list(selected_gene_groups = selected_gene_groups,
              mean_signal = mean_signal,
              dispersion_signal = dispersion_signal,
              mean_noise = mean_noise,
              dispersion_noise = dispersion_noise))
}

Select_Core_cells <- function(srat,reduc = "umap", npc = 2, min_knn = 50, n_neighbors, seed){
  set.seed(seed)
  cell_coords = Embeddings(srat[[reduc]])[,1:npc]
  cell_names = rownames(cell_coords)
  
  # Find candidates with FPS (oversampling)
  num_candidate <- n_neighbors * 10
  dist_matrix <- rdist::rdist(cell_coords)
  fps_order <- rdist::farthest_point_sampling(dist_matrix, k = num_candidate)
  
  candidate_indices <- fps_order
  candidate_names <- cell_names[candidate_indices]
  
  p1 = DimPlot(srat,cells.highlight = candidate_names)
  knn_buffer_list <- FNN::get.knn(cell_coords,
                           k = min_knn, algorithm = "kd_tree")
  exclusion_map <- lapply(1:nrow(knn_buffer_list$nn.index), function(i) {
    neighbor_indices <- knn_buffer_list$nn.index[i, ]
    return(neighbor_indices)
  })
  
  # Greedy search
  final_core_indices <- c()
  is_excluded <- rep(FALSE, nrow(cell_coords))
  names(is_excluded) <- cell_names
  
  # Loop through our FPS-prioritized candidates
  for (candidate_index in candidate_indices) {
    candidate_name <- cell_names[candidate_index]
    # 1. Check if this candidate is already "occupied"
    if (is_excluded[candidate_name]) {
      next
    }
    
    # 2. If it's not occupied, this is our new core!
    final_core_indices <- c(final_core_indices, candidate_index)
    
    # 3. Find its *exclusion buffer* neighbors
    buffer_neighbor_names <- cell_names[exclusion_map[[candidate_name]]]
    
    # 4. Add the core *and* its buffer neighbors to the excluded set
    cells_to_exclude <- c(candidate_name, buffer_neighbor_names)
    is_excluded[cells_to_exclude] <- TRUE
    
    # 5. Check if we have enough cores
    if (length(final_core_indices) >= n_neighbors) {
      break # We're done
    }
  }
  final_core_names <- cell_names[final_core_indices]
  p2 = DimPlot(srat,cells.highlight = final_core_names)
  return(final_core_names)
  # return(p1+p2)
}

Inject_local_signals_v2 <- function(srat, condition_vec = NULL, assay = "RNA",
                                    layer = "counts",
                                    reduc = "pca", npc = 10, min_knn = 50,
                                    ncell_per_neighbors_ls,
                                    zinb_prob = 0.2, diff.pct = 0.8, 
                                    gene_params, seed,
                                    cov_strength = 0,
                                    cov_K = 3,
                                    cov_loading_sd = 1,
                                    diffuse = FALSE,
                                    overlap_frac = NULL, cap_background = TRUE,
                                    batch_downsample = 1){
  
  merged_sample = FALSE
  set.seed(seed)
  srat$'signal_cells' = NA
  srat[[assay]]@meta.features$'signal_genes' = NA
  n_neighbors = length(ncell_per_neighbors_ls)
  if(!is.null(condition_vec) & !is.null(names(ncell_per_neighbors_ls))){
    merged_sample = TRUE
    # print("merged_sample!")
    condition_neighbors = do.call(rbind,strsplit(names(ncell_per_neighbors_ls),"_"))[,1]
  }
  
  count_mtx = as.matrix(GetAssayData(srat, assay = assay, layer = layer))
  core_cell_ls = Select_Core_cells(srat,reduc = "umap", npc=2, min_knn, n_neighbors, seed)
  
  for(id in seq_len(n_neighbors)){
    core_cell = core_cell_ls[id]
    # Define Cell neighborhood
    lowdim_embedding <- Embeddings(srat[[reduc]])[,1:npc]
    core_cell_coords <- lowdim_embedding[core_cell, ]
    distances <- apply(lowdim_embedding, 1, function(x) sqrt(sum((x - core_cell_coords)^2)))
    if(merged_sample){
      if(condition_neighbors[id] == "c"){distances = distances}
      if(condition_neighbors[id] == "1"){distances = distances[condition_vec == levels(condition_vec)[1]]}
      if(condition_neighbors[id] == "2"){distances = distances[condition_vec == levels(condition_vec)[2]]}
    }
    incore_cells <- names(sort(distances)[1:ncell_per_neighbors_ls[id]]) # Include the core cell itself
    
    # --- Overlap logic: store neighborhood 1 and override neighborhood 2 ---
    if (!is.null(overlap_frac) && id == 1) {
      prev_incore_cells <- incore_cells
      prev_core_cell <- core_cell
    }
    
    if (!is.null(overlap_frac) && id == 2) {
      # Pick core 2 at edge of neighborhood 1 (farthest from core 1)
      prev_core_coords <- lowdim_embedding[prev_core_cell, ]
      dists_from_core1 <- apply(lowdim_embedding[prev_incore_cells, , drop = FALSE], 1, 
                                 function(x) sqrt(sum((x - prev_core_coords)^2)))
      core_cell <- names(sort(dists_from_core1, decreasing = TRUE))[1]
      
      # Recompute distances from core 2
      core_cell_coords <- lowdim_embedding[core_cell, ]
      distances <- apply(lowdim_embedding, 1, function(x) sqrt(sum((x - core_cell_coords)^2)))
      
      # Strict overlap control
      n_shared <- floor(ncell_per_neighbors_ls[id] * overlap_frac)
      n_unique <- ncell_per_neighbors_ls[id] - n_shared
      
      # Shared: n_shared cells from neighborhood 1, closest to core 2
      in_prev_sorted <- names(sort(distances[prev_incore_cells]))
      shared_cells <- in_prev_sorted[1:n_shared]
      
      # Unique: n_unique cells NOT in neighborhood 1, closest to core 2
      outside_prev <- setdiff(colnames(count_mtx), prev_incore_cells)
      outside_sorted <- names(sort(distances[outside_prev]))
      unique_cells <- outside_sorted[1:n_unique]
      
      incore_cells <- c(shared_cells, unique_cells)
    }
    
    # Perturb gene expression
    selected_genes = gene_params$'selected_gene_groups'[[id]]
    n_genes = length(selected_genes)
    n_incore = length(incore_cells)
    
    # --- Diffuse logic: compute cell-specific mu for incore cells ---
    if (diffuse) {
      incore_dists <- distances[incore_cells]
      max_dist <- max(incore_dists)
      min_dist <- min(incore_dists)
      if (max_dist > min_dist) {
        alpha <- 1 - (incore_dists - min_dist) / (max_dist - min_dist)
      } else {
        alpha <- rep(1, n_incore)
      }
      mu_incore <- alpha * gene_params$'mean_signal' + (1 - alpha) * gene_params$'mean_noise'
    }
    
    if(cov_strength > 0){
      # -----------------------------------------------------------------
      # K-factor latent covariance model (replaces the single shared
      # scalar z_j). Each cell draws cov_K shared latent factors; each
      # gene has its own signed loadings on those factors. The log-scale
      # multiplier for gene g in cell j is  exp( w_g . f_j ), so the
      # covariance between two genes is sigma^2 * (w_a . w_b):
      #   - it VARIES across gene pairs (structured), because loadings
      #     differ per gene, and
      #   - it can be POSITIVE OR NEGATIVE, because loadings are mean-0.
      # cov_strength keeps its meaning (overall covariance strength) via
      # sigma(cov_strength); cov_strength = 0 routes to the independent
      # branch below, so "0 = independent" is preserved.
      #
      # sigma is chosen so the per-gene multiplier variance matches the
      # old Gamma model's CV^2 = cov_strength/(1-cov_strength) on average,
      # keeping the cov_strength axis comparable to the original S7-D.
      # Each gene's multiplier is mean-1 centered so the SIGNAL LEVEL is
      # unchanged and only the covariance structure is added.
      # -----------------------------------------------------------------
      sigma = sqrt(-log(1 - cov_strength) / cov_K)  # 0 -> 0 (independent), monotone in cov_strength
      
      # Loadings: fixed per neighborhood, signed (mean 0) => +/- structured covariance
      W = matrix(rnorm(n_genes * cov_K, mean = 0, sd = cov_loading_sd),
                 nrow = n_genes, ncol = cov_K)
      logvar_g = (sigma^2) * rowSums(W^2)           # per-gene log-variance, for mean-1 centering
      
      rownames(W) <- selected_genes
      srat@misc[[paste0("cov_loadings_", names(gene_params$selected_gene_groups)[id])]] <- list(W = W, sigma = sigma)
      
      repeat {
        # Shared latent factors per cell: f_j ~ N(0, sigma^2)
        Fmat = matrix(rnorm(n_incore * cov_K, mean = 0, sd = sigma),
                      nrow = n_incore, ncol = cov_K)
        U = W %*% t(Fmat)                           # n_genes x n_incore log-scale effect
        M = exp(U - 0.5 * logvar_g)                 # multiplier, E_f[M] = 1 per gene
        
        # For each gene, sample NB with cell-specific mean = base * M[g, ]
        signal_mat = matrix(0, nrow = n_genes, ncol = n_incore)
        for (i in seq_along(selected_genes)) {
          base_i = if (diffuse) mu_incore else gene_params$'mean_signal'
          cell_means = base_i * M[i, ]
          signal_mat[i, ] = rnbinom(n_incore, 
                                    mu = cell_means, 
                                    size = gene_params$'dispersion_signal')
        }
        
        # Background and dropout (same as original, per gene)
        excore_cells = setdiff(colnames(count_mtx), incore_cells)
        p_base <- max(min(1 - diff.pct, 1), 0)
        if(cap_background){
          max_excore <- n_incore
          p_expr <- if (length(excore_cells) > 0) min(p_base, max_excore / length(excore_cells)) else 0
        }else{
          p_expr <- p_base
        }
                
        all_ok = TRUE
        for (i in seq_along(selected_genes)) {
          gene = selected_genes[i]
          count_mtx[gene, ] = 0
          
          background <- sapply(excore_cells, function(x) {
            if (runif(1) < p_expr) {
              rnbinom(1, mu = gene_params$'mean_noise', size = gene_params$'dispersion_noise')
            } else { 0 }
          })
          
          zero_mask <- rbinom(ncol(count_mtx), size = 1, prob = 1 - zinb_prob)
          expr_vec <- c(signal_mat[i, ], background)[match(colnames(count_mtx),
                                                           c(incore_cells, excore_cells))] * zero_mask
          names(expr_vec) = colnames(count_mtx)
          
          if (!any(expr_vec[incore_cells] != 0)) {
            all_ok = FALSE
            break
          }
          count_mtx[gene, ] = expr_vec
        }
        if (all_ok) break
      }
      
    } else {
      for (i in seq_along(selected_genes)) {
        gene = selected_genes[i]
        count_mtx_changed = count_mtx[,names(distances)]
        count_mtx_changed[gene,] = 0 # initialize
        
        repeat{
          ## Signal for incore cells
          if (diffuse) {
            signal_vec = sapply(seq_along(incore_cells), function(j) {
              rnbinom(1, mu = mu_incore[j], size = gene_params$'dispersion_signal')
            })
          } else {
            signal_vec = rnbinom(length(incore_cells), 
                                 mu = gene_params$'mean_signal', 
                                 size = gene_params$'dispersion_signal')
          }
          
          ## Background noise
          # Cap the excore expression probability by the pct-diff, 
          # with a hard limit ensuring excore expressers never exceed 
          # the incore size.
          excore_cells = setdiff(colnames(count_mtx_changed), incore_cells)
          p_base <- max(min(1 - diff.pct, 1), 0)
          if(cap_background){
            max_excore <- n_incore
            p_expr <- if (length(excore_cells) > 0) min(p_base, max_excore / length(excore_cells)) else 0
          }else{
            p_expr <- p_base
          }
        
          background <- sapply(excore_cells, function(x) {
            if (runif(1) < p_expr) {
              rnbinom(1, mu = gene_params$'mean_noise', 
                      size = gene_params$'dispersion_noise')
            } else {
              0  # No expression
            }
          })
          
          # Dropout
          zero_mask <- rbinom(ncol(count_mtx_changed), size = 1, prob = 1 - zinb_prob)
          expr_vec <- (c(signal_vec,background))[match(colnames(count_mtx_changed),
                                                       c(incore_cells,excore_cells))] * zero_mask
          names(expr_vec) = colnames(count_mtx_changed)
          if(any(expr_vec[incore_cells] !=0)){
            count_mtx[gene,colnames(count_mtx_changed)] = expr_vec
            break
          }
        }
      }
    }
    
    
    
    srat$'signal_cells'[incore_cells] = names(gene_params$'selected_gene_groups')[id]
    srat[[assay]]@meta.features[selected_genes,'signal_genes'] = names(gene_params$'selected_gene_groups')[id]
  
  }
  
  # batch dropout for signal cells
  if(batch_downsample < 1){
    cat("downsample UMI by: ",batch_downsample)
    signal_cells <- names(which(!is.na(srat$signal_cells)))
    count_mtx[, signal_cells] <- as.matrix(scuttle::downsampleMatrix(
      count_mtx[, signal_cells], prop = batch_downsample, bycol = TRUE))
  }
  
  # Update srat obj
  srat[[assay]] <- SetAssayData(srat[[assay]], layer = layer, new.data = as(count_mtx, "dgCMatrix"))
  srat <- DietSeurat(
    srat,
    layers = "counts",    
    assays = "RNA",    
    dimreducs = NULL, 
    graphs = NULL
  )
  srat[["RNA"]]$scale.data = NULL; srat[["RNA"]]$data = NULL
  return(srat)
}


# ---------------------------------------------------------------------------
# Summarize_injected_covariance
#
# Diagnostic for the K-factor covariance model. For each injected module it
# computes the empirical gene-gene correlation among the injected genes,
# measured across that module's injected cells, and reports how much of the
# covariance is negative and its range. This is the evidence for R3-Q1: it
# shows the design produces STRUCTURED (pair-varying) and SIGNED (+/-)
# covariance, unlike the old single-scalar model (which was uniform and
# strictly positive).
#
# Returns, per module and pooled: the correlation matrix, the fraction of
# off-diagonal gene pairs that are negatively correlated, and the min / mean /
# max correlation. Feed $cor_matrices into a heatmap (e.g. pheatmap) or
# $pooled$offdiag into a histogram for the figure / point-to-point reply.
# ---------------------------------------------------------------------------
Summarize_injected_covariance <- function(srat) {
  cov_info <- srat@misc[grep("^cov_loadings_", names(srat@misc))]
  
  if (length(cov_info) == 0) {
    # cov_strength = 0, no loadings stored => all zeros
    return(list(cov_matrices = list(), per_module = list(),
                pooled = list(offdiag = 0, frac_negative = 0,
                              min = 0, mean = 0, max = 0)))
  }
  
  offdiag_of <- function(m) m[upper.tri(m)]
  cor_matrices <- list()
  per_module   <- list()
  
  for (nm in names(cov_info)) {
    mod   <- sub("^cov_loadings_", "", nm)
    W     <- cov_info[[nm]]$W
    sigma <- cov_info[[nm]]$sigma
    
    # Theoretical pairwise covariance of log-scale multipliers
    cov_mat <- (sigma^2) * (W %*% t(W))
    rownames(cov_mat) <- colnames(cov_mat) <- rownames(W)
    
    od <- offdiag_of(cov_mat)
    cor_matrices[[mod]] <- cov_mat
    per_module[[mod]] <- list(
      n_genes       = nrow(cov_mat),
      frac_negative = mean(od < 0),
      min = min(od), mean = mean(od), max = max(od)
    )
  }
  
  pooled_od <- unlist(lapply(cor_matrices, offdiag_of))
  list(cor_matrices = cor_matrices, per_module = per_module,
       pooled = list(offdiag = pooled_od,
                     frac_negative = mean(pooled_od < 0),
                     min = min(pooled_od), mean = mean(pooled_od),
                     max = max(pooled_od)))
}

AddRandomNoise <- function(srat, n_shuffle_genes = 50, candidate_genes,
                           assay = "RNA", layer = "count", seed = 233){
  # Add noise signals by randomly select some genes and shuffle
  set.seed(seed)
  selected_genes <- sample(candidate_genes, n_shuffle_genes)
  
  count_mtx = as.matrix(GetAssayData(srat, assay = assay, layer = layer))
  set.seed(seed + 1)
  count_mtx[selected_genes,] = t(apply(count_mtx[selected_genes,],1,sample))
  Counts(srat[[assay]]) = count_mtx
  srat[[assay]]@meta.features[selected_genes,'signal_genes'] = "noise"
  
  return(srat)
}
