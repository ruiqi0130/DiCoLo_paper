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
                                    diffuse = FALSE,
                                    overlap_frac = NULL){
  
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
      # Generate shared latent variable per cell (Gamma-distributed)
      # This induces positive correlation across genes within same cell
      # cov_strength controls how much variance comes from shared vs independent
      # Higher shape = less variance in latent = less correlation
      shape_param = (1 - cov_strength) / cov_strength  # maps [0,1] -> [Inf, 0]
      
      repeat {
        # Shared latent: z_j for each cell in neighborhood
        latent = rgamma(n_incore, shape = shape_param, rate = shape_param)  # mean = 1
        
        # For each gene, sample NB with cell-specific mean = mu_signal * z_j
        signal_mat = matrix(0, nrow = n_genes, ncol = n_incore)
        for (i in seq_along(selected_genes)) {
          if (diffuse) {
            cell_means = mu_incore * latent
          } else {
            cell_means = gene_params$'mean_signal' * latent
          }
          signal_mat[i, ] = rnbinom(n_incore, 
                                    mu = cell_means, 
                                    size = gene_params$'dispersion_signal')
        }
        
        # Background and dropout (same as original, per gene)
        excore_cells = setdiff(colnames(count_mtx), incore_cells)
        p_base <- max(min(1 - diff.pct, 1), 0)
        max_excore <- n_incore
        p_expr <- if (length(excore_cells) > 0) min(p_base, max_excore / length(excore_cells)) else 0
        
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
          # max_excore <- length(incore_cells)
          # p_expr <- if (length(excore_cells)>0) min(p_base, max_excore / length(excore_cells)) else 0
          p_expr <- p_base
          
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
