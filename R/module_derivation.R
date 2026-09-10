# =============================================================================
# module_derivation.R
#
# One function implementing the DiCoLo module-derivation workflow for a single
# direction of the contrast (see the package vignette,
# https://klugerlab.github.io/DiCoLo/articles/DiCoLo_demo.html):
#
#   Step 2  ComputeGraphOperator -> ComputeDifferentialOperator
#   Step 3  RunSVD -> FindKneePoint -> SelectSignificantGenes
#   Down    ClusterGenes -> RankGeneModules
#
# It returns the modules plus the four derivation panels (A-D) that make up
# Supplemental Fig S2 / S3 (SmoM2) and Fig S6 (Wls). Both direction-specific
# scripts call it, so the two contrasts are derived by identical code.
# =============================================================================

# -----------------------------------------------------------------------------
# DeriveDiCoLoModules
#
#   emd_list   named list of two gene-gene OT distance matrices, already
#              restricted to a common gene set and in the same gene order.
#              Names are the short condition labels, e.g. c("MUT", "CTL").
#   tsne_list  named list of two gene t-SNE embeddings, same names/order.
#   comp       name of the condition the operator is computed FOR (the
#              "compared" condition). The other one is the projection.
#   plot.fig   whether FindKneePoint should draw the spectrum (panel A).
#
# Value: list with
#   modules     factor of module assignments, names = genes, levels prefixed
#               with `comp` (e.g. MUT1, MUT2, ...)
#   eigen       the RunSVD output
#   n_eigvec    number of leading eigenvectors kept at the knee point
#   indicators  gene x eigenvector 0/1 significance matrix
#   panels      list(A, B, C, D) of the derivation figure panels
# -----------------------------------------------------------------------------
DeriveDiCoLoModules <- function(emd_list, tsne_list, comp,
                                min_gene = 5, deepSplit = 0, lfdr_thresh = 0.2) {

  stopifnot(length(emd_list) == 2, length(tsne_list) == 2)
  stopifnot(identical(names(emd_list), names(tsne_list)))
  stopifnot(comp %in% names(emd_list))

  labels <- names(emd_list)
  proj   <- setdiff(labels, comp)

  # --- Step 2: graph operators and the differential operator -----------------
  graph_list <- lapply(emd_list, ComputeGraphOperator)
  diff.op    <- ComputeDifferentialOperator(graph_list[[comp]],
                                            graph_list[[proj]])

  # --- Step 3: spectral decomposition ----------------------------------------
  eigen_list <- RunSVD(diff.op, eig_keep = nrow(diff.op))

  # Panel A: eigenvalue spectrum, top eigenvalues selected at the knee point.
  # FindKneePoint draws the spectrum itself, so it is captured as a base plot.
  panelA <- function() {
    FindKneePoint(
      eigen_list$values[1:floor(sqrt(length(eigen_list$values)))],
      plot.fig = TRUE)
  }
  n_eigvec <- FindKneePoint(
    eigen_list$values[1:floor(sqrt(length(eigen_list$values)))],
    plot.fig = FALSE) - 1

  # Panel B: gene embeddings of both conditions, coloured by the loadings of
  # each retained eigenvector. Row 1 = comp, row 2 = proj.
  panelB <- wrap_plots(lapply(seq_len(n_eigvec), function(i) {
    gene_partition <- setNames(eigen_list$vectors[, i],
                               rownames(eigen_list$vectors))
    p_comp <- VisualizeGeneTSNE(gene_embedding = tsne_list[[comp]],
                                gene_partition = gene_partition) +
      ggtitle(sprintf("Eigenvector %d", i)) + labs(color = "loadings")
    p_proj <- VisualizeGeneTSNE(gene_embedding = tsne_list[[proj]],
                                gene_partition = gene_partition) +
      ggtitle(NULL) + labs(color = "loadings")
    p_comp / p_proj
  }), nrow = 1)

  # --- Significant genes per eigenvector (localFDR) --------------------------
  gene_loadings <- eigen_list$vectors[, seq_len(n_eigvec), drop = FALSE]
  gene_indictors <- SelectSignificantGenes(gene_loadings,lfdr_thresh = lfdr_thresh)

  # Panel C: candidate genes in red on the comp gene embedding.
  panelC <- wrap_plots(lapply(seq_len(ncol(gene_indictors)), function(i) {
    VisualizeGeneTSNE(gene_embedding = tsne_list[[comp]],
                      gene_partition = as.factor(gene_indictors[, i]),
                      text = FALSE,
                      module_color = c("1" = "red", "0" = "lightgrey")) +
      labs(color = "Signif. genes") + ggtitle(sprintf("Eigenvector %d", i))
  }), nrow = 1) + plot_layout(guides = "collect")

  # --- Cluster the candidate genes into modules ------------------------------
  signf_genes <- rownames(gene_indictors)[rowSums(gene_indictors) > 0]
  gene_modules <- ClusterGenes(
    as.dist(emd_list[[comp]][signf_genes, signf_genes]),
    min_gene = min_gene, deepSplit = deepSplit)
  gene_modules <- RankGeneModules(gene_modules, gene_indictors,
                                  eigen_list$values)
  levels(gene_modules) <- paste0(comp, levels(gene_modules))

  # Panel D: the modules overlaid on both gene embeddings (comp left).
  panelD <- wrap_plots(lapply(c(comp, proj), function(lab) {
    VisualizeGeneTSNE(gene_embedding = tsne_list[[lab]],
                      gene_partition = gene_modules,
                      text = FALSE) + ggtitle(lab)
  }), ncol = 2) + plot_layout(guides = "collect")

  list(modules    = gene_modules,
       eigen      = eigen_list,
       n_eigvec   = n_eigvec,
       indicators = gene_indictors,
       panels     = list(A = panelA, B = panelB, C = panelC, D = panelD))
}


# -----------------------------------------------------------------------------
# LoadAlignedGeneEMD
#
# Load the per-condition OT distance matrices and restrict them to the genes
# present in both, in a common order. Returns a named list.
# -----------------------------------------------------------------------------
LoadAlignedGeneEMD <- function(emd_paths) {
  emd_list <- lapply(emd_paths, function(p) LoadGeneEMD(file.path(p, "")))
  names(emd_list) <- names(emd_paths)
  g <- Reduce(intersect, lapply(emd_list, rownames))
  lapply(emd_list, function(x) x[g, g])
}
