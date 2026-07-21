library(ggplotify)
VisualizeGeneHeatmap = function (mat, gene_partition = NULL, 
                                 gene_order = NULL, clustering_method = "complete",
                                 highlighted_genes = NULL,module_color = NULL) 
{
  if(is.null(gene_order)){
    gene_order = 1:nrow(mat); cluster_rows = TRUE; cluster_cols = TRUE
  }else{
    cluster_rows = FALSE; cluster_cols = FALSE; clustering_method = NULL
  }
  if(!is.null(gene_partition)){
    gene_partition = as.factor(gene_partition)
    if(is.null(module_color)){
      module_color = setNames(colorRampPalette(brewer.pal(12, "Paired"))(nlevels(gene_partition)), 
                              levels(gene_partition))
      module_color[is.na(names(module_color))] = "lightgrey"
      module_color[names(module_color) == "Other"] = "lightgrey"
    }
    annotation_col = data.frame(Module = gene_partition)
    annotation_colors = list(Module = module_color)
  }else{
    annotation_col = NULL; annotation_colors = NULL
  }
  if(!is.null(highlighted_genes)){
    rownames(mat) = ifelse(rownames(mat) %in% highlighted_genes, rownames(mat),"")
    show_rownames = TRUE
  }else{
    show_rownames = FALSE
  }
  
  nice_ticks <- seq(0, max(mat[gene_order,gene_order]), length.out = 2)
  
  p = pheatmap::pheatmap(mat[gene_order,gene_order],
                         annotation_col = annotation_col, 
                         annotation_colors = annotation_colors,
                         cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                         treeheight_row = 0, treeheight_col = 0, 
                         legend_breaks = nice_ticks,
                         legend_labels = round(nice_ticks,3),
                         clustering_method = clustering_method,
                         # color = colorRampPalette(rev(brewer.pal(n = 7, name = "Reds")))(100),
                         color = colorRampPalette(c("white", "#ffcccc", "#ff6666", "#cc0000"))(100),
                         show_colnames = FALSE, show_rownames = show_rownames,
                         silent = TRUE)
  return(p)
}


Visualize_Eigenvec <- function(eig_vec,gene_partition){
  df = data.frame(x = 1:length(eig_vec), y = eig_vec)
  sizes = table(gene_partition)
  edges <- c(0, cumsum(sizes))
  bands <- data.frame(
    label = names(sizes),
    xmin  = head(edges, -1),
    xmax  = tail(edges, -1)
  )
  cols <- color_m[names(sizes)]
  cols <- sapply(cols, function(clr) {
    colorRampPalette(c("white", clr))(3)[3]   # pick a lighter tone
  })
  
  p = ggplot(df, aes(x, y)) +
    geom_rect(
      data = bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = label),
      inherit.aes = FALSE,
      alpha = 0.3
    ) +
    geom_line(color = "black", linewidth = 0.6) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    labs(x = NULL, y = NULL) + 
    theme(
      axis.text.x  = element_text(size = 15),
      axis.text.y  = element_text(size = 15),
      # axis.line    = element_line(linewidth = 0.8),
      # axis.ticks   = element_line(linewidth = 0.7)
    )
  return(p)
}
