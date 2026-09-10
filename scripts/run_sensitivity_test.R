# =============================================================================
# run_sensitivity_test.R   (Fig S7)
#
# Robustness of DiCoLo to cell-graph construction parameters, on the SmoM2
# dataset (differential operator of SmoM2 relative to CTL). Stability is the
# Jaccard index of the selected genes against the default setting
# (npc = 10, K = 10, reduction = "dm"), coloured by the number of leading
# eigenvectors used for gene selection.
#
#   Left   (Ndim)       varying diffusion-map components, K = 10
#   Middle (KNN)        varying cell-graph neighborhood size, npc = 10
#   Right  (Embeddings) top-10 PCs vs the default top-10 dm
# =============================================================================

# Load library ----
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

# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()
dir.path      <- DATA_DIR
figure.path   <- FIGURE_DIR
script_dir    <- OT_SCRIPT_DIR
dir.path.data <- file.path(DATA_DIR, "smom2")

# Load SmoM2 data -----


data_S_MUT <- readRDS(file = file.path(dir.path.data,"data_S_smom2_dermal_E13.5_MUT.rds"))
data_S_CTL <- readRDS(file = file.path(dir.path.data,"data_S_smom2_dermal_E13.5_CTL.rds"))
common_genes = SelectCommonGenes(data_S_MUT,data_S_CTL,
                                 ngenes = 500)

# Test PC
for(npc in c(5,10,15,20,25)){
  dir.path.MUT = file.path(dir.path,"sensitivity_test",sprintf("GeneTrajectory_MUT_%spc",npc))
  dir.path.CTL = file.path(dir.path,"sensitivity_test",sprintf("GeneTrajectory_CTL_%spc",npc))
  if(!dir.exists(dir.path.MUT) & !file.exists(file.path(dir.path.MUT,"emd.csv"))){
    ComputeGeneEMD(data_S_MUT, common_genes, dir.path = dir.path.MUT, 
                   script_dir = script_dir, npc = npc, K = 10)
  }
  if(!dir.exists(dir.path.CTL) & !file.exists(file.path(dir.path.CTL,"emd.csv"))){
    ComputeGeneEMD(data_S_CTL, common_genes, dir.path = dir.path.CTL, 
                   script_dir = script_dir, npc = npc, K = 10)
  }
}

# Test kNN
for(K in c(5,10,15,20,25)){
  dir.path.MUT = file.path(dir.path,"sensitivity_test",sprintf("GeneTrajectory_MUT_%sknn",K))
  dir.path.CTL = file.path(dir.path,"sensitivity_test",sprintf("GeneTrajectory_CTL_%sknn",K))
  if(!dir.exists(dir.path.MUT) & !file.exists(file.path(dir.path.MUT,"emd.csv"))){
    ComputeGeneEMD(data_S_MUT, common_genes, dir.path = dir.path.MUT, 
                   script_dir = script_dir, npc = 10, K = K)
  }
  if(!dir.exists(dir.path.CTL) & !file.exists(file.path(dir.path.CTL,"emd.csv"))){
    ComputeGeneEMD(data_S_CTL, common_genes, dir.path = dir.path.CTL, 
                   script_dir = script_dir, npc = 10, K = K)
  }
}

# Test Different embeddings (top10 PC vs. top10 DM)
dir.path.MUT = file.path(dir.path,"sensitivity_test","GeneTrajectory_MUT_10PC_only")
dir.path.CTL = file.path(dir.path,"sensitivity_test","GeneTrajectory_CTL_10PC_only")
if(!dir.exists(dir.path.MUT) & !file.exists(file.path(dir.path.MUT,"emd.csv"))){
  ComputeGeneEMD(data_S_MUT, common_genes, dir.path = dir.path.MUT, 
                 script_dir = script_dir, reduction = "pca")
}
if(!dir.exists(dir.path.CTL) & !file.exists(file.path(dir.path.CTL,"emd.csv"))){
  ComputeGeneEMD(data_S_CTL, common_genes, dir.path = dir.path.CTL, 
                 script_dir = script_dir, reduction = "pca")
}

# Sensitivity Test
df_sensitivity = do.call(rbind,lapply(c("pc","knn"),function(test_id){
  gene_ls = lapply(c(5,10,15,20,25), function(id){
    dir.path.MUT = file.path(dir.path,"sensitivity_test",sprintf("GeneTrajectory_MUT_%s%s",id,test_id))
    dir.path.CTL = file.path(dir.path,"sensitivity_test",sprintf("GeneTrajectory_CTL_%s%s",id,test_id))
    
    gene.dist.mat_MUT = LoadGeneEMD(file.path(dir.path.MUT,""))
    gene.dist.mat_CTL = LoadGeneEMD(file.path(dir.path.CTL,""))
    
    gene.graph_MUT = ComputeGraphOperator(gene.dist.mat_MUT)
    gene.graph_CTL = ComputeGraphOperator(gene.dist.mat_CTL)
    
    diff.op_MUT = ComputeDifferentialOperator(gene.graph_MUT,gene.graph_CTL)
    # diff.op_CTL = ComputeDifferentialOperator(gene.graph_CTL,gene.graph_MUT)
    
    diff_op = diff.op_MUT
    eigen_list <- RunSVD(diff_op, eig_keep = nrow(diff_op))
    gene_ls = lapply(c(2:8),function(n_eigvec){
      # n_eigvec = FindKneePoint(eigen_list$values[1:floor(sqrt(length(eigen_list$values)))],plot.fig = TRUE) - 1
      gene_loadings = eigen_list$vectors[,1:n_eigvec]
      gene_indictors = SelectSignificantGenes(gene_loadings, lfdr_thresh = )
      signf_genes = rownames(gene_indictors)[rowSums(gene_indictors) > 0]
      signf_genes
    })
    names(gene_ls) = paste0("top_eigvec",c(2:8))
    
    gene_ls
  })
  names(gene_ls) = paste0(test_id,c(5,10,15,20,25))
  baseline = paste0(test_id,10)
  df_sensitivity = do.call(rbind,lapply(names(gene_ls),function(para){
    gene_set =gene_ls[[para]]
    jaccard_idx = unlist(lapply(names(gene_set),function(eig_id){
      genes = gene_set[[eig_id]]
      genes_baseline = gene_ls[[baseline]][[eig_id]]
      length(intersect(genes,genes_baseline))/
        length(union(genes,genes_baseline))
    }))
    data.frame(Variable = para,
               TopN_eigvec = names(gene_set),
               Jaccard = jaccard_idx)  
    
  }))
  df_sensitivity
}))
df_sensitivity$'Type' = gsub("[0-9]", "", df_sensitivity$Variable)
df_sensitivity$'Value' = gsub("[^0-9]", "", df_sensitivity$Variable)
df_sensitivity$TopN_eigvec = gsub("[^0-9]", "", df_sensitivity$TopN_eigvec)

write.csv(df_sensitivity, file = file.path(dir.path,"sensitivity_test","jaccard_result_ndim_knn.csv"), 
          row.names = FALSE)

# Sensitivity Test to DM vs. PC
gene_ls = lapply(c("10PC_only","10pc"),function(reduc){
  dir.path.MUT = file.path(dir.path,"sensitivity_test",sprintf("GeneTrajectory_MUT_%s",reduc))
  dir.path.CTL = file.path(dir.path,"sensitivity_test",sprintf("GeneTrajectory_CTL_%s",reduc))
  
  gene.dist.mat_MUT = LoadGeneEMD(file.path(dir.path.MUT,""))
  gene.dist.mat_CTL = LoadGeneEMD(file.path(dir.path.CTL,""))
  
  gene.graph_MUT = ComputeGraphOperator(gene.dist.mat_MUT)
  gene.graph_CTL = ComputeGraphOperator(gene.dist.mat_CTL)
  
  diff.op_MUT = ComputeDifferentialOperator(gene.graph_MUT,gene.graph_CTL)
  # diff.op_CTL = ComputeDifferentialOperator(gene.graph_CTL,gene.graph_MUT)
  
  diff_op = diff.op_MUT
  eigen_list <- RunSVD(diff_op, eig_keep = nrow(diff_op))
  gene_ls = lapply(c(2:8),function(n_eigvec){
    # n_eigvec = FindKneePoint(eigen_list$values[1:floor(sqrt(length(eigen_list$values)))],plot.fig = TRUE) - 1
    gene_loadings = eigen_list$vectors[,1:n_eigvec]
    gene_indictors = SelectSignificantGenes(gene_loadings, lfdr_thresh = )
    signf_genes = rownames(gene_indictors)[rowSums(gene_indictors) > 0]
    signf_genes
  })
  names(gene_ls) = paste0("top_eigvec",c(2:8))
  gene_ls
})
names(gene_ls) = c("10PC_only","10pc")
baseline = "10pc"
df_sensitivity = do.call(rbind,lapply(names(gene_ls),function(para){
  gene_set =gene_ls[[para]]
  jaccard_idx = unlist(lapply(names(gene_set),function(eig_id){
    genes = gene_set[[eig_id]]
    genes_baseline = gene_ls[[baseline]][[eig_id]]
    length(intersect(genes,genes_baseline))/
      length(union(genes,genes_baseline))
  }))
  data.frame(Variable = para,
             TopN_eigvec = names(gene_set),
             Jaccard = jaccard_idx)  
  
}))
df_sensitivity$TopN_eigvec = gsub("[^0-9]", "", df_sensitivity$TopN_eigvec)
df_sensitivity = df_sensitivity %>% filter(Variable == "10PC_only") %>% select(TopN_eigvec,Jaccard)
write.csv(df_sensitivity, file = file.path(dir.path,"sensitivity_test","jaccard_result_dm_vs_pc.csv"), 
          row.names = FALSE)

# Visualize
df_sensitivity = read.csv(file.path(dir.path,"sensitivity_test","jaccard_result_ndim_knn.csv"))
df_sensitivity$Type = toupper(df_sensitivity$Type); df_sensitivity$Type = gsub("PC","Ndim",df_sensitivity$Type)
df_sensitivity$Value = as.factor(as.numeric(df_sensitivity$Value))
df_sensitivity1 = read.csv(file.path(dir.path,"sensitivity_test","jaccard_result_dm_vs_pc.csv"))
df_sensitivity = rbind(df_sensitivity, 
                       data.frame(Variable = NA, df_sensitivity1, Type = "Embeddings", Value = "Top10 PCs"))
df_sensitivity$Type = factor(df_sensitivity$Type,levels = c("Ndim","KNN","Embeddings"))
df_sensitivity$TopN_eigvec = as.factor(df_sensitivity$TopN_eigvec)
df_sensitivity = df_sensitivity %>% filter(TopN_eigvec %in% c(2,3,4,5))
p1 = ggplot(df_sensitivity, aes(x = Value, y = Jaccard, 
                           color = TopN_eigvec, group = TopN_eigvec)) +
  facet_wrap(~Type, scales = "free_x") +
  ylim(0,1) +
  geom_line() +
  geom_point() +
  labs(
    # title = "Robustness of Mutant-Specific Gene Detection",
    x = "Parameter Value",
    y = "Jaccard Index (vs. Default)",
    color = "TopN Eigenvectors"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
png(file.path(figure.path,"figS7_sensitivity.png"),width = 300*8, height = 300*4, res = 300)
p1
dev.off()

