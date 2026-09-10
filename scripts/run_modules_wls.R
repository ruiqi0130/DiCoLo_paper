# =============================================================================
# run_modules_wls.R   (Fig S5, Fig S6)
#
# Wls-KO vs CTL analysis:
#   condition UMAP, cell types, module activity scores -> Fig S5A / S5B
#   MsigDB pathway enrichment of the modules           -> Fig S5C
#   eigenvalue spectrum, gene embeddings, locFDR,
#   hierarchical clustering                            -> Fig S6A-D
#
# -----------------------------------------------------------------------------
# STATUS: pending migration to the DiCoLo package API.
#
# This script still calls helper functions from an earlier internal codebase
# rather than their DiCoLo equivalents, and is not runnable as written. Before
# release:
#
#   1. Replace each call tagged TODO(api) with its DiCoLo equivalent
#      (ComputeGeneEMD, LoadGeneEMD, SelectCommonGenes, ComputeGraphOperator,
#      ComputeDifferentialOperator, RunSVD, SelectSignificantGenes,
#      Find_kneepoint, VisualizeGeneHeatmap, ...).
#   2. Point the paths tagged TODO(path) at DATA_DIR.
#   3. Give the outputs tagged TODO(name) their manuscript figure names.
#   4. Run end to end and confirm the saved *_dicolo_modules.rds matches the
#      copy distributed in DiCoLo_data.
# -----------------------------------------------------------------------------
# =============================================================================

# --- Paths, python backend, helper functions (see config.R at repo root) -----
source("config.R")
source_helpers()
figure.path <- FIGURE_DIR

require(GeneTrajectory)
require(LocalizedMarkerDetector)
require(plot3D)
require(plotly)
require(ggplot2)
require(patchwork)
require(viridis)
require(scales)
require(Seurat)
require(SeuratWrappers)
require(RColorBrewer)

# TODO(api): replace with the DiCoLo equivalent of differential_graph_function.R
# TODO(api): replace with the DiCoLo equivalent of simulation_function.R
# TODO(api): replace with the DiCoLo equivalent of benchmarking_function.R
# TODO(api): replace with the DiCoLo equivalent of load_paths.R

# Load samples -------
sample_ls = paste0("E14.5_",c("MUT","CTL"))
sample_ls = paste0("wls_dermal_",sample_ls)
sample_ls_short = c("MUT","CTL")
dir.path.data = file.path(DATA_DIR, "wlsko")   # TODO(path): confirm subdirectory
dir.path = file.path(DATA_DIR, "wlsko")   # TODO(path): confirm subdirectory
data_S_ls = lapply(sample_ls,function(sample_name){
  readRDS(file.path(dir.path.data,"process_data",sprintf("data_S_%s.rds",sample_name)))
})
names(data_S_ls) = sample_ls

# Visualize Merged embedding -----
data_S_merge = merge(data_S_ls[[1]],data_S_ls[-1])
data_S_merge = data_S_merge %>% NormalizeData() %>% 
  FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE) %>% 
  RunUMAP(dims = 1:15,seed = 42)
png(file.path(figure.path,"TODO_name_f5-0.png"),width = 300*5, height = 300*3, res = 300)
DimPlot(data_S_merge, group.by = 'orig.ident', shuffle = TRUE) + NoAxes() + 
  ggtitle(NULL)
dev.off()

# RunEMD -----
## Select genes ----
common_genes = Select_common_genes(data_S_ls[[1]],data_S_ls[[2]],
                                   ngenes = 500)
## Gene-gene EMD ----
for(sample_name in sample_ls){
  folder.path = file.path(dir.path,"GeneTrajectory",paste0("GeneTrajectory_",sample_name),"")
  # folder.path = file.path(dir.path,"GeneTrajectory","trajectory_compare_denoise",paste0("GeneTrajectory_",sample_name),"")
  Obtain_GeneEMD(data_S_ls[[sample_name]], common_genes, backend = TRUE, 
                 dir.path = folder.path)
}

# Load gene EMD distance
gene_emd_ls = lapply(sample_ls,function(sample_name){
  # folder.path = file.path(dir.path,"GeneTrajectory","trajectory_compare")
  # dir.path = file.path(dir.path,"GeneTrajectory","trajectory_compare_denoise")
  folder.path <- file.path(dir.path,"GeneTrajectory",paste0("GeneTrajectory_",sample_name,"/"))
  Load_GeneEMD(folder.path)
})
# Align gene name
g = Reduce(intersect, lapply(gene_emd_ls,function(x) rownames(x)))
gene_emd_ls = lapply(gene_emd_ls, function(x) x[g,g])

# Load gene TSNE embedding
gene_tsne_ls = lapply(gene_emd_ls,function(x){
  ObtainGeneTSNE(as.dist(x))
})

# DUFS ----------
gene_graph_ls = lapply(gene_emd_ls, Obtain_gene_laplacian)

E.list_1 = Run_DUFS(gene_graph_ls[[1]],gene_graph_ls[[2]], eig_keep = nrow(gene_graph_ls[[1]]))
E.list_2 = Run_DUFS(gene_graph_ls[[2]],gene_graph_ls[[1]], eig_keep = nrow(gene_graph_ls[[1]]))

sample_diff = 2
E.list = get(paste0("E.list_",sample_diff))

## knee-point cutoff --------
png(file.path(figure.path,"TODO_name_f5-sup1.png"),width = 300*5, height = 300*5, res = 300)
n_eigvec_knee = Find_kneepoint(E.list$values[1:floor(sqrt(length(E.list$values)))],plot.fig = TRUE) - 1
dev.off()

## Eigen vecs ----------
n_eigvec = n_eigvec_knee
pl = lapply(1:n_eigvec,function(i){
  gene_partition = setNames(E.list$vectors[,i],rownames(E.list$vectors))
  p1 = VisualizeGeneTSNE(gene_embedding = gene_tsne_ls[[1]],
                         gene_partition = gene_partition, filtered = FALSE) + 
    ggtitle(sample_ls_short[1]) +
    labs(color = paste0("Eig_vec",i))
  p2 = VisualizeGeneTSNE(gene_embedding = gene_tsne_ls[[2]],
                         gene_partition = gene_partition, filtered = FALSE) + 
    ggtitle(sample_ls_short[2]) +
    labs(color = paste0("Eig_vec",i))
  p1 / p2
})
png(file.path(figure.path,"TODO_name_f5-sup2.png"),width = 300*4*n_eigvec, height = 300*6, res = 300)
wrap_plots(pl, nrow = 1)
dev.off()

## locfdr ----------
gene_loadings = E.list$vectors[,1:n_eigvec]
res = select_modules_locfdr(gene_loadings,lfdr_thresh = 0.2)
modules = res$'modules'; gene_indictors = res$gene_indicator
pl = lapply(1:ncol(gene_indictors),function(i){
  VisualizeGeneTSNE(gene_embedding = gene_tsne_ls[[sample_diff]],
                    gene_partition = as.factor(gene_indictors[,i]), 
                    filtered = FALSE, text = FALSE,
                    module_color = c("1"="red","0"="lightgrey")) + 
    labs(color = "Signif. genes") + ggtitle(NULL)
})
png(file.path(figure.path,"TODO_name_f5-sup3.png"),width = 300*4*n_eigvec, height = 300*3, res = 300)
wrap_plots(pl, nrow = 1) + plot_layout(guides = "collect")
dev.off()

## Obtain modules -------
signf_genes = unique(unlist(modules, use.names = FALSE))
gene_dist = as.dist(gene_emd_ls[[sample_diff]][signf_genes,signf_genes])
merged_modules = ClusterGenes_simplify(gene_dist,min_gene = 5, deepSplit = 0)
merged_modules = rank_modules(merged_modules,res$gene_indicator,E.list$values)
levels(merged_modules) = paste0(sample_ls_short[sample_diff],levels(merged_modules))
saveRDS(merged_modules,file = file.path(dir.path,paste0(sample_ls_short[sample_diff],"_dicolo_modules.rds")))

pl = lapply(1:length(gene_tsne_ls),function(i){
  VisualizeGeneTSNE(gene_embedding = gene_tsne_ls[[i]],
                    gene_partition = merged_modules, 
                    filtered = FALSE, text = FALSE) + 
    ggtitle(NULL)
  # ggtitle(sample_ls_short[i])
})
png(file.path(figure.path,"TODO_name_f5-sup4.png"),width = 300*10, height = 300*4, res = 300)
wrap_plots(pl,ncol = 2) + plot_layout(guides = "collect")
dev.off()

## Gene Modules on CellUmap ----------
merged_modules_ls = split(names(merged_modules), merged_modules)
lapply(names(merged_modules_ls),function(i){
  message(i,": ",paste0(merged_modules_ls[[i]],collapse = ", "),"\n")
})
module_prefix = "Module_"
data_S_ls = lapply(data_S_ls, function(data_S){
  DefaultAssay(data_S) = "RNA"
  # data_S = AddModuleActivityScore(data_S, gene_partition = merged_modules)
  data_S <- AddModuleScore(
    object = data_S,
    features = merged_modules_ls,
    ctrl = 100,
    name = module_prefix
  )
  module_score = data_S@meta.data[,grepl(module_prefix,colnames(data_S@meta.data))]
  module_score = apply(module_score,2,function(x) scale(x)[,1])
  module_score = apply(module_score,2,function(x) (x - min(x)) / (max(x) - min(x)))
  data_S = AddMetaData(data_S,module_score)
})
pl = lapply(data_S_ls,function(data_S){
  FeaturePlot(data_S, features = paste0(module_prefix,seq_len(length(merged_modules_ls))), 
              # cols = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(100),
              ncol = length(merged_modules_ls)) & 
    scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds")) & NoAxes()
})
png(file.path(figure.path,"TODO_name_f5-1.png"),width = 300*4*length(merged_modules_ls), height = 300*8, res = 300)
(wrap_plots(pl,nrow = 2) & NoAxes() & ggtitle(NULL)) + plot_layout(guides = "collect")
dev.off()

color_ct = setNames(c("darkred","darkblue","#AC9362"),c("DC","LD","UD"))
pl = lapply(data_S_ls,function(data_S){
  DimPlot(data_S, group.by = "celltype", cols = color_ct) + NoAxes() + ggtitle(NULL)
})
png(file.path(figure.path,"TODO_name_f4-1-1.png"),width = 300*4, height = 300*6, res = 300)
(wrap_plots(pl,nrow = 2) & NoAxes() & ggtitle(NULL)) + plot_layout(guides = "collect")
dev.off()

## Pathway analysis --------
merged_modules = do.call(c,lapply(sample_ls_short,function(sample_name){
  readRDS(file.path(dir.path,paste0(sample_name,"_dicolo_modules.rds")))
}))
merged_modules_ls = split(names(merged_modules), merged_modules)
library("ReactomePA")
library("org.Mm.eg.db")
library("clusterProfiler")
library("msigdbr")
m_df = msigdbr(species = "Mus musculus", category = "H")
m_df$gs_name = gsub("HALLMARK_","",m_df$gs_name)
msig_result <- do.call(rbind,lapply(names(merged_modules_ls)[5:6], function(i){
  df = enricher(gene=merged_modules_ls[[i]],
                TERM2GENE = m_df[,c("gs_name","gene_symbol")],
                pAdjustMethod = "none",pvalueCutoff = 0.1)@result
  df$'modules' = i
  df = df %>% filter(p.adjust < 0.2)
  # df = df[1:5,]
  df <- df[complete.cases(df), ]
  df
}))
msig_result$GeneRatio = unlist(lapply(msig_result$GeneRatio,function(x) eval(parse(text = x))))
msig_result$logP = -log(msig_result$pvalue,base = 10)
msig_result$Description <- factor(msig_result$Description,
                                  levels = rev(unique(msig_result$Description[order(msig_result$pvalue)])))

msig_result$modules = factor(msig_result$modules,levels = names(merged_modules_ls)[5:6])
levels(msig_result$modules) = gsub("MUT","Wlsko_M",levels(msig_result$modules))
levels(msig_result$modules) = gsub("CTL","CTL_M",levels(msig_result$modules))
# Dotplot
p = ggplot(msig_result, aes(x=modules, y=Description)) +
  geom_point(aes(size=GeneRatio, color=logP), alpha=0.7) +
  scale_color_gradientn(colors=colorRampPalette(c("blue","red"))(99), name="-log10_pvalue") +
  scale_size_continuous(name = "Gene Ratio", limits = c(0, 1)) + 
  theme_minimal() +
  labs(size="Gene Ratio", x=NULL,y=NULL) + 
  theme(legend.position = "bottom",
        legend.box = "vertical",
        axis.text.x = element_text(size = 10, angle = 45,hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10)
  )
ggsave(filename = file.path(figure.path,"f5-3.png"), 
       plot = p, width = 8, height = 5)
