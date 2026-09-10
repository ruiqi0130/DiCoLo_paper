# =============================================================================
# run_modules_smom2.R   (Fig S2, Fig S3, Supplemental Table 2)
#
# Derivation of the SmoM2-specific and CTL-specific co-localized gene modules:
#   eigenvalue spectrum with knee-point selection -> Fig S2A / S3A
#   gene t-SNE coloured by eigenvector loadings   -> Fig S2B / S3B
#   locFDR candidate genes                        -> Fig S2C / S3C
#   hierarchical clustering into modules          -> Fig S2D / S3D
#
# Also writes <MUT|CTL>_dicolo_modules.rds, the upstream dependency of
# Fig 4B, 4C, 4D, Fig S4 and run_concordance_analysis.R, and the source of
# the gene lists in Supplemental Table 2.
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

# Construct gene EMD distance by GeneTrajectory --------
# remotes::install_github("KlugerLab/GeneTrajectory")
# remotes::install_github("KlugerLab/LocalizedMarkerDetector")
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
sample_ls = paste0("E13.5_",c("MUT","CTL"))
sample_ls = paste0("smom2_dermal_",sample_ls)
sample_ls_short = c("MUT","CTL")
dir.path.data = file.path(DATA_DIR, "smom2")   # TODO(path): confirm subdirectory
dir.path = file.path(DATA_DIR, "smom2")   # TODO(path): confirm subdirectory
sample_name = sample_ls[1]
data_S_ls = lapply(sample_ls,function(sample_name){
  readRDS(file.path(dir.path.data,"process_data",sprintf("data_S_%s.rds",sample_name)))
})
names(data_S_ls) = sample_ls
# data_S <- data_S_ls[[sample_name]]
# DimPlot(data_S, reduction = "umap", group.by = "Phase")

# Visualize Merged embedding -----
data_S_merge = merge(data_S_ls[[1]],data_S_ls[-1])
data_S_merge = data_S_merge %>% NormalizeData() %>% 
  FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE) %>% 
  RunUMAP(dims = 1:15,seed = 42)
png(file.path(figure.path,"TODO_name_f4-0.png"),width = 300*5, height = 300*3, res = 300)
DimPlot(data_S_merge, group.by = 'orig.ident', shuffle = TRUE) + NoAxes() + 
  ggtitle(NULL)
dev.off()

# Load reference gene modules ---------
dc_genes = read.csv(file.path(DATA_DIR, "smom2", "gene_list_dc.csv"))   # TODO(path)[,2]
arrested_genes = read.csv(file.path(DATA_DIR, "smom2", "gene_list_gli3KO.csv"))   # TODO(path)[,2]
shh_genes = c("Ptch1","Dkk1","Gli1")
# Add module score
data_S_ls = lapply(data_S_ls, function(data_S){
  data_S <- AddModuleScore(
    object = data_S,
    features = list(dc_genes, arrested_genes),
    ctrl = 100,
    name = 'RefModule'
  )
  data_S
})
pl = lapply(data_S_ls,function(data_S){
  FeaturePlot(data_S, paste0("RefModule",1:2), 
              cols = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(100)) & NoAxes()
})
wrap_plots(pl, nrow = 2)

# RunDUFS -----
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

# Load gene DM embedding
gene_dm_ls = lapply(gene_emd_ls,function(gene.dist.mat){
  GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
})

# Load gene trajectories
gene_trajectory_ls = lapply(1:2, function(i){
  if(i == 1){
    # E13.5 MUT
    # gene_trajectory <- ExtractGeneTrajectory(gene_dm_ls[[i]], gene_emd_ls[[i]], N = 4, t.list = c(1,3,2,2), K = 5) # before alra
    gene_trajectory <- ExtractGeneTrajectory(gene_dm_ls[[i]], gene_emd_ls[[i]], N = 3, t.list = c(2,1,1), K = 5) # after alra
  }else{
    # E13.5 CTL
    # gene_trajectory <- ExtractGeneTrajectory(gene_dm_ls[[i]], gene_emd_ls[[i]], N = 3, t.list = c(3,2,2), K = 5) # before alra
    gene_trajectory <- ExtractGeneTrajectory(gene_dm_ls[[i]], gene_emd_ls[[i]], N = 3, t.list = c(1,2,1), K = 5) # after alra
  }
  gene_trajectory
})

## DUFS ----------
gene_graph_ls = lapply(gene_emd_ls, Obtain_gene_laplacian)
# Eigen decomposition
eigen_ls = lapply(gene_graph_ls, Run_SVD)

## Differential operator --------
E.list_1 = Run_DUFS(gene_graph_ls[[1]],gene_graph_ls[[2]], eig_keep = nrow(gene_graph_ls[[1]]))
E.list_2 = Run_DUFS(gene_graph_ls[[2]],gene_graph_ls[[1]], eig_keep = nrow(gene_graph_ls[[1]]))
E.list_common = Run_Pcommon(gene_graph_ls[[1]],gene_graph_ls[[2]])


# Visualize Genes-------
E.list = E.list_1

##  Eigen values -----
regulator = 0.047
png("~/1.png",width = 300*5, height = 300*5, res = 300)
plot(E.list$values, 
     main = paste0("regulator ",regulator))
dev.off()

## Gene trajectory on Gene TSNE -------
pl = lapply(1:2,function(i){
  gene_partition = as.factor(setNames(gene_trajectory_ls[[i]]$selected,rownames(gene_trajectory_ls[[i]])))
  levels(gene_partition)
  levels(gene_partition)[1] = NA
  if(i == 1){
    module_color = setNames(colorRampPalette(brewer.pal(8, "Set1"))(nlevels(gene_partition)), 
                            levels(gene_partition))
  }else{
    module_color = setNames(colorRampPalette(brewer.pal(12, "Paired"))(nlevels(gene_partition)), 
                            levels(gene_partition))
  }
  
  pl = lapply(1:2,function(j){
    VisualizeGeneTSNE(gene_embedding = gene_tsne_ls[[j]],
                      gene_partition = gene_partition, filtered = FALSE) + 
      scale_color_manual(values = module_color,
                         na.value = "lightgrey") + 
      ggtitle(sprintf("%s (color by %s)",sample_ls_short[j],sample_ls_short[i]))
  })
  pl
})
pl = unlist(pl, recursive = FALSE)
wrap_plots(pl[c(1,4,3,2)], nrow = 2) + plot_layout(guides = "collect")


## Eigen vecs of gene graph -------
pl = lapply(1:2,function(i){
  E.list = eigen_ls[[i]]
  pl = lapply(1:2,function(j){
    pl = lapply(2:6,function(eig_id){
      gene_partition = setNames(E.list$vectors[,eig_id],rownames(E.list$vectors))
      VisualizeGeneTSNE(gene_embedding = gene_tsne_ls[[j]],
                        gene_partition = gene_partition, 
                        filtered = FALSE) + 
        ggtitle(sprintf("%s (color by eigvec of %s)",
                        sample_ls_short[j],
                        sample_ls_short[i])) + 
        labs(color = paste0("Eig_vec",eig_id))
    })
    wrap_plots(pl, ncol = 1)
  })
})
pl = unlist(pl, recursive = FALSE)
wrap_plots(pl[1:2], nrow = 1) + plot_layout(guides = "collect")
wrap_plots(pl[3:4], nrow = 1) + plot_layout(guides = "collect")


## P_diff Eigen vecs ----------
pl = lapply(1:5,function(i){
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
png(file.path(figure.path,"TODO_name_f4-sup2.png"),width = 300*20, height = 300*6, res = 300)
wrap_plots(pl, nrow = 1)
dev.off()

## Partition an eigen vector to gene modules ------
sample_id = 1
eig_vec_id = 1
res = plot_discrete_gene_embed(gene_dist = gene_emd_ls[[sample_id]], 
                               gene_embedding = gene_tsne_ls[[sample_id]], 
                               eig_vec = E.list$vectors[,eig_vec_id],
                               lower_bound = -0.05, upper_bound = 0.05, 
                               min_gene = 3)
res$plot
pl = CustomModulePlot(data_S_ls[[sample_id]], 
                      feature_partition = res$gene_discrete$gene_partition, 
                      assays = "alra")
wrap_plots(pl, nrow = 2)
FeaturePlot(data_S_ls[[i]], features = names(res$gene_discrete$gene_partition)[res$gene_discrete$gene_partition %in% "neg_2"],
            order = TRUE, ncol = 4) & NoAxes()


## Meta-features --------
eigvec_selected = 1:5
pl = lapply(data_S_ls,function(srat){
  cos_sim = data.frame(compute_cell_scores(srat, eigen_vec = E.list$vectors[,eigvec_selected], method = "cosine"))
  srat = AddMetaData(srat, metadata = cos_sim)
  FeaturePlot(srat, features = colnames(cos_sim), ncol = length(eigvec_selected)) & NoAxes() & 
    scale_color_gradient2(low = "blue",high = "red",mid = "white",midpoint = 0)
})
png("~/1.png",width = 300*20, height = 300*6, res = 300)
wrap_plots(pl, nrow = 2)
dev.off()

# Obtain Gene modules ----------
sample_diff = 1
E.list = get(paste0("E.list_",sample_diff))

png("~/1.png",width = 300*5, height = 300*5, res = 300)
plot(E.list$values)
dev.off()

## Define cutoff ----------
### kneepoint cutoff -------
png(file.path(figure.path,"TODO_name_f4-sup1.png"),width = 300*5, height = 300*5, res = 300)
n_eigvec_knee = Find_kneepoint(E.list$values[1:floor(sqrt(length(E.list$values)))],plot.fig = TRUE) - 1
dev.off()

### Noise floor cutoff ------
common_genes = do.call(union,lapply(gene_graph_ls,function(x) rownames(x)))
data_S_ls_test = lapply(data_S_ls,function(data_S){
  data_S = Denoise_seurat(data_S, method = "alra")
  DefaultAssay(data_S) = "alra"
  data_S
})
names(data_S_ls_test) = sample_ls
first_eigval_ls = lapply(1:2,function(i){
  permute.path = file.path(dir.path,"dropout_ds",paste0("dropout_ds",i))
  dir.create(permute.path,recursive = TRUE)
  permute_E_list = Permutation_test(data_S_ls_test[[i]],
                                    input_genes = common_genes,
                                    res.path = permute.path, seed = 233, 
                                    n_downsample = 5,nrsvd = FALSE)
  first_eigval_ls = unlist(lapply(permute_E_list,function(x) x$values[1]))
  first_eigval_ls
})

png("~/1.png",width = 300*5, height = 300*5, res = 300)
plot(seq_along(E.list$values[1:30]), E.list$values[1:30], pch = 16,
     xlab = "Eigenvalue index", ylab = "Eigenvalue",
     ylim = range(c(E.list$values[1:30], first_eigval_ls[[1]],
                    first_eigval_ls[[2]])))
abline(h = first_eigval_ls[[1]], lty = 3, col = "red")
abline(h = first_eigval_ls[[2]], lty = 3, col = "blue")
abline(h = mean(c(first_eigval_ls[[1]],
                  first_eigval_ls[[2]])), lty = 1, col = "green")
dev.off()

n_eigvec_noise_floor = sum(E.list$values > mean(c(first_eigval_ls[[1]],first_eigval_ls[[2]])))

## locfdr ----------
# n_eigvec = min(n_eigvec_noise_floor,n_eigvec_knee)
n_eigvec = n_eigvec_knee
gene_loadings = E.list$vectors[,1:n_eigvec]
res = select_modules_locfdr(gene_loadings,lfdr_thresh = 0.2)
modules = res$'modules'; gene_indictors = res$gene_indicator
pl = lapply(names(modules),function(module){
  g = modules[[module]]
  VisualizeGeneTSNE(gene_embedding = gene_tsne_ls[[sample_diff]],
                    gene_partition = setNames(rep(module,length(g)),g), 
                    filtered = FALSE)
})
png("~/1.png",width = 300*4*9, height = 300*3*2, res = 300)
wrap_plots(pl,nrow = 2)
dev.off()

pl = lapply(1:ncol(gene_indictors),function(i){
  VisualizeGeneTSNE(gene_embedding = gene_tsne_ls[[sample_diff]],
                    gene_partition = as.factor(gene_indictors[,i]), 
                    filtered = FALSE, text = FALSE,
                    module_color = c("1"="red","0"="lightgrey")) + 
    labs(color = "Signif. genes") + ggtitle(NULL)
})
png(file.path(figure.path,"TODO_name_f4-sup3.png"),width = 300*4*8, height = 300*3, res = 300)
wrap_plots(pl, nrow = 1) + plot_layout(guides = "collect")
dev.off()

# merged_modules = merge_modules_by_overlap(modules,overlap_thresh = 0.5)
# merged_mod_vec = setNames(rep(names(merged_modules), lengths(merged_modules)),
#                           unlist(merged_modules))
# names(merged_modules) = paste0(sample_ls_short[sample_diff],"_",names(merged_modules))
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
png(file.path(figure.path,"TODO_name_f4-sup4.png"),width = 300*10, height = 300*4, res = 300)
wrap_plots(pl,ncol = 2) + plot_layout(guides = "collect")
dev.off()

## Gene Modules on CellUmap ----------
merged_modules_ls = split(names(merged_modules), merged_modules)
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
png(file.path(figure.path,"TODO_name_f4-1.png"),width = 300*4*length(merged_modules_ls), height = 300*8, res = 300)
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
msig_result <- do.call(rbind,lapply(names(merged_modules_ls)[1:5], function(i){
  df = enricher(gene=merged_modules_ls[[i]],
           TERM2GENE = m_df[,c("gs_name","gene_symbol")],
           pAdjustMethod = "none",pvalueCutoff = 0.1)@result
  df$'modules' = i
  df = df %>% filter(p.adjust < 0.2)
  df = df[1:5,]
  df <- df[complete.cases(df), ]
  df
}))
msig_result$GeneRatio = unlist(lapply(msig_result$GeneRatio,function(x) eval(parse(text = x))))
msig_result$logP = -log(msig_result$pvalue,base = 10)
msig_result$Description <- factor(msig_result$Description,
                                 levels = rev(unique(msig_result$Description[order(msig_result$pvalue)])))

msig_result$modules = factor(msig_result$modules,levels = names(merged_modules_ls)[1:5])
levels(msig_result$modules) = gsub("MUT","SmoM2_M",levels(msig_result$modules))
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
ggsave(filename = file.path(figure.path,"f4-2.png"), 
       plot = p, width = 8, height = 5)

# Benchmark w/ other methods -------
# TODO(api): replace with the DiCoLo equivalent of benchmarking_function.R
srat1 = data_S_ls[[1]]
srat2 = data_S_ls[[2]]
res_memento_path = RunMemento(srat1, srat2, input_genes = common_genes,
           backend = TRUE, dir.path = dir.path)
res_memento_path = file.path(dir.path,"Memento","memento_result.csv")
res_memento = read.csv(res_memento_path)
res_DGCA = RunDGCA(srat1, srat2, input_genes = common_genes)
res_lemur = RunLEMUR(srat1, srat2, input_genes = common_genes)
res_milode = RunMiloDE(srat1, srat2, input_genes = common_genes)
res_emdDUFS = E.list_1$vectors

saveRDS(list(res_memento = res_memento,
             res_milode = res_milode,
             res_lemur = res_lemur,
             res_DGCA = res_DGCA), file = file.path(dir.path,"method_res.rds"))

method_ls = c("memento","milode","lemur","DGCA","emdDUFS")
res = lapply(method_ls,function(method){
  df = Generate_rank_table(de_res = get(paste0("res_",method)), 
                           method = method, input_genes = common_genes, 
                           direc = "condition1")
  get_auc(real_score = setNames(df$score,df$gene),
          gt_gene_ls = c(dc_genes,arrested_genes),
          metric = "auprc", plot = TRUE)
})
names(res) = method_ls

## Visualize AUC curve ----------
df = do.call(rbind,lapply(method_ls,function(method){
  df = res[[method]]$curve
  df$'method' = method
  df
})) %>% arrange(method, recall)

ggplot(df, aes(x = recall, y = precision, color = method)) +
  geom_line()
