# =============================================================================
# run_concordance_analysis.R   (Supplemental Tables 3 and 4)
#
# Concordance between the DiCoLo-defined SmoM2 modules and the outputs of the
# competing methods on the same dataset. Note that the modules are defined by
# DiCoLo, so this describes differences in OUTPUT STRUCTURE, not method
# accuracy.
#
#   Supplemental Table 3 - DGCA: significant gene pairs classified as
#       within_module / cross_module / one_in_module / neither.
#   Supplemental Table 4 - LEMUR and miloDE: how many module genes are
#       recovered as significant DE genes, and what fraction of each method's
#       total significant list those genes represent.
#
# Requires <MUT|CTL>_dicolo_modules.rds, produced by run_modules_smom2.R and
# also distributed in DiCoLo_data/ (see README - "Derived data provenance").
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
dir.path    <- DATA_DIR
figure.path <- FIGURE_DIR
data.path <- file.path(DATA_DIR, "smom2")
dir.path.data <- data.path

# Load data ------
sample_ls = paste0("smom2_dermal_",paste0("E13.5_",c("MUT","CTL")))
data_S_ls = lapply(sample_ls,function(sample_name){
  readRDS(file.path(dir.path.data,sprintf("data_S_%s.rds",sample_name)))
})
names(data_S_ls) = sample_ls
                       
# Load results -----
common_genes = SelectCommonGenes(data_S_ls[[1]],data_S_ls[[2]],
                                 ngenes = 500)
gene.dist.mat_ls = lapply(sample_ls,function(sample_name){
  LoadGeneEMD(file.path(data.path,"parameters",sprintf("GeneTrajectory_%s",sample_name),""))
})
gene_modules = readRDS(file.path(data.path,"MUT_dicolo_modules.rds"))
module_gene_list <- split(names(gene_modules), gene_modules)

### Benchmarking ------
res_DGCA <- RunDGCA(data_S_ls[[1]], data_S_ls[[2]], input_genes = common_genes)
res_lemur <- RunLEMUR(data_S_ls[[1]], data_S_ls[[2]], input_genes = common_genes)
res_milode <- RunMiloDE(data_S_ls[[1]], data_S_ls[[2]], input_genes = common_genes)
res_ls = list(res_milode = res_milode,
              res_lemur = res_lemur,
              res_DGCA = res_DGCA)
saveRDS(res_ls, file = file.path(data.path,"parameters","method_res_concordance.rds"))

method_ls = c("milode","lemur", "DGCA")
rank_df = do.call(rbind, lapply(method_ls, function(method){
  df = Generate_rank_table(de_res = res_ls[[paste0("res_",method)]], 
                           method = method, input_genes = common_genes, 
                           direc = "condition1")
  df$'method' = method
  df
}))
pl = lapply(names(module_gene_list),function(modules){
  pl = lapply(method_ls,function(method){
    plot_module_enrichment(rank_df, module_genes = module_gene_list[[modules]],
                           method_name = method)
  })
})
pl = unlist(pl, recursive = FALSE)

# Test DGCA
direc = "condition1"
sig_pairs <- res_DGCA %>%
  filter(direction == direc) %>%
  filter(Classes != "NonSig") %>%
  mutate(
    condition1 = unlist(lapply(strsplit(as.character(Classes), "/"), function(x) x[[1]])),
    condition2 = unlist(lapply(strsplit(as.character(Classes), "/"), function(x) x[[2]]))
  ) %>% filter(.data[[direc]] == "+") %>% 
  filter(pValDiff_adj < 0.05)
table(sig_pairs$pair_type)
sig_pairs <- sig_pairs %>%
  mutate(
    mod1 = gene_modules[Gene1],
    mod2 = gene_modules[Gene2],
    pair_type = case_when(
      is.na(mod1) & is.na(mod2) ~ "neither",
      is.na(mod1) | is.na(mod2) ~ "one_in_module",
      mod1 == mod2 ~ "within_module",
      TRUE ~ "cross_module"))
sig_pairs %>%
  filter(pair_type == "within_module") %>%
  mutate(module_pair = paste(
    pmin(as.character(mod1), as.character(mod2), na.rm = T), 
    pmax(as.character(mod1), as.character(mod2), na.rm = T), 
    sep = " × "
  )) %>% count(module_pair)
sig_pairs %>%
  filter(pair_type == "cross_module") %>%
  mutate(module_pair = paste(
    pmin(as.character(mod1), as.character(mod2)), 
    pmax(as.character(mod1), as.character(mod2)), 
    sep = " × "
  )) %>% count(module_pair)
sig_pairs %>%
  filter(pair_type == "one_in_module") %>%
  mutate(module_pair = paste(
    pmin(as.character(mod1), as.character(mod2), na.rm = T), 
    pmax(as.character(mod1), as.character(mod2), na.rm = T), 
    sep = " × "
  )) %>% count(module_pair)

all_module_genes <- names(gene_modules)
module_recovery <- list()

# Test lemur
sig_genes = res_lemur %>% mutate(padj = adj_pval) %>% 
  filter(direction == direc) %>% filter(padj < 0.05) %>%.$name
n_module <- length(all_module_genes)
n_recovered <- sum(all_module_genes %in% sig_genes)
module_recovery$lemur <- data.frame(
  method = "LEMUR",
  total_sig = length(sig_genes),
  module_recovered = sprintf("%d/%d (%0.1f%%)", n_recovered, n_module,
                             100 * n_recovered / n_module),
  module_frac_of_sig = round(100 * n_recovered / length(sig_genes), 1),
  not_in_module = length(sig_genes) - n_recovered
)

# Test milode
sig_genes = res_milode %>% mutate(padj = pval_corrected_across_genes) %>% 
  filter(direction == direc) %>% filter(padj < 0.05) %>% 
  arrange(.,padj,desc(abs(logFC))) %>% distinct(.,gene, .keep_all = TRUE) %>%.$gene
n_module <- length(all_module_genes)
n_recovered <- sum(all_module_genes %in% sig_genes)
module_recovery$milode <- data.frame(
  method = "miloDE",
  total_sig = length(sig_genes),
  module_recovered = sprintf("%d/%d (%0.1f%%)", n_recovered, n_module,
                             100 * n_recovered / n_module),
  module_frac_of_sig = round(100 * n_recovered / length(sig_genes), 1),
  not_in_module = length(sig_genes) - n_recovered
)
module_recovery <- do.call(rbind, module_recovery)

# =============================================================================
# Write the supplemental tables
# =============================================================================
table.path <- file.path(FIGURE_DIR, "tables")
if (!dir.exists(table.path)) dir.create(table.path, recursive = TRUE)

write.csv(sig_pairs %>% count(pair_type),
          file.path(table.path, "supplemental_table3_dgca_pair_classes.csv"),
          row.names = FALSE)

write.csv(module_recovery,
          file.path(table.path, "supplemental_table4_module_recovery.csv"),
          row.names = FALSE)

message("Wrote Supplemental Tables 3 and 4 to ", table.path)
