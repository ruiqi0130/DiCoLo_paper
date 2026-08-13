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
setwd("/data/ruiqi/DiCoLo_paper")
source("./code/simulation_functions.R")
source("./code/benchmarking_functions.R")
# figure.path = "/banach1/ruiqi/bi_gene_graphs/figures"
figure.path = "./figures"

# Load data ------
dir.path = "DiCoLo_data"
# pbmc10k
data.path = file.path(dir.path,"pbmc10k")
sample_ls = c("monocyte1","monocyte2")
data_S_merge = readRDS(file.path(data.path,"pbmc10k.rds"))
data_S_ls = SplitObject(data_S_merge, split.by = "batch")
names(data_S_ls) = sample_ls

# pbmcsca
data.path = file.path(dir.path,"pbmcsca")
sample_ls = c("Drop-seq","inDrops")
data_S_merge = readRDS(file.path(data.path,"pbmcsca.rds"))
data_S_ls = SplitObject(data_S_merge, split.by = "batch")
names(data_S_ls) = sample_ls

# pancreas
data.path = file.path(dir.path,"panc8")
sample_ls = c("human2","human3")
data_S_merge = readRDS(file.path(data.path,"panc8.rds"))
data_S_ls = SplitObject(data_S_merge, split.by = "batch")
names(data_S_ls) = sample_ls


data_S_ls <- lapply(data_S_ls,function(data_S){
  data_S <- data_S %>% NormalizeData() %>% 
    FindVariableFeatures() %>% ScaleData() %>%
    RunPCA(npcs = 50, verbose = FALSE) %>% 
    RunUMAP(dims = 1:10,seed = 42)
  data_S
})

# Set parameters
zinb_prob_range = seq(0.1,0.7,0.1)
ncell_range = c(seq(0.03,0.04,0.01),seq(0.05,0.3,0.05))
ngene_range = seq(5,50,5)
cov_strength_range = seq(0,0.8,0.2)
diff.pct_range = seq(0.9,0.1,-0.1)
diffuse_range = c(FALSE, TRUE)
overlap_range = seq(0.5, 1.0, 0.1)
batch_downsample_range = c(1, 0.75, 0.5, 0.3, 0.1)
para_test_ls = c("ncell","ngene","zinb_prob",
                 "cov_strength","diff.pct",
                 "diffuse","overlap","batch_downsample")

## Default parameters ------
ncell = 0.3 # default 0.1
ngene = 15
zinb_prob = 0.6 # default 0.4
diff.pct = 0.7
cov_strength = 0
mean_quantile = 0.75 # default 0.95
diffuse = FALSE
method_ls = c("milode","DGCA","lemur","DiCoLo","memento")
n_neighbors_all_ls = list(c(1,0),c(0,1))
batch_downsample = 1

## Define tested parameters
for(para_id in 1:length(para_test_ls)){
  for(n_neighbors_ls in n_neighbors_all_ls){
    
para_test = para_test_ls[para_id]
test_id = paste0(n_neighbors_ls,collapse = "vs")
res.path = file.path(data.path,"parameters",para_test,test_id)
if(!dir.exists(res.path)){
  dir.create(res.path,recursive = TRUE)
}

n_reps = 10
seed_grid = make_seed_grid(n_reps,2, base_seed = para_id)
# n_reps = 5
## Simulation params for reproducibility -----
for(rep_id in seq_len(n_reps)){
  rep_path = file.path(res.path,paste0("rep",rep_id))
  if(!dir.exists(rep_path)){
    dir.create(rep_path,recursive = TRUE)
  }
  
  lapply(get(paste0(para_test,"_range")),function(para_value){
    assign(para_test,para_value)
    
    params_ls = lapply(1:2, function(i){
      n_neighbors = n_neighbors_ls[i]
      if(n_neighbors==0){
        return(NULL)
      }else{
        # Override for overlap: need 2 neighborhoods
        if(para_test == "overlap"){
          n_neighbors = 2
          ncell_per_neighbors_ls = rep(floor(ncell * ncol(data_S_ls[[i]])), 2)
          ngene_per_neighbors_ls = rep(ngene, 2)
          overlap_val = para_value
        } else {
          ncell_per_neighbors_ls = floor(ncell * ncol(data_S_ls[[i]]))
          ngene_per_neighbors_ls = ngene
          overlap_val = NULL
        }
        seed = seed_grid %>% filter(sample == i, rep == rep_id) %>%.$seed
        gene_params = Generate_gene_parameters(data_S_ls[i][[1]], data_S_ls[-i][[1]], 
                                               assay = "RNA", layer = "counts",
                                               n_neighbors = n_neighbors,
                                               ngene_per_neighbors_ls = ngene_per_neighbors_ls, 
                                               seed = seed, mean_quantile = mean_quantile)
        return(list(ncell_per_neighbors_ls = ncell_per_neighbors_ls,
                    seed = seed,
                    zinb_prob = zinb_prob,
                    diff.pct = diff.pct,
                    gene_params = gene_params,
                    cov_strength = cov_strength,
                    diffuse = diffuse,
                    overlap_frac = overlap_val,
                    batch_downsample = batch_downsample))
      }
    })
    # save params
    lapply(1:2,function(i){
      saveRDS(params_ls[[i]],file = file.path(rep_path,sprintf("%s_params%d.rds",para_value,i)))
    })
  })
}


## Generate simulated sample list ------
for(rep_id in 1:n_reps){
  cat("rep",rep_id)
  rep_path = file.path(res.path,paste0("rep",rep_id))
  simulation_ls <- lapply(get(paste0(para_test,"_range")), function(para_value){
    params_ls = lapply(1:2,function(i){
      readRDS(file = file.path(rep_path,sprintf("%s_params%d.rds",para_value,i)))
    })
    injected_genes_ls = lapply(params_ls,function(params){
      if(is.null(params)){return(NULL)}
      with(as.list(params),{
        return(setNames(unlist(mapply(function(genes, val) rep(val, length(genes)), 
                                      gene_params$selected_gene_groups, names(gene_params$selected_gene_groups), SIMPLIFY = FALSE)),
                        unlist(gene_params$selected_gene_groups)))
      })
    })
    injected_genes = lapply(1:2,function(i){
      c_g = intersect(names(injected_genes_ls[i][[1]]),names(injected_genes_ls[-i][[1]]))
      a = injected_genes_ls[i][[1]]
      a1 = ifelse(names(a) %in% c_g, paste0("c_",a), paste0(i,"_",a))
      setNames(a1,names(a))
    })
    injected_genes = unlist(injected_genes)
    injected_genes = injected_genes[!duplicated(names(injected_genes))]
    
    # inject signals
    data_S_ls_test = lapply(1:2,function(i){
      params = params_ls[[i]]
      if(is.null(params)){
        return(data_S_ls[[i]])
      }else{
        with(as.list(params), {
          cat(sprintf("  Injecting sample %d, ncell=%s, ngene = %s, dropout = %s, cov_strength = %s\n ", 
                      i, ncell_per_neighbors_ls[1],length(gene_params$'selected_gene_groups'[[1]]),zinb_prob,cov_strength))
          srat_injected <- Inject_local_signals_v2(data_S_ls[[i]], assay = "RNA",
                                                   layer = "counts",
                                                   reduc = "pca", npc = 10, min_knn = max(ncell_per_neighbors_ls), 
                                                   ncell_per_neighbors_ls = ncell_per_neighbors_ls,
                                                   zinb_prob = zinb_prob, diff.pct = diff.pct, 
                                                   gene_params = gene_params, seed = seed, 
                                                   cov_strength = cov_strength,
                                                   diffuse = diffuse,
                                                   overlap_frac = overlap_frac, 
                                                   cap_background = (para_test != "diff.pct"),
                                                   batch_downsample = batch_downsample)
          cat(sprintf("  Done injecting sample %d\n", i))
          srat_injected
        })
      }
    })
    
    common_genes = union(SelectCommonGenes(data_S_ls_test[[1]],data_S_ls_test[[2]],
                                             ngenes = 500),names(injected_genes))
    
    return(list(data_S_ls_test=data_S_ls_test,
                common_genes=common_genes,
                injected_genes=injected_genes,
                params_ls = params_ls))
  })
  names(simulation_ls) = get(paste0(para_test,"_range"))
  
  ## Results for each method -----
  sample_diff = which(n_neighbors_ls!=0)
  ### DiCoLo------
  # Run EMD
  if("DiCoLo" %in% method_ls){
    all_common_genes = lapply(simulation_ls,function(x) sort(unique(x[["common_genes"]])))
    # test_common = lapply(all_common_genes,function(x) setdiff(x,Reduce(intersect, all_common_genes)))
    # flag = all(lengths(test_common) == 0)
    common_genes_overall = Reduce(union, all_common_genes)
    
    emd_paths = lapply(get(paste0(para_test,"_range")),function(para_value){
      simu_obj = simulation_ls[[as.character(para_value)]]
      with(as.list(simu_obj), {
        lapply(1:length(data_S_ls_test),function(i){
          if(is.null(params_ls[[i]])){
            tmp_path = file.path(rep_path,paste0("GeneTrajectory",i))
            common_genes = common_genes_overall
          }else{
            tmp_path = file.path(rep_path,paste0(para_value,"_GeneTrajectory",i))
          }
          if(!dir.exists(tmp_path) & !file.exists(file.path(tmp_path,"emd.csv"))){
            data_S = data_S_ls_test[[i]]
            data_S = data_S %>% NormalizeData() %>%
              FindVariableFeatures() %>%
              ScaleData(verbose = FALSE) %>% 
              RunPCA(npcs = 50, verbose = FALSE)
            tryCatch({
              ComputeGeneEMD(data_S, common_genes, dir.path = tmp_path)
              message("Run on backend\n")
            }, error = function(e){
              message(sprintf("GeneEMD failed @ %s: %s", para_value, conditionMessage(e)))
              dir.create(tmp_path, showWarnings = FALSE)
              write.csv(NULL, file.path(tmp_path, "emd.csv"))  # placeholder
            })
            return(tmp_path)
          }
        })
      })
    })
    emd_paths = unlist(emd_paths)
    
    # Wait for all jobs
    repeat {
      still_running <- !sapply(emd_paths, function(p) file.exists(file.path(p, "emd.csv")))
      if (!any(still_running)) break
      Sys.sleep(10)
    }
    message("All jobs finished. Proceeding...")
    
    # Compute Differential Operator
    res_DiCoLo = lapply(get(paste0(para_test,"_range")),function(para_value){
      simu_obj = simulation_ls[[as.character(para_value)]]
      with(as.list(simu_obj), {
        # Load gene EMD distance
        gene_emd_ls = lapply(1:length(data_S_ls_test),function(i){
          if(is.null(params_ls[[i]])){
            tmp_path = file.path(rep_path,paste0("GeneTrajectory",i))
          }else{
            tmp_path = file.path(rep_path,paste0(para_value,"_GeneTrajectory",i))
          }
          LoadGeneEMD(file.path(tmp_path,""))
        })
        if(any(sapply(gene_emd_ls, is.null))){
          return(NULL)
        }
        # Align gene name
        g = Reduce(intersect, lapply(gene_emd_ls,function(x) rownames(x)))
        gene_emd_ls = lapply(gene_emd_ls, function(x) x[g,g])
        
        gene_graph_ls = lapply(gene_emd_ls, ComputeGraphOperator)
        diff.op = ComputeDifferentialOperator(gene_graph_ls[sample_diff][[1]], gene_graph_ls[-sample_diff][[1]])
        E.list = RunSVD(diff.op, eig_keep = nrow(diff.op))
        return(E.list)
      })
    })
    names(res_DiCoLo) = get(paste0(para_test,"_range"))
    res_DiCoLo = lapply(res_DiCoLo,function(x){
      x$vectors
    })
  }else{
    res_DiCoLo = NULL
  }

  
  ### MiloDE ------
  if("milode" %in% method_ls){
    res_milode = lapply(get(paste0(para_test,"_range")),function(para_value){
      simu_obj = simulation_ls[[as.character(para_value)]]
      with(as.list(simu_obj), {
        res = RunMiloDE(data_S_ls_test[[1]], data_S_ls_test[[2]], input_genes = common_genes,
                        query_id = sample_diff)
        return(res)
      })
    })
    names(res_milode) = get(paste0(para_test,"_range"))
  }else{
    res_milode = NULL
  }
  
  ### LEMUR ------
  if("lemur" %in% method_ls){
    res_lemur = lapply(get(paste0(para_test,"_range")),function(para_value){
      simu_obj = simulation_ls[[as.character(para_value)]]
      with(as.list(simu_obj), {
        res = RunLEMUR(data_S_ls_test[[1]], data_S_ls_test[[2]], input_genes = common_genes)
        return(res)
      })
    })
    names(res_lemur) = get(paste0(para_test,"_range"))
  }else{
    res_lemur = NULL
  }
  
  ### DGCA ------
  if("DGCA" %in% method_ls){
    res_DGCA = lapply(get(paste0(para_test,"_range")),function(para_value){
      simu_obj = simulation_ls[[as.character(para_value)]]
      with(as.list(simu_obj), {
        res = RunDGCA(data_S_ls_test[[1]], data_S_ls_test[[2]], input_genes = common_genes)
        return(res)
      })
    })
    names(res_DGCA) = get(paste0(para_test,"_range"))
  }else{
    res_DGCA = NULL
  }
  
  ### Memento ------
  if("memento" %in% method_ls){
    res_memento = lapply(get(paste0(para_test,"_range")), function(para_value){
      simu_obj = simulation_ls[[as.character(para_value)]]
      with(as.list(simu_obj), {
        tryCatch({
          res = RunMemento(data_S_ls_test[[1]], data_S_ls_test[[2]],
                           input_genes = common_genes,
                           num_boot = 1000L, num_cpus = 4L)
          return(res)
        }, error = function(e){
          message(sprintf("memento failed @ %s=%s: %s", para_test, para_value, conditionMessage(e)))
          return(NULL)
        })
      })
    })
    names(res_memento) = get(paste0(para_test,"_range"))
  }else{
    res_memento = NULL
  }
  
  common_genes_list = lapply(get(paste0(para_test,"_range")), function(para_value){
    simu_obj = simulation_ls[[as.character(para_value)]]
    with(as.list(simu_obj), {
      return(common_genes)
    })
  })
  names(common_genes_list) = get(paste0(para_test,"_range"))
  saveRDS(list(res_DiCoLo = res_DiCoLo,
               res_milode = res_milode,
               res_lemur = res_lemur,
               res_DGCA = res_DGCA,
               res_memento = res_memento,
               common_genes = common_genes_list), 
          file = file.path(rep_path,"method_res.rds"))
  
}

  }# n_neighborhood_id
} # para_id

### Benchmarking result ------
method_ls = c("milode","DGCA","lemur","DiCoLo","memento")

zinb_prob_range = seq(0.1,0.7,0.1)
ncell_range = c(seq(0.03,0.04,0.01),seq(0.05,0.3,0.05))
ngene_range = seq(5,50,5)
cov_strength_range = seq(0,0.8,0.2)
diff.pct_range = seq(0.9,0.1,-0.1)
diffuse_range = c(FALSE, TRUE)
overlap_range = seq(0.5, 1.0, 0.1)
para_test_ls = c("ncell","ngene","zinb_prob",
                 "cov_strength","diff.pct",
                 "diffuse","overlap")
n_neighbors_all_ls = list(c(1,0),c(0,1))

df = do.call(rbind,lapply(para_test_ls,function(para_test){
  lapply(n_neighbors_all_ls,function(n_neighbors_ls){
    test_id = paste0(n_neighbors_ls,collapse = "vs")
    res.path = file.path(data.path,"parameters",para_test,test_id)
    if(!dir.exists(res.path)) return(NULL)
    n_reps = length(list.files(res.path, pattern = "^rep"))
    if(n_reps == 0) return(NULL)
    res = do.call(rbind,lapply(1:n_reps,function(rep_id){
      rep_path = file.path(res.path,paste0("rep",rep_id))
      if(!file.exists(file.path(rep_path,"method_res.rds"))) return(NULL)
      method_res = readRDS(file.path(rep_path,"method_res.rds"))
      if(is.null(method_res$'common_genes')){method_res$'common_genes' = lapply(method_res$'res_DiCoLo',function(x) rownames(x))}
      names(method_res$'common_genes') = get(paste0(para_test,"_range"))
      simulation_ls <- lapply(get(paste0(para_test,"_range")), function(para_value){
        params_ls = lapply(1:2,function(i){
          readRDS(file = file.path(rep_path,sprintf("%s_params%d.rds",para_value,i)))
        })
        injected_genes_ls = lapply(params_ls,function(params){
          if(is.null(params)){return(NULL)}
          with(as.list(params),{
            return(setNames(unlist(mapply(function(genes, val) rep(val, length(genes)), 
                                          gene_params$selected_gene_groups, names(gene_params$selected_gene_groups), SIMPLIFY = FALSE)),
                            unlist(gene_params$selected_gene_groups)))
          })
        })
        injected_genes = lapply(1:2,function(i){
          c_g = intersect(names(injected_genes_ls[i][[1]]),names(injected_genes_ls[-i][[1]]))
          a = injected_genes_ls[i][[1]]
          a1 = ifelse(names(a) %in% c_g, paste0("c_",a), paste0(i,"_",a))
          setNames(a1,names(a))
        })
        injected_genes = unlist(injected_genes)
        injected_genes = injected_genes[!duplicated(names(injected_genes))]
        return(list(injected_genes=injected_genes,
                    params_ls = params_ls))
      })
      names(simulation_ls) = get(paste0(para_test,"_range"))
      
      ## Results for each method
      sample_diff = which(n_neighbors_ls!=0)
      res = lapply(get(paste0(para_test,"_range")),function(para_value){
        simu_obj = simulation_ls[[as.character(para_value)]]
        with(c(as.list(simu_obj), as.list(method_res)), {
          score = do.call(c,lapply(method_ls,function(method){
            if(!as.character(para_value) %in% names(get(paste0("res_",method))) | is.null(get(paste0("res_",method))[[as.character(para_value)]]) ){
              return(NA)
            }
            common_genes = common_genes[[as.character(para_value)]]
            if(is.null(common_genes)){
              common_genes = rownames(res_DiCoLo[[as.character(para_value)]])
            }
            df = Generate_rank_table(de_res = get(paste0("res_",method))[[as.character(para_value)]], 
                                     method = method, input_genes = common_genes, 
                                     direc = paste0("condition",sample_diff))
            get_auc(real_score = setNames(-df$rank,df$gene),
                    gt_gene_ls = names(injected_genes),
                    metric = "auprc", plot = FALSE)
          }))
          data.frame(method = method_ls, score = score)
          
        })
      })
      res = do.call(rbind,res)
      res[,para_test] = rep(get(paste0(para_test,"_range")),each = length(method_ls))
      # if(para_test=="ncell"){
      #   res[,para_test] = floor(ncol(data_S_ls[[sample_diff]]) * res[,para_test])
      # }
      res[,para_test] = as.factor(res[,para_test])
      res$'replicate' = rep_id
      res
    }))
    colnames(res)[colnames(res) == para_test] = "para_grid"
    res[,"para_test"] = para_test
    res[,"direc"] = sample_ls[n_neighbors_ls!=0]
    res
  })
}))
df = do.call(rbind,df)
write.csv(df,file = file.path(data.path,"parameters",para_test,"benchmarking_result.csv"),row.names = FALSE)
df <- read.csv(file.path(data.path,"parameters",para_test,"benchmarking_result.csv"))

# write.csv(df,file = file.path(data.path,"benchmarking_result.csv"),row.names = FALSE)
# df <- read.csv(file.path(data.path,"benchmarking_result.csv"))
# Plot
df$method = factor(df$method, levels = c("DiCoLo","DGCA","lemur","milode"))
p = ggplot(data = df, 
       aes(x = para_grid,
           y = score,color = method
       )) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0.5,1)) +
  labs(
       x = "Background expression fraction",
       # x = "cov strength",
       # x = "UMI Downsampling Proportion",
       y = "Normalized AUPRC") + 
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"))
ggsave(file.path(figure.path,"figS7D.png"), p, width = 7, height = 5)


p = ggplot(data = df %>% filter(method %in% c("memento","DGCA")) %>%
             mutate(para_grid = factor(para_grid, levels = sort(unique(as.numeric(as.character(para_grid))))),
                    para_test = recode(para_test,
                                       ncell = "neighborhood size",
                                       ngene = "number of genes",
                                       zinb_prob = "dropout rate")), 
           aes(x = para_grid,
               y = score,color = method
           )) +
  geom_boxplot() + facet_wrap(~para_test, scales = "free") + 
  labs(
    x = "",
    # x = "cov strength",
    y = "Normalized AUPRC") + 
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    # axis.text.x = element_text(size = 15),
    # axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"))
ggsave(file.path(figure.path,"figS8.png"), p, width = 10, height = 5)

df_summary <- df %>%
  group_by(method, para_grid) %>%
  summarise(
    median_score = median(score, na.rm = TRUE),
    lower = pmax(0,median_score - sd(score, na.rm = TRUE) / sqrt(n())),  # Standard error (SE)
    upper = pmin(1,median_score + sd(score, na.rm = TRUE) / sqrt(n()))
  )
ggplot(data = df_summary, 
       aes(x = para_grid, y = median_score, 
           group = method, color = method,
       )) +
  geom_line() +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.4,linewidth = 0.8) + 
  labs(x = para_test, 
       y = "Normalized AUPRC", 
       color = "Method") + 
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"))
