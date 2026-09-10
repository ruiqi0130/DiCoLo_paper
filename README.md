# DiCoLo — code for the paper

Code accompanying *Integration-free and cluster-free detection of localized
differential gene co-expression in single-cell data with DiCoLo*
(Li, Yang, Su, Jaffe, Lindenbaum & Kluger).

The DiCoLo method itself is an R package in a separate repository:
<https://github.com/KlugerLab/DiCoLo>. This repository contains only the
analysis and figure code for the paper.

---

## Quick start

```bash
git clone https://github.com/ruiqi0130/DiCoLo_paper.git
cd DiCoLo_paper

# Download and unpack the data archive (see "Data" below)
#   -> creates ./DiCoLo_data

export DICOLO_PYTHON=$(which python3)     # needs `memento` and `pot`
Rscript -e 'rmarkdown::render("DiCoLo_figure_code.Rmd")'
```

All paths live in `config.R` and are overridable with environment variables;
no script contains a machine-specific path. Run everything from the repository
root.

| Variable | Default | Meaning |
| --- | --- | --- |
| `DICOLO_ROOT` | auto-detected | Repository root |
| `DICOLO_DATA` | `<root>/DiCoLo_data` | Unpacked Zenodo archive |
| `DICOLO_FIGS` | `<root>/figures` | Where figures and tables are written |
| `DICOLO_PYTHON` | `which python3` | Interpreter for Memento and the OT backend |
| `DICOLO_OT_SCRIPTS` | package default | `script_dir` passed to `ComputeGeneEMD()` |

---

## Repository layout

```
config.R                      paths, python backend, helper loading
DiCoLo_figure_code.Rmd        Fig 2, 3, 4, S4
R/
  simulation_functions.R      signal injection, covariance model
  benchmarking_functions.R    method wrappers, ranking, AUPRC/AUROC
  plotting_functions.R        shared plotting helpers
  module_derivation.R         DiCoLo module-derivation workflow (Fig S2/S3/S6)
scripts/
  run_benchmark_on_simulations.R   Fig 3, Fig S8A-E, Fig S9
  aggregate_benchmark_results.R    per-dataset CSV consumed by Fig 3
  run_benchmark_on_realdata.R      Fig 4E
  run_concordance_analysis.R       Supplemental Tables 3, 4
  run_modules_smom2.R              Fig S2, S3, Supplemental Table 2
  run_modules_wls.R                Fig S5, S6
  run_sensitivity_test.R           Fig S7
  run_covariance_diagnosis.R       Fig S8F, S8G
  run_null_experiment.R            Fig S10
  run_cellprop_test.R              Fig S11
  run_graphical_abstract.R         graphical abstract
```

## Intermediate Data

`DiCoLo_data/` is distributed via Zenodo (DOI: `<DOI>`) and is not tracked in
git. Expected layout:

```
DiCoLo_data/
  pbmc10k/   pbmc10k.rds,  benchmarking_result.csv,  null_experiments/, simulation_parameters/…
  pbmcsca/   pbmcsca.rds,  benchmarking_result.csv,  null_experiments/…
  panc8/     panc8.rds,    benchmarking_result.csv,  null_experiments/…
  smom2/     {MUT,CTL}_dicolo_modules.rds
  					 benchmarking_result.csv
             RunTime
             Sennett_gene_list.xlsx
             parameters/GeneTrajectory_{MUT,CTL}/emd.csv
             parameters/method_res.rds
             cellprop_test/stat_table.csv
  wlsko/     CTL_dicolo_modules.rds
             parameters/GeneTrajectory_{MUT,CTL}/emd.csv
  sensitivity_test/   jaccard_result_{ndim_knn,dm_vs_pc}.csv
```

Preprocessed mouse dermal Seurat objects are also on figshare
([SmoM2 mutant](https://doi.org/10.6084/m9.figshare.26507098), [wlsko mutant](https://doi.org/10.6084/m9.figshare.25243225)).
The Sennett *et al.* (2015) DC signature is Table S2 of that paper
(`1-s2.0-S153458071500430X-mmc2.xlsx`), saved as `Sennett_gene_list.xlsx`.

---

## Dependencies

The full `sessionInfo()` from the machine used for the reported analyses:

```r
> sessionInfo()
R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] DGCA_1.0.3             lemur_1.3.6            miloDE_0.0.0.9000     
 [4] reticulate_1.37.0      DiCoLo_1.0.0           ggrepel_0.9.5         
 [7] AUCell_1.28.0          tidyr_1.3.1            harmony_1.2.0         
[10] Rcpp_1.0.12            SeuratWrappers_0.4.0   org.Hs.eg.db_3.20.0   
[13] AnnotationDbi_1.68.0   IRanges_2.40.1         S4Vectors_0.44.0      
[16] Biobase_2.66.0         BiocGenerics_0.52.0    clusterProfiler_4.14.6
[19] pheatmap_1.0.12        RColorBrewer_1.1-3     patchwork_1.2.0       
[22] ggplotify_0.1.2        ggplot2_3.5.1          CSCORE_1.1.0          
[25] dplyr_1.1.4            Seurat_5.1.0           SeuratObject_5.0.2    
[28] sp_2.1-4              

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2           dichromat_2.0-0.1          
  [3] GSEABase_1.68.0             rARPACK_0.11-0             
  [5] nnet_7.3-19                 locfdr_1.1-8               
  [7] goftest_1.2-3               Biostrings_2.74.1          
  [9] vctrs_0.6.5                 ggtangle_0.1.2.001         
 [11] spatstat.random_3.2-3       digest_0.6.35              
 [13] png_0.1-8                   Augur_1.0.3                
 [15] deldir_2.0-4                parallelly_1.37.1          
 [17] MASS_7.3-60.2               reshape2_1.4.4             
 [19] foreach_1.5.2               httpuv_1.6.15              
 [21] qvalue_2.38.0               withr_3.0.0                
 [23] xfun_0.44                   ggfun_0.2.1                
 [25] ggpubr_0.6.0                survival_3.5-8             
 [27] memoise_2.0.1               ggbeeswarm_0.7.2           
 [29] parsnip_1.2.1               gson_0.1.0                 
 [31] gtools_3.9.5                tidytree_0.4.6             
 [33] zoo_1.8-12                  pbapply_1.7-2              
 [35] R.oo_1.26.0                 Formula_1.2-5              
 [37] KEGGREST_1.46.0             promises_1.3.0             
 [39] httr_1.4.7                  rstatix_0.7.2              
 [41] globals_0.16.3              fitdistrplus_1.1-11        
 [43] rstudioapi_0.16.0           UCSC.utils_1.2.0           
 [45] miniUI_0.1.1.1              generics_0.1.3             
 [47] DOSE_4.0.1                  base64enc_0.1-3            
 [49] zlibbioc_1.52.0             ScaledMatrix_1.14.0        
 [51] ggraph_2.2.1                polyclip_1.10-6            
 [53] randomForest_4.7-1.1        GenomeInfoDbData_1.2.13    
 [55] SparseArray_1.6.2           doParallel_1.0.17          
 [57] pracma_2.4.4                xtable_1.8-4               
 [59] stringr_1.5.1               evaluate_0.24.0            
 [61] S4Arrays_1.6.0              preprocessCore_1.68.0      
 [63] GenomicRanges_1.58.0        irlba_2.3.5.1              
 [65] colorspace_2.1-0            ROCR_1.0-11                
 [67] readxl_1.4.3                spatstat.data_3.0-4        
 [69] magrittr_2.0.3              lmtest_0.9-40              
 [71] glmGamPoi_1.14.3            viridis_0.6.5              
 [73] later_1.3.2                 ggtree_3.14.0              
 [75] lattice_0.22-6              mapproj_1.2.11             
 [77] spatstat.geom_3.2-9         future.apply_1.11.2        
 [79] scattermore_1.2             XML_3.99-0.16.1            
 [81] scuttle_1.16.0              cowplot_1.1.3              
 [83] matrixStats_1.3.0           RcppAnnoy_0.0.22           
 [85] Hmisc_5.1-3                 class_7.3-22               
 [87] pillar_1.9.0                nlme_3.1-164               
 [89] iterators_1.0.14            compiler_4.4.0             
 [91] beachmat_2.22.0             RSpectra_0.16-1            
 [93] stringi_1.8.4               gower_1.0.1                
 [95] tensor_1.5                  SummarizedExperiment_1.36.0
 [97] lubridate_1.9.3             plyr_1.8.9                 
 [99] crayon_1.5.2                abind_1.4-5                
[101] gridGraphics_0.5-1          locfit_1.5-9.9             
[103] pals_1.8                    graphlayouts_1.1.1         
[105] bit_4.0.5                   RcppGreedySetCover_0.1.0   
[107] fastmatch_1.1-4             fastcluster_1.2.6          
[109] codetools_0.2-20            recipes_1.0.10             
[111] BiocSingular_1.22.0         plotly_4.10.4              
[113] mime_0.12                   rsample_1.2.1              
[115] splines_4.4.0               fastDummies_1.7.3          
[117] sparseMatrixStats_1.14.0    cellranger_1.1.0           
[119] knitr_1.47                  blob_1.2.4                 
[121] utf8_1.2.4                  fs_1.6.4                   
[123] checkmate_2.3.1             listenv_0.9.1              
[125] DelayedMatrixStats_1.24.0   ggsignif_0.6.4             
[127] tibble_3.2.1                Matrix_1.6-4               
[129] statmod_1.5.0               tweenr_2.0.3               
[131] pkgconfig_2.0.3             tools_4.4.0                
[133] cachem_1.1.0                RSQLite_2.3.7              
[135] numDeriv_2016.8-1.1         viridisLite_0.4.2          
[137] DBI_1.2.3                   impute_1.80.0              
[139] fastmap_1.2.0               rmarkdown_2.27             
[141] scales_1.3.0                grid_4.4.0                 
[143] pbmcapply_1.5.1             ica_1.0-3                  
[145] broom_1.0.6                 FNN_1.1.4                  
[147] BiocManager_1.30.23         dotCall64_1.1-1            
[149] carData_3.0-5               graph_1.84.1               
[151] RANN_2.6.1                  rpart_4.1.23               
[153] farver_2.1.2                tidygraph_1.3.1            
[155] yaml_2.3.8                  foreign_0.8-86             
[157] MatrixGenerics_1.14.0       cli_3.6.2                  
[159] purrr_1.0.2                 tester_0.2.0               
[161] leiden_0.4.3.1              lifecycle_1.0.4            
[163] uwot_0.2.2                  lava_1.8.0                 
[165] backports_1.5.0             BiocParallel_1.36.0        
[167] annotate_1.84.0             timechange_0.3.0           
[169] gtable_0.3.5                ggridges_0.5.6             
[171] yardstick_1.3.1             progressr_0.14.0           
[173] parallel_4.4.0              ape_5.8                    
[175] limma_3.62.2                jsonlite_1.8.8             
[177] edgeR_4.4.2                 RcppHNSW_0.6.0             
[179] bit64_4.0.5                 Rtsne_0.17                 
[181] yulab.utils_0.2.5           BiocNeighbors_1.20.2       
[183] spatstat.utils_3.0-4        GOSemSim_2.32.0            
[185] R.utils_2.12.3              timeDate_4032.109          
[187] lazyeval_0.2.2              shiny_1.8.1.1              
[189] dynamicTreeCut_1.63-1       htmltools_0.5.8.1          
[191] enrichplot_1.26.6           GO.db_3.20.0               
[193] sctransform_0.4.1           rappdirs_0.3.3             
[195] glue_1.7.0                  spam_2.10-0                
[197] XVector_0.46.0              treeio_1.37.0.001          
[199] GeneTrajectory_1.0.0        gridExtra_2.3              
[201] igraph_2.0.3                R6_2.5.1                   
[203] SingleCellExperiment_1.28.1 labeling_0.4.3             
[205] cluster_2.1.6               aplot_0.2.2                
[207] GenomeInfoDb_1.42.3         ipred_0.9-14               
[209] WGCNA_1.72-5                DelayedArray_0.32.0        
[211] tidyselect_1.2.1            vipor_0.4.7                
[213] htmlTable_2.4.2             maps_3.4.2                 
[215] ggforce_0.4.2               car_3.1-2                  
[217] future_1.33.2               rsvd_1.0.5                 
[219] munsell_0.5.1               KernSmooth_2.23-22         
[221] furrr_0.3.1                 miloR_2.2.0                
[223] PRROC_1.3.1                 data.table_1.15.4          
[225] htmlwidgets_1.6.4           fgsea_1.32.4               
[227] rlang_1.1.4                 spatstat.sparse_3.0-3      
[229] spatstat.explore_3.2-7      remotes_2.5.0              
[231] fansi_1.0.6                 hardhat_1.4.0              
[233] beeswarm_0.4.0              prodlim_2023.08.28      
```

R (≥ 4.4) packages from GitHub, pinned by commit:

```r
# Commit hashes used for the reported analyses are given in parentheses.
remotes::install_github("KlugerLab/DiCoLo")
remotes::install_github("KlugerLab/GeneTrajectory")
remotes::install_github("KlugerLab/LocalizedMarkerDetector")
remotes::install_github("MarioniLab/miloDE")
remotes::install_github("const-ae/lemur")
devtools::install_github("andymckenzie/DGCA")
```

From CRAN/Bioconductor: `Seurat`, `SingleCellExperiment`, `dplyr`, `tibble`,
`tidyr`, `igraph`, `pdist`, `reshape2`, `ggplot2`, `ggplotify`, `patchwork`,
`viridis`, `scales`, `RColorBrewer`, `plotly`, `pROC`, `PRROC`, `locfdr`,
`readxl`, `msigdbr` (v7.5.1), `clusterProfiler` (v4.14.6), `ReactomePA`,
`org.Mm.eg.db`, `scuttle`, `SeuratWrappers`, `reticulate`.

Python (via `reticulate`, selected by `DICOLO_PYTHON`):
`memento`, `pot`, `numpy`, `scipy`, `anndata` .

---

## Code and data availability

> All custom code is available at <https://github.com/ruiqi0130/DiCoLo_paper> and is provided as Supplemental
> Code. The DiCoLo R package is available at
> <https://github.com/KlugerLab/DiCoLo>. Processed data required to reproduce all figures are deposited at Zenodo (DOI: `<data DOI>`).
