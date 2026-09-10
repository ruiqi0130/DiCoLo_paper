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
DiCoLo_figure_code.Rmd        Fig 2, 3, 4A-D, S4
R/
  simulation_functions.R      signal injection, covariance model
  benchmarking_functions.R    method wrappers, ranking, AUPRC/AUROC
  plotting_functions.R        shared plotting helpers
scripts/
  run_benchmark_on_simulations.R   Fig 3 (data), Fig S8A-E, Fig S9
  aggregate_benchmark_results.R    per-dataset CSV consumed by Fig 3
  run_benchmark_on_realdata.R      Fig 4E
  run_concordance_analysis.R       Supplemental Tables 3, 4
  run_modules_smom2.R              Fig S2, S3, Supplemental Table 2  [INCOMPLETE]
  run_modules_wls.R                Fig S5, S6                        [INCOMPLETE]
  run_sensitivity_test.R           Fig S7
  run_covariance_diagnosis.R       Fig S8F, S8G
  run_null_experiment.R            Fig S10
  run_cellprop_test.R              Fig S11
  run_graphical_abstract.R         graphical abstract
```

> **`run_modules_smom2.R` and `run_modules_wls.R` are marked `[INCOMPLETE]`.**
> Both still call helper functions from an earlier internal codebase rather than
> their DiCoLo equivalents, and are not runnable as written. Every affected line
> is tagged `TODO` in the file. See "Known gaps".

---

## Figure map

Figure numbers below are the **manuscript** numbers, and output filenames match
them.

| Manuscript | Produced by | Output |
| --- | --- | --- |
| Fig 1 | schematic, no code | — |
| Fig 2A, 2B, 2C | `DiCoLo_figure_code.Rmd` | inline |
| Fig 3 | `run_benchmark_on_simulations.R` → `aggregate_benchmark_results.R` → Rmd | inline |
| Fig 4A | `DiCoLo_figure_code.Rmd` | inline |
| Fig 4B | `DiCoLo_figure_code.Rmd` | `fig4B_graph_operator_heatmap.png` |
| Fig 4C | `DiCoLo_figure_code.Rmd` | inline |
| Fig 4D | `DiCoLo_figure_code.Rmd` | `fig4D_pathway_enrichment.png` |
| Fig 4E | `scripts/run_benchmark_on_realdata.R` | `fig4E_signature_recovery.png` |
| Fig S1 | Rmd Fig 2A chunk, re-run per dataset — see "Known gaps" | inline |
| Fig S2 A–D | `scripts/run_modules_smom2.R` | *(pending, see "Known gaps")* |
| Fig S3 | `scripts/run_modules_smom2.R` | *(pending)* |
| Fig S4 | `DiCoLo_figure_code.Rmd` | `figS4_morans_I.png` |
| Fig S5 | `scripts/run_modules_wls.R` | *(pending)* |
| Fig S6 | `scripts/run_modules_wls.R` | *(pending)* |
| Fig S7 | `scripts/run_sensitivity_test.R` | `figS7_sensitivity.png` |
| Fig S8 A | `run_benchmark_on_simulations.R` (`para_test = diff.pct`) | `figS8A_diff.pct.png` |
| Fig S8 B | `run_benchmark_on_simulations.R` (`diffuse`) | `figS8B_diffuse.png` |
| Fig S8 C | `run_benchmark_on_simulations.R` (`overlap`) | `figS8C_overlap.png` |
| Fig S8 D | `run_benchmark_on_simulations.R` (`batch_downsample`) | `figS8D_batch_downsample.png` |
| Fig S8 E | `run_benchmark_on_simulations.R` (`cov_strength`) | `figS8E_cov_strength.png` |
| Fig S8 F | `scripts/run_covariance_diagnosis.R` | `figS8F_covariance_violin.png` |
| Fig S8 G | `scripts/run_covariance_diagnosis.R` | `figS8G_covariance_heatmap.png` |
| Fig S9 | `run_benchmark_on_simulations.R` | `figS9_memento_dgca.png` |
| Fig S10 | `scripts/run_null_experiment.R` | `figS10_null_and_negative_control.png` |
| Fig S11 | `scripts/run_cellprop_test.R` | `figS11_celltype_abundance.png` |
| Suppl. Table 1 | dataset summary statistics, no code | — |
| Suppl. Table 2 | `run_modules_smom2.R` (`*_dicolo_modules.rds`) | *(pending)* |
| Suppl. Table 3 | `scripts/run_concordance_analysis.R` | `tables/supplemental_table3_dgca_pair_classes.csv` |
| Suppl. Table 4 | `scripts/run_concordance_analysis.R` | `tables/supplemental_table4_module_recovery.csv` |

---

## Data

`DiCoLo_data/` is distributed via Zenodo (DOI: **TODO**) and is not tracked in
git. Expected layout:

```
DiCoLo_data/
  pbmc10k/   pbmc10k.rds,  benchmarking_result.csv,  parameters/…
  pbmcsca/   pbmcsca.rds,  benchmarking_result.csv,  parameters/…
  panc8/     panc8.rds,    benchmarking_result.csv,  parameters/…
  smom2/     data_S_smom2_dermal_E13.5_{MUT,CTL}.rds
             {MUT,CTL}_dicolo_modules.rds
             parameters/GeneTrajectory_{MUT,CTL}/emd.csv
  wlsko/     …
  Sennett_gene_list.xlsx
```

Preprocessed mouse dermal Seurat objects are also on figshare
([SmoM2 mutant](https://figshare.com/ndownloader/files/55008218),
[wild-type](https://figshare.com/ndownloader/files/55008224)).
The Sennett *et al.* (2015) DC signature is Table S2 of that paper
(`1-s2.0-S153458071500430X-mmc2.xlsx`), saved as `Sennett_gene_list.xlsx`.

### Derived data provenance

The archive ships **computed intermediates alongside raw input**, so that
figures can be regenerated in minutes rather than re-running the full
optimal-transport pipeline. Some files are therefore read from a different
location than the script that produced them writes to; the table below gives
the mapping.

| File in `DiCoLo_data/` | Written by | Read by |
| --- | --- | --- |
| `<dataset>/parameters/<para_test>/rep*/params*.rds` | `run_benchmark_on_simulations.R` | same script on resume; `run_covariance_diagnosis.R` |
| `<dataset>/parameters/<para_test>/rep*/method_res.rds` | `run_benchmark_on_simulations.R` | same script (aggregation stage) |
| `<dataset>/parameters/<para_test>/benchmarking_result.csv` | `run_benchmark_on_simulations.R` | `aggregate_benchmark_results.R`; Fig S8 panels |
| `<dataset>/benchmarking_result.csv` | `aggregate_benchmark_results.R` | Rmd, Fig 3 |
| `smom2/<MUT\|CTL>_dicolo_modules.rds` | `run_modules_smom2.R` | Rmd Fig 4B/4C/4D and S4; `run_concordance_analysis.R` |
| `wlsko/<CTL>_dicolo_modules.rds` | `run_modules_wls.R` | `run_modules_wls.R` (plotting stage) |
| `*/parameters/GeneTrajectory_*/emd.csv` | `DiCoLo::ComputeGeneEMD()` | every script that builds a graph operator |

Every script checks whether its intermediates already exist and skips
recomputation if so, which is what makes the two modes below possible.

---

## Reproducing

### Mode 1 — from archived intermediates (minutes)

Unpack `DiCoLo_data/` as distributed and run the figure code. All gene–gene OT
distances, simulation parameters, per-method results and module assignments are
already present, so nothing expensive is recomputed.

```bash
Rscript -e 'rmarkdown::render("DiCoLo_figure_code.Rmd")'   # Fig 2, 3, 4A-D, S4
for s in run_benchmark_on_realdata run_concordance_analysis \
         run_sensitivity_test run_covariance_diagnosis \
         run_null_experiment run_cellprop_test; do
  Rscript scripts/$s.R
done
Rscript scripts/run_benchmark_on_simulations.R    # plotting stage only
```

### Mode 2 — full recomputation (days)

Delete the `parameters/` subdirectories and the derived files listed above,
then run in this order:

1. `scripts/run_modules_smom2.R`, `scripts/run_modules_wls.R` — module derivation
2. `scripts/run_benchmark_on_simulations.R` — per-`para_test` sweeps, all datasets
3. `scripts/aggregate_benchmark_results.R`
4. the remaining `scripts/run_*.R`
5. `DiCoLo_figure_code.Rmd`

Approximate cost: gene–gene OT distances dominate (hours per dataset per
condition, parallelised over genes); the simulation benchmark is 10 replicates
× 8 parameter sweeps × 3 datasets; Memento adds roughly 10 minutes per
simulation replicate, which is why Fig S9 is restricted to pbmc10k with 5
replicates.

---

## Dependencies

**TODO — paste `sessionInfo()` from the machine used for the final run.**

R (≥ 4.4) packages from GitHub, pinned by commit:

```r
remotes::install_github("KlugerLab/DiCoLo")            # TODO: pin tag/commit
remotes::install_github("KlugerLab/GeneTrajectory")    # TODO: pin
remotes::install_github("KlugerLab/LocalizedMarkerDetector")  # TODO: pin
remotes::install_github("MarioniLab/miloDE")           # TODO: pin
remotes::install_github("const-ae/lemur")              # TODO: pin
devtools::install_github("andymckenzie/DGCA")          # TODO: pin
```

From CRAN/Bioconductor: `Seurat`, `SingleCellExperiment`, `dplyr`, `tibble`,
`tidyr`, `igraph`, `pdist`, `reshape2`, `ggplot2`, `ggplotify`, `patchwork`,
`viridis`, `scales`, `RColorBrewer`, `plotly`, `pROC`, `PRROC`, `locfdr`,
`readxl`, `msigdbr` (v7.5.1), `clusterProfiler` (v4.14.6), `ReactomePA`,
`org.Mm.eg.db`, `scuttle`, `SeuratWrappers`, `reticulate`.

Python (via `reticulate`, selected by `DICOLO_PYTHON`):
`memento`, `pot`, `numpy`, `scipy`, `anndata`. **TODO — pin versions.**

---

## Known gaps

1. **`run_modules_smom2.R` / `run_modules_wls.R` are pending migration to the
   DiCoLo package API** and are not runnable as written; every affected line is
   tagged `TODO`. They are the sole source of Fig S2, S3, S5, S6, Supplemental
   Table 2 and `*_dicolo_modules.rds` — which Fig 4B/4C/4D, Fig S4 and
   `run_concordance_analysis.R` all consume.
2. **Fig S1 has no dedicated script.** It is the merged-batch UMAP from the
   Fig 2A chunk applied to pbmc10k, pancreas and pbmcsca. Either add a small
   loop over the three datasets or note in the caption that it reuses that code.
3. **Supplemental Table 1** (per-dataset cell counts, median nUMI) is reported
   in the text with no script behind it. Worth a five-line summary script.

---

## Code and data availability

> All custom code is available at <https://github.com/ruiqi0130/DiCoLo_paper>
> (archived at Zenodo, DOI: TODO) and is provided as Supplemental Code. The
> DiCoLo R package is available at <https://github.com/KlugerLab/DiCoLo>
> (vTODO, archived at Zenodo, DOI: TODO). Processed data required to reproduce
> all figures are deposited at Zenodo (DOI: TODO).

## Citation

TODO
