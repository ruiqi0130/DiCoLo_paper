# =============================================================================
# aggregate_benchmark_results.R
#
# run_benchmark_on_simulations.R writes one CSV per stress-test parameter:
#
#   DiCoLo_data/<dataset>/parameters/<para_test>/benchmarking_result.csv
#
# Figure 3 in DiCoLo_figure_code.Rmd reads a single concatenated table per
# dataset:
#
#   DiCoLo_data/<dataset>/benchmarking_result.csv
#
# This script performs that concatenation. It is the missing link between the
# two paths: the per-dataset file shipped in the Zenodo DiCoLo_data archive is
# exactly this script's output, so Figure 3 can be reproduced without re-running
# the benchmark. Running the script over freshly computed per-para_test CSVs
# regenerates it.
#
# Usage (from the repository root):
#   Rscript scripts/aggregate_benchmark_results.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

source("config.R")

DATASETS <- c("pbmc10k", "pbmcsca", "panc8")

# Only the three parameters that appear in main Figure 3. The remaining sweeps
# (diff.pct, diffuse, overlap, batch_downsample, cov_strength) feed Fig S8 and
# are plotted directly from their per-para_test CSVs.
MAIN_PARA_TESTS <- c("ncell", "ngene", "zinb_prob")

for (dataset in DATASETS) {
  param_root <- file.path(DATA_DIR, dataset, "parameters")
  if (!dir.exists(param_root)) {
    message(sprintf("[%s] skipped, no parameters/ directory", dataset)); next
  }

  parts <- lapply(MAIN_PARA_TESTS, function(para_test) {
    csv <- file.path(param_root, para_test, "benchmarking_result.csv")
    if (!file.exists(csv)) {
      message(sprintf("[%s] missing %s", dataset, csv)); return(NULL)
    }
    read.csv(csv)
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]

  if (length(parts) == 0) {
    message(sprintf("[%s] nothing to aggregate", dataset)); next
  }

  df <- bind_rows(parts)
  out <- file.path(DATA_DIR, dataset, "benchmarking_result.csv")
  write.csv(df, file = out, row.names = FALSE)
  message(sprintf("[%s] wrote %s (%d rows, para_test: %s)",
                  dataset, out, nrow(df),
                  paste(sort(unique(df$para_test)), collapse = ", ")))
}
