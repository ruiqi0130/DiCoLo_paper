# =============================================================================
# config.R - central paths and environment settings
#
# Every script in scripts/ and the figure notebook source this file. Nothing
# here is machine-specific: override any value with an environment variable
# before starting R, e.g.
#
#   export DICOLO_DATA=/scratch/me/DiCoLo_data
#   export DICOLO_PYTHON=$(which python3)
#
# Source it from the repository root:  source("config.R")
# =============================================================================

# --- Repository root ---------------------------------------------------------
# Resolved once, so scripts can be launched from anywhere (Rscript, RStudio,
# knitr) without setwd().
PROJECT_ROOT <- local({
  p <- Sys.getenv("DICOLO_ROOT", "")
  if (nzchar(p)) return(normalizePath(p, mustWork = TRUE))
  # walk up from the working directory until config.R is found
  d <- normalizePath(getwd())
  for (i in 1:5) {
    if (file.exists(file.path(d, "config.R"))) return(d)
    d <- dirname(d)
  }
  stop("Could not locate the repository root. Set DICOLO_ROOT.")
})

# --- Data / output locations -------------------------------------------------
# DATA_DIR must point at the unpacked Zenodo archive `DiCoLo_data/`.
# See README.md ("Data") for the expected directory layout and for which files
# inside it are archived derived products rather than raw input.
DATA_DIR   <- Sys.getenv("DICOLO_DATA",
                         file.path(PROJECT_ROOT, "DiCoLo_data"))
FIGURE_DIR <- Sys.getenv("DICOLO_FIGS",
                         file.path(PROJECT_ROOT, "figures"))

if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

# --- Shared code -------------------------------------------------------------
R_DIR <- file.path(PROJECT_ROOT, "R")

source_helpers <- function() {
  source(file.path(R_DIR, "simulation_functions.R"))
  source(file.path(R_DIR, "benchmarking_functions.R"))
  source(file.path(R_DIR, "plotting_functions.R"))
  invisible(TRUE)
}

# --- Python backend ----------------------------------------------------------
# Used by (a) reticulate for Memento and (b) the POT optimal-transport backend
# behind DiCoLo::ComputeGeneEMD(). Set DICOLO_PYTHON to a interpreter that has
# `memento` and `pot` installed; see README.md ("Dependencies").
PYTHON_BIN <- Sys.getenv("DICOLO_PYTHON", Sys.which("python3"))

# Directory holding the DiCoLo package's OT helper scripts; passed to
# ComputeGeneEMD(script_dir = ). NULL uses the package default.
OT_SCRIPT_DIR <- local({
  p <- Sys.getenv("DICOLO_OT_SCRIPTS", "")
  if (nzchar(p)) p else NULL
})

use_dicolo_python <- function() {
  if (!nzchar(PYTHON_BIN)) {
    stop("No python3 found. Set DICOLO_PYTHON to a suitable interpreter.")
  }
  reticulate::use_python(PYTHON_BIN, required = TRUE)
  invisible(PYTHON_BIN)
}

# --- Reproducibility ---------------------------------------------------------
GLOBAL_SEED <- 42L

message(sprintf("[config] root   : %s", PROJECT_ROOT))
message(sprintf("[config] data   : %s", DATA_DIR))
message(sprintf("[config] figures: %s", FIGURE_DIR))
