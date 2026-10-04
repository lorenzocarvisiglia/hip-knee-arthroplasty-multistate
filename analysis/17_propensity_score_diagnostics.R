source("config.R")

# propensity-score diagnostics for the overlap-weighted sensitivity analysis
# the implementation is split into small modules to keep the public repository readable

source(file.path("R", "ps_diagnostics_01_setup.R"))
source(file.path("R", "ps_diagnostics_02_analysis.R"))
source(file.path("R", "ps_diagnostics_03_figures_qc.R"))
source(file.path("R", "ps_diagnostics_04_export.R"))
