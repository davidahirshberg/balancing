# Source all package modules in dependency order.
# Usage: source("R/all.R") from the balancing/ root.
#
# Under devtools::load_all(), files are sourced via Collate order
# in DESCRIPTION, so this file is a no-op.

if (!identical(Sys.getenv("DEVTOOLS_LOAD"), "true") &&
    !environmentName(parent.env(environment())) %in% c("balancing", "imports:balancing")) {
  source("R/kernel.R")       # kernel.R sources utils.R
  source("R/dispersion.R")
  source("R/solver.R")
  source("R/grid.R")
  source("R/estimands.R")
  source("R/survival.R")
  source("R/survival_estimate.R")
  source("R/estimate.R")
}
