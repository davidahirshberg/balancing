#!/usr/bin/env Rscript
## Simple coverage test: does surv_dr + bootstrap work?
##
## One DGP (paper_cts), one kernel (smooth), one eta (1/n), two estimands.
## No grid sweep, no tuning selection, no CV.
##
## Cluster script: runs one rep, saves raw dr_result + boot_result objects.
## Postprocessing (coverage, summaries) happens locally.
##
## Usage:
##   Rscript examples/test-simple.R --rep 1 --n 200 --boot 200
##   # or via SLURM array (SLURM_ARRAY_TASK_ID sets rep)

# ============================================================
# Parse arguments
# ============================================================
args = commandArgs(trailingOnly = TRUE)
n_obs = 200; boot_reps = 200; rep_id = NULL

i = 1
while (i <= length(args)) {
  if (args[i] == "--n")         { n_obs = as.integer(args[i + 1]); i = i + 2 }
  else if (args[i] == "--boot") { boot_reps = as.integer(args[i + 1]); i = i + 2 }
  else if (args[i] == "--rep")  { rep_id = as.integer(args[i + 1]); i = i + 2 }
  else { i = i + 1 }
}

if (is.null(rep_id)) {
  rep_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))
  if (is.na(rep_id)) stop("Need --rep or SLURM_ARRAY_TASK_ID")
}

# ============================================================
# Setup
# ============================================================
source("R/kernel.R")
source("R/dispersions.R")
source("R/estimands.R")
source("R/survival.R")
source("R/dgp.R")
source("R/surv_estimate.R")

horizon = 1
p = 2
kern = direct_product(matern_kernel(sigma = 2, nu = 3/2), iw = 2, levels = c(0, 1))
lambda_disp = entropy_dispersion()
gamma_disp_fn = target_scaled_entropy
estimand_names = c("surv_prob", "rmst")

make_estimand = function(name) {
  switch(name, surv_prob = survival_probability_ate(), rmst = rmst_ate())
}

dgp = paper_cts_dgp(horizon = horizon)

# ============================================================
# Run one rep
# ============================================================
cat(sprintf("rep %d: n=%d, boot=%d, eta=1/n\n", rep_id, n_obs, boot_reps))
t0 = proc.time()
set.seed(rep_id)
dat = dgp$generate(n_obs, p = p)
n = n_obs
eta = 1 / n

rep_results = list()
for (en in estimand_names) {
  estimand = make_estimand(en)
  t1 = proc.time()

  log = make_logger()
  dr = survival_effect(list(T_obs = dat$T_obs, D = dat$D, Z = dat$Z),
                       kern, estimand, lambda_disp,
                       eta_lam = eta, eta_gam = eta, horizon = horizon,
                       gamma_disp_fn = gamma_disp_fn, log = log)
  boot = bootstrap(dr, boot_reps, log = log)

  elapsed_en = (proc.time() - t1)[3]
  cat(sprintf("  %s: est=%.5f, se=%.5f (%.1fs)\n",
              en, dr$est, dr$se, elapsed_en))
  rep_results[[en]] = list(dr = dr, boot = boot, log_data = log$all())
}

elapsed = (proc.time() - t0)[3]
cat(sprintf("rep %d done: %.1fs\n", rep_id, elapsed))

# Save raw objects
out_file = sprintf("examples/simple-n%d-rep%03d.rds", n_obs, rep_id)
saveRDS(list(rep_id = rep_id, n = n_obs, eta = eta,
             boot_reps = boot_reps, results = rep_results),
        out_file)
cat(sprintf("Saved to %s\n", out_file))
