#!/usr/bin/env Rscript
## Regression test for .phi product kernel path.
## Run BEFORE and AFTER optimizing outer() in .phi.kernel_bregman.
##
## Usage:
##   Rscript test-phi-regression.R save    # save reference output
##   Rscript test-phi-regression.R check   # compare against reference

args = commandArgs(trailingOnly = TRUE)
mode = if (length(args) > 0) args[1] else "check"
ref_file = "test-phi-reference.rds"

setwd(path.expand("~/work/balancing"))
source("R/all.R")
source("R/tuning.R")
source("R/dgp.R")

set.seed(42)
n = 50
res = 4
horiz = 1
M = res
du = horiz / M

# Covariates: treatment A and continuous X
A = rbinom(n, 1, 0.5)
X = matrix(rnorm(n * 2), n, 2)
T_obs = sample(1:M, n, replace = TRUE) * du
D = rbinom(n, 1, 0.7)
Z = cbind(A, X)

base = matern_kernel(sigma = 2, nu = 3/2)
kern = direct_product_kernel(base, iw = 2, levels = c(0, 1))

# Build mesh projection (creates stacked data with time column)
observations = list(T_obs = T_obs, D = D, Z = Z)
mps = mesh_project(observations, horiz, M, n_folds = 2, kern, discrete = TRUE)

# Fit lambda on fold 1
mp = mps[[1]]
lambda_disp = softplus_dispersion()
eta_lam = 1.0
lam_fit = fit_lambda(mp, kern, eta_lam, lambda_disp)

# Test .phi on different Z inputs:
# 1. eval_Z bound path (the hot path we're optimizing)
lam_bound = .bind(lam_fit, mp$eval_Z)
test_Zs = list()

# Predict at each mesh point for all eval subjects
for (k in 1:M) {
  Z_k = cbind(mp$mesh[k], mp$eval_Z)
  test_Zs[[sprintf("mesh_k%d", k)]] = Z_k
}

# Predict with counterfactual arms
eval_X = mp$eval_Z[, -1, drop = FALSE]
test_Zs[["arm1"]] = cbind(mp$mesh[1], 1, eval_X)
test_Zs[["arm0"]] = cbind(mp$mesh[1], 0, eval_X)

# Compute .phi for each test case
results = lapply(test_Zs, function(Z_test) .phi(lam_bound, Z_test))

# Also test fit_gamma end-to-end (uses .phi internally)
eta_gam = 1.0
estimand = survival_probability_ate()
gam_res = fit_gamma(mp, lam_fit, kern, eta_gam, estimand)

# Estimate path
grid = discrete_grid(mp$mesh, du = mp$du_bin)
est = estimate(mp, lam_fit, gam_res, estimand, grid)

results$gamma_fit_alpha = gam_res$gamma_fit$alpha
results$estimate_terms = est$terms
results$estimate_direct = est$direct

if (mode == "save") {
  saveRDS(results, ref_file)
  cat(sprintf("Reference saved to %s (%d test cases)\n", ref_file, length(results)))
  for (nm in names(results))
    cat(sprintf("  %s: length=%d, range=[%.6g, %.6g]\n",
                nm, length(results[[nm]]),
                min(results[[nm]]), max(results[[nm]])))
} else {
  if (!file.exists(ref_file)) stop("Reference file not found. Run with 'save' first.")
  ref = readRDS(ref_file)
  all_pass = TRUE
  for (nm in names(ref)) {
    if (!nm %in% names(results)) {
      cat(sprintf("FAIL: %s missing from results\n", nm))
      all_pass = FALSE
      next
    }
    diff = max(abs(results[[nm]] - ref[[nm]]))
    status = if (diff < 1e-6) "PASS" else "FAIL"
    if (status == "FAIL") all_pass = FALSE
    cat(sprintf("  %s: %s (max diff = %.2e)\n", nm, status, diff))
  }
  if (all_pass) cat("\nAll regression tests PASSED.\n")
  else { cat("\nSome tests FAILED.\n"); quit(status = 1) }
}
