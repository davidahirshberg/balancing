## Reimplement the spinoff3 survival simulation on top of the balancing package.
## This stress-tests bregbal_surv: training mesh, sign-flip gamma, compensator,
## and the full DR pipeline.

# Load package
pkg_dir = "/Users/skip/work/balancing/R"
for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R",
            "weights.R", "estimands.R", "tuning.R", "crossfit.R",
            "survival.R", "estimate.R", "dgp.R"))
  source(file.path(pkg_dir, f))

# ============================================================
# DGP: paper continuous-time (Section 6.2)
# ============================================================
dgp = paper_cts_dgp(horizon = 1)

set.seed(0)
big = dgp$generate(5000, p = 2)
truth_surv = big$psi1_true - big$psi0_true
truth_rmst = big$rmst1_true - big$rmst0_true
cat(sprintf("True surv prob ATE = %.5f\n", truth_surv))
cat(sprintf("True RMST ATE = %.5f\n", truth_rmst))

# ============================================================
# Simulation (small scale — API stress test, not power study)
# ============================================================
kern = matern_kernel(sigma = 2, nu = 3/2)
est_surv = surv_prob_estimand()

n = 200; n_reps = 3
cat(sprintf("n=%d, n_reps=%d\n\n", n, n_reps))

results = data.frame(rep = integer(), estimand = character(),
                     est = numeric(), se = numeric(), truth = numeric(),
                     stringsAsFactors = FALSE)

for (rep_i in 1:n_reps) {
  set.seed(rep_i)
  dat = dgp$generate(n, p = 2)

  res = bregbal_surv(dat$T_obs, dat$D, dat$A, dat$X, kern,
                     estimand = est_surv,
                     horizon = dgp$horizon,
                     lam_dispersion = entropy_dispersion(),
                     gam_dispersion = entropy_dispersion(),
                     eta_lam = 0.5,
                     eta_gam = 5,
                     M_train = 15,
                     Q_comp = 50,
                     n_folds = 2)

  results = rbind(results, data.frame(
    rep = rep_i, estimand = "surv_prob",
    est = res$est, se = res$se, truth = truth_surv))

  if (rep_i %% 5 == 0) cat(sprintf("  rep %d/%d done\n", rep_i, n_reps))
}

# Summary
cat("\n=== Results ===\n")
bias = mean(results$est - results$truth)
rmse = sqrt(mean((results$est - results$truth)^2))
cov = mean(abs(results$est - results$truth) < 1.96 * results$se)
cat(sprintf("surv_prob  bias=%.4f  rmse=%.4f  coverage=%.2f\n", bias, rmse, cov))
