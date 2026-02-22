## Reimplement the spinoff4 binary-outcome simulation on top of the balancing package.
## This is the API stress test: if this is clunky to write, the package needs work.

# Load package
pkg_dir = "/Users/skip/work/balancing/R"
for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R",
            "weights.R", "estimands.R", "tuning.R", "crossfit.R", "estimate.R"))
  source(file.path(pkg_dir, f))

# ============================================================
# DGP: Hainmueller design 3 (binary outcome)
# ============================================================
Sigma_hain = matrix(c(2, 1, -1, 1, 1, -0.5, -1, -0.5, 1), nrow = 3)
R_hain = cov2cor(Sigma_hain)

generate_binary_hain = function(n, zeta = sqrt(10)) {
  X = matrix(rnorm(n * 3), n, 3) %*% chol(R_hain)
  linpred = X %*% c(1, 2, -1)
  pi_x = pnorm(as.numeric(linpred) / zeta)
  W = rbinom(n, 1, pi_x)
  s = X[, 1] + X[, 2] + X[, 3]
  mu1 = plogis((s^2 - 2) / 3)
  mu0 = plogis(s / 3)
  Y = rbinom(n, 1, W * mu1 + (1 - W) * mu0)
  list(X = X, W = W, Y = Y, pi = pi_x, mu1 = mu1, mu0 = mu0,
       tau = mu1 - mu0, tsm1 = mean(mu1), tsm0 = mean(mu0),
       ate = mean(mu1 - mu0), rr = mean(mu1) / mean(mu0), vcate = var(mu1 - mu0))
}

# ============================================================
# Simulation
# ============================================================
set.seed(0)
big = generate_binary_hain(1e6)
truth_tsm = big$tsm1
truth_ate = big$ate
truth_rr = big$rr
truth_vcate = big$vcate
cat(sprintf("True TSM1=%.4f, ATE=%.4f, RR=%.4f, VCATE=%.6f\n",
            truth_tsm, truth_ate, truth_rr, truth_vcate))

kern_x = matern_kernel(sigma = 2, nu = 3/2)
kern = product_kernel(kern_x, iw = 1)
eta_grid = 4^seq(-3, 6, by = 0.5)

n = 500; n_reps = 20
cat(sprintf("n=%d, n_reps=%d\n\n", n, n_reps))

est_tsm = tsm_estimand(arm = 1)
est_rr = rr_estimand()

results = data.frame(rep = integer(), estimand = character(), method = character(),
                     est = numeric(), se = numeric(), truth = numeric(),
                     stringsAsFactors = FALSE)

for (rep_i in 1:n_reps) {
  set.seed(rep_i)
  dat = generate_binary_hain(n)

  # --- Method 1: Lepski tuning ---
  res_lepski = bregbal(dat$Y, dat$W, dat$X, kern,
                       estimand = est_tsm,
                       dispersion = softplus_dispersion(),
                       eta_grid = eta_grid,
                       tuning = lepski_tuning(C = 3),
                       eta_gam = 0.1,
                       n_folds = 3)

  results = rbind(results, data.frame(
    rep = rep_i, estimand = "tsm1", method = "lepski",
    est = res_lepski$est, se = res_lepski$se, truth = truth_tsm))

  # --- Method 2: Oracle tuning ---
  res_oracle = bregbal(dat$Y, dat$W, dat$X, kern,
                       estimand = est_tsm,
                       dispersion = softplus_dispersion(),
                       eta_grid = eta_grid,
                       tuning = oracle_tuning(truth_tsm),
                       eta_gam = 0.1,
                       n_folds = 3)

  results = rbind(results, data.frame(
    rep = rep_i, estimand = "tsm1", method = "oracle",
    est = res_oracle$est, se = res_oracle$se, truth = truth_tsm))

  # --- Method 3: Fixed tuning ---
  res_fixed = bregbal(dat$Y, dat$W, dat$X, kern,
                      estimand = est_tsm,
                      dispersion = softplus_dispersion(),
                      eta_grid = eta_grid,
                      tuning = 1.0,  # bare number -> fixed
                      eta_gam = 0.1,
                      n_folds = 3)

  results = rbind(results, data.frame(
    rep = rep_i, estimand = "tsm1", method = "fixed_1.0",
    est = res_fixed$est, se = res_fixed$se, truth = truth_tsm))

  if (rep_i %% 5 == 0) cat(sprintf("  rep %d/%d done\n", rep_i, n_reps))
}

# Summary
cat("\n=== Results ===\n")
for (m in unique(results$method)) {
  sub = results[results$method == m, ]
  bias = mean(sub$est - sub$truth)
  rmse = sqrt(mean((sub$est - sub$truth)^2))
  cov = mean(abs(sub$est - sub$truth) < 1.96 * sub$se)
  cat(sprintf("%-12s  bias=%.4f  rmse=%.4f  coverage=%.2f\n", m, bias, rmse, cov))
}
