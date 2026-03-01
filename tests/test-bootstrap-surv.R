## Smoke test: bregbal_surv with bootstrap.
## Verify bootstrap machinery runs without error and produces sensible output.
##
## Flow: grid sweep → pick a point → bootstrap via init.

pkg_dir = "/Users/skip/work/balancing/R"
for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R",
            "weights.R", "estimands.R", "tuning.R", "crossfit.R",
            "survival.R", "estimate.R", "dgp.R"))
  source(file.path(pkg_dir, f))

dgp = paper_cts_dgp(horizon = 1)
set.seed(42)
big = dgp$generate(5000, p = 2)
truth = big$psi1_true - big$psi0_true
cat(sprintf("True ATE = %.5f\n", truth))

kern = matern_kernel(sigma = 2, nu = 3/2)

set.seed(1)
dat = dgp$generate(200, p = 2)

# --- 1. Grid sweep (no bootstrap yet) ---
cat("Running grid sweep...\n")
t0 = proc.time()
grid = bregbal_surv(dat$T_obs, dat$D, dat$A, dat$X, kern,
                    estimand = surv_prob_estimand(),
                    horizon = dgp$horizon,
                    lam_dispersion = entropy_dispersion(),
                    gam_dispersion = entropy_dispersion(),
                    sigma2_lam_grid = 2^(-6:6),
                    sigma2_gam_grid = 2^(-6:6),
                    M_train = 15,
                    Q_comp = 50,
                    n_folds = 2,
                    bootstrap_reps = 20)
t_grid = (proc.time() - t0)[3]
cat(sprintf("Grid done in %.1f seconds\n", t_grid))
cat(sprintf("Path rows: %d\n", nrow(grid$path)))

# --- 2. Select a point via Lepski ---
n = nrow(dat$X)
path = grid$path
lam_vals = sort(unique(path$sigma2_lam))
gam_vals = sort(unique(path$sigma2_gam))

jg_mid = ceiling(length(gam_vals) / 2)
lam_slice = path[path$sigma2_gam == gam_vals[jg_mid], ]
lam_slice = lam_slice[order(lam_slice$sigma2_lam), ]

sel = lepski_tuning()$select(tuning_context(
  lam_vals / n, n = n,
  dir_val = function() lam_slice$dir_val,
  dir_se = function() lam_slice$dr_se,
  loo_scores = function() lam_slice$loo_lam))

sel_sigma2_lam = sel$eta * n
gam_slice = path[path$sigma2_lam == sel_sigma2_lam, ]
gam_slice = gam_slice[order(gam_slice$sigma2_gam), ]
sel_sigma2_gam = gam_vals[which.min(gam_slice$loo_gam)]

cat(sprintf("Selected: sigma2_lam=%.3g, sigma2_gam=%.3g\n",
            sel_sigma2_lam, sel_sigma2_gam))

# --- 3. Bootstrap at selected point via init ---
cat("Running bootstrap (20 reps) via init...\n")
t0 = proc.time()
res = bregbal_surv(dat$T_obs, dat$D, dat$A, dat$X, kern,
                   estimand = surv_prob_estimand(),
                   horizon = dgp$horizon,
                   sigma2_lam = sel_sigma2_lam,
                   sigma2_gam = sel_sigma2_gam,
                   bootstrap_reps = 20,
                   init = grid$init)
t_boot = (proc.time() - t0)[3]
cat(sprintf("Bootstrap done in %.1f seconds\n\n", t_boot))

cat(sprintf("est = %.5f  (truth = %.5f)\n", res$est, truth))
cat(sprintf("se_eif = %.5f\n", res$se))
cat(sprintf("se_amle = %.5f\n", res$se_amle))
cat(sprintf("se_boot = %.5f\n", res$se_boot))
cat(sprintf("boot_ates range: [%.4f, %.4f]\n", min(res$boot_ates), max(res$boot_ates)))
cat(sprintf("boot_ses range: [%.5f, %.5f]\n", min(res$boot_ses), max(res$boot_ses)))

# Bootstrap-t CI
t_stats = (res$boot_ates - res$est) / res$boot_ses
q = quantile(t_stats, c(0.025, 0.975))
ci_boot_t = res$est - q[2] * res$se_amle
ci_boot_t_hi = res$est - q[1] * res$se_amle
cat(sprintf("Bootstrap-t CI: [%.4f, %.4f]\n", ci_boot_t, ci_boot_t_hi))

# Percentile CI
ci_pct = quantile(res$boot_ates, c(0.025, 0.975))
cat(sprintf("Percentile CI:  [%.4f, %.4f]\n", ci_pct[1], ci_pct[2]))

# Normal CI
ci_norm = res$est + c(-1, 1) * 1.96 * res$se_amle
cat(sprintf("Normal CI:      [%.4f, %.4f]\n", ci_norm[1], ci_norm[2]))

cat(sprintf("\nTruth covered by boot-t: %s\n", truth >= ci_boot_t & truth <= ci_boot_t_hi))
cat(sprintf("Truth covered by pct:    %s\n", truth >= ci_pct[1] & truth <= ci_pct[2]))
cat(sprintf("Truth covered by normal: %s\n", truth >= ci_norm[1] & truth <= ci_norm[2]))
