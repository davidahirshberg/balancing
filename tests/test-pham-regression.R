## Regression test: our discrete survival estimator vs reference Pham code.
##
## Reference: /Users/skip/work/survival/kernel-amle.qmd
##
## What matches exactly:
## 1. Lambda predictions: pooled kernel_bregman with product_kernel + quadratic_dispersion
##    = per-timestep kernel.ridge from reference. Regularization: sigma^2 = 2*eta.
## 2. Direct term (product-formula survival): both use S(tau) = prod(1-h_k).
## 3. Gamma within a single arm: when all subjects share the same treatment,
##    h_K = K*Y and our kernel regression gamma = AMLE minimax gamma.
##
## What differs (by design):
## 4. Per-arm gamma with mixed treatment: the AMLE functional dpsi(K(Z_j,.))
##    evaluates K(Za, Z_j) with counterfactual treatment, creating cross-arm
##    information flow. Our product kernel K(Z_i, Z_j) blocks this.
##    The AMLE gamma for arm a is zero for wrong-arm subjects; ours isn't.
##    Both converge to the same Riesz representer asymptotically.

# ============================================================
# Source our package
# ============================================================
pkg_dir = "/Users/skip/work/balancing/R"
for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R",
            "estimands.R", "survival.R", "estimate.R"))
  source(file.path(pkg_dir, f))

# ============================================================
# Reference functions (extracted from kernel-amle.qmd)
# ============================================================

ref_outerproduct = function(A, B, K) {
  A = atleast_2d(A); B = atleast_2d(B)
  O = array(dim = c(nrow(A), nrow(B)))
  for (ib in 1:nrow(B)) O[, ib] = K(A, B[ib, , drop = FALSE])
  O
}

ref_K.matern = function(z, zp, sigma = 1, nu = 1/2) {
  difference.norm.sq = colSums((t(z) - c(zp))^2)
  d = (sqrt(2 * nu) / sigma) * sqrt(difference.norm.sq)
  ifelse(d == 0, 1, (2^(1 - nu) / gamma(nu)) * d^nu * besselK(d, nu))
}

ref_direct.product.kernel = function(K, iw = 1) {
  function(wx, wxp) {
    niw = setdiff(1:ncol(wx), iw)
    w  = wx[, iw];   x  = wx[, niw, drop = FALSE]
    wp = wxp[, iw];  xp = wxp[, niw, drop = FALSE]
    ifelse(w == wp, K(x, xp), 0)
  }
}

ref_kernel.ridge = function(Y, Z, Kfun, sigma) {
  K = ref_outerproduct(Z, Z, Kfun)
  beta = solve(K + sigma^2 * diag(nrow(K)), Y)
  structure(list(beta = beta, K = Kfun, Z = Z), class = "ref_kernel.regression")
}

predict.ref_kernel.regression = function(obj, newdata) {
  z = newdata$Z
  KzZ = ref_outerproduct(z, obj$Z, obj$K)
  KzZ %*% obj$beta
}

# Reference survival function (product formula)
ref_S = function(t) {
  stopifnot(t == round(t), t >= 0)
  function(lambda) {
    function(Z) {
      if (t == 0) rep(1, nrow(atleast_2d(Z)))
      else {
        haz_mat = sapply(1:t, function(u) lambda(u)(atleast_2d(Z)))
        if (is.null(dim(haz_mat))) haz_mat = matrix(haz_mat, nrow = 1)
        apply(1 - haz_mat, 1, prod)
      }
    }
  }
}

# ============================================================
# DGP (same as kernel-amle.qmd)
# ============================================================

ref_ht = function(u) {
  if (u == 1) {
    function(Z) rep(0, nrow(atleast_2d(Z)))
  } else {
    function(Z) {
      Z = atleast_2d(Z)
      a = Z[, 1]; x = Z[, 2:ncol(Z), drop = FALSE]
      plogis(-2 + 0.5 * a + 0.2 * rowSums(x) + 0.1 * u)
    }
  }
}

ref_hC = function(u) {
  if (u == 1) {
    function(Z) rep(0, nrow(atleast_2d(Z)))
  } else {
    function(Z) {
      Z = atleast_2d(Z)
      a = Z[, 1]; x = Z[, 2:ncol(Z), drop = FALSE]
      plogis(-3 - 0.3 * a + 0.2 * rowSums(a * x^2) + 0.05 * u)
    }
  }
}

ref_pi = function(X) plogis(X %*% rep(0.3, ncol(X)))

ref_generate_sim_data = function(n, t_max, p, ht, hC) {
  Sigma = 0.8 * diag(p) + 0.2 * matrix(1, p, p)
  N = matrix(rnorm(n * p), n, p)
  X = N %*% chol(Sigma)
  A = rbinom(n, 1, ref_pi(X))
  Z = cbind(A, X)

  T.gt = matrix(NA, n, t_max)
  C.gt = matrix(NA, n, t_max)
  T.gt[, 1] = rbinom(n, 1, 1 - ht(1)(Z))
  C.gt[, 1] = rbinom(n, 1, 1 - hC(1)(Z))
  for (u in 2:t_max) {
    T.gt[, u] = T.gt[, u - 1] * rbinom(n, 1, 1 - ht(u)(Z))
    C.gt[, u] = C.gt[, u - 1] * rbinom(n, 1, 1 - hC(u)(Z))
  }

  T.tilde.gt = T.gt * C.gt
  E = rowSums(T.gt) <= rowSums(C.gt)
  Y = diag(E) %*% cbind(1 - T.tilde.gt[, 1], -t(apply(T.tilde.gt, 1, diff)))
  G = cbind(rep(1, n), T.tilde.gt[, 1:(t_max - 1)])

  list(X = X, A = A, Z = Z, T.gt = T.gt, C.gt = C.gt,
       Y = Y, G = G, E = E)
}

# ============================================================
# Test setup
# ============================================================

set.seed(42)
n = 100
t_max = 4
p = 2
arm = 1

sigma_ref = 1
eta_ours = sigma_ref^2 / 2  # sigma^2 = 2*eta

nu = 3/2; sigma_kern = 5
ref_Kfun.x = function(x, xp) ref_K.matern(x, xp, sigma = sigma_kern, nu = nu)
ref_Kfun = ref_direct.product.kernel(ref_Kfun.x, iw = 1)

our_kern_x = matern_kernel(sigma = sigma_kern, nu = nu)
our_kern = product_kernel(our_kern_x, iw = c(1, 2))

sim = ref_generate_sim_data(n, t_max, p, ref_ht, ref_hC)
Y_mat = sim$Y; G_mat = sim$G; Z = sim$Z

fold = rep(1:2, length.out = n)
train_idx = which(fold == 1)
test_idx  = which(fold == 2)

Z_train = Z[train_idx, ]; Z_test = Z[test_idx, ]
Y_train = Y_mat[train_idx, ]; G_train = G_mat[train_idx, ]
Y_test  = Y_mat[test_idx, ];  G_test  = G_mat[test_idx, ]

cat("=== Pham Regression Test ===\n")
cat("n =", n, ", t_max =", t_max, ", p =", p, "\n\n")

# ============================================================
# Test 1: Lambda predictions
# ============================================================
cat("--- Test 1: Lambda predictions ---\n")

ref_lambda_models = list()
for (u in 1:t_max) {
  ar = G_train[, u] == 1
  if (!any(ar)) next
  ref_lambda_models[[u]] = ref_kernel.ridge(
    Y = Y_train[ar, u], Z = Z_train[ar, , drop = FALSE],
    Kfun = ref_Kfun, sigma = sigma_ref)
}
ref_lambdahat = function(u) {
  function(Znew) {
    if (is.null(ref_lambda_models[[u]])) return(rep(0, nrow(atleast_2d(Znew))))
    as.vector(predict(ref_lambda_models[[u]], newdata = list(Z = Znew)))
  }
}

# Ours: pooled product-kernel ridge
our_mesh = 1:t_max
pool_rows = which(G_train == 1, arr.ind = TRUE)
Z_pool = cbind(u = our_mesh[pool_rows[, "col"]],
               Z_train[pool_rows[, "row"], , drop = FALSE])
Y_pool = Y_train[cbind(pool_rows[, "row"], pool_rows[, "col"])]

our_lam_fit = kernel_bregman(Y_pool, Z_pool, our_kern,
                              eta = eta_ours,
                              dispersion = quadratic_dispersion(),
                              intercept = FALSE)

max_err = 0
for (u in 1:t_max) {
  ref_pred = ref_lambdahat(u)(Z_test)
  our_pred = predict_phi(our_lam_fit, cbind(u, Z_test))
  err = max(abs(ref_pred - our_pred))
  max_err = max(max_err, err)
  cat(sprintf("  u=%d: max|ref-ours| = %.2e\n", u, err))
}
cat(sprintf("  PASS (%.2e < 1e-8)\n\n", max_err))
stopifnot(max_err < 1e-8)

# ============================================================
# Test 2: Direct term (survival probability)
# ============================================================
cat("--- Test 2: Direct term ---\n")

Za = with.treatment(Z_test, arm)
ref_direct = ref_S(t_max)(ref_lambdahat)(Za)

estimand = surv_prob_estimand()
our_h_fn = function(k) predict_phi(our_lam_fit, cbind(our_mesh[k], Za))
our_direct = estimand$psi_discrete(our_h_fn, t_max)

direct_err = max(abs(ref_direct - our_direct))
cat(sprintf("  max|S_ref - S_ours| = %.2e\n", direct_err))
cat(sprintf("  PASS (%.2e < 1e-8)\n\n", direct_err))
stopifnot(direct_err < 1e-8)

# ============================================================
# Test 3: Gamma within single arm (no cross-arm issue)
# ============================================================
cat("--- Test 3: Gamma within single arm ---\n")

# Restrict to test subjects with A = arm at timestep u=2
u_test = 2
ar_all = G_test[, u_test] == 1
ar_arm = ar_all & (Z_test[, 1] == arm)
idx_arm = which(ar_arm)
n_arm = length(idx_arm)

Za_arm = with.treatment(Z_test[idx_arm, , drop = FALSE], arm)
r_u_arm = -(ref_S(t_max)(ref_lambdahat)(Za_arm)) *
            (ref_S(u_test - 1)(ref_lambdahat)(Za_arm)) /
            (ref_S(u_test)(ref_lambdahat)(Za_arm))

# Reference gamma: (K + sigma^2 I)^{-1} h_K on arm-a subjects only
dpsi_arm = function(m) mean(r_u_arm * m(Za_arm))
K_arm_ref = ref_outerproduct(Z_test[idx_arm, ], Z_test[idx_arm, ], ref_Kfun)
KZjs = lapply(1:n_arm, function(j)
  function(z) ref_Kfun(z, Z_test[idx_arm[j], , drop = FALSE]))
dpsiK = sapply(KZjs, dpsi_arm)
gamma_ref = solve(K_arm_ref + sigma_ref^2 * diag(n_arm), n_arm * dpsiK)

# Our gamma: kernel_bregman on same arm-a subjects with same targets
Z_pool_arm = cbind(u = u_test, Z_test[idx_arm, , drop = FALSE])
our_gam_fit_arm = kernel_bregman(r_u_arm, Z_pool_arm, our_kern,
                                  eta = eta_ours,
                                  dispersion = quadratic_dispersion(),
                                  intercept = FALSE)
gamma_ours = predict.kernel_bregman(our_gam_fit_arm, Z_pool_arm)

gam_err = max(abs(gamma_ref - gamma_ours))
cat(sprintf("  u=%d, n_arm=%d: max|gamma_ref - gamma_ours| = %.2e\n",
            u_test, n_arm, gam_err))
cat(sprintf("  PASS (%.2e < 1e-8)\n\n", gam_err))
stopifnot(gam_err < 1e-8)

# ============================================================
# Test 4: Full DR with shared nuisances
# ============================================================
cat("--- Test 4: Full DR (shared lambda, reference gamma) ---\n")

# Use reference's per-arm gamma (the AMLE formulation).
# Run the reference DR and our DR formula with the SAME gamma.
# This tests that our correction formula (dirac - compensator) matches
# the reference's sum_u G * gamma * (Y - lambda).

# Compute reference gamma at all test points for each timestep
ref_gammahat_list = list()
for (u in 1:t_max) {
  ar = G_test[, u] == 1
  if (!any(ar)) { ref_gammahat_list[[u]] = rep(0, nrow(Z_test)); next }

  Za_all = with.treatment(Z_test, arm)
  r_u = -(ref_S(t_max)(ref_lambdahat)(Za_all)) *
          (ref_S(u - 1)(ref_lambdahat)(Za_all)) /
          (ref_S(u)(ref_lambdahat)(Za_all))

  dpsi_u = function(m) mean(r_u * m(Za_all))
  gamma_vals = rep(0, nrow(Z_test))
  Z_ar = Z_test[ar, , drop = FALSE]
  K_ar = ref_outerproduct(Z_ar, Z_ar, ref_Kfun)
  KZjs_ar = lapply(1:sum(ar), function(j)
    function(z) ref_Kfun(z, Z_ar[j, , drop = FALSE]))
  dpsiK_ar = sapply(KZjs_ar, dpsi_u)
  gamma_vals[ar] = solve(K_ar + sigma_ref^2 * diag(sum(ar)), sum(ar) * dpsiK_ar)
  ref_gammahat_list[[u]] = gamma_vals
}

# Reference DR: direct + sum_u G * gamma * (Y - lambda)
ref_correction = rep(0, nrow(Z_test))
for (u in 1:t_max) {
  resid_u = Y_test[, u] - ref_lambdahat(u)(Z_test)
  ref_correction = ref_correction + G_test[, u] * ref_gammahat_list[[u]] * resid_u
}
ref_dr_terms = ref_direct + ref_correction

# Our DR using same gamma: dirac - compensator
# Dirac: sum_u 1{event at u} * gamma(u)
our_dirac = rep(0, nrow(Z_test))
for (u in 1:t_max) {
  events_u = which(G_test[, u] == 1 & Y_test[, u] == 1)
  if (length(events_u) > 0)
    our_dirac[events_u] = our_dirac[events_u] + ref_gammahat_list[[u]][events_u]
}

# Compensator: sum_u G * gamma * lambda
our_comp = rep(0, nrow(Z_test))
for (u in 1:t_max) {
  ar = which(G_test[, u] == 1)
  if (length(ar) == 0) next
  lam_u = predict_phi(our_lam_fit, cbind(u, Z_test[ar, , drop = FALSE]))
  our_comp[ar] = our_comp[ar] + ref_gammahat_list[[u]][ar] * lam_u
}

our_dr_terms = our_direct + our_dirac - our_comp

# Compare per-subject DR terms
dr_err = max(abs(ref_dr_terms - our_dr_terms))
cat(sprintf("  Per-subject max|DR_ref - DR_ours| = %.2e\n", dr_err))
cat(sprintf("  ref mean: %.6f, our mean: %.6f\n",
            mean(ref_dr_terms), mean(our_dr_terms)))
cat(sprintf("  PASS (%.2e < 1e-8)\n\n", dr_err))
stopifnot(dr_err < 1e-8)

# ============================================================
# Sanity check vs truth
# ============================================================
cat("--- Sanity check ---\n")
truth = mean(ref_S(t_max)(ref_ht)(with.treatment(Z, arm)))
cat(sprintf("  True S(%d|A=%d) = %.4f\n", t_max, arm, truth))
cat(sprintf("  DR estimate    = %.4f (|bias| = %.4f)\n",
            mean(ref_dr_terms), abs(mean(ref_dr_terms) - truth)))

cat("\n=== All tests passed! ===\n")
