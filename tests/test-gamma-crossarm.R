## Verify: gamma matches when h_K = K*Y (no cross-arm issue).
## The discrepancy arises because the reference functional uses counterfactual
## evaluation K(Za, Z_j) while our product kernel uses K(Z_i, Z_j).
##
## Test: within a SINGLE arm block (all subjects treated), no cross-arm issue.
## Gamma should match exactly.

pkg_dir = "/Users/skip/work/balancing/R"
for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R"))
  source(file.path(pkg_dir, f))

ref_outerproduct = function(A, B, K) {
  A = atleast_2d(A); B = atleast_2d(B)
  O = array(dim = c(nrow(A), nrow(B)))
  for (ib in 1:nrow(B)) O[, ib] = K(A, B[ib, , drop = FALSE])
  O
}
ref_K.matern = function(z, zp, sigma = 5, nu = 3/2) {
  difference.norm.sq = colSums((t(z) - c(zp))^2)
  d = (sqrt(2 * nu) / sigma) * sqrt(difference.norm.sq)
  (1 + d) * exp(-d)
}

no_null = function(Z) matrix(nrow = nrow(atleast_2d(Z)), ncol = 0)

set.seed(42)
n = 20
p = 2
X = matrix(rnorm(n * p), n, p)
A = rbinom(n, 1, 0.5)
Z = cbind(A, X)

sigma_ref = 1
# Paper-scale eta: solver ridge = n*eta*I, reference = sigma^2*I → eta = sigma^2/n
eta_from_sigma = function(sigma, n) sigma^2 / n
arm = 1
Za = Z; Za[, 1] = arm
r_u = rnorm(n)

ref_Kfun.x = function(x, xp) ref_K.matern(x, xp, sigma = 5, nu = 3/2)
ref_Kfun = function(wx, wxp) {
  niw = 2:ncol(wx)
  w  = wx[, 1];   x  = wx[, niw, drop = FALSE]
  wp = wxp[, 1];  xp = wxp[, niw, drop = FALSE]
  ifelse(w == wp, ref_Kfun.x(x, xp), 0)
}
our_kern_x = matern_kernel(sigma = 5, nu = 3/2, null_basis = no_null)

# ============================================================
# Test 1: Single arm only (all subjects have A=arm)
# ============================================================
cat("=== Test 1: Single arm (no cross-arm) ===\n")
idx_arm = which(A == arm)
n_arm = length(idx_arm)
Z_arm = Z[idx_arm, , drop = FALSE]
r_arm = r_u[idx_arm]

# Reference minimax weights on arm-a subjects only
dpsi_arm = function(m) mean(r_arm * m(Za[idx_arm, , drop = FALSE]))
K_ref = ref_outerproduct(Z_arm, Z_arm, ref_Kfun)
KZjs = lapply(1:n_arm, function(j) function(z) ref_Kfun(z, Z_arm[j, , drop = FALSE]))
dpsiK = sapply(KZjs, dpsi_arm)
h_K = n_arm * dpsiK
gamma_ref = solve(K_ref + sigma_ref^2 * diag(n_arm), h_K)

# Our kernel regression on arm-a subjects only (no product kernel needed)
# Pooled with single timestep u=1
Z_pool = cbind(u = 1, Z_arm)
our_kern_pooled = product_kernel(our_kern_x, iw = c(1, 2))
K_pool = kernel_matrix(Z_pool, Z_pool, our_kern_pooled)
B_pool = null_basis(Z_pool, our_kern_pooled)
tgt = dpsi_from_response(r_arm, K_pool, B_pool)
our_fit = kernel_bregman(Z_pool, our_kern_pooled,
                          eta = eta_from_sigma(sigma_ref, nrow(Z_pool)),
                          dispersion = quadratic_dispersion(),
                          target = tgt$target, K = K_pool)
gamma_ours = predict.kernel_bregman(our_fit, Z_pool)

# Check K*Y = h_K within this single arm
K_x_arm = ref_outerproduct(X[idx_arm, , drop = FALSE],
                            X[idx_arm, , drop = FALSE], ref_Kfun.x)
KY_arm = K_x_arm %*% r_arm
cat("max|K*Y - h_K|:", max(abs(KY_arm - h_K)), "\n")
cat("gamma_ref:", round(gamma_ref, 4), "\n")
cat("gamma_ours:", round(gamma_ours, 4), "\n")
cat("max|gamma_ref - gamma_ours|:", max(abs(gamma_ref - gamma_ours)), "\n\n")

# ============================================================
# Test 2: Both arms, but reference uses OBSERVED treatment (not counterfactual)
# ============================================================
cat("=== Test 2: Both arms, observed treatment functional ===\n")

# Modified reference: dpsi uses Z (observed), not Za (counterfactual)
dpsi_obs = function(m) mean(r_u * m(Z))  # observed Z, not Za
K_full = ref_outerproduct(Z, Z, ref_Kfun)
KZjs_full = lapply(1:n, function(j) function(z) ref_Kfun(z, Z[j, , drop = FALSE]))
dpsiK_full = sapply(KZjs_full, dpsi_obs)
h_K_full = n * dpsiK_full
gamma_ref_obs = solve(K_full + sigma_ref^2 * diag(n), h_K_full)

# Our kernel regression
Z_pool_full = cbind(u = 1, Z)
K_pool_full = kernel_matrix(Z_pool_full, Z_pool_full, our_kern_pooled)

B_pool_full = null_basis(Z_pool_full, our_kern_pooled)
tgt_full = dpsi_from_response(r_u, K_pool_full, B_pool_full)
our_fit_full = kernel_bregman(Z_pool_full, our_kern_pooled,
                               eta = eta_from_sigma(sigma_ref, nrow(Z_pool_full)),
                               dispersion = quadratic_dispersion(),
                               target = tgt_full$target, K = K_pool_full)
gamma_ours_full = predict.kernel_bregman(our_fit_full, Z_pool_full)

# Check K*Y = h_K
K_ZZ = ref_outerproduct(Z, Z, ref_Kfun)
KY_full = K_ZZ %*% r_u
cat("max|K*Y - h_K|:", max(abs(KY_full - h_K_full)), "\n")
cat("max|gamma_ref - gamma_ours|:", max(abs(gamma_ref_obs - gamma_ours_full)), "\n\n")

# ============================================================
# Test 3: Both arms, reference uses COUNTERFACTUAL (the actual case)
# ============================================================
cat("=== Test 3: Both arms, counterfactual functional (actual reference) ===\n")

dpsi_cfl = function(m) mean(r_u * m(Za))  # counterfactual Za
KZjs_full2 = lapply(1:n, function(j) function(z) ref_Kfun(z, Z[j, , drop = FALSE]))
dpsiK_full2 = sapply(KZjs_full2, dpsi_cfl)
h_K_cfl = n * dpsiK_full2
gamma_ref_cfl = solve(K_full + sigma_ref^2 * diag(n), h_K_cfl)

# Check K*Y vs h_K — these should DIFFER due to cross-arm
KY_same = K_ZZ %*% r_u
cat("max|K*Y - h_K|:", max(abs(KY_same - h_K_cfl)), "\n")
cat("  (This is the cross-arm contribution)\n")
cat("gamma_ref_cfl for A=0:", round(gamma_ref_cfl[A == 0], 4), "\n")
cat("gamma_ours for A=0:   ", round(gamma_ours_full[A == 0], 4), "\n")
cat("max|gamma_ref_cfl - gamma_ours|:", max(abs(gamma_ref_cfl - gamma_ours_full)), "\n")
