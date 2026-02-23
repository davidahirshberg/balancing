## Check the regularization scaling: is sigma^2 = 2*eta or 2*n*eta?
## Lambda test already passed with sigma^2 = 2*eta.
## Test both to be sure.

pkg_dir = "/Users/skip/work/balancing/R"
for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R"))
  source(file.path(pkg_dir, f))

ref_outerproduct = function(A, B, K) {
  A = atleast_2d(A); B = atleast_2d(B)
  O = array(dim = c(nrow(A), nrow(B)))
  for (ib in 1:nrow(B)) O[, ib] = K(A, B[ib, , drop = FALSE])
  O
}
ref_K.matern = function(z, zp, sigma = 1, nu = 3/2) {
  difference.norm.sq = colSums((t(z) - c(zp))^2)
  d = (sqrt(2 * nu) / sigma) * sqrt(difference.norm.sq)
  (1 + d) * exp(-d)
}

set.seed(1)
n = 30
p = 2
X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)

sigma_kern = 5
ref_Kfun = function(x, xp) ref_K.matern(x, xp, sigma = sigma_kern)
our_kern = matern_kernel(sigma = sigma_kern, nu = 3/2)

# Reference ridge: (K + sigma^2 I) beta = Y
# Prediction at new points: K_new %*% beta
sigma_ref = 1
K_ref = ref_outerproduct(X, X, ref_Kfun)
beta_ref = solve(K_ref + sigma_ref^2 * diag(n), Y)

X_test = matrix(rnorm(10 * p), 10, p)
K_test_ref = ref_outerproduct(X_test, X, ref_Kfun)
pred_ref = K_test_ref %*% beta_ref

# Our solver: (K + 2*eta*I) alpha = Y
# Try eta = sigma^2 / 2
eta1 = sigma_ref^2 / 2
fit1 = kernel_bregman(Y, X, our_kern, eta = eta1,
                       dispersion = quadratic_dispersion(), intercept = FALSE)
pred1 = predict_phi(fit1, X_test)

# Try eta = sigma^2 / (2*n)
eta2 = sigma_ref^2 / (2 * n)
fit2 = kernel_bregman(Y, X, our_kern, eta = eta2,
                       dispersion = quadratic_dispersion(), intercept = FALSE)
pred2 = predict_phi(fit2, X_test)

cat("=== Regularization scaling test ===\n")
cat("n =", n, ", sigma_ref =", sigma_ref, "\n")
cat("eta1 = sigma^2/2 =", eta1, "\n")
cat("eta2 = sigma^2/(2n) =", eta2, "\n\n")

cat("Predictions at test points:\n")
cat("ref:                ", round(pred_ref[1:5], 6), "\n")
cat("ours (sigma^2=2eta):", round(pred1[1:5], 6), "\n")
cat("ours (sigma^2=2n*eta):", round(pred2[1:5], 6), "\n\n")

cat("Max error (sigma^2=2eta):  ", max(abs(pred_ref - pred1)), "\n")
cat("Max error (sigma^2=2n*eta):", max(abs(pred_ref - pred2)), "\n")

# Also check: reference minimax weights vs our gamma at TRAINING points
# Reference: gamma = (K + sigma^2 I)^{-1} (n * dpsiK)
# For a simple functional: dpsi(m) = mean(m(X)) = (1/n) sum m(X_i)
# dpsiK[j] = dpsi(K(X_j, .)) = (1/n) sum_i K(X_i, X_j) = (1/n) * (K %*% rep(1,n))[j]
# h_K = n * dpsiK = K %*% rep(1, n)

dpsiK = K_ref %*% rep(1, n) / n
h_K = n * dpsiK

gamma_ref = solve(K_ref + sigma_ref^2 * diag(n), h_K)

# Our regression with Y = rep(1, n) (the per-subject functional: dpsi_i(m) = m(X_i))
Y_simple = rep(1, n)

# Predicted gamma at training points: K * alpha
gam_fit1 = kernel_bregman(Y_simple, X, our_kern, eta = eta1,
                           dispersion = quadratic_dispersion(), intercept = FALSE)
gamma_ours1 = predict.kernel_bregman(gam_fit1, X)

gam_fit2 = kernel_bregman(Y_simple, X, our_kern, eta = eta2,
                           dispersion = quadratic_dispersion(), intercept = FALSE)
gamma_ours2 = predict.kernel_bregman(gam_fit2, X)

cat("\n=== Minimax weights comparison ===\n")
cat("Simple functional: dpsi(m) = mean(m(X))\n")
cat("gamma_ref (first 5):           ", round(gamma_ref[1:5], 6), "\n")
cat("gamma_ours (sigma^2=2eta):     ", round(gamma_ours1[1:5], 6), "\n")
cat("gamma_ours (sigma^2=2n*eta):   ", round(gamma_ours2[1:5], 6), "\n\n")

cat("Max error (sigma^2=2eta):  ", max(abs(gamma_ref - gamma_ours1)), "\n")
cat("Max error (sigma^2=2n*eta):", max(abs(gamma_ref - gamma_ours2)), "\n")

# Key identity check: does K*Y = h_K?
# h_K = K * rep(1,n). K*Y = K * rep(1,n). So h_K = K*Y trivially! ✓
# (no cross-arm issue since there's no treatment grouping here)
cat("\nmax|K*Y - h_K| =", max(abs(K_ref %*% Y_simple - h_K)), "\n")
