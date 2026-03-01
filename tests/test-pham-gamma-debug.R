## Debug: trace the gamma discrepancy between reference and ours.
##
## Focus on a single timestep u, single arm a, to see exactly where they diverge.

pkg_dir = "/Users/skip/work/balancing/R"
for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R",
            "estimands.R", "survival.R", "estimate.R"))
  source(file.path(pkg_dir, f))

# Reference functions
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
ref_kernel.minimax.weights = function(dpsi, Z, Kfun, sigma) {
  K = ref_outerproduct(Z, Z, Kfun)
  KZjs = lapply(1:nrow(Z), function(j) function(z) Kfun(z, Z[j, , drop = FALSE]))
  dpsiK = sapply(KZjs, dpsi)
  gamma = solve(K + sigma^2 * diag(nrow(K)), nrow(Z) * dpsiK)
  list(gamma = gamma, K = K, dpsiK = dpsiK, h_K = nrow(Z) * dpsiK)
}

# Small example: 20 subjects, timestep u=2, arm a=1
set.seed(42)
n = 20
p = 2
X = matrix(rnorm(n * p), n, p)
A = rbinom(n, 1, 0.5)
Z = cbind(A, X)

sigma_ref = 1
eta_ours = sigma_ref^2 / 2
nu = 3/2; sigma_kern = 5

ref_Kfun.x = function(x, xp) ref_K.matern(x, xp, sigma = sigma_kern, nu = nu)
ref_Kfun = ref_direct.product.kernel(ref_Kfun.x, iw = 1)
our_kern_x = matern_kernel(sigma = sigma_kern, nu = nu)

# Fake dot_psi values: r_u[i] for each subject, evaluated at counterfactual Za
# (these would normally come from the estimand, but let's use simple values)
arm = 1
Za = with.treatment(Z, arm)
r_u = rnorm(n)  # pretend these are -S(tau)/(1-h_u) values

cat("=== Setup ===\n")
cat("n =", n, ", arm =", arm, "\n")
cat("n_treated =", sum(A == 1), ", n_control =", sum(A == 0), "\n\n")

# ============================================================
# Reference gamma
# ============================================================
# dpsi_u(m) = mean(r_u * m(Za))
dpsi_u = function(m) mean(r_u * m(Za))

ref = ref_kernel.minimax.weights(dpsi_u, Z, ref_Kfun, sigma = sigma_ref)

cat("=== Reference ===\n")
cat("h_K (first 5):", round(ref$h_K[1:5], 4), "\n")
cat("gamma (first 5):", round(ref$gamma[1:5], 4), "\n")
cat("gamma for A=0:", round(ref$gamma[A == 0], 4), "\n")
cat("gamma for A=1:", round(ref$gamma[A == 1], 4), "\n\n")

# ============================================================
# Our gamma via kernel_bregman with quadratic dispersion
# ============================================================
# Pooled data: just one timestep u=1 for simplicity
u_val = 1
Z_pool = cbind(u = rep(u_val, n), Z)

# Our targets: Y_gam[i] = r_u[i]
# (In the actual code, Y_gam = W_pool * deriv_vals, but for this test
# we're comparing per-arm gamma, so Y_gam = r_u directly)
Y_gam = r_u

our_kern_pooled = product_kernel(our_kern_x, iw = c(1, 2))  # group on (u, W)

our_fit = kernel_bregman(Y_gam, Z_pool, our_kern_pooled,
                          n_eta = eta_ours,
                          dispersion = quadratic_dispersion())

# Gamma values at training points = predict at same points
our_gamma = predict.kernel_bregman(our_fit, Z_pool)

cat("=== Ours (kernel regression of r_u targets) ===\n")
cat("gamma (first 5):", round(our_gamma[1:5], 4), "\n")
cat("gamma for A=0:", round(our_gamma[A == 0], 4), "\n")
cat("gamma for A=1:", round(our_gamma[A == 1], 4), "\n\n")

# ============================================================
# Trace the algebra
# ============================================================
# Reference: (K_ZZ + sigma^2 I) gamma = h_K
# where h_K[j] = n * dpsi_u(K(Z_j, .)) = sum_i r_u[i] * K(Za_i, Z_j)
#
# K(Za_i, Z_j) = product_kernel((a, X_i), (A_j, X_j))
#              = K_x(X_i, X_j) if a == A_j, else 0
#
# So for A_j = a=1: h_K[j] = sum_{ALL i} r_u[i] * K_x(X_i, X_j)
# For A_j = 0:      h_K[j] = 0

# Let's verify:
K_x_full = ref_outerproduct(X, X, ref_Kfun.x)  # n x n
h_K_manual = numeric(n)
for (j in 1:n) {
  if (A[j] == arm) {
    h_K_manual[j] = sum(r_u * K_x_full[, j])  # sum over ALL i
  } else {
    h_K_manual[j] = 0
  }
}
cat("=== Manual h_K verification ===\n")
cat("max|h_K_ref - h_K_manual| =", max(abs(ref$h_K - h_K_manual)), "\n\n")

# Ours: (K_pool + 2*eta*I) alpha = Y_gam
# K_pool is block-diagonal: K_x for A=1 block, K_x for A=0 block
# gamma at training points = K_pool * alpha
#
# For A=1 block: (K_x[11] + 2*eta*I) alpha_1 = r_u[A=1]
#   gamma_1 = K_x[11] * alpha_1 = K_x[11] * (K_x[11] + 2*eta*I)^{-1} * r_u[A=1]
#
# Reference for A=1 block: (K_ZZ[11] + sigma^2*I) gamma_1 = h_K[A=1]
#   K_ZZ[11] = K_x[11] (product kernel restricted to same arm)
#   gamma_1 = (K_x[11] + sigma^2*I)^{-1} * h_K[A=1]
#
# KEY DIFFERENCE:
# Ours: gamma_1 = K_x[11] * (K_x[11]+sigma^2*I)^{-1} * r_u[A=1]
#                = (K_x[11]+sigma^2*I)^{-1} * K_x[11] * r_u[A=1]
#                = (K_x[11]+sigma^2*I)^{-1} * {sum of K_x * r_u WITHIN arm 1}
#
# Reference: gamma_1 = (K_x[11]+sigma^2*I)^{-1} * h_K[A=1]
#                     = (K_x[11]+sigma^2*I)^{-1} * {sum of K_x * r_u over ALL subjects}

# These differ because: K_x[11] * r_u[A=1] sums K_x(X_i, X_j) * r_u[i] only for i in arm 1
# while h_K[A=1] sums K_x(X_i, X_j) * r_u[i] for ALL i

# Let's verify this is the source of the discrepancy:
idx1 = which(A == 1)
idx0 = which(A == 0)

K_x_11 = K_x_full[idx1, idx1]  # within-arm kernel
r_u_1 = r_u[idx1]

# RHS in our solve: K_x[11] * r_u[A=1]
our_rhs_1 = as.vector(K_x_11 %*% r_u_1)

# RHS in reference solve: h_K[A=1] = sum over ALL i of r_u[i] * K_x(X_i, X_j) for j in arm 1
ref_rhs_1 = as.vector(K_x_full[idx1, ] %*% r_u)  # full column sum

cat("=== RHS comparison (arm 1 block) ===\n")
cat("our_rhs (K_x[11]*r_u[1]):", round(our_rhs_1, 3), "\n")
cat("ref_rhs (K_x[1,]*r_u)  :", round(ref_rhs_1, 3), "\n")
cat("difference:", round(ref_rhs_1 - our_rhs_1, 3), "\n")
cat("The difference = K_x[10]*r_u[0] (cross-arm contribution):",
    round(as.vector(K_x_full[idx1, idx0] %*% r_u[idx0]), 3), "\n\n")

# So the question is: should our targets include the cross-arm contribution?
# The Riesz representer gamma_psi for E[psi(a, X)] satisfies:
# E[gamma_psi(Z) * m(Z)] = E[dot_psi(a, X) * m(a, X)]
#
# In the RKHS with product kernel:
# gamma_psi(a, x) = dot_psi(a, x) for the correct arm
# gamma_psi(1-a, x) = 0 for the wrong arm
#
# The minimax weights approximate this by solving in the RKHS.
# h_K encodes the target functional: h_K[j] = E_n[dot_psi * K(Za, Z_j)]
# This uses counterfactual Za, so the cross-arm kernel K(Za, Z_j) is zero
# when A_j != a. But the sum over i includes ALL subjects evaluating at Za.
#
# In our regression: Y[i] = dot_psi[i] for all i, and the kernel separates
# by arm. So within the arm-a block, we're regressing dot_psi on K_x using
# only arm-a subjects. This misses the cross-arm contribution.

cat("=== Conclusion ===\n")
cat("The discrepancy comes from the cross-arm kernel contribution in h_K.\n")
cat("Reference h_K[j] = sum_{ALL i} r_u[i] * K_x(X_i, X_j) for j in arm a.\n")
cat("Our K_x[aa] * r_u[a] = sum_{i in arm a} r_u[i] * K_x(X_i, X_j).\n")
cat("Difference = sum_{i in arm 1-a} r_u[i] * K_x(X_i, X_j) (cross-arm).\n")
