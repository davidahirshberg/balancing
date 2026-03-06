#' ---
#' title: "Regression Tests: External Method Correspondences"
#' ---
#'
#' This file documents the mathematical relationships between the
#' `balancing` package and three external implementations:
#' Pham's survival weights, the retargeted-mean estimator, and synthdid.
#' Each section derives the correspondence and verifies it numerically.
#'
#' Run from `balancing/` root: `Rscript tests/test-external.R`

source("R/all.R")

tol = 1e-6
pass_all = TRUE

report = function(name, maxdiff, tol) {
  status = if (maxdiff < tol) "PASS" else "FAIL"
  cat(sprintf("  [%s] %s: max |diff| = %.2e (tol = %.0e)\n", status, name, maxdiff, tol))
  if (status == "FAIL") pass_all <<- FALSE
}

#' # Test 1: Pham Survival Code (Norm Penalty)
#'
#' ## Mathematical Correspondence
#'
#' **Paper** (Pham et al.): The kernel balancing weights solve
#' $$(K + \sigma^2 I)\gamma = n\hat{P}\dot\psi_K$$
#' where $\hat{P}\dot\psi_K$ is the Riesz representer of the
#' linear functional $\dot\psi$ projected onto the span of
#' kernel basis functions, and $\sigma^2$ is the regularization
#' parameter. No intercept (norm penalty, not seminorm).
#'
#' For the mean functional $\dot\psi(f) = n^{-1}\sum_i f(Z_i)$:
#' $$[\hat{P}\dot\psi_K]_j = n^{-1}\sum_i K(Z_i, Z_j)$$
#' so $n\hat{P}\dot\psi_K = K\mathbf{1}$.
#'
#' **Code**: Our `kernel_bregman` with quadratic dispersion,
#' no null space, $\eta = \sigma^2/n$ solves
#' $$K(K + n\eta I)\alpha = c$$
#' where $c = K r$ is the projected target ($r = \mathbf{1}$).
#' Since $n\eta = \sigma^2$:
#' $$\alpha = (K + \sigma^2 I)^{-1} r, \quad
#'   \phi = K\alpha = K(K + \sigma^2 I)^{-1}\mathbf{1}$$
#'
#' $K$ and $(K + \sigma^2 I)^{-1}$ share eigenvectors (both are
#' functions of $K$), so they commute:
#' $$\phi = (K + \sigma^2 I)^{-1} K\mathbf{1} = \gamma_{\text{Pham}}$$

cat("\n=== Test 1: Pham survival code (norm penalty) ===\n")

set.seed(42)
n = 80; p = 3
Z = matrix(rnorm(n * p), n, p)
sigma = 2.0
eta = sigma^2 / n

no_null = function(Z) matrix(nrow = nrow(atleast_2d(Z)), ncol = 0)
kern = matern_kernel(sigma = 1.5, nu = 3/2, null_basis = no_null)
K = kernel_matrix(Z, Z, kern)

#' **Reference**: $\gamma = (K + \sigma^2 I)^{-1} K\mathbf{1}$
r = rep(1, n)
gamma_pham = solve(K + sigma^2 * diag(n), K %*% r)

#' **Ours**: `kernel_bregman(quadratic, no null)` with `target = K %*% r`
target = K %*% r
fit = kernel_bregman(Z, kern, eta = eta,
                     dispersion = quadratic_dispersion(),
                     target = target, K = K, tol = 1e-10)
phi_ours = as.vector(K %*% fit$alpha)

report("gamma (phi vs closed-form)", max(abs(phi_ours - gamma_pham)), tol)
report("dchis(phi) = phi (quadratic)", max(abs(.dchis(fit$dispersion, phi_ours) - phi_ours)), 1e-15)


#' # Test 2: Retargeted-Mean Estimator (Seminorm Penalty)
#'
#' ## Mathematical Correspondence
#'
#' **Paper** (Hirshberg & Wager): The seminorm version adds an
#' unpenalized intercept. The augmented KKT system is
#' $$\begin{pmatrix} K + \sigma^2 I & \mathbf{1} \\
#'   \mathbf{1}' & 0 \end{pmatrix}
#'   \begin{pmatrix} \gamma \\ \nu \end{pmatrix} =
#'   \begin{pmatrix} K\mathbf{1} \\ n \end{pmatrix}$$
#' where $\nu$ is the Lagrange multiplier for the constraint
#' $\sum_i \gamma_i = n$.
#'
#' **Code**: Our `kernel_bregman` with quadratic dispersion and
#' `null_basis = matrix(1, n, 1)` (intercept in $\ker(\rho)$):
#' $$\phi = K\alpha + \beta\mathbf{1}, \quad
#'   \gamma = \dot\chi^*(\phi) = \phi$$
#'
#' The solver's FOC for $(\alpha, \beta)$ with $c = K\mathbf{1}$,
#' $c_0 = n$, and $n\eta = \sigma^2$:
#' $$K\phi + \sigma^2 K\alpha = c, \quad \mathbf{1}'\phi = c_0$$
#'
#' Substituting $K\alpha = \phi - \beta\mathbf{1}$ and left-multiplying
#' the first equation by $K^{-1}$:
#' $$(K + \sigma^2 I)\alpha + \beta\mathbf{1} = \mathbf{1}$$
#'
#' With $\gamma = \phi = K\alpha + \beta\mathbf{1}$ and
#' $\sum\gamma = n$, this is the same system as above (with $\nu$
#' absorbing a function of $\beta$ and $\sigma^2$).

cat("\n=== Test 2: retargeted-mean (seminorm penalty) ===\n")

set.seed(123)
n = 60; p = 2
Z = matrix(rnorm(n * p), n, p)
sigma = 1.5
eta = sigma^2 / n

kern_s = matern_kernel(sigma = 1.0, nu = 3/2)  # null_basis = intercept
K = kernel_matrix(Z, Z, kern_s)

#' **Reference**: augmented KKT solve
ones = rep(1, n)
M = rbind(cbind(K + sigma^2 * diag(n), ones), c(ones, 0))
rhs = c(K %*% ones, n)
sol = solve(M, rhs)
gamma_rm = sol[1:n]

#' **Ours**: `kernel_bregman(quadratic, intercept)`
B = null_basis(Z, kern_s)
r = rep(1, n)
target = K %*% r
target_null = as.vector(crossprod(B, r))  # = n

fit_s = kernel_bregman(Z, kern_s, eta = eta,
                       dispersion = quadratic_dispersion(),
                       target = target, target_null = target_null,
                       K = K, tol = 1e-10)
gamma_ours = as.vector(K %*% fit_s$alpha) + as.vector(B %*% fit_s$beta)

report("gamma (seminorm)", max(abs(gamma_ours - gamma_rm)), tol)
report("sum(gamma) = n (ours)", abs(sum(gamma_ours) - n), tol)
report("sum(gamma) = n (ref)", abs(sum(gamma_rm) - n), tol)


#' # Test 3: Synthetic Diff-in-Diff (Entropy + Seminorm)
#'
#' ## Mathematical Correspondence
#'
#' **Paper** (Arkhangelsky et al.): synthdid constructs unit weights
#' $\omega$ on the simplex $\{\omega \geq 0, \sum\omega_i = 1\}$
#' using a linear kernel $K_{ij} = Y_{\text{pre},i}' Y_{\text{pre},j}$
#' on pre-treatment outcomes.
#'
#' The simplex constraint decomposes in the Bregman framework as:
#'
#' 1. **Positivity** $\omega \geq 0$: entropy dispersion
#'    ($\chi = $ negative entropy $\Rightarrow \gamma = \exp(\phi) > 0$).
#' 2. **Sum to 1**: constants in $\ker(\rho)$ (unpenalized intercept
#'    enforces $\sum \gamma_i = c$ for a constant $c$ determined by
#'    the target).
#'
#' The kernel itself is linear on pre-treatment outcome vectors:
#' $k(i, j) = Y_{\text{pre},i}' Y_{\text{pre},j}$.
#' synthdid's regularization $\zeta^2$ maps to $\eta = \zeta^2 / N_0$.
#'
#' **Note**: The objectives differ — synthdid minimizes a constrained
#' quadratic while our framework minimizes $\sum\exp(\phi_i) +
#' (\eta/2)\|\alpha\|_K^2 - c'\alpha$. The weights are structurally
#' analogous (positive, sum-constrained, ridge-regularized on the
#' same kernel) but not numerically identical.

cat("\n=== Test 3: synthdid-style (entropy + seminorm, linear kernel) ===\n")

set.seed(99)
N0 = 30; T0 = 20
Y0 = matrix(rnorm(N0 * T0), N0, T0)

#' **Kernel**: $K_{ij} = Y_{\text{pre},i}' Y_{\text{pre},j}$
K_lin = tcrossprod(Y0)

#' **Target**: $k_i = \sum_t Y_{it} \bar{Y}_t$ where
#' $\bar{Y}_t = N_0^{-1}\sum_i Y_{it}$. This targets the
#' average pre-treatment pattern.
k_target = Y0 %*% colMeans(Y0)

zeta_sq = N0
eta = zeta_sq / N0  # = 1.0

#' **Ours**: `kernel_bregman(entropy, intercept)` on linear kernel
r = k_target / N0
target_ent = K_lin %*% r
B_null = matrix(1, N0, 1)
target_null_ent = as.vector(crossprod(B_null, r))

kern_int = matern_kernel(sigma = 1)  # placeholder; K passed directly
fit_ent = kernel_bregman(matrix(rnorm(N0), N0, 1), kern_int, eta = eta,
                         dispersion = entropy_dispersion(),
                         target = target_ent, target_null = target_null_ent,
                         K = K_lin, tol = 1e-10)
phi_ent = as.vector(K_lin %*% fit_ent$alpha) + fit_ent$beta[1]
gamma_ent = .dchis(fit_ent$dispersion, phi_ent)  # = exp(phi)

#' **Reference**: unconstrained ridge (for context)
omega_ridge = solve(K_lin + zeta_sq * diag(N0), k_target)

cat(sprintf("  Entropy weights: min=%.4f, max=%.4f, sum=%.4f\n",
            min(gamma_ent), max(gamma_ent), sum(gamma_ent)))
cat(sprintf("  Ridge weights:   min=%.4f, max=%.4f, sum=%.4f\n",
            min(omega_ridge), max(omega_ridge), sum(omega_ridge)))

#' Structural checks: entropy weights are positive by construction.
#' Sum constraint is enforced by the null space.
report("entropy: all positive", max(0, -min(gamma_ent)), 1e-15)
# Ridge can go negative (no simplex constraint); entropy cannot
has_negative_ridge = any(omega_ridge < -1e-10)
cat(sprintf("  Ridge has negative weights: %s (expected TRUE)\n", has_negative_ridge))


#' # Summary

cat("\n")
if (pass_all) {
  cat("All tests PASSED.\n")
} else {
  cat("Some tests FAILED.\n")
  quit(status = 1)
}
