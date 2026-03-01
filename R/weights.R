## Balancing weights: augmented system with sum constraint.
##
## Solves: (K + n*eta*I) gamma = n hP dot_psi, subject to sum constraint.
## eta is paper-scale: penalty = (eta/2)||f||^2. Ridge diagonal = n*eta*I.
## The constraint is enforced via Lagrange multiplier (augmented system).
##
## For ATE with product kernel: separate gamma1, gamma0 with arm-specific
## sum constraints. For TSM: single gamma with sum-to-n constraint.
##
## This wraps kernel_bregman with quadratic_dispersion (or any other).
## The augmented system is about the constraint, not the dispersion.

#' Compute balancing weights for binary treatment.
#'
#' Solves the augmented system:
#'   [K + n*eta*I,  c] [gamma]   [n hP dot_psi]
#'   [c',           0] [nu   ] = [n           ]
#' where c is the constraint vector (arm indicator).
#'
#' @param Z Covariate matrix with treatment in column iw. (n x p)
#' @param kern Kernel object.
#' @param eta Regularization (paper scale). Penalty = (eta/2)||f||^2.
#' @param iw Treatment column index (default 1).
#' @return List with gamma1, gamma0 (n-vectors of weights).
balancing_weights = function(Z, kern, eta, iw = 1) {
  Z = atleast_2d(Z)
  n = nrow(Z)
  W = Z[, iw]
  K = kernel_matrix(Z, Z, kern)
  A = K + n * eta * diag(n)

  # gamma1: constraint on W==1 units
  c1 = as.numeric(W == 1)
  M1 = rbind(cbind(A, c1), c(c1, 0))
  Z_w1 = Z; Z_w1[, iw] = 1
  n_hP_dpsi_1 = colSums(kernel_matrix(Z_w1, Z, kern))
  sol1 = solve(M1, c(n_hP_dpsi_1, n))
  gamma1 = sol1[1:n]

  # gamma0: constraint on W==0 units
  c0 = as.numeric(W == 0)
  M0 = rbind(cbind(A, c0), c(c0, 0))
  Z_w0 = Z; Z_w0[, iw] = 0
  n_hP_dpsi_0 = colSums(kernel_matrix(Z_w0, Z, kern))
  sol0 = solve(M0, c(n_hP_dpsi_0, n))
  gamma0 = sol0[1:n]

  list(gamma1 = gamma1, gamma0 = gamma0, K = K)
}
