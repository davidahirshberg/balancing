# ============================================================
# Grid: time-axis measure for unified discrete/continuous survival
# ============================================================
#
#' A grid represents the support of the cumulative hazard measure Lambda.
#' All survival operations (survival curve, compensator, RMST) reduce to
#' integration against this measure.
#'
#' The grid carries Lebesgue quadrature weights for int f du:
#'   Discrete: rectangle (exact for piecewise-constant functions).
#'   Continuous: composite Simpson (O(h^4), Q must be even).
#'
#' surv_curve uses trapezoid cumulative (O(h^2)) because cumulative
#' integrals don't reduce to a weight vector.

# ============================================================
# Constructors
# ============================================================

#' Discrete grid: atoms at mesh points.
#' Model output = hazard probability = dLambda directly.
discrete_grid = function(mesh, horizon = max(mesh), du = NULL) {
  M = length(mesh)
  if (is.null(du)) du = mesh[1]
  list(points = mesh, M = M, du = du, horizon = horizon, discrete = TRUE,
       weights = rep(du, M))
}

#' Continuous grid: fine evaluation mesh with Simpson quadrature.
#' Q must be even.
continuous_grid = function(horizon, Q) {
  if (Q %% 2 != 0) stop("Q must be even for Simpson quadrature")
  du = horizon / Q
  M = Q + 1
  w = numeric(M)
  w[1] = du / 3; w[M] = du / 3
  for (j in 2:Q) w[j] = if (j %% 2 == 0) 4 * du / 3 else 2 * du / 3
  list(points = (0:Q) * du, M = M, du = du, horizon = horizon, discrete = FALSE,
       weights = w)
}

# ============================================================
# Survival: S(t) = exp(-int lambda du) or prod(1 - h_k)
# ============================================================

#' Survival curve from hazard increments.
#'
#' Discrete: S(u_k) = prod_{j<=k} (1 - dLambda_j). Exact.
#' Continuous: S(u_k) = exp(-trapezoid_cumsum(dLambda)). O(h^2).
#'
#' @param dLambda Matrix (n x M) of hazard increments.
#' @param grid Grid object. If NULL, uses discrete (product integral).
#' @return Matrix (n x M) of survival probabilities.
surv_curve = function(dLambda, grid = NULL) {
  n = nrow(dLambda); M = ncol(dLambda)
  if (!is.null(grid) && !grid$discrete) {
    # Exp-trapezoid: S(u_k) = exp(-int_0^{u_k} lambda du).
    # Pairwise cumulative: cumhaz(u_1) = 0,
    #   cumhaz(u_k) = cumhaz(u_{k-1}) + (dLambda_{k-1} + dLambda_k) / 2.
    # (dLambda_k = lambda_k * du, so the /2 is the trapezoid average.)
    cum = matrix(0, n, M)
    for (k in 2:M)
      cum[, k] = cum[, k - 1] + (dLambda[, k - 1] + dLambda[, k]) / 2
    return(exp(-cum))
  }
  S = matrix(1, n, M)
  S[, 1] = 1 - dLambda[, 1]
  for (k in 2:M) S[, k] = S[, k - 1] * (1 - dLambda[, k])
  S
}

#' Survival probability at horizon: S(tau).
#'
#' Discrete: prod(1 - dLambda_k). Continuous: exp(-sum w_k lambda_k)
#' where w_k are the grid's quadrature weights and lambda_k = dLambda_k / du.
#'
#' @param dLambda Matrix (n x M).
#' @param grid Grid object. If NULL, uses discrete.
#' @return Vector (n).
surv_prob = function(dLambda, grid = NULL) {
  n = nrow(dLambda); M = ncol(dLambda)
  if (!is.null(grid) && !grid$discrete) {
    # int lambda du = sum w_k * lambda_k = sum (w_k / du) * dLambda_k
    return(as.vector(exp(-(dLambda %*% (grid$weights / grid$du)))))
  }
  n_over = sum(dLambda >= 1)
  if (n_over > 0)
    stop(structure(
      list(message  = sprintf("dLambda >= 1 at %d entries; dispersion unbounded (use softplus or entropy)", n_over),
           n_over   = n_over),
      class = c("dLambda_overflow", "error", "condition")))
  S = rep(1, n)
  for (k in 1:M) S = S * (1 - dLambda[, k])
  S
}

# ============================================================
# Derivatives of survival functionals
# ============================================================

#' Gateaux derivative of S(tau) w.r.t. dLambda_k.
#'
#' Discrete: -S(tau) / (1 - dLambda_k). Chain rule of product.
#' Continuous: -S(tau). Derivative of exp(-integral) is just -S.
#'
#' @param dLambda Matrix (n x M).
#' @param grid Grid object. If NULL, uses discrete.
#' @return Function(k) returning n-vector.
surv_prob_dot = function(dLambda, grid = NULL) {
  S_tau = surv_prob(dLambda, grid)
  if (!is.null(grid) && !grid$discrete) {
    return(function(k) -S_tau)
  }
  function(k) -S_tau / (1 - dLambda[, k])
}

# ============================================================
# Integration against Lebesgue measure
# ============================================================

#' Integral against Lebesgue measure: sum_k f_k * w_k.
#' Uses the grid's quadrature weights.
#'
#' @param f Matrix (n x M).
#' @param grid Grid object.
#' @return Vector (n).
integrate_du = function(f, grid) {
  as.vector(f %*% grid$weights)
}

# ============================================================
# Integration against dLambda (compensator, dpsi targets)
# ============================================================

#' Truncated integral: sum_{k: u_k <= T_i} f_k * dLambda_k.
#'
#' @param f Matrix (n x M).
#' @param dLambda Matrix (n x M).
#' @param grid Grid object.
#' @param upper Vector (n) of upper limits.
#' @return Vector (n).
integrate_dLambda = function(f, dLambda, grid, upper) {
  n = nrow(f); M = ncol(f)
  result = numeric(n)
  for (k in 1:M) {
    ar = which(upper >= grid$points[k])
    if (length(ar) == 0) break
    result[ar] = result[ar] + f[ar, k] * dLambda[ar, k]
  }
  result
}

# ============================================================
# RMST: int_0^tau S(u) du
# ============================================================

#' RMST from dLambda: int_0^tau S(u) du.
#' Uses grid's quadrature weights for the outer integral.
#'
#' @param dLambda Matrix (n x M).
#' @param grid Grid object.
#' @return Vector (n).
rmst = function(dLambda, grid) {
  S = surv_curve(dLambda, grid)
  integrate_du(S, grid)
}

#' Gateaux derivative of RMST w.r.t. dLambda_k.
#'
#' Discrete: -tail_k / (1 - dLambda_k) where tail_k = sum_{j>=k} S_j * du.
#' Continuous: -(total - cum_k) via trapezoid cumulative integral of S.
#'
#' @param dLambda Matrix (n x M).
#' @param grid Grid object.
#' @return Function(k) returning n-vector.
rmst_dot = function(dLambda, grid) {
  S = surv_curve(dLambda, grid)
  M = ncol(dLambda); du = grid$du
  if (!grid$discrete) {
    # Trapezoid cumulative integral of S
    cum = matrix(0, nrow(dLambda), M)
    for (k in 2:M)
      cum[, k] = cum[, k - 1] + (S[, k - 1] + S[, k]) / 2 * du
    total = cum[, M]
    return(function(k) -(total - cum[, k]))
  }
  # Discrete: rectangle tail integral
  tail = matrix(0, nrow(dLambda), M)
  tail[, M] = S[, M] * du
  for (k in (M - 1):1) tail[, k] = tail[, k + 1] + S[, k] * du
  function(k) -tail[, k] / (1 - dLambda[, k])
}

# ============================================================
# Materialize dLambda matrix from hazard function
# ============================================================

#' Build n x M matrix of hazard increments from h_fn(k).
#'
#' @param h_fn Function(k) returning n-vector of dLambda at grid point k.
#' @param M Number of grid points.
#' @return Matrix (n x M).
materialize_dLambda = function(h_fn, M) {
  n = length(h_fn(1))
  dL = matrix(0, n, M)
  for (k in 1:M) dL[, k] = h_fn(k)
  dL
}
