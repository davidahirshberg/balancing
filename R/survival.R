## Survival-specific machinery.
##
## Training mesh (Rao-Blackwellization), compensator integral,
## and survival-specific DR estimation.

# ============================================================
# Training mesh
# ============================================================

#' Build training mesh for hazard estimation.
#' Pools (time, subject) pairs into a single regression problem.
#'
#' @param T_obs Observed times (n).
#' @param D Event indicators (n).
#' @param Z Covariate matrix with treatment in column 1 (n x p).
#' @param horizon Maximum follow-up time.
#' @param M Number of time bins.
#' @return List with Y (n x M event indicators), G (n x M at-risk indicators),
#'   Z, mesh (midpoints), du (bin width), M.
build_training_mesh = function(T_obs, D, Z, horizon, M, discrete = FALSE) {
  n = length(T_obs)
  du = horizon / M
  mesh = if (discrete) (1:M) * du else ((1:M) - 0.5) * du

  Y = matrix(0, n, M)
  G = matrix(0, n, M)
  for (k in 1:M) {
    u_lo = (k - 1) * du; u_hi = k * du
    G[, k] = as.numeric(T_obs > u_lo)
    Y[, k] = as.numeric(T_obs > u_lo & T_obs <= u_hi & D == 1)
  }

  list(Y = Y, G = G, Z = Z, mesh = mesh, du = du, M = M)
}

#' Pool training mesh into (Z_pool, Y_pool) for kernel regression.
#' Extracts at-risk (u, Z) pairs and corresponding outcomes.
at_risk_pairs = function(train) {
  rows = which(train$G == 1, arr.ind = TRUE)
  Z_pool = cbind(u = train$mesh[rows[, "col"]], train$Z[rows[, "row"], , drop = FALSE])
  Y_pool = train$Y[cbind(rows[, "row"], rows[, "col"])]
  list(Z_pool = Z_pool, Y_pool = Y_pool, i_pool = rows[, "row"],
       u_pool = train$mesh[rows[, "col"]], du_bin = train$du)
}

# ============================================================
# Compensator: int_0^{T_i} gamma(u) dLambda(u)
# ============================================================

#' Compensator integral against the cumulative hazard measure.
#'
#' Discrete: sum_{k: T_i >= u_k} gamma_k * h_k (hazard probability).
#' Continuous: Simpson's rule on gamma * lambda with per-person
#'   endpoint corrections and stub at T_i.
#'
#' @param lambda_fit_ Bound subjects for lambda model.
#' @param gamma_fit_ Bound subjects for gamma model.
#' @param T_obs Observed times.
#' @param Z Covariate matrix (treatment in col 1).
#' @param comp_grid Grid object for compensator evaluation.
#' @param lambda_dispersion Dispersion for lambda model.
#' @param dchis_gamma Function(Z, phi) -> gamma. Z = (u, W, X) matrix.
#' @param du_bin Training mesh bin width.
#' @return List with comp (compensator values) and noise_var.
compensator = function(lambda_fit_, gamma_fit_, T_obs, Z,
                       comp_grid, lambda_dispersion, dchis_gamma, du_bin) {
  dchis_lam = \(p) .dchis(lambda_dispersion, p)
  n = length(T_obs)
  T_eff = pmin(T_obs, comp_grid$horizon)

  if (comp_grid$discrete) {
    # Discrete: simple sum of gamma * dLambda (hazard probability)
    comp = numeric(n)
    noise_var = numeric(n)
    for (k in 1:comp_grid$M) {
      u = comp_grid$points[k]
      ar = which(T_eff >= u)
      if (length(ar) == 0) break
      Z_ar = cbind(u, Z[ar, , drop = FALSE])
      phi_lam = \(Z) .phi(lambda_fit_[ar], Z)
      phi_gam = \(Z) .phi(gamma_fit_[ar], Z)
      dLambda_ar = dchis_lam(phi_lam(Z_ar))
      gamma_ar   = dchis_gamma(Z_ar, phi_gam(Z_ar))
      comp[ar] = comp[ar] + gamma_ar * dLambda_ar
      noise_var[ar] = noise_var[ar] + gamma_ar^2 * dLambda_ar
    }
    return(list(comp = comp, noise_var = noise_var))
  }

  # Continuous: composite Simpson's rule with per-person endpoint corrections.
  # Integrates gamma(u) * lambda(u) du over [0, min(T_i, horizon)].
  Q = comp_grid$M - 1
  h = comp_grid$du
  if (Q %% 2 != 0) stop("Q_comp must be even for Simpson's rule")

  fi = findInterval(T_eff, comp_grid$points)
  k_max = pmin(fi - 1, Q)
  frac = T_eff - comp_grid$points[fi]
  n_simp = (k_max %/% 2) * 2

  # Full Simpson weights
  sw = numeric(Q + 1)
  sw[1] = h / 3; sw[Q + 1] = h / 3
  if (Q >= 2) for (j in 2:Q) sw[j] = if (j %% 2 == 0) 4 * h / 3 else 2 * h / 3

  comp = numeric(n)
  noise_var = numeric(n)

  for (k in 0:Q) {
    u = comp_grid$points[k + 1]
    ar = which(k_max >= k)
    if (length(ar) == 0) break

    Z_ar = cbind(u, Z[ar, , drop = FALSE])
    phi_lam = \(Z) .phi(lambda_fit_[ar], Z)
    phi_gam = \(Z) .phi(gamma_fit_[ar], Z)
    lambdaZ_ar = pmax(dchis_lam(phi_lam(Z_ar)) / du_bin, 0)
    gammaZ_ar  = dchis_gamma(Z_ar, phi_gam(Z_ar))

    # Per-person weight at this grid point
    w = rep(sw[k + 1], length(ar))
    km_ar = k_max[ar]; ns_ar = n_simp[ar]; fr_ar = frac[ar]

    # Simpson endpoint adjustments
    at_se = (k == ns_ar) & (ns_ar >= 2) & (k > 0)
    w[at_se] = h / 3
    at_se_lo = at_se & (km_ar > ns_ar)
    w[at_se_lo] = w[at_se_lo] + h / 2
    at_se_stub = at_se & (km_ar == ns_ar) & (fr_ar > 1e-12)
    w[at_se_stub] = w[at_se_stub] + fr_ar[at_se_stub] / 2

    if (k == 0) {
      no_simp = ns_ar < 2
      has_leftover = no_simp & (km_ar > 0)
      w[has_leftover] = h / 2
      only_stub = no_simp & (km_ar == 0) & (fr_ar > 1e-12)
      w[only_stub] = fr_ar[only_stub] / 2
      dead = no_simp & (km_ar == 0) & (fr_ar <= 1e-12)
      w[dead] = 0
    }

    at_lo = (k == km_ar) & (k == ns_ar + 1) & (k > 0)
    w[at_lo] = h / 2
    at_lo_stub = at_lo & (fr_ar > 1e-12)
    w[at_lo_stub] = w[at_lo_stub] + fr_ar[at_lo_stub] / 2

    comp[ar] = comp[ar] + gammaZ_ar * lambdaZ_ar * w
    noise_var[ar] = noise_var[ar] + gammaZ_ar^2 * lambdaZ_ar * w
  }

  # Stub endpoints at T_eff
  has_stub = which(frac > 1e-12 & T_eff > 0)
  if (length(has_stub) > 0) {
    Z_stub = cbind(T_eff[has_stub], Z[has_stub, , drop = FALSE])
    phi_lam_stub = \(Z) .phi(lambda_fit_[has_stub], Z)
    phi_gam_stub = \(Z) .phi(gamma_fit_[has_stub], Z)
    lambdaZ_stub = pmax(dchis_lam(phi_lam_stub(Z_stub)) / du_bin, 0)
    gammaZ_stub  = dchis_gamma(Z_stub, phi_gam_stub(Z_stub))
    comp[has_stub] = comp[has_stub] + gammaZ_stub * lambdaZ_stub * frac[has_stub] / 2
    noise_var[has_stub] = noise_var[has_stub] + gammaZ_stub^2 * lambdaZ_stub * frac[has_stub] / 2
  }

  list(comp = comp, noise_var = noise_var)
}
