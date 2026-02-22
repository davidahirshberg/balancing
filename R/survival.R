## Survival-specific machinery.
##
## Training mesh (Rao-Blackwellization), compensator integral (Simpson's rule),
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
build_training_mesh = function(T_obs, D, Z, horizon, M) {
  n = length(T_obs)
  du = horizon / M
  mesh = ((1:M) - 0.5) * du

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
pool_mesh = function(train) {
  rows = which(train$G == 1, arr.ind = TRUE)
  Z_pool = cbind(u = train$mesh[rows[, "col"]], train$Z[rows[, "row"], , drop = FALSE])
  Y_pool = train$Y[cbind(rows[, "row"], rows[, "col"])]
  list(Z_pool = Z_pool, Y_pool = Y_pool, i_pool = rows[, "row"],
       u_pool = train$mesh[rows[, "col"]], du_bin = train$du)
}

# ============================================================
# Hazard function helpers
# ============================================================

#' Build a hazard evaluation function from a fitted model.
#' Returns haz_fn(k) giving n-vector of hazard rates at eval_grid[k].
#' Also returns du (grid spacing).
make_haz_fn = function(model, Z, horizon, Q, dispersion, lam_bound = NULL, du_bin = NULL) {
  z = atleast_2d(Z)
  du = horizon / Q
  eval_grid = (0:Q) * du

  if (!is.null(lam_bound)) {
    haz_fn = function(k) pmax(dispersion$dchis(eval_phi(lam_bound, eval_grid[k])) / du_bin, 0)
  } else {
    haz_fn = function(k) {
      phi = predict_phi(model, cbind(eval_grid[k], z))
      pmax(dispersion$dchis(phi) / du_bin, 0)
    }
  }
  list(haz_fn = haz_fn, du = du)
}

# ============================================================
# Compensator integral (Simpson's rule)
# ============================================================

#' Compensator integral: int_0^{min(T_i, horizon)} gamma(u, Z_i) lambda(u, Z_i) du.
#' Composite Simpson's rule with endpoint correction.
#'
#' Uses bind_subjects for fast evaluation. Gamma is computed as
#' W_i * base_gam_dispersion$dchis(W_i * phi) for sign-flipped weights,
#' applied per-subset at each quadrature point.
#'
#' @param lam_bound Bound subjects for lambda model.
#' @param gam_bound Bound subjects for gamma model.
#' @param T_obs Observed times.
#' @param D Event indicators.
#' @param Z Covariate matrix (treatment in col 1).
#' @param horizon Maximum follow-up.
#' @param Q Number of quadrature points (even).
#' @param lam_dispersion Dispersion for lambda model.
#' @param gam_dispersion Base dispersion for gamma model (before sign-flip).
#' @param du_bin Bin width for hazard conversion.
#' @param W_eval Per-subject sign vector (2*W - 1), length n. If NULL, no sign-flip.
#' @return List with comp (compensator values) and noise_var (int gamma^2 lambda du).
compensator_simpson = function(lam_bound, gam_bound, T_obs, D, Z, horizon, Q,
                               lam_dispersion, gam_dispersion, du_bin,
                               W_eval = NULL) {
  n = length(T_obs)
  if (Q %% 2 != 0) Q = Q + 1
  h = horizon / Q
  grid = (0:Q) * h

  T_eff = pmin(T_obs, horizon)
  fi = findInterval(T_eff, grid)
  k_max = pmin(fi - 1, Q)
  frac = T_eff - grid[fi]
  n_simp = (k_max %/% 2) * 2

  # Simpson weights for full grid
  sw = numeric(Q + 1)
  sw[1] = h / 3; sw[Q + 1] = h / 3
  if (Q >= 2) for (j in 2:Q) sw[j] = if (j %% 2 == 0) 4 * h / 3 else 2 * h / 3

  comp = numeric(n)
  noise_var = numeric(n)

  for (k in 0:Q) {
    u = grid[k + 1]
    ar = which(k_max >= k)
    if (length(ar) == 0) break

    # Lambda and gamma at (u, Z[ar])
    lam_phi = eval_phi(lam_bound, u, ar)
    lam_vals = pmax(lam_dispersion$dchis(lam_phi) / du_bin, 0)

    gam_phi = eval_phi(gam_bound, u, ar)
    if (!is.null(W_eval)) {
      W_ar = W_eval[ar]
      gam_vals = W_ar * gam_dispersion$dchis(W_ar * gam_phi)
    } else {
      gam_vals = gam_dispersion$dchis(gam_phi)
    }

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

    comp[ar] = comp[ar] + gam_vals * lam_vals * w
    noise_var[ar] = noise_var[ar] + gam_vals^2 * lam_vals * w
  }

  # Stub endpoints at T_eff
  has_stub = which(frac > 1e-12 & T_eff > 0)
  if (length(has_stub) > 0) {
    lam_phi_stub = eval_phi_vec(lam_bound, T_eff[has_stub], has_stub)
    lam_stub = pmax(lam_dispersion$dchis(lam_phi_stub) / du_bin, 0)

    gam_phi_stub = eval_phi_vec(gam_bound, T_eff[has_stub], has_stub)
    if (!is.null(W_eval)) {
      W_stub = W_eval[has_stub]
      gam_stub = W_stub * gam_dispersion$dchis(W_stub * gam_phi_stub)
    } else {
      gam_stub = gam_dispersion$dchis(gam_phi_stub)
    }

    comp[has_stub] = comp[has_stub] + gam_stub * lam_stub * frac[has_stub] / 2
    noise_var[has_stub] = noise_var[has_stub] + gam_stub^2 * lam_stub * frac[has_stub] / 2
  }

  list(comp = comp, noise_var = noise_var)
}
