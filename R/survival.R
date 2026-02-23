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
#' Uses bind_subjects for fast evaluation. Gamma = gam_dispersion$dchis(phi).
#' The caller passes the appropriate dispersion (already sign-flipped if needed).
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
#' @param W_eval Per-subject sign vector (2*W - 1). If NULL, no sign-flip.
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

# ============================================================
# Discrete compensator
# ============================================================

#' Discrete compensator: sum_{k: T_i >= u_k} gamma(u_k, Z_i) * h_k(Z_i).
#'
#' For discrete-time data, the compensator is a sum over mesh points,
#' not an integral. h_k are hazard probabilities (not rates), so no du factor.
#'
#' @param lam_bound Bound subjects for lambda model.
#' @param gam_bound Bound subjects for gamma model.
#' @param T_obs Observed times (discrete, on the mesh grid).
#' @param D Event indicators.
#' @param mesh Training mesh midpoints (length M).
#' @param lam_dispersion Dispersion for lambda model.
#' @param gam_dispersion Base dispersion for gamma model (before sign-flip).
#' @param W_eval Per-subject sign vector (2*W - 1). If NULL, no sign-flip.
#' @return List with comp (compensator values) and noise_var.
compensator_discrete = function(lam_bound, gam_bound, T_obs, D, mesh,
                                lam_dispersion, gam_dispersion,
                                W_eval = NULL) {
  n = length(T_obs)
  M = length(mesh)
  comp = rep(0, n)
  noise_var = rep(0, n)

  for (k in 1:M) {
    u = mesh[k]
    ar = which(T_obs >= u)
    if (length(ar) == 0) next

    # Hazard probability (not rate — no du division)
    lam_phi = eval_phi(lam_bound, u, ar)
    h_vals = lam_dispersion$dchis(lam_phi)

    # Gamma
    gam_phi = eval_phi(gam_bound, u, ar)
    if (!is.null(W_eval)) {
      W_ar = W_eval[ar]
      gam_vals = W_ar * gam_dispersion$dchis(W_ar * gam_phi)
    } else {
      gam_vals = gam_dispersion$dchis(gam_phi)
    }

    comp[ar] = comp[ar] + gam_vals * h_vals
    noise_var[ar] = noise_var[ar] + gam_vals^2 * h_vals
  }

  list(comp = comp, noise_var = noise_var)
}

# ============================================================
# Bootstrap precomputation
# ============================================================

#' Stack a list of kernel matrices for batched BLAS multiplies.
.stack_K_list = function(K_list) do.call(rbind, K_list)

#' Build haz_fn from (alpha, beta0) + precomputed stacked kernel matrices.
#' One BLAS multiply for all Q+1 grid points. Returns haz_fn(k) compatible
#' with estimand$psi / estimand$dot_psi.
#'
#' @param alpha Kernel coefficients.
#' @param beta0 Intercept.
#' @param K_list List of Q+1 kernel matrices (n_eval x n_centers).
#' @param du_bin Training mesh bin width (for hazard conversion).
#' @param dispersion Dispersion object (dchis maps phi -> hazard increment).
#' @param K_stacked Pre-stacked matrices (optional, computed if NULL).
#' @return List with haz_fn(k) and du = du_comp.
make_haz_fn_alpha = function(alpha, beta0, K_list, du_comp, du_bin,
                             dispersion, K_stacked = NULL) {
  if (is.null(K_stacked)) K_stacked = .stack_K_list(K_list)
  n = nrow(K_list[[1]]); Q1 = length(K_list)
  phi_all = as.vector(K_stacked %*% alpha) + beta0
  haz_mat = matrix(pmax(dispersion$dchis(phi_all) / du_bin, 0), nrow = n, ncol = Q1)
  list(haz_fn = function(k) haz_mat[, k], du = du_comp)
}

#' Precompute kernel matrices for bootstrap refitting.
#'
#' Avoids recomputing O(n^2) kernel matrices per bootstrap rep.
#' Returns everything needed for one-step Newton refits of lambda and gamma,
#' plus prediction infrastructure for DR terms.
#'
#' @param lam_fit Fitted lambda model (kernel_bregman object).
#' @param gam_fit Fitted gamma model (kernel_bregman object).
#' @param eval_Z Evaluation subjects' covariate matrix (treatment in col 1).
#' @param T_obs Observed times for eval subjects.
#' @param D Event indicators for eval subjects.
#' @param horizon Maximum follow-up time.
#' @param Q_comp Quadrature points for compensator.
#' @param kern Kernel object (same for lambda and gamma).
#' @param du_bin Training mesh bin width.
#' @param W_eval Per-subject sign vector (2*W - 1).
#' @param eval_train Training mesh for eval fold (for gamma pool mapping).
#' @param lam_train Training mesh for lambda fold (for lambda pool mapping).
#' @return Precomputed structure for bootstrap_surv_fold.
precompute_boot_surv = function(lam_fit, gam_fit, eval_Z, T_obs, D,
                                horizon, Q_comp, kern, du_bin,
                                W_eval, eval_train, lam_train) {
  n_eval = nrow(atleast_2d(eval_Z))
  du_comp = horizon / Q_comp
  comp_grid = (0:Q_comp) * du_comp

  # === Lambda precomputation ===
  # Covariate-only squared distances between eval subjects and lambda centers.
  # Lambda centers are (u, Z) where Z includes treatment; split off u.
  Z_pool_lam = lam_fit$Z
  u_centers_lam = Z_pool_lam[, 1]
  cov_centers_lam = Z_pool_lam[, -1, drop = FALSE]

  Z1 = with.treatment(eval_Z, 1)
  Z0 = with.treatment(eval_Z, 0)
  D2_cov_surv1 = .fast_sqdist(Z1, cov_centers_lam)
  D2_cov_surv0 = .fast_sqdist(Z0, cov_centers_lam)
  D2_cov_obs = .fast_sqdist(eval_Z, cov_centers_lam)

  # Build K_list for each setting at each compensator grid point.
  # K[i,j] = k((u_k, Z_i), (u_j_center, Z_j_center)).
  K_lam_surv1 = K_lam_surv0 = K_lam_obs = vector("list", Q_comp + 1)
  for (k in 0:Q_comp) {
    delta_u2 = (comp_grid[k + 1] - u_centers_lam)^2
    D2_1 = sweep(D2_cov_surv1, 2, delta_u2, "+"); D2_1[D2_1 < 0] = 0
    D2_0 = sweep(D2_cov_surv0, 2, delta_u2, "+"); D2_0[D2_0 < 0] = 0
    D2_o = sweep(D2_cov_obs, 2, delta_u2, "+"); D2_o[D2_o < 0] = 0
    K_lam_surv1[[k + 1]] = .kernel_from_D2(D2_1, kern)
    K_lam_surv0[[k + 1]] = .kernel_from_D2(D2_0, kern)
    K_lam_obs[[k + 1]] = .kernel_from_D2(D2_o, kern)
  }
  K_stacked_surv1 = .stack_K_list(K_lam_surv1)
  K_stacked_surv0 = .stack_K_list(K_lam_surv0)
  K_stacked_obs = .stack_K_list(K_lam_obs)

  # === Gamma precomputation ===
  Z_pool_gam = gam_fit$Z
  T_eff = pmin(T_obs, horizon)

  # Gamma prediction kernel matrices at compensator grid
  K_pred_list = list(); at_risk_list = list(); W_ar_list = list()
  for (k in 0:Q_comp) {
    u = comp_grid[k + 1]
    ar = which(T_eff >= u)
    at_risk_list[[k + 1]] = ar
    if (length(ar) == 0) {
      K_pred_list[[k + 1]] = NULL; W_ar_list[[k + 1]] = NULL; next
    }
    z_aug = cbind(u, eval_Z[ar, , drop = FALSE])
    K_pred_list[[k + 1]] = kernel_matrix(z_aug, Z_pool_gam, kern)
    W_ar_list[[k + 1]] = W_eval[ar]
  }

  # Dirac prediction at event times
  event_idx = which(D == 1 & T_obs <= horizon)
  if (length(event_idx) > 0) {
    z_ev = cbind(T_obs[event_idx], eval_Z[event_idx, , drop = FALSE])
    K_pred_dirac = kernel_matrix(z_ev, Z_pool_gam, kern)
    W_dirac = W_eval[event_idx]
  } else {
    K_pred_dirac = NULL; W_dirac = NULL
  }

  # === Training data for onestep refits ===
  # Lambda pool mapping: which person each pooled observation came from
  lam_rows = which(lam_train$G == 1, arr.ind = TRUE)
  Y_pool_lam = lam_train$Y[cbind(lam_rows[, "row"], lam_rows[, "col"])]
  i_pool_lam = lam_rows[, "row"]

  # Gamma pool mapping: person index + mesh column for derivative lookup
  gam_rows = which(eval_train$G == 1, arr.ind = TRUE)
  i_pool_gam = gam_rows[, "row"]
  gam_mesh_u = eval_train$mesh
  gam_mesh_col = gam_rows[, "col"]
  W_pool_gam = 2 * eval_train$Z[gam_rows[, "row"], 1] - 1

  list(
    # Lambda prediction
    K_lam_surv1 = K_lam_surv1, K_stacked_surv1 = K_stacked_surv1,
    K_lam_surv0 = K_lam_surv0, K_stacked_surv0 = K_stacked_surv0,
    K_lam_obs = K_lam_obs, K_stacked_obs = K_stacked_obs,
    # Lambda training
    Y_pool_lam = Y_pool_lam, i_pool_lam = i_pool_lam,
    # Gamma prediction
    K_pred_list = K_pred_list, at_risk_list = at_risk_list, W_ar_list = W_ar_list,
    K_pred_dirac = K_pred_dirac, W_dirac = W_dirac, event_idx = event_idx,
    # Gamma training
    i_pool_gam = i_pool_gam, gam_mesh_u = gam_mesh_u,
    gam_mesh_col = gam_mesh_col, W_pool_gam = W_pool_gam,
    # Shared
    n_eval = n_eval, du_bin = du_bin, du_comp = du_comp, comp_grid = comp_grid
  )
}

#' DR terms from precomputed gamma alpha + kernel matrices.
#'
#' Bootstrap-fast alternative to compensator_simpson: uses precomputed
#' K_pred_list for gamma, precomputed lam_vals for lambda.
#' Rectangle rule (comp_grid spacing is fine enough).
#'
#' @param gam_alpha Gamma kernel coefficients (from onestep_bregman).
#' @param gam_beta0 Gamma intercept.
#' @param direct Per-subject direct term (from bootstrap lambda).
#' @param precomp Precomputed structure from precompute_boot_surv.
#' @param lam_vals List of at-risk lambda values at each grid point
#'   (lam_vals[[k]][j] = hazard for at-risk subject j at grid point k).
#' @param gam_dispersion Base dispersion for gamma (before sign-flip).
#' @return List with terms (DR estimate per subject) and noise_var.
dr_from_alpha = function(gam_alpha, gam_beta0, direct, precomp,
                         lam_vals, gam_dispersion) {
  n = precomp$n_eval
  du = precomp$du_comp
  compensator = rep(0, n)
  noise_var = rep(0, n)

  for (k in seq_along(precomp$K_pred_list)) {
    ar = precomp$at_risk_list[[k]]
    if (length(ar) == 0) next
    phi_k = as.vector(precomp$K_pred_list[[k]] %*% gam_alpha) + gam_beta0
    W_k = precomp$W_ar_list[[k]]
    if (!is.null(W_k)) {
      gamma_k = W_k * gam_dispersion$dchis(W_k * phi_k)
    } else {
      gamma_k = gam_dispersion$dchis(phi_k)
    }
    compensator[ar] = compensator[ar] + gamma_k * lam_vals[[k]] * du
    noise_var[ar] = noise_var[ar] + gamma_k^2 * lam_vals[[k]] * du
  }

  dirac = rep(0, n)
  if (length(precomp$event_idx) > 0) {
    phi_ev = as.vector(precomp$K_pred_dirac %*% gam_alpha) + gam_beta0
    W_ev = precomp$W_dirac
    if (!is.null(W_ev)) {
      dirac[precomp$event_idx] = W_ev * gam_dispersion$dchis(W_ev * phi_ev)
    } else {
      dirac[precomp$event_idx] = gam_dispersion$dchis(phi_ev)
    }
  }

  terms = direct + dirac - compensator
  list(terms = terms, noise_var = noise_var)
}
