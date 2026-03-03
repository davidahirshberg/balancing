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
# Hazard function helpers
# ============================================================

#' Build a hazard evaluation function from a fitted model.
#' Returns haz_fn(k) giving n-vector of hazard rates at eval_grid[k].
#' Also returns du (grid spacing).
make_haz_fn = function(model, Z, horizon, Q, dispersion, lambda_fit_ = NULL, du_bin = NULL) {
  z = atleast_2d(Z)
  du = horizon / Q
  eval_grid = (0:Q) * du

  dchis = \(p) .dchis(dispersion, p)
  if (!is.null(lambda_fit_)) {
    phi_lam = \(Z) .phi(lambda_fit_, Z)
    haz_fn = function(k) pmax(dchis(phi_lam(cbind(eval_grid[k], z))) / du_bin, 0)
  } else {
    phi_mod = \(Z) .phi(model, Z)
    haz_fn = function(k) pmax(dchis(phi_mod(cbind(eval_grid[k], z))) / du_bin, 0)
  }
  list(haz_fn = haz_fn, du = du)
}

# ============================================================
# Compensator integral (Simpson's rule)
# ============================================================

#' Compensator integral: int_0^{min(T_i, horizon)} gamma(u, Z_i) lambda(u, Z_i) du.
#' Composite Simpson's rule with endpoint correction.
#'
#' Gamma = .dchis(lambda_dispersion, phi).
#' The caller passes the appropriate dispersion (already sign-flipped if needed).
#'
#' @param lambda_fit_ Bound subjects for lambda model.
#' @param gamma_fit_ Bound subjects for gamma model.
#' @param T_obs Observed times.
#' @param D Event indicators.
#' @param Z Covariate matrix (treatment in col 1).
#' @param horizon Maximum follow-up.
#' @param Q Number of quadrature points (even).
#' @param lambda_dispersion Dispersion for lambda model.
#' @param dchis_gamma Function(Z, phi) -> gamma. Z = (u, W, X) matrix.
#' @param du_bin Bin width for hazard conversion.
#' @return List with comp (compensator values) and noise_var (int gamma^2 lambda du).
compensator_simpson = function(lambda_fit_, gamma_fit_, T_obs, D, Z, horizon, Q,
                               lambda_dispersion, dchis_gamma, du_bin) {
  dchis_lam = \(p) .dchis(lambda_dispersion, p)
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

# ============================================================
# Discrete compensator
# ============================================================

#' Discrete compensator: sum_{k: T_i >= u_k} gamma(u_k, Z_i) * h_k(Z_i).
#'
#' For discrete-time data, the compensator is a sum over mesh points,
#' not an integral. h_k are hazard probabilities (not rates), so no du factor.
#'
#' @param lambda_fit_ Bound subjects for lambda model.
#' @param gamma_fit_ Bound subjects for gamma model.
#' @param T_obs Observed times (discrete, on the mesh grid).
#' @param D Event indicators.
#' @param mesh Training mesh midpoints (length M).
#' @param lambda_dispersion Dispersion for lambda model.
#' @param dchis_gamma Function(Z, phi) -> gamma. Z = (u, W, X) matrix.
#' @param Z Covariate matrix (W, X) — no time column.
#' @return List with comp (compensator values) and noise_var.
compensator_discrete = function(lambda_fit_, gamma_fit_, T_obs, D, mesh,
                                lambda_dispersion, dchis_gamma, Z) {
  dchis_lam = \(p) .dchis(lambda_dispersion, p)
  n = length(T_obs)
  M = length(mesh)
  comp = rep(0, n)
  noise_var = rep(0, n)

  for (k in 1:M) {
    u = mesh[k]
    ar = which(T_obs >= u)
    if (length(ar) == 0) next

    Z_ar = cbind(u, Z[ar, , drop = FALSE])
    phi_lam = \(Z) .phi(lambda_fit_[ar], Z)
    phi_gam = \(Z) .phi(gamma_fit_[ar], Z)
    lambdaZ_ar = dchis_lam(phi_lam(Z_ar))
    gammaZ_ar  = dchis_gamma(Z_ar, phi_gam(Z_ar))

    comp[ar] = comp[ar] + gammaZ_ar * lambdaZ_ar
    noise_var[ar] = noise_var[ar] + gammaZ_ar^2 * lambdaZ_ar
  }

  list(comp = comp, noise_var = noise_var)
}

# ============================================================
# Bootstrap precomputation
# ============================================================

#' Stack a list of kernel matrices for batched BLAS multiplies.
.stack_K_list = function(K_list) do.call(rbind, K_list)

#' Build haz_fn from (alpha, beta) + precomputed stacked kernel matrices.
#' One BLAS multiply for all Q+1 grid points. Returns haz_fn(k) compatible
#' with estimand$psi / estimand$dot_psi.
#'
#' @param alpha Kernel coefficients.
#' @param beta Null-space coefficients (d-vector).
#' @param K_list List of Q+1 kernel matrices (n_eval x n_centers).
#' @param B_stacked Stacked null-space basis matrices (n_eval*(Q+1) x d).
#' @param du_bin Training mesh bin width (for hazard conversion).
#' @param dispersion Dispersion object (dchis maps phi -> hazard increment).
#' @param K_stacked Pre-stacked matrices (optional, computed if NULL).
#' @return List with haz_fn(k) and du = du_comp.
make_haz_fn_alpha = function(alpha, beta, K_list, B_stacked, du_comp, du_bin,
                             dispersion, K_stacked = NULL) {
  if (is.null(K_stacked)) K_stacked = .stack_K_list(K_list)
  n = nrow(K_list[[1]]); Q1 = length(K_list)
  phi_all = as.vector(K_stacked %*% alpha)
  if (length(beta) > 0) phi_all = phi_all + as.vector(B_stacked %*% beta)
  dchis = \(p) .dchis(dispersion, p)
  haz_mat = matrix(pmax(dchis(phi_all) / du_bin, 0), nrow = n, ncol = Q1)
  list(haz_fn = function(k) haz_mat[, k], du = du_comp)
}

#' Precompute kernel matrices for bootstrap refitting.
#'
#' Avoids recomputing O(n^2) kernel matrices per bootstrap rep.
#' Returns everything needed for one-step Newton refits of lambda and gamma,
#' plus prediction infrastructure for DR terms.
#'
#' @param lambda_fit Fitted lambda model (kernel_bregman object).
#' @param gamma_fit Fitted gamma model (kernel_bregman object).
#' @param eval_Z Evaluation subjects' covariate matrix (treatment in col 1).
#' @param T_obs Observed times for eval subjects.
#' @param D Event indicators for eval subjects.
#' @param horizon Maximum follow-up time.
#' @param Q_comp Quadrature points for compensator.
#' @param kern Kernel object (same for lambda and gamma).
#' @param du_bin Training mesh bin width.
#' @param eval_train Training mesh for eval fold (for gamma pool mapping).
#' @param lambda_train Training mesh for lambda fold (for lambda pool mapping).
#' @return Precomputed structure for bootstrap_surv_fold.
precompute_boot_surv = function(lambda_fit, gamma_fit, eval_Z, T_obs, D,
                                horizon, Q_comp, kern, du_bin,
                                eval_train, lambda_train) {
  n_eval = nrow(atleast_2d(eval_Z))
  du_comp = horizon / Q_comp
  comp_grid = (0:Q_comp) * du_comp

  # === Lambda precomputation ===
  # Covariate-only squared distances between eval subjects and lambda centers.
  # Lambda centers are (u, Z) where Z includes treatment; split off u.
  Z_pool_lambda = lambda_fit$Z
  u_centers_lambda = Z_pool_lambda[, 1]
  cov_centers_lambda = Z_pool_lambda[, -1, drop = FALSE]

  Z1 = with.treatment(eval_Z, 1)
  Z0 = with.treatment(eval_Z, 0)
  D2_cov_surv1 = .fast_sqdist(Z1, cov_centers_lambda)
  D2_cov_surv0 = .fast_sqdist(Z0, cov_centers_lambda)
  D2_cov_obs = .fast_sqdist(eval_Z, cov_centers_lambda)

  # Build K_list for each setting at each compensator grid point.
  # K[i,j] = k((u_k, Z_i), (u_j_center, Z_j_center)).
  K_lambda_surv1 = K_lambda_surv0 = K_lambda_obs = vector("list", Q_comp + 1)
  for (k in 0:Q_comp) {
    delta_u2 = (comp_grid[k + 1] - u_centers_lambda)^2
    D2_1 = sweep(D2_cov_surv1, 2, delta_u2, "+"); D2_1[D2_1 < 0] = 0
    D2_0 = sweep(D2_cov_surv0, 2, delta_u2, "+"); D2_0[D2_0 < 0] = 0
    D2_o = sweep(D2_cov_obs, 2, delta_u2, "+"); D2_o[D2_o < 0] = 0
    K_lambda_surv1[[k + 1]] = .kernel_from_D2(D2_1, kern)
    K_lambda_surv0[[k + 1]] = .kernel_from_D2(D2_0, kern)
    K_lambda_obs[[k + 1]] = .kernel_from_D2(D2_o, kern)
  }
  K_stacked_surv1 = .stack_K_list(K_lambda_surv1)
  K_stacked_surv0 = .stack_K_list(K_lambda_surv0)
  K_stacked_obs = .stack_K_list(K_lambda_obs)

  # Lambda null-space basis at prediction points (same for all grid points —
  # null_basis depends on covariates, not time, for product kernels).
  B_surv1 = null_basis(cbind(0, Z1), kern)
  B_surv0 = null_basis(cbind(0, Z0), kern)
  B_obs = null_basis(cbind(0, eval_Z), kern)
  B_stacked_surv1 = do.call(rbind, rep(list(B_surv1), Q_comp + 1))
  B_stacked_surv0 = do.call(rbind, rep(list(B_surv0), Q_comp + 1))
  B_stacked_obs = do.call(rbind, rep(list(B_obs), Q_comp + 1))

  # === Gamma precomputation ===
  Z_pool_gamma = gamma_fit$Z
  T_eff = pmin(T_obs, horizon)

  # Gamma prediction kernel matrices and null-space bases at compensator grid
  K_pred_list = list(); B_pred_list = list()
  at_risk_list = list(); W_ar_list = list()
  for (k in 0:Q_comp) {
    u = comp_grid[k + 1]
    ar = which(T_eff >= u)
    at_risk_list[[k + 1]] = ar
    if (length(ar) == 0) {
      K_pred_list[[k + 1]] = NULL; B_pred_list[[k + 1]] = NULL
      W_ar_list[[k + 1]] = NULL; next
    }
    z_aug = cbind(u, eval_Z[ar, , drop = FALSE])
    K_pred_list[[k + 1]] = kernel_matrix(z_aug, Z_pool_gamma, kern)
    B_pred_list[[k + 1]] = null_basis(z_aug, kern)
    # Z at risk stored for dchis_gamma lookup (no longer need per-subject W)
  }

  # Dirac prediction at event times
  event_idx = which(D == 1 & T_obs <= horizon)
  if (length(event_idx) > 0) {
    z_ev = cbind(T_obs[event_idx], eval_Z[event_idx, , drop = FALSE])
    K_pred_dirac = kernel_matrix(z_ev, Z_pool_gamma, kern)
    B_pred_dirac = null_basis(z_ev, kern)
    # (W_dirac removed: dchis_gamma handles sign internally)
  } else {
    K_pred_dirac = NULL; B_pred_dirac = NULL; W_dirac = NULL
  }

  # === Training data for onestep refits ===
  # Lambda pool mapping: which person each pooled observation came from
  lambda_rows = which(lambda_train$G == 1, arr.ind = TRUE)
  Y_pool_lambda = lambda_train$Y[cbind(lambda_rows[, "row"], lambda_rows[, "col"])]
  i_pool_lambda = lambda_rows[, "row"]

  # Gamma pool mapping: person index + mesh column for derivative lookup
  gamma_rows = which(eval_train$G == 1, arr.ind = TRUE)
  i_pool_gamma = gamma_rows[, "row"]
  gamma_mesh_u = eval_train$mesh
  gamma_mesh_col = gamma_rows[, "col"]
  W_pool_gamma = 2 * eval_train$Z[gamma_rows[, "row"], 1] - 1

  list(
    # Lambda prediction
    K_lambda_surv1 = K_lambda_surv1, K_stacked_surv1 = K_stacked_surv1,
    K_lambda_surv0 = K_lambda_surv0, K_stacked_surv0 = K_stacked_surv0,
    K_lambda_obs = K_lambda_obs, K_stacked_obs = K_stacked_obs,
    B_stacked_surv1 = B_stacked_surv1, B_stacked_surv0 = B_stacked_surv0,
    B_stacked_obs = B_stacked_obs,
    # Lambda training
    Y_pool_lambda = Y_pool_lambda, i_pool_lambda = i_pool_lambda,
    # Gamma prediction
    K_pred_list = K_pred_list, B_pred_list = B_pred_list,
    at_risk_list = at_risk_list,
    K_pred_dirac = K_pred_dirac, B_pred_dirac = B_pred_dirac,
    event_idx = event_idx,
    # Gamma training
    i_pool_gamma = i_pool_gamma, gamma_mesh_u = gamma_mesh_u,
    gamma_mesh_col = gamma_mesh_col, W_pool_gamma = W_pool_gamma,
    # Shared
    n_eval = n_eval, Z = eval_Z, T_obs = T_obs,
    du_bin = du_bin, du_comp = du_comp, comp_grid = comp_grid
  )
}

#' DR terms from precomputed gamma alpha + kernel matrices.
#'
#' Bootstrap-fast alternative to compensator_simpson: uses precomputed
#' K_pred_list for gamma, precomputed lambdaZ_list for lambda.
#' Rectangle rule (comp_grid spacing is fine enough).
#'
#' @param gamma_alpha Gamma kernel coefficients (from onestep_bregman).
#' @param gamma_beta Gamma null-space coefficients (d-vector).
#' @param direct Per-subject direct term (from bootstrap lambda).
#' @param precomp Precomputed structure from precompute_boot_surv.
#' @param lambdaZ_list List of at-risk lambda values at each grid point
#'   (lambdaZ_list[[k]][j] = hazard for at-risk subject j at grid point k).
#' @param dchis_gamma Function(Z, phi) -> gamma. Z = (u, W, X) matrix.
#' @return List with terms (DR estimate per subject) and noise_var.
dr_from_alpha = function(gamma_alpha, gamma_beta, direct, precomp,
                         lambdaZ_list, dchis_gamma) {
  n = precomp$n_eval
  du = precomp$du_comp
  Z = precomp$Z
  compensator = rep(0, n)
  noise_var = rep(0, n)

  for (k in seq_along(precomp$K_pred_list)) {
    ar = precomp$at_risk_list[[k]]
    if (length(ar) == 0) next
    phi_k = as.vector(precomp$K_pred_list[[k]] %*% gamma_alpha)
    B_k = precomp$B_pred_list[[k]]
    if (!is.null(B_k) && length(gamma_beta) > 0) phi_k = phi_k + as.vector(B_k %*% gamma_beta)
    Z_ar = cbind(precomp$comp_grid[k], Z[ar, , drop = FALSE])
    gamma_k = dchis_gamma(Z_ar, phi_k)
    compensator[ar] = compensator[ar] + gamma_k * lambdaZ_list[[k]] * du
    noise_var[ar] = noise_var[ar] + gamma_k^2 * lambdaZ_list[[k]] * du
  }

  dirac = rep(0, n)
  if (length(precomp$event_idx) > 0) {
    phi_ev = as.vector(precomp$K_pred_dirac %*% gamma_alpha)
    if (!is.null(precomp$B_pred_dirac) && length(gamma_beta) > 0)
      phi_ev = phi_ev + as.vector(precomp$B_pred_dirac %*% gamma_beta)
    Z_ev = cbind(precomp$T_obs[precomp$event_idx],
                 Z[precomp$event_idx, , drop = FALSE])
    dirac[precomp$event_idx] = dchis_gamma(Z_ev, phi_ev)
  }

  terms = direct + dirac - compensator
  list(terms = terms, noise_var = noise_var)
}
