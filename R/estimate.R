## Top-level estimation pipelines.

# ============================================================
# Single-outcome pipeline
# ============================================================

#' Bregman balancing weights estimator for single-outcome data.
#'
#' @param Y Outcome vector.
#' @param W Treatment vector (binary 0/1).
#' @param X Covariate matrix.
#' @param kern Kernel object.
#' @param estimand Estimand (default: tsm_estimand(1)).
#' @param dispersion Dispersion for outcome model (default: softplus for binary).
#' @param eta_grid Regularization grid for outcome model.
#' @param tuning Tuning strategy for outcome eta (default: lepski).
#' @param eta_gam Regularization for weight model.
#' @param n_folds Number of cross-fitting folds.
#' @return List with $est, $se, $eif, $selected_eta.
bregbal = function(Y, W, X, kern,
                   estimand = tsm_estimand(1),
                   dispersion = softplus_dispersion(),
                   eta_grid = 4^seq(-3, 6, by = 0.5),
                   tuning = lepski_tuning(),
                   eta_gam = 0.1,
                   n_folds = 3,
                   seed = NULL) {
  crossfit_single(Y, W, X, kern, estimand, dispersion,
                  eta_grid, tuning, eta_gam, n_folds, seed)
}

# ============================================================
# Survival pipeline
# ============================================================

#' Bregman balancing weights estimator for survival data.
#'
#' @param T_obs Observed times.
#' @param D Event indicators.
#' @param W Treatment (binary 0/1).
#' @param X Covariate matrix.
#' @param kern Kernel object.
#' @param estimand Survival estimand (default: surv_prob_estimand()).
#' @param horizon Maximum follow-up time.
#' @param lam_dispersion Dispersion for hazard model (default: entropy = Poisson).
#' @param gam_dispersion Base dispersion for weight model (default: entropy).
#'   Will be sign-flipped for ATE.
#' @param eta_lam Regularization for hazard model (or tuning strategy).
#' @param eta_gam Regularization for weight model.
#' @param M_train Number of time bins for training mesh.
#' @param Q_comp Quadrature points for compensator.
#' @param n_folds Number of cross-fitting folds.
#' @return List with $est, $se, $terms, $direct, $correction.
bregbal_surv = function(T_obs, D, W, X, kern,
                        estimand = surv_prob_estimand(),
                        horizon,
                        lam_dispersion = entropy_dispersion(),
                        gam_dispersion = entropy_dispersion(),
                        eta_lam = 0.5,
                        eta_gam = 5,
                        M_train = 15,
                        Q_comp = 50,
                        n_folds = 2) {
  n = length(T_obs)
  Z = cbind(W, X)
  fold_id = make_folds(n, n_folds)

  ate_terms = rep(NA, n)

  for (ff in 1:n_folds) {
    lam_idx = fold_id != ff
    eval_idx = which(fold_id == ff)

    # Build training mesh and fit hazard
    train = build_training_mesh(T_obs[lam_idx], D[lam_idx],
                                Z[lam_idx, , drop = FALSE], horizon, M_train)
    pooled = pool_mesh(train)

    # Resolve eta_lam (could be a number or tuning strategy)
    eta_lam_val = if (is.numeric(eta_lam)) eta_lam else stop("Tuning for survival eta_lam not yet implemented")

    lam_fit = kernel_bregman(pooled$Y_pool, pooled$Z_pool, kern,
                             eta = eta_lam_val, dispersion = lam_dispersion)

    # Lambda function closure
    du_bin = pooled$du_bin
    lambda_fn = function(u, z) {
      phi = predict_phi(lam_fit, cbind(u, atleast_2d(z)))
      pmax(lam_dispersion$dchis(phi) / du_bin, 0)
    }

    # Fit weights (sign-flipped for ATE)
    eval_Z = Z[eval_idx, , drop = FALSE]
    W_eval = 2 * eval_Z[, 1] - 1  # {-1, +1}
    gam_disp = signflip(gam_dispersion, W_eval)

    # dot_psi targets for weight fitting
    Q_surv = 100
    h_eval = make_haz_fn(lam_fit, eval_Z, horizon, Q_surv, lam_dispersion, du_bin = du_bin)
    deriv = estimand$dot_psi(h_eval$haz_fn, Q_surv, h_eval$du)

    # Pool weight targets onto training mesh
    train_eval = build_training_mesh(T_obs[eval_idx], D[eval_idx], eval_Z, horizon, M_train)
    pooled_eval = pool_mesh(train_eval)
    W_pool = 2 * pooled_eval$Z_pool[, 2] - 1  # treatment is col 2 after u
    mesh_u = unique(pooled_eval$u_pool)
    deriv_vals = matrix(NA, nrow(eval_Z), length(mesh_u))
    for (j in seq_along(mesh_u)) deriv_vals[, j] = deriv(mesh_u[j])
    mesh_idx = match(pooled_eval$u_pool, mesh_u)
    Y_gam = W_pool * deriv_vals[cbind(pooled_eval$i_pool, mesh_idx)]

    gam_disp_pool = signflip(gam_dispersion, W_pool)
    gam_fit = kernel_bregman(Y_gam, pooled_eval$Z_pool, kern,
                             eta = eta_gam, dispersion = gam_disp_pool)

    # DR estimate: direct + dirac - compensator
    Z1 = with.treatment(eval_Z, 1)
    Z0 = with.treatment(eval_Z, 0)
    h1 = make_haz_fn(lam_fit, Z1, horizon, Q_surv, lam_dispersion, du_bin = du_bin)
    h0 = make_haz_fn(lam_fit, Z0, horizon, Q_surv, lam_dispersion, du_bin = du_bin)
    direct = estimand$psi(h1$haz_fn, Q_surv, h1$du) -
             estimand$psi(h0$haz_fn, Q_surv, h0$du)

    # Bind subjects for compensator
    lam_bound = bind_subjects(lam_fit, eval_Z)
    gam_bound = bind_subjects(gam_fit, eval_Z)

    comp = compensator_simpson(lam_bound, gam_bound,
                               T_obs[eval_idx], D[eval_idx], eval_Z, horizon, Q_comp,
                               lam_dispersion, gam_disp_pool, du_bin)

    # Dirac
    n_eval = length(eval_idx)
    dirac = rep(0, n_eval)
    events = which(D[eval_idx] == 1 & T_obs[eval_idx] <= horizon)
    if (length(events) > 0) {
      phi_ev = eval_phi_vec(gam_bound, T_obs[eval_idx[events]], events)
      dirac[events] = gam_disp_pool$dchis(phi_ev)
      # Wait — gam_disp_pool has W_pool, but dirac needs W from eval subjects.
      # Need to use eval-level dispersion with per-event W.
      W_ev = 2 * eval_Z[events, 1] - 1
      gam_disp_ev = signflip(gam_dispersion, W_ev)
      dirac[events] = gam_disp_ev$dchis(phi_ev)
    }

    correction = dirac - comp$comp
    ate_terms[eval_idx] = direct + correction
  }

  est = mean(ate_terms)
  se = sd(ate_terms) / sqrt(n)
  list(est = est, se = se, terms = ate_terms)
}
