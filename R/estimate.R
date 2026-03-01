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
                   n_folds = 3) {
  crossfit_single(Y, W, X, kern, estimand, dispersion,
                  eta_grid, tuning, eta_gam, n_folds)
}

# ============================================================
# Gamma target: kernel mean embedding of dot_psi
# ============================================================

gam_embedding = function(stacked_gam, eval_X, mesh_u, mesh,
                         kern, estimand, predict_haz, M_train,
                         is_rmst = FALSE, du_bin = NULL) {
  # Compute kernel mean embedding of dot_psi:
  #   tilde_Y_j = (1 / n_eval) * sum_{i,m} r_{im} K(Z_{im}^ctf, Z_j)
  #
  # The sum is over ALL eval-fold subjects (not just at-risk), because
  # dot_psi is a population-level functional. The at-risk set only
  # restricts where gamma is trained (which Z_j appear).
  #
  # For TSM: Z_i^ctf = (u_m, arm, X_i).
  # For ATE: two-arm sum, Z_i^ctf_a = (u_m, a, X_i).

  Z_obs = stacked_gam$Z_pool
  n_pool = nrow(Z_obs)
  n_eval = nrow(eval_X)

  # Build the counterfactual evaluation grid: all subjects x all mesh times.
  # This is n_eval x n_mesh_u points.
  n_mu = length(mesh_u)

  if (!isTRUE(estimand$contrast)) {
    arm = estimand$arm

    # r_subj[i, j] = dot_psi(mesh_u[j]) for subject i
    mesh_k = match(mesh_u, mesh)
    r_subj = array(NA_real_, c(n_eval, n_mu))
    if (is_rmst) {
      dp = estimand$dot_psi_Z_discrete(predict_haz, eval_X, M_train, du_bin)
    } else {
      dp = estimand$dot_psi_Z_discrete(predict_haz, eval_X, M_train)
    }
    for (j in seq_along(mesh_u)) r_subj[, j] = dp(mesh_k[j])

    # Counterfactual grid: (u_m, arm, X_i) for all i, m
    Z_ctf_list = lapply(seq_along(mesh_u), function(j) {
      cbind(mesh_u[j], arm, eval_X)
    })
    Z_ctf = do.call(rbind, Z_ctf_list)  # (n_eval * n_mu) x p
    r_ctf = as.vector(r_subj)  # length n_eval * n_mu

    # Cross-kernel: K(Z_obs, Z_ctf) -- n_pool x (n_eval * n_mu)
    K_cross = kernel_matrix(Z_obs, Z_ctf, kern)
    return(as.vector(K_cross %*% r_ctf) / n_eval)
  }

  # ATE: two-arm sum.
  if (is_rmst) {
    dp1 = estimand$dot_psi_arm_discrete(predict_haz, eval_X, M_train, du_bin, arm = 1)
    dp0 = estimand$dot_psi_arm_discrete(predict_haz, eval_X, M_train, du_bin, arm = 0)
  } else {
    dp1 = estimand$dot_psi_arm_discrete(predict_haz, eval_X, M_train, arm = 1)
    dp0 = estimand$dot_psi_arm_discrete(predict_haz, eval_X, M_train, arm = 0)
  }

  mesh_k = match(mesh_u, mesh)
  r1_subj = r0_subj = array(NA_real_, c(n_eval, n_mu))
  for (j in seq_along(mesh_u)) {
    r1_subj[, j] = dp1(mesh_k[j])
    r0_subj[, j] = dp0(mesh_k[j])
  }

  Z_ctf1_list = lapply(seq_along(mesh_u), function(j) cbind(mesh_u[j], 1, eval_X))
  Z_ctf0_list = lapply(seq_along(mesh_u), function(j) cbind(mesh_u[j], 0, eval_X))
  Z_ctf1 = do.call(rbind, Z_ctf1_list)
  Z_ctf0 = do.call(rbind, Z_ctf0_list)

  K1 = kernel_matrix(Z_obs, Z_ctf1, kern)
  K0 = kernel_matrix(Z_obs, Z_ctf0, kern)

  as.vector(K1 %*% as.vector(r1_subj) + K0 %*% as.vector(r0_subj)) / n_eval
}

# ============================================================
# Survival pipeline
# ============================================================

#' Bregman balancing weights estimator for survival data.
#'
#' Computes the DR estimator over a (sigma2_lam × sigma2_gam) grid.
#' Returns path data for external tuning/selection, and optionally
#' runs bootstrap at a single (sigma2_lam, sigma2_gam) point.
#'
#' When called with init (from a previous grid call), skips setup
#' and grid computation, reuses cached fits and bootstrap weights.
#'
#' @param T_obs Observed event/censoring times.
#' @param D Event indicator (1 = event, 0 = censored).
#' @param W Treatment vector (binary 0/1).
#' @param X Covariate matrix.
#' @param kern Kernel object.
#' @param estimand Estimand object.
#' @param horizon Time horizon for survival estimand.
#' @param sigma2_lam Scalar lambda regularization (sigma^2 = n*eta).
#' @param sigma2_lam_grid Grid of sigma2_lam values for sweep.
#' @param sigma2_gam Scalar gamma regularization.
#' @param sigma2_gam_grid Grid of sigma2_gam values for sweep.
#' @param bootstrap_reps Number of bootstrap replications (0 = none).
#' @param init Initialization object from a previous call (reuses fits,
#'   meshes, kernel matrices, bootstrap weights).
#' @param log Logging function: log(fmt, ...) for progress messages.
#' @return List with $est, $se, $path, $init, and optionally $boot_*.
bregbal_surv = function(T_obs, D, W, X, kern,
                        estimand = surv_prob_ate(),
                        horizon,
                        lam_dispersion = entropy_dispersion(),
                        gam_dispersion = entropy_dispersion(),
                        sigma2_lam = 1,
                        sigma2_lam_grid = NULL,
                        sigma2_gam = 5,
                        sigma2_gam_grid = NULL,
                        M_train = 15,
                        Q_comp = 50,
                        n_folds = 2,
                        time = "continuous",
                        bootstrap_reps = 0,
                        bootstrap_gamma_onestep = FALSE,
                        init = NULL,
                        log = function(...) invisible()) {
  n = length(T_obs)
  Z = cbind(W, X)
  fold_id = make_folds(n, n_folds)
  discrete = (time == "discrete")
  contrast = isTRUE(estimand$contrast)
  is_rmst = grepl("rmst", estimand$name)

  # Convert sigma^2 to eta = sigma^2 / n
  eta_lam = sigma2_lam / n
  eta_lam_grid = if (!is.null(sigma2_lam_grid)) sigma2_lam_grid / n else NULL
  eta_gam = sigma2_gam / n
  eta_gam_grid = if (!is.null(sigma2_gam_grid)) sigma2_gam_grid / n else NULL

  do_sweep = !is.null(eta_lam_grid)
  m_lam = if (do_sweep) length(eta_lam_grid) else 1L
  m_gam = if (!is.null(eta_gam_grid)) length(eta_gam_grid) else 1L

  eta_lam_vec = if (do_sweep) eta_lam_grid else eta_lam
  eta_gam_vec = if (!is.null(eta_gam_grid)) eta_gam_grid else eta_gam

  # --- If init provided, skip setup and grid; go straight to bootstrap ---
  if (!is.null(init)) {
    return(bregbal_surv_bootstrap(init, eta_lam, eta_gam, bootstrap_reps, log))
  }

  # --- Fresh computation: setup + grid sweep ---
  ate_terms = rep(NA, n)
  ate_direct = rep(NA, n)
  ate_noise_var = rep(NA, n)

  m_total = m_lam * m_gam
  if (m_total > 1) {
    path_dir_val = numeric(m_total)
    path_dr_est = numeric(m_total)
    path_dr_se = numeric(m_total)
    path_loo_lam = rep(NA_real_, m_total)
    path_loo_gam = rep(NA_real_, m_total)
    path_counts = numeric(m_total)
  }

  # Bootstrap weights: generate once, reuse across init calls
  w_boot = if (bootstrap_reps > 0) {
    array(rpois(n * bootstrap_reps, 1), dim = c(n, bootstrap_reps))
  } else NULL

  if (bootstrap_reps > 0) {
    all_boot_ates = array(0, bootstrap_reps)
    all_boot_ses = array(0, bootstrap_reps)
  }

  # Per-fold state, saved for init
  fold_state = vector("list", n_folds)

  for (ff in 1:n_folds) {
    log("fold %d/%d", ff, n_folds)

    if (n_folds == 2) {
      lam_fold = fold_id != ff
      gam_fold = fold_id == ff
      eval_fold = gam_fold
    } else {
      lam_fold = fold_id == ff
      gam_fold = fold_id == (ff %% n_folds) + 1
      eval_fold = fold_id == ((ff + 1) %% n_folds) + 1
    }

    n_eval = sum(eval_fold)

    # --- Fold data ---
    lam_train = build_training_mesh(T_obs[lam_fold], D[lam_fold],
                                    Z[lam_fold, , drop = FALSE], horizon, M_train)
    stacked_lam = pool_mesh(lam_train)
    mesh = lam_train$mesh
    du_bin = stacked_lam$du_bin

    eval_Z = Z[eval_fold, , drop = FALSE]
    eval_X = X[eval_fold, , drop = FALSE]
    W_eval = if (contrast) 2 * eval_Z[, 1] - 1 else NULL

    gam_Z = Z[gam_fold, , drop = FALSE]
    gam_train = build_training_mesh(T_obs[gam_fold], D[gam_fold], gam_Z, horizon, M_train)
    stacked_gam = pool_mesh(gam_train)
    mesh_u = unique(stacked_gam$u_pool)

    if (contrast) {
      W_pool_gam = 2 * stacked_gam$Z_pool[, 2] - 1
      gam_disp_train = signflip(gam_dispersion, W_pool_gam)
    } else {
      gam_disp_train = gam_dispersion
    }

    K_lam = kernel_matrix(stacked_lam$Z_pool, stacked_lam$Z_pool, kern)
    B_lam = null_basis(stacked_lam$Z_pool, kern)
    tgt_lam = dpsi_from_response(stacked_lam$Y_pool, K_lam, B_lam)

    i_pool_lam = stacked_lam$i_pool
    i_pool_gam = stacked_gam$i_pool

    # --- fit_lambda ---
    fit_lambda = function(eta_val, alpha0 = NULL, beta0 = NULL) {
      lam_fit = kernel_bregman(stacked_lam$Z_pool, kern,
                               eta = eta_val, dispersion = lam_dispersion,
                               target = tgt_lam$target,
                               target_null = tgt_lam$target_null,
                               K = K_lam,
                               alpha0 = alpha0, beta0 = beta0)

      predict_haz = function(k, Z_ev) {
        lam_dispersion$dchis(predict_phi(lam_fit, cbind(mesh[k], Z_ev)))
      }

      if (discrete) {
        direct = if (is_rmst) estimand$direct_discrete(predict_haz, eval_X, M_train, du_bin)
                 else estimand$direct_discrete(predict_haz, eval_X, M_train)
      } else {
        Q_surv = 100
        predict_haz_cts = function(k, Z_ev) {
          phi = predict_phi(lam_fit, cbind((k - 1) * horizon / Q_surv, Z_ev))
          pmax(lam_dispersion$dchis(phi) / du_bin, 0)
        }
        direct = estimand$direct_cts(predict_haz_cts, eval_X, horizon, Q_surv)

        predict_haz = function(k, Z_ev) {
          lam_dispersion$dchis(predict_phi(lam_fit, cbind(mesh[k], Z_ev)))
        }
      }

      c_gam = gam_embedding(stacked_gam, eval_X, mesh_u, mesh,
                            kern, estimand, predict_haz, M_train,
                            is_rmst = is_rmst, du_bin = du_bin)
      u_col = stacked_gam$Z_pool[, 1]
      n_u = as.numeric(table(u_col)[as.character(u_col)])
      c_gam = c_gam * n_u

      lam_bound = bind_subjects(lam_fit, eval_Z)

      list(lam_fit = lam_fit, direct = direct, c_gam = c_gam,
           lam_bound = lam_bound)
    }

    # --- fit_gamma ---
    fit_gamma = function(lam_res, eta_gam_val, alpha0 = NULL, beta0 = NULL) {
      gam_fit = kernel_bregman(stacked_gam$Z_pool, kern,
                               eta = eta_gam_val, dispersion = gam_disp_train,
                               target = lam_res$c_gam,
                               alpha0 = alpha0, beta0 = beta0)

      gam_bound = bind_subjects(gam_fit, eval_Z)

      if (discrete) {
        comp = compensator_discrete(lam_res$lam_bound, gam_bound,
                                    T_obs[eval_fold], D[eval_fold], mesh,
                                    lam_dispersion, gam_dispersion,
                                    W_eval = W_eval)
      } else {
        comp = compensator_simpson(lam_res$lam_bound, gam_bound,
                                   T_obs[eval_fold], D[eval_fold], eval_Z, horizon, Q_comp,
                                   lam_dispersion, gam_dispersion, du_bin,
                                   W_eval = W_eval)
      }

      dirac = rep(0, n_eval)
      T_obs_eval = T_obs[eval_fold]
      D_eval = D[eval_fold]
      events = which(D_eval == 1 & T_obs_eval <= horizon)
      if (length(events) > 0) {
        t_ev = T_obs_eval[events]
        bin_idx = pmin(ceiling(t_ev / du_bin), M_train)
        t_mesh = mesh[bin_idx]
        phi_ev = eval_phi_vec(gam_bound, t_mesh, events)
        if (contrast) {
          W_ev = W_eval[events]
          dirac[events] = W_ev * gam_dispersion$dchis(W_ev * phi_ev)
        } else {
          dirac[events] = gam_dispersion$dchis(phi_ev)
        }
      }

      correction = dirac - comp$comp
      list(direct = lam_res$direct, correction = correction,
           noise_var = comp$noise_var,
           lam_fit = lam_res$lam_fit, gam_fit = gam_fit,
           lam_bound = lam_res$lam_bound, gam_bound = gam_bound,
           events = events)
    }

    # --- Grid sweep ---
    log("  lambda sweep: %d points", m_lam)
    lam_results = vector("list", m_lam)
    for (jl in 1:m_lam) {
      prev = if (jl > 1 && !is.null(lam_results[[jl - 1]])) lam_results[[jl - 1]]$lam_fit else NULL
      lam_results[jl] = list(tryCatch(
        fit_lambda(eta_lam_vec[jl],
                   alpha0 = prev$alpha, beta0 = prev$beta),
        error = function(e) { log("    fit_lambda error jl=%d: %s", jl, e$message); NULL }))
    }

    log("  gamma sweep: %d×%d grid", m_lam, m_gam)
    grid_results = vector("list", m_lam * m_gam)
    dim(grid_results) = c(m_lam, m_gam)
    for (jl in 1:m_lam) {
      if (is.null(lam_results[[jl]])) next
      prev_gam = NULL
      for (jg in 1:m_gam) {
        grid_results[jl, jg] = list(tryCatch(
          fit_gamma(lam_results[[jl]], eta_gam_vec[jg],
                    alpha0 = prev_gam$alpha, beta0 = prev_gam$beta),
          error = function(e) NULL))
        if (!is.null(grid_results[[jl, jg]]))
          prev_gam = grid_results[[jl, jg]]$gam_fit
      }
    }

    # Accumulate path data
    if (m_total > 1) {
      for (jl in 1:m_lam) {
        for (jg in 1:m_gam) {
          r = grid_results[[jl, jg]]
          if (is.null(r)) next
          idx = (jl - 1) * m_gam + jg
          dr_terms = r$direct + r$correction
          path_dir_val[idx] = path_dir_val[idx] + mean(r$direct)
          path_dr_est[idx] = path_dr_est[idx] + mean(dr_terms)
          path_dr_se[idx] = path_dr_se[idx] + sd(dr_terms) / sqrt(n_eval)
          if (!is.null(lam_results[[jl]]))
            path_loo_lam[idx] = loo_loss(lam_results[[jl]]$lam_fit)
          path_loo_gam[idx] = loo_loss(r$gam_fit)
          path_counts[idx] = path_counts[idx] + 1
        }
      }
    }

    # Save fold state for init reuse
    fold_state[[ff]] = list(
      lam_fold = lam_fold, gam_fold = gam_fold, eval_fold = eval_fold,
      n_eval = n_eval,
      T_obs_eval = T_obs[eval_fold], D_eval = D[eval_fold],
      stacked_lam = stacked_lam, stacked_gam = stacked_gam,
      mesh = mesh, du_bin = du_bin, mesh_u = mesh_u,
      K_lam = K_lam, B_lam = B_lam, tgt_lam = tgt_lam,
      i_pool_lam = i_pool_lam, i_pool_gam = i_pool_gam,
      gam_disp_train = gam_disp_train,
      eval_Z = eval_Z, eval_X = eval_X, W_eval = W_eval,
      lam_results = lam_results, grid_results = grid_results)

    # --- Bootstrap at single point (if no grid, just sigma2_lam × sigma2_gam) ---
    if (bootstrap_reps > 0 && m_total == 1) {
      res_sel = grid_results[[1, 1]]
      if (is.null(res_sel)) stop("Fit at (sigma2_lam, sigma2_gam) failed")

      ate_terms[eval_fold] = res_sel$direct + res_sel$correction
      ate_direct[eval_fold] = res_sel$direct
      ate_noise_var[eval_fold] = res_sel$noise_var

      log("  bootstrap: %d reps", bootstrap_reps)
      boot_res = bootstrap_fold(res_sel$lam_fit, res_sel$gam_fit,
                                 res_sel$lam_bound, res_sel$gam_bound,
                                 res_sel$events, bootstrap_reps,
                                 lam_fold, gam_fold, eval_fold,
                                 n_eval, n,
                                 T_obs[eval_fold], D[eval_fold],
                                 stacked_lam, stacked_gam,
                                 mesh, mesh_u, du_bin,
                                 i_pool_lam, i_pool_gam,
                                 gam_disp_train,
                                 eval_Z, eval_X, W_eval,
                                 kern, lam_dispersion, gam_dispersion,
                                 estimand, horizon, Q_comp,
                                 discrete, contrast, is_rmst,
                                 eta_gam, bootstrap_gamma_onestep,
                                 w_boot)
      all_boot_ates = all_boot_ates + boot_res$ates
      all_boot_ses = all_boot_ses + boot_res$ses
    } else if (m_total == 1) {
      res_sel = grid_results[[1, 1]]
      if (!is.null(res_sel)) {
        ate_terms[eval_fold] = res_sel$direct + res_sel$correction
        ate_direct[eval_fold] = res_sel$direct
        ate_noise_var[eval_fold] = res_sel$noise_var
      }
    }
  }

  # --- Assemble output ---
  out = list()

  if (m_total == 1) {
    out$est = mean(ate_terms)
    out$se = sd(ate_terms) / sqrt(n)
    out$se_amle = sqrt((var(ate_direct) + mean(ate_noise_var)) / n)
    out$terms = ate_terms
    out$direct = ate_direct
    out$noise_var = ate_noise_var

    if (bootstrap_reps > 0) {
      out$boot_ates = all_boot_ates / n_folds
      out$boot_ses = all_boot_ses / n_folds
      out$se_boot = sd(out$boot_ates)
    }
  }

  if (m_total > 1) {
    ok = path_counts > 0
    path_dir_val[ok] = path_dir_val[ok] / path_counts[ok]
    path_dr_est[ok] = path_dr_est[ok] / path_counts[ok]
    path_dr_se[ok] = path_dr_se[ok] / path_counts[ok]
    path_dir_val[!ok] = NA
    path_dr_est[!ok] = NA
    path_dr_se[!ok] = NA

    out$path = data.frame(
      sigma2_lam = rep(eta_lam_vec * n, each = m_gam),
      sigma2_gam = rep(eta_gam_vec * n, times = m_lam),
      dir_val = path_dir_val,
      dr_est = path_dr_est,
      dr_se = path_dr_se,
      loo_lam = path_loo_lam,
      loo_gam = path_loo_gam)
  }

  # Init object for reuse
  out$init = list(
    fold_state = fold_state,
    n = n, n_folds = n_folds, fold_id = fold_id,
    eta_lam_vec = eta_lam_vec, eta_gam_vec = eta_gam_vec,
    m_lam = m_lam, m_gam = m_gam,
    kern = kern,
    lam_dispersion = lam_dispersion, gam_dispersion = gam_dispersion,
    estimand = estimand, horizon = horizon, Q_comp = Q_comp,
    M_train = M_train,
    discrete = discrete, contrast = contrast, is_rmst = is_rmst,
    bootstrap_gamma_onestep = bootstrap_gamma_onestep,
    w_boot = w_boot)

  out
}

# ============================================================
# Bootstrap using init object
# ============================================================

#' Run bootstrap at a specific (sigma2_lam, sigma2_gam) using cached fits.
#'
#' Looks up the fit at the requested grid point, verifies the FOC,
#' and runs the bootstrap loop reusing all cached per-fold state.
bregbal_surv_bootstrap = function(init, eta_lam, eta_gam, bootstrap_reps, log) {
  n = init$n
  n_folds = init$n_folds

  # Find grid indices for requested eta
  jl = which.min(abs(init$eta_lam_vec - eta_lam))
  jg = which.min(abs(init$eta_gam_vec - eta_gam))

  ate_terms = rep(NA, n)
  ate_direct = rep(NA, n)
  ate_noise_var = rep(NA, n)

  all_boot_ates = if (bootstrap_reps > 0) array(0, bootstrap_reps) else NULL
  all_boot_ses = if (bootstrap_reps > 0) array(0, bootstrap_reps) else NULL

  for (ff in 1:n_folds) {
    fs = init$fold_state[[ff]]
    res = fs$grid_results[[jl, jg]]
    if (is.null(res)) stop(sprintf("No fit at grid[%d,%d] for fold %d", jl, jg, ff))

    # FOC check: verify cached fit is valid
    check_foc(res$lam_fit, "lambda", ff, log)
    check_foc(res$gam_fit, "gamma", ff, log)

    ate_terms[fs$eval_fold] = res$direct + res$correction
    ate_direct[fs$eval_fold] = res$direct
    ate_noise_var[fs$eval_fold] = res$noise_var

    if (bootstrap_reps > 0) {
      log("  fold %d/%d bootstrap: %d reps", ff, n_folds, bootstrap_reps)
      boot_res = bootstrap_fold(
        res$lam_fit, res$gam_fit,
        res$lam_bound, res$gam_bound,
        res$events, bootstrap_reps,
        fs$lam_fold, fs$gam_fold, fs$eval_fold,
        fs$n_eval, n,
        fs$T_obs_eval, fs$D_eval,
        fs$stacked_lam, fs$stacked_gam,
        fs$mesh, fs$mesh_u, fs$du_bin,
        fs$i_pool_lam, fs$i_pool_gam,
        fs$gam_disp_train,
        fs$eval_Z, fs$eval_X, fs$W_eval,
        init$kern, init$lam_dispersion, init$gam_dispersion,
        init$estimand, init$horizon, init$Q_comp,
        init$discrete, init$contrast, init$is_rmst,
        eta_gam, init$bootstrap_gamma_onestep,
        init$w_boot)
      all_boot_ates = all_boot_ates + boot_res$ates
      all_boot_ses = all_boot_ses + boot_res$ses
    }
  }

  out = list(
    est = mean(ate_terms),
    se = sd(ate_terms) / sqrt(n),
    se_amle = sqrt((var(ate_direct) + mean(ate_noise_var)) / n),
    terms = ate_terms, direct = ate_direct, noise_var = ate_noise_var,
    init = init)

  if (bootstrap_reps > 0) {
    out$boot_ates = all_boot_ates / n_folds
    out$boot_ses = all_boot_ses / n_folds
    out$se_boot = sd(out$boot_ates)
  }

  out
}

# ============================================================
# FOC check
# ============================================================

#' Verify the first-order condition of a cached fit.
#' At convergence: K(w*mu + n*eta*alpha) = target, so
#' grad_alpha = K(w*mu) + n*eta*K*alpha - target should be ~0.
check_foc = function(fit, label, fold, log) {
  K = fit$K
  nn = nrow(K)
  phi = as.vector(K %*% fit$alpha)
  if (length(fit$beta) > 0) phi = phi + as.vector(fit$B %*% fit$beta)

  w = if (!is.null(fit$w)) fit$w else rep(1, nn)
  mu = fit$dispersion$dchis(phi)

  grad_a = as.vector(K %*% (w * mu)) + nn * fit$eta * as.vector(K %*% fit$alpha) - fit$target
  rel_err = sqrt(sum(grad_a^2)) / (sqrt(sum(fit$target^2)) + 1e-12)

  if (rel_err > 1e-3) {
    log("  WARNING: %s FOC check fold %d: rel_err=%.2e", label, fold, rel_err)
  }
}

# ============================================================
# Bootstrap loop (extracted, shared by direct and init paths)
# ============================================================

bootstrap_fold = function(lam_fit, gam_fit, lam_bound, gam_bound,
                          events, n_boot,
                          lam_fold, gam_fold, eval_fold,
                          n_eval, n,
                          T_obs_eval, D_eval,
                          stacked_lam, stacked_gam,
                          mesh, mesh_u, du_bin,
                          i_pool_lam, i_pool_gam,
                          gam_disp_train,
                          eval_Z, eval_X, W_eval,
                          kern, lam_dispersion, gam_dispersion,
                          estimand, horizon, Q_comp,
                          discrete, contrast, is_rmst,
                          eta_gam, bootstrap_gamma_onestep,
                          w_boot) {
  boot_ates = array(NA_real_, n_boot)
  boot_ses = array(NA_real_, n_boot)
  M_train = length(mesh)

  for (b in 1:n_boot) {
    w_subj = if (!is.null(w_boot)) w_boot[, b] else rpois(n, 1)
    w_pool_lam = w_subj[lam_fold][i_pool_lam]
    w_pool_gam = w_subj[gam_fold][i_pool_gam]
    w_eval = w_subj[eval_fold]

    # One-step Newton for lambda
    tgt_b = dpsi_from_response(stacked_lam$Y_pool, lam_fit$K, lam_fit$B,
                                w_pool_lam)
    lam_boot = onestep_bregman(lam_fit,
                                tgt_b$target, tgt_b$target_null,
                                w_pool_lam)

    # Rebuild predict_haz from bootstrap lambda
    predict_haz_boot = function(k, Z_ev) {
      lam_dispersion$dchis(predict_phi_alpha(
        lam_boot$alpha, lam_boot$beta,
        lam_fit, cbind(mesh[k], Z_ev)))
    }

    # Recompute direct term
    if (discrete) {
      if (is_rmst) {
        direct_b = estimand$direct_discrete(predict_haz_boot, eval_X,
                                             M_train, du_bin)
      } else {
        direct_b = estimand$direct_discrete(predict_haz_boot, eval_X, M_train)
      }
    } else {
      Q_surv = 100
      predict_haz_boot_cts = function(k, Z_ev) {
        phi = predict_phi_alpha(lam_boot$alpha, lam_boot$beta,
                                lam_fit,
                                cbind((k - 1) * horizon / Q_surv, Z_ev))
        pmax(lam_dispersion$dchis(phi) / du_bin, 0)
      }
      direct_b = estimand$direct_cts(predict_haz_boot_cts, eval_X,
                                      horizon, Q_surv)

      predict_haz_boot = function(k, Z_ev) {
        lam_dispersion$dchis(predict_phi_alpha(
          lam_boot$alpha, lam_boot$beta,
          lam_fit, cbind(mesh[k], Z_ev)))
      }
    }

    # Recompute gamma target from bootstrap lambda
    c_gam_b = gam_embedding(stacked_gam, eval_X, mesh_u, mesh,
                             kern, estimand, predict_haz_boot, M_train,
                             is_rmst = is_rmst, du_bin = du_bin)
    u_col = stacked_gam$Z_pool[, 1]
    n_u = as.numeric(table(u_col)[as.character(u_col)])
    c_gam_b = c_gam_b * n_u

    # Re-solve gamma
    if (bootstrap_gamma_onestep) {
      gam_boot = onestep_bregman(gam_fit,
                                  target_new = c_gam_b,
                                  w_new = w_pool_gam)
    } else {
      gam_boot = kernel_bregman(stacked_gam$Z_pool, kern,
                                 eta = eta_gam,
                                 dispersion = gam_disp_train,
                                 target = c_gam_b,
                                 w = w_pool_gam,
                                 K = gam_fit$K,
                                 alpha0 = gam_fit$alpha,
                                 beta0 = gam_fit$beta)
    }

    # Patch bound objects with bootstrap alphas
    lam_bound_b = lam_bound
    lam_bound_b$alpha = lam_boot$alpha
    lam_bound_b$beta = lam_boot$beta

    gam_bound_b = gam_bound
    gam_bound_b$alpha = gam_boot$alpha
    gam_bound_b$beta = gam_boot$beta

    # Compensator from bootstrap nuisances
    if (discrete) {
      comp_b = compensator_discrete(lam_bound_b, gam_bound_b,
                                     T_obs_eval, D_eval,
                                     mesh, lam_dispersion, gam_dispersion,
                                     W_eval = W_eval)
    } else {
      comp_b = compensator_simpson(lam_bound_b, gam_bound_b,
                                    T_obs_eval, D_eval,
                                    eval_Z, horizon, Q_comp,
                                    lam_dispersion, gam_dispersion, du_bin,
                                    W_eval = W_eval)
    }

    # Dirac from bootstrap gamma
    dirac_b = rep(0, n_eval)
    if (length(events) > 0) {
      t_ev = T_obs_eval[events]
      bin_idx = pmin(ceiling(t_ev / du_bin), M_train)
      t_mesh = mesh[bin_idx]
      phi_ev = eval_phi_vec(gam_bound_b, t_mesh, events)
      if (contrast) {
        W_ev = W_eval[events]
        dirac_b[events] = W_ev * gam_dispersion$dchis(W_ev * phi_ev)
      } else {
        dirac_b[events] = gam_dispersion$dchis(phi_ev)
      }
    }

    correction_b = dirac_b - comp_b$comp
    dr_terms_b = direct_b + correction_b

    sw = sum(w_eval)
    boot_ates[b] = sum(w_eval * dr_terms_b) / sw
    boot_ses[b] = sd(dr_terms_b) / sqrt(n_eval)
  }

  list(ates = boot_ates, ses = boot_ses)
}
