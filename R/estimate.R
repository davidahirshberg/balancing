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

# ---- Gamma embedding ----
#
# Kernel mean embedding of dot_psi for the weight loss.
# tilde_Y_j = (1/n_eval) sum_{i,m} dot_psi(u_m, a, X_i) K(Z_{im}^ctf, Z_j)
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

#' Bregman balancing weights estimator for survival data.
#'
#' @param T_obs Observed times.
#' @param D Event indicators.
#' @param W Treatment (binary 0/1).
#' @param X Covariate matrix.
#' @param kern Kernel object.
#' @param estimand Survival estimand. Must provide high-level interface:
#'   direct_discrete/direct_cts and dot_psi_Z_discrete/dot_psi_Z_cts.
#'   Use surv_prob_ate(), surv_prob_tsm(arm), rmst_ate(), rmst_tsm(arm).
#' @param horizon Maximum follow-up time.
#' @param lam_dispersion Dispersion for hazard model (default: entropy = Poisson).
#' @param gam_dispersion Dispersion for weight model. The user passes the
#'   appropriate dispersion directly (e.g. signflip(entropy, W) for ATE,
#'   plain entropy for TSM). The pipeline does not modify it.
#' @param eta_lam Regularization for hazard model (or tuning strategy).
#' @param eta_gam Regularization for weight model.
#' @param M_train Number of time bins for training mesh.
#' @param Q_comp Quadrature points for compensator (continuous only).
#' @param n_folds Number of cross-fitting folds.
#' @param time "continuous" or "discrete". Controls which DR estimator is used.
#' @param bootstrap_reps Number of bootstrap replications (0 = none).
#' @return List with $est, $se, $se_amle, $terms, $direct, $noise_var.
#'   If bootstrap_reps > 0: also $boot_ates, $boot_ses (per-rep AMLE SEs).
bregbal_surv = function(T_obs, D, W, X, kern,
                        estimand = surv_prob_ate(),
                        horizon,
                        lam_dispersion = entropy_dispersion(),
                        gam_dispersion = entropy_dispersion(),
                        sigma2_lam = 1,
                        sigma2_lam_grid = NULL,
                        tuning = lepski_tuning(),
                        sigma2_gam = 5,
                        sigma2_gam_grid = NULL,
                        M_train = 15,
                        Q_comp = 50,
                        n_folds = 2,
                        time = "continuous",
                        bootstrap_reps = 0,
                        bootstrap_gamma_onestep = FALSE) {
  n = length(T_obs)
  Z = cbind(W, X)
  fold_id = make_folds(n, n_folds)
  discrete = (time == "discrete")
  contrast = isTRUE(estimand$contrast)
  is_rmst = grepl("rmst", estimand$name)

  tuning = as_tuning(tuning)

  # Convert sigma^2 grids to eta = sigma^2 / n
  eta_lam = sigma2_lam / n
  eta_lam_grid = if (!is.null(sigma2_lam_grid)) sigma2_lam_grid / n else NULL
  eta_gam = sigma2_gam / n
  eta_gam_grid = if (!is.null(sigma2_gam_grid)) sigma2_gam_grid / n else NULL

  do_sweep = !is.null(eta_lam_grid)
  m_lam = if (do_sweep) length(eta_lam_grid) else 1L
  m_gam = if (!is.null(eta_gam_grid)) length(eta_gam_grid) else 1L

  ate_terms = rep(NA, n)
  ate_direct = rep(NA, n)
  ate_noise_var = rep(NA, n)

  # Path data: one entry per (eta_lam, eta_gam) grid point, averaged over folds.
  # Stored as flat vectors indexed by (jl - 1) * m_gam + jg.
  m_total = m_lam * m_gam
  if (m_total > 1) {
    path_dir_val = numeric(m_total)
    path_dr_est = numeric(m_total)
    path_dr_se = numeric(m_total)
    path_loo_lam = rep(NA_real_, m_total)
    path_loo_gam = rep(NA_real_, m_total)
    path_counts = numeric(m_total)
  }

  # Pre-allocate bootstrap accumulators
  if (bootstrap_reps > 0) {
    all_boot_ates = array(0, bootstrap_reps)
    all_boot_ses = array(0, bootstrap_reps)
  }

  for (ff in 1:n_folds) {
    # Fold assignment:
    #   2-fold (balancing): lambda on ~ff, gamma+eval on ff.
    #     Gamma is trained on the people it weights.
    #   3-fold (plug-in): lambda on ff, gamma on (ff%%3)+1, eval on ((ff+1)%%3)+1.
    #     Standard honest sample splitting.
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

    # --- Fold data (computed once, shared across eta grid) ---
    lam_train = build_training_mesh(T_obs[lam_fold], D[lam_fold],
                                    Z[lam_fold, , drop = FALSE], horizon, M_train)
    stacked_lam = pool_mesh(lam_train)
    mesh = lam_train$mesh
    du_bin = stacked_lam$du_bin

    eval_Z = Z[eval_fold, , drop = FALSE]
    eval_X = X[eval_fold, , drop = FALSE]
    W_eval = if (contrast) 2 * eval_Z[, 1] - 1 else NULL

    # Gamma fold data: same as eval in 2-fold, separate in 3-fold
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

    # Pre-compute lambda kernel matrix and targets (O(n^2), once per fold)
    K_lam = kernel_matrix(stacked_lam$Z_pool, stacked_lam$Z_pool, kern)
    B_lam = null_basis(stacked_lam$Z_pool, kern)
    tgt_lam = dpsi_from_response(stacked_lam$Y_pool, K_lam, B_lam)

    # Pool indices for bootstrap resampling
    i_pool_lam = stacked_lam$i_pool
    i_pool_gam = stacked_gam$i_pool

    # --- fit_lambda: closure over fold data, returns lambda-side results ---
    # Everything that depends on eta_lam but not eta_gam.
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

        # Gamma needs mesh-indexed predict_haz (not the fine-grid version)
        predict_haz = function(k, Z_ev) {
          lam_dispersion$dchis(predict_phi(lam_fit, cbind(mesh[k], Z_ev)))
        }
      }

      # Gamma target: kernel mean embedding of dot_psi, scaled by n_u.
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

    # --- fit_gamma: given lambda results, fit gamma and compute DR terms ---
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

      # Dirac: evaluate gamma at event times
      dirac = rep(0, n_eval)
      events = which(D[eval_fold] == 1 & T_obs[eval_fold] <= horizon)
      if (length(events) > 0) {
        t_ev = T_obs[eval_fold][events]
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

    # --- bootstrap_fold: closure over fold data ---
    bootstrap_fold = function(lam_fit, gam_fit, lam_bound, gam_bound,
                              events, n_boot) {
      boot_ates = array(NA_real_, n_boot)
      boot_ses = array(NA_real_, n_boot)

      for (b in 1:n_boot) {
        # One weight per subject, sliced by fold membership
        w_subj = rpois(n, 1)
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
                                         T_obs[eval_fold], D[eval_fold],
                                         mesh, lam_dispersion, gam_dispersion,
                                         W_eval = W_eval)
        } else {
          comp_b = compensator_simpson(lam_bound_b, gam_bound_b,
                                        T_obs[eval_fold], D[eval_fold],
                                        eval_Z, horizon, Q_comp,
                                        lam_dispersion, gam_dispersion, du_bin,
                                        W_eval = W_eval)
        }

        # Dirac from bootstrap gamma
        dirac_b = rep(0, n_eval)
        if (length(events) > 0) {
          t_ev = T_obs[eval_fold][events]
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

    # --- Execute: sweep (eta_lam × eta_gam) grid ---
    eta_lam_vec = if (do_sweep) eta_lam_grid else eta_lam
    eta_gam_vec = if (!is.null(eta_gam_grid)) eta_gam_grid else eta_gam

    # Sweep lambda (with warm starts)
    lam_results = vector("list", m_lam)
    for (jl in 1:m_lam) {
      prev = if (jl > 1 && !is.null(lam_results[[jl - 1]])) lam_results[[jl - 1]]$lam_fit else NULL
      lam_results[[jl]] = tryCatch(
        fit_lambda(eta_lam_vec[jl],
                   alpha0 = prev$alpha, beta0 = prev$beta),
        error = function(e) NULL)
    }

    # Sweep gamma at each lambda (with warm starts along gamma axis)
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

    # Tuning selection: Lepski/CV along lambda axis (at fixed gamma for now)
    # Use middle gamma if grid, otherwise the only gamma
    jg_sel = if (m_gam > 1) ceiling(m_gam / 2) else 1L

    dir_vals = vapply(1:m_lam, function(jl) {
      r = grid_results[[jl, jg_sel]]
      if (!is.null(r)) mean(r$direct) else NA_real_
    }, numeric(1))
    dir_ses = vapply(1:m_lam, function(jl) {
      r = grid_results[[jl, jg_sel]]
      if (!is.null(r)) sd(r$direct) / sqrt(n_eval) else NA_real_
    }, numeric(1))

    ctx = tuning_context(eta_lam_vec, n = n_eval,
      dir_val = function() dir_vals,
      dir_se = function() dir_ses,
      loo_scores = function() {
        vapply(1:m_lam, function(jl) {
          if (!is.null(lam_results[[jl]])) loo_loss(lam_results[[jl]]$lam_fit)
          else Inf
        }, numeric(1))
      })
    sel_lam = tuning$select(ctx)
    jl_sel = sel_lam$idx

    # Gamma selection: CV along gamma axis at selected lambda (if grid)
    if (m_gam > 1) {
      gam_loo = vapply(1:m_gam, function(jg) {
        r = grid_results[[jl_sel, jg]]
        if (!is.null(r)) loo_loss(r$gam_fit) else Inf
      }, numeric(1))
      jg_sel = which.min(gam_loo)
    }

    res_sel = grid_results[[jl_sel, jg_sel]]
    if (is.null(res_sel)) stop("Selected (eta_lam, eta_gam) failed")

    ate_terms[eval_fold] = res_sel$direct + res_sel$correction
    ate_direct[eval_fold] = res_sel$direct
    ate_noise_var[eval_fold] = res_sel$noise_var

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

    # --- Bootstrap (at selected point only) ---
    if (bootstrap_reps > 0) {
      boot_res = bootstrap_fold(res_sel$lam_fit, res_sel$gam_fit,
                                 res_sel$lam_bound, res_sel$gam_bound,
                                 res_sel$events, bootstrap_reps)
      all_boot_ates = all_boot_ates + boot_res$ates
      all_boot_ses = all_boot_ses + boot_res$ses
    }
  }

  est = mean(ate_terms)
  se = sd(ate_terms) / sqrt(n)
  se_amle = sqrt((var(ate_direct) + mean(ate_noise_var)) / n)

  out = list(est = est, se = se, se_amle = se_amle,
             terms = ate_terms, direct = ate_direct, noise_var = ate_noise_var)

  if (bootstrap_reps > 0) {
    out$boot_ates = all_boot_ates / n_folds
    out$boot_ses = all_boot_ses / n_folds
    out$se_boot = sd(out$boot_ates)
  }

  if (m_total > 1) {
    # Average path data over folds
    ok = path_counts > 0
    path_dir_val[ok] = path_dir_val[ok] / path_counts[ok]
    path_dr_est[ok] = path_dr_est[ok] / path_counts[ok]
    path_dr_se[ok] = path_dr_se[ok] / path_counts[ok]
    path_dir_val[!ok] = NA
    path_dr_est[!ok] = NA
    path_dr_se[!ok] = NA

    # Build path data frame (one row per grid point)
    eta_lam_vec = if (do_sweep) eta_lam_grid else eta_lam
    eta_gam_vec = if (!is.null(eta_gam_grid)) eta_gam_grid else eta_gam
    out$path = data.frame(
      sigma2_lam = rep(eta_lam_vec * n, each = m_gam),
      sigma2_gam = rep(eta_gam_vec * n, times = m_lam),
      dir_val = path_dir_val,
      dr_est = path_dr_est,
      dr_se = path_dr_se,
      loo_lam = path_loo_lam,
      loo_gam = path_loo_gam)
    out$tuning = tuning$name

    # Final selection on fold-averaged path (for reporting)
    eta_lam_vec_final = if (do_sweep) eta_lam_grid else eta_lam
    dir_vals_avg = path_dir_val[seq(1, m_total, by = m_gam)]  # lambda axis, first gamma
    dir_ses_avg = path_dr_se[seq(1, m_total, by = m_gam)]
    ctx_final = tuning_context(eta_lam_vec_final, n = n,
      dir_val = function() dir_vals_avg,
      dir_se = function() dir_ses_avg,
      loo_scores = function() {
        path_loo_lam[seq(1, m_total, by = m_gam)]
      })
    final_lam = tuning$select(ctx_final)
    out$selected_sigma2_lam = final_lam$eta * n
    if (m_gam > 1) {
      gam_slice = path_loo_gam[((final_lam$idx - 1) * m_gam + 1):(final_lam$idx * m_gam)]
      out$selected_sigma2_gam = eta_gam_vec[which.min(gam_slice)] * n
    }
  }

  out
}
