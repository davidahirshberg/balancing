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
gam_embedding = function(pooled_eval, eval_X, mesh_u, mesh,
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

  Z_obs = pooled_eval$Z_pool
  n_pool = nrow(Z_obs)
  n_eval = nrow(eval_X)

  # Build the counterfactual evaluation grid: all subjects × all mesh times.
  # This is n_eval × n_mesh_u points.
  n_mu = length(mesh_u)

  if (!isTRUE(estimand$contrast)) {
    arm = estimand$arm

    # r_subj[i, j] = dot_psi(mesh_u[j]) for subject i
    mesh_k = match(mesh_u, mesh)
    r_subj = matrix(NA, n_eval, n_mu)
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

    # Cross-kernel: K(Z_obs, Z_ctf) — n_pool x (n_eval * n_mu)
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
  r1_subj = r0_subj = matrix(NA, n_eval, n_mu)
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
                        eta_lam = 0.5,
                        eta_lam_grid = NULL,
                        tuning = lepski_tuning(),
                        eta_gam = 5,
                        M_train = 15,
                        Q_comp = 50,
                        n_folds = 2,
                        time = "continuous",
                        intercept = TRUE,
                        bootstrap_reps = 0) {
  n = length(T_obs)
  Z = cbind(W, X)
  fold_id = make_folds(n, n_folds)
  discrete = (time == "discrete")
  contrast = isTRUE(estimand$contrast)
  is_rmst = grepl("rmst", estimand$name)

  tuning = as_tuning(tuning)
  do_sweep = !is.null(eta_lam_grid)
  m_grid = if (do_sweep) length(eta_lam_grid) else 1L

  ate_terms = rep(NA, n)
  ate_direct = rep(NA, n)
  ate_noise_var = rep(NA, n)

  # Path data (when sweeping grid): one entry per grid point, averaged over folds
  if (do_sweep) {
    path_dir_val = numeric(m_grid)   # mean direct term at each eta_lam
    path_dir_se = numeric(m_grid)    # SE of direct term
    path_dr_est = numeric(m_grid)    # mean DR estimate
    path_dr_se = numeric(m_grid)     # SE of DR terms
    path_counts = numeric(m_grid)    # fold count (for averaging)
  }

  # --- Helper: run the full pipeline at one eta_lam value ---
  # Returns list(direct, correction, noise_var) — per-eval-subject vectors
  run_at_eta = function(eta_val, pooled, lam_train, eval_idx,
                        eval_Z, eval_X, W_eval,
                        pooled_eval, mesh_u, gam_disp_train) {
    lam_fit = kernel_bregman(pooled$Y_pool, pooled$Z_pool, kern,
                             eta = eta_val, dispersion = lam_dispersion,
                             intercept = intercept)
    du_bin = pooled$du_bin
    mesh = lam_train$mesh

    predict_haz = function(k, Z_ev) {
      lam_dispersion$dchis(predict_phi(lam_fit, cbind(mesh[k], Z_ev)))
    }

    if (discrete) {
      if (is_rmst) {
        direct = estimand$direct_discrete(predict_haz, eval_X, M_train, du_bin)
      } else {
        direct = estimand$direct_discrete(predict_haz, eval_X, M_train)
      }

      # Weight model: kernel mean embedding of dot_psi as weight_target.
      # The weight loss has linear term c'α (paired with representer coefficients),
      # not Y'Kα (paired with function values). weight_target mode handles this.
      # Scale by n_u (at-risk count per time slot): the χ loss sums over training
      # points while dpsiK averages over eval subjects, so the FOC has
      # K(K+2ηI)α = n_u · dpsiK at each time block.
      c_gam = gam_embedding(pooled_eval, eval_X, mesh_u, mesh,
                            kern, estimand, predict_haz, M_train,
                            is_rmst = is_rmst, du_bin = du_bin)
      u_col = pooled_eval$Z_pool[, 1]
      n_u = as.numeric(table(u_col)[as.character(u_col)])
      c_gam = c_gam * n_u

      gam_fit = kernel_bregman(Y = NULL, pooled_eval$Z_pool, kern,
                               eta = eta_gam, dispersion = gam_disp_train,
                               intercept = intercept,
                               weight_target = c_gam)

      lam_bound = bind_subjects(lam_fit, eval_Z)
      gam_bound = bind_subjects(gam_fit, eval_Z)

      comp = compensator_discrete(lam_bound, gam_bound,
                                  T_obs[eval_idx], D[eval_idx], mesh,
                                  lam_dispersion, gam_dispersion,
                                  W_eval = W_eval)
    } else {
      Q_surv = 100
      predict_haz_cts = function(k, Z_ev) {
        phi = predict_phi(lam_fit, cbind((k - 1) * horizon / Q_surv, Z_ev))
        pmax(lam_dispersion$dchis(phi) / du_bin, 0)
      }

      direct = estimand$direct_cts(predict_haz_cts, eval_X, horizon, Q_surv)

      # Weight model: kernel mean embedding of dot_psi as weight_target.
      # Need mesh-indexed predict_haz (not the fine-grid predict_haz_cts).
      predict_haz_mesh = function(k, Z_ev) {
        lam_dispersion$dchis(predict_phi(lam_fit, cbind(mesh[k], Z_ev)))
      }
      c_gam = gam_embedding(pooled_eval, eval_X, mesh_u, mesh,
                            kern, estimand, predict_haz_mesh, M_train,
                            is_rmst = is_rmst, du_bin = du_bin)
      u_col = pooled_eval$Z_pool[, 1]
      n_u = as.numeric(table(u_col)[as.character(u_col)])
      c_gam = c_gam * n_u

      gam_fit = kernel_bregman(Y = NULL, pooled_eval$Z_pool, kern,
                               eta = eta_gam, dispersion = gam_disp_train,
                               intercept = intercept,
                               weight_target = c_gam)

      lam_bound = bind_subjects(lam_fit, eval_Z)
      gam_bound = bind_subjects(gam_fit, eval_Z)

      comp = compensator_simpson(lam_bound, gam_bound,
                                 T_obs[eval_idx], D[eval_idx], eval_Z, horizon, Q_comp,
                                 lam_dispersion, gam_dispersion, du_bin,
                                 W_eval = W_eval)
    }

    # Dirac: evaluate gamma at event times.
    # Map T_obs to mesh midpoints (gamma model centers live at midpoints,
    # and product kernels require exact time match).
    n_eval = length(eval_idx)
    dirac = rep(0, n_eval)
    events = which(D[eval_idx] == 1 & T_obs[eval_idx] <= horizon)
    if (length(events) > 0) {
      t_ev = T_obs[eval_idx[events]]
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
    list(direct = direct, correction = correction, noise_var = comp$noise_var)
  }

  for (ff in 1:n_folds) {
    lam_idx = fold_id != ff
    eval_idx = which(fold_id == ff)

    # Build training mesh (shared across eta_lam grid)
    lam_train = build_training_mesh(T_obs[lam_idx], D[lam_idx],
                                    Z[lam_idx, , drop = FALSE], horizon, M_train)
    pooled = pool_mesh(lam_train)

    eval_Z = Z[eval_idx, , drop = FALSE]
    eval_X = X[eval_idx, , drop = FALSE]
    W_eval = if (contrast) 2 * eval_Z[, 1] - 1 else NULL

    # Build eval fold mesh (shared across eta_lam grid)
    eval_train = build_training_mesh(T_obs[eval_idx], D[eval_idx], eval_Z, horizon, M_train)
    pooled_eval = pool_mesh(eval_train)
    mesh_u = unique(pooled_eval$u_pool)

    if (contrast) {
      W_pool = 2 * pooled_eval$Z_pool[, 2] - 1
      gam_disp_train = signflip(gam_dispersion, W_pool)
    } else {
      gam_disp_train = gam_dispersion
    }

    if (do_sweep) {
      # --- Grid sweep: run pipeline at each eta_lam ---
      grid_results = vector("list", m_grid)
      dir_vals = numeric(m_grid)
      dir_ses = numeric(m_grid)

      for (jj in 1:m_grid) {
        grid_results[[jj]] = tryCatch(
          run_at_eta(eta_lam_grid[jj], pooled, lam_train, eval_idx,
                     eval_Z, eval_X, W_eval,
                     pooled_eval, mesh_u, gam_disp_train),
          error = function(e) NULL)

        if (!is.null(grid_results[[jj]])) {
          d = grid_results[[jj]]$direct
          dir_vals[jj] = mean(d)
          dir_ses[jj] = sd(d) / sqrt(length(d))
        } else {
          dir_vals[jj] = NA
          dir_ses[jj] = NA
        }
      }

      # Tuning selection (pluggable: Lepski, oracle, fixed, rule, ...)
      sel = tuning$select(dir_vals, dir_ses, eta_lam_grid, n = length(eval_idx))
      selected = sel$idx

      # Use selected grid point for the main estimate
      res_sel = grid_results[[selected]]
      if (is.null(res_sel)) stop("Tuning-selected eta_lam failed at grid point ", selected)

      ate_terms[eval_idx] = res_sel$direct + res_sel$correction
      ate_direct[eval_idx] = res_sel$direct
      ate_noise_var[eval_idx] = res_sel$noise_var

      # Accumulate path data (always, regardless of selection strategy)
      n_eval = length(eval_idx)
      for (jj in 1:m_grid) {
        if (!is.null(grid_results[[jj]])) {
          r = grid_results[[jj]]
          dr_terms = r$direct + r$correction
          path_dir_val[jj] = path_dir_val[jj] + mean(r$direct)
          path_dir_se[jj] = path_dir_se[jj] + sd(r$direct) / sqrt(n_eval)
          path_dr_est[jj] = path_dr_est[jj] + mean(dr_terms)
          path_dr_se[jj] = path_dr_se[jj] + sd(dr_terms) / sqrt(n_eval)
          path_counts[jj] = path_counts[jj] + 1
        }
      }

    } else {
      # --- Fixed eta_lam (no grid) ---
      res = run_at_eta(eta_lam, pooled, lam_train, eval_idx,
                       eval_Z, eval_X, W_eval,
                       pooled_eval, mesh_u, gam_disp_train)

      ate_terms[eval_idx] = res$direct + res$correction
      ate_direct[eval_idx] = res$direct
      ate_noise_var[eval_idx] = res$noise_var
    }

    if (bootstrap_reps > 0) {
      stop("Bootstrap not yet updated for new estimand interface. ",
           "Use bootstrap_reps = 0 for now.")
    }
  }

  est = mean(ate_terms)
  se = sd(ate_terms) / sqrt(n)
  se_amle = sqrt((var(ate_direct) + mean(ate_noise_var)) / n)

  out = list(est = est, se = se, se_amle = se_amle,
             terms = ate_terms, direct = ate_direct, noise_var = ate_noise_var)

  if (do_sweep) {
    # Average path data over folds
    ok = path_counts > 0
    path_dir_val[ok] = path_dir_val[ok] / path_counts[ok]
    path_dir_se[ok] = path_dir_se[ok] / path_counts[ok]
    path_dr_est[ok] = path_dr_est[ok] / path_counts[ok]
    path_dr_se[ok] = path_dr_se[ok] / path_counts[ok]
    path_dir_val[!ok] = NA
    path_dir_se[!ok] = NA
    path_dr_est[!ok] = NA
    path_dr_se[!ok] = NA

    # Final selection on fold-averaged path (for reporting)
    final_sel = tuning$select(path_dir_val, path_dir_se, eta_lam_grid, n = n)
    out$selected_eta_lam = final_sel$eta
    out$tuning = tuning$name
    out$path = data.frame(
      eta_lam = eta_lam_grid,
      dir_val = path_dir_val,
      dir_se = path_dir_se,
      dr_est = path_dr_est,
      dr_se = path_dr_se,
      selected = (1:m_grid == final_sel$idx))
  }

  out
}
