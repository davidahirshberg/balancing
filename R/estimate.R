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
# Gamma target: project dpsi measure onto representer basis
# ============================================================

#' Build the gamma target measure from the estimand's dpsi_grid.
#'
#' Returns a measure list(Z, r) suitable for project_target.
#' w_eval: per-subject bootstrap weights (NULL = all ones).
gam_measure = function(estimand, predict_haz, eval_Z, mesh_u, mesh,
                       M_train, du_bin, w_eval = NULL) {
  emb = estimand$dpsi_grid(predict_haz, eval_Z, mesh_u, mesh, M_train, du_bin)
  if (!is.null(w_eval)) {
    n_eval = nrow(atleast_2d(eval_Z))
    n_per_subj = length(emb$r) %/% n_eval
    emb$r = emb$r * rep(w_eval, times = n_per_subj)
  }
  emb
}

# ============================================================
# Survival pipeline
# ============================================================

#' Bregman balancing weights estimator for survival data.
#'
#' Computes the DR estimator at a single (sigma2_lam, sigma2_gam).
#' Returns $est, $se, $init. Pass init to subsequent calls to reuse
#' fold setup (meshes, kernel matrices, lambda fits).
#'
#' For grid sweeps, the caller loops over sigma2 values and passes init.
#' For bootstrap, the caller passes init and bootstrap_reps > 0.
#'
#' @param T_obs Observed event/censoring times.
#' @param D Event indicator (1 = event, 0 = censored).
#' @param W Treatment vector (binary 0/1).
#' @param X Covariate matrix.
#' @param kern Kernel object.
#' @param estimand Estimand object.
#' @param horizon Time horizon for survival estimand.
#' @param sigma2_lam Lambda regularization (sigma^2 = n*eta).
#' @param sigma2_gam Gamma regularization.
#' @param bootstrap_reps Number of bootstrap replications (0 = none).
#' @param init Initialization object from a previous call.
#' @param log Logging function: log(fmt, ...) for progress messages.
#' @return List with $est, $se, $se_amle, $init,
#'   and optionally $boot_ates, $boot_ses, $se_boot.
bregbal_surv = function(T_obs, D, W, X, kern,
                        estimand = surv_prob_ate(),
                        horizon,
                        lam_dispersion = entropy_dispersion(),
                        gam_dispersion = shifted_entropy_dispersion(),
                        sigma2_lam = 1,
                        sigma2_gam = 5,
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

  eta_lam = sigma2_lam / n
  eta_gam = sigma2_gam / n

  # --- If init provided, reuse fold setup ---
  if (!is.null(init)) {
    return(bregbal_surv_init(init, eta_lam, eta_gam, bootstrap_reps, log))
  }

  # --- Fresh computation ---
  ate_terms = rep(NA, n)
  ate_direct = rep(NA, n)
  ate_noise_var = rep(NA, n)

  # Bootstrap weights: generate once, stored in init for reuse
  w_boot = if (bootstrap_reps > 0) {
    array(rpois(n * bootstrap_reps, 1), dim = c(n, bootstrap_reps))
  } else NULL

  if (bootstrap_reps > 0) {
    all_boot_ates = array(0, bootstrap_reps)
    all_boot_ses = array(0, bootstrap_reps)
  }

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
    stacked_lam = at_risk_pairs(lam_train)
    mesh = lam_train$mesh
    du_bin = stacked_lam$du_bin

    eval_Z = Z[eval_fold, , drop = FALSE]

    gam_Z = Z[gam_fold, , drop = FALSE]
    gam_train = build_training_mesh(T_obs[gam_fold], D[gam_fold], gam_Z, horizon, M_train)
    stacked_gam = at_risk_pairs(gam_train)
    mesh_u = unique(stacked_gam$u_pool)

    W_pool_gam = estimand$W_fn(stacked_gam$Z_pool[, 2])
    gam_disp_train = signflip(gam_dispersion, W_pool_gam)
    # dchis_gamma(Z, phi): Z = (u, treatment, X...), phi = dual variable
    # Default: sign-flip of base dispersion using estimand's W_fn
    dchis_gamma = function(Z, phi) {
      W = estimand$W_fn(Z[, 2])
      W * gam_dispersion$dchis(W * phi)
    }

    K_lam = kernel_matrix(stacked_lam$Z_pool, stacked_lam$Z_pool, kern)
    Y_pool = stacked_lam$Y_pool
    # Lambda target: mu_lam = sum Y_j delta_{Z_j}  (event indicators at D pairs)
    c_lam = project_target(stacked_lam$Z_pool, kern,
                           list(Z = stacked_lam$Z_pool, r = Y_pool),
                           K_cross = K_lam)
    tgt_lam = split_target(c_lam, nrow(stacked_lam$Z_pool))

    i_pool_lam = stacked_lam$i_pool
    i_pool_gam = stacked_gam$i_pool

    # --- Lambda fit ---
    lam_res = fit_lambda_fold(stacked_lam, kern, eta_lam, lam_dispersion,
                               tgt_lam, K_lam, mesh, du_bin,
                               eval_Z, stacked_gam, mesh_u,
                               estimand, horizon, discrete, M_train)

    # --- Gamma fit ---
    res = fit_gamma_fold(lam_res, stacked_gam, kern, eta_gam, gam_disp_train,
                          eval_Z, T_obs[eval_fold], D[eval_fold], mesh, du_bin,
                          lam_dispersion, dchis_gamma,
                          horizon, Q_comp, n_eval, discrete)

    ate_terms[eval_fold] = res$direct + res$correction
    ate_direct[eval_fold] = res$direct
    ate_noise_var[eval_fold] = res$noise_var

    # --- Bootstrap ---
    if (bootstrap_reps > 0) {
      log("  bootstrap: %d reps", bootstrap_reps)
      boot_res = bootstrap_fold(res$lam_fit, res$gam_fit,
                                 res$lam_bound, res$gam_bound,
                                 res$events, bootstrap_reps,
                                 lam_fold, gam_fold, eval_fold,
                                 n_eval, n,
                                 T_obs[eval_fold], D[eval_fold],
                                 stacked_lam, stacked_gam,
                                 mesh, mesh_u, du_bin,
                                 i_pool_lam, i_pool_gam,
                                 gam_disp_train, dchis_gamma,
                                 eval_Z,
                                 kern, lam_dispersion,
                                 estimand, horizon, Q_comp,
                                 discrete,
                                 eta_gam, bootstrap_gamma_onestep,
                                 w_boot)
      all_boot_ates = all_boot_ates + boot_res$ates
      all_boot_ses = all_boot_ses + boot_res$ses
    }

    # Save fold state for init reuse
    fold_state[[ff]] = list(
      lam_fold = lam_fold, gam_fold = gam_fold, eval_fold = eval_fold,
      n_eval = n_eval,
      T_obs_eval = T_obs[eval_fold], D_eval = D[eval_fold],
      stacked_lam = stacked_lam, stacked_gam = stacked_gam,
      mesh = mesh, du_bin = du_bin, mesh_u = mesh_u,
      K_lam = K_lam, tgt_lam = tgt_lam,
      i_pool_lam = i_pool_lam, i_pool_gam = i_pool_gam,
      gam_disp_train = gam_disp_train, dchis_gamma = dchis_gamma,
      eval_Z = eval_Z)
  }

  # --- Assemble output ---
  out = list(
    est = mean(ate_terms),
    se = sd(ate_terms) / sqrt(n),
    # noise_var can be negative when quadratic dispersion predicts negative hazard;
    # clamp per-subject noise_var at 0 since it estimates a conditional variance
    se_amle = sqrt((var(ate_direct) + mean(pmax(ate_noise_var, 0))) / n),
    terms = ate_terms,
    direct = ate_direct,
    noise_var = ate_noise_var)

  if (bootstrap_reps > 0) {
    out$boot_ates = all_boot_ates / n_folds
    out$boot_ses = all_boot_ses / n_folds
    out$se_boot = sd(out$boot_ates)
  }

  out$init = list(
    fold_state = fold_state,
    n = n, n_folds = n_folds, fold_id = fold_id,
    kern = kern,
    lam_dispersion = lam_dispersion, gam_dispersion = gam_dispersion,
    estimand = estimand, horizon = horizon, Q_comp = Q_comp,
    M_train = M_train,
    discrete = discrete,
    bootstrap_gamma_onestep = bootstrap_gamma_onestep,
    w_boot = w_boot)

  out
}

# ============================================================
# Shared helpers: compute_direct and correction_terms (DR corrections)
# ============================================================

#' Compute direct estimand from a predict_phi function.
#' Returns list(direct, predict_haz) — predict_haz is the discrete-mesh
#' version used by dpsi_grid.
compute_direct = function(predict_phi_fn, lam_dispersion, estimand, eval_Z,
                           M_train, du_bin, mesh, horizon, discrete) {
  predict_haz = function(k, Z_ev) {
    lam_dispersion$dchis(predict_phi_fn(cbind(mesh[k], Z_ev)))
  }

  if (discrete) {
    direct = estimand$direct_discrete(predict_haz, eval_Z, M_train, du_bin)
  } else {
    Q_surv = 100
    predict_haz_cts = function(k, Z_ev) {
      phi = predict_phi_fn(cbind((k - 1) * horizon / Q_surv, Z_ev))
      pmax(lam_dispersion$dchis(phi) / du_bin, 0)
    }
    direct = estimand$direct_cts(predict_haz_cts, eval_Z, horizon, Q_surv)
  }

  list(direct = direct, predict_haz = predict_haz)
}

#' Compute DR correction terms (compensator + dirac) and noise variance.
#'
#' @param dchis_gamma Function(Z, phi) -> gamma. Evaluates the weight
#'   dispersion dchis at Z = (u, W, X) points with dual index phi.
correction_terms = function(lam_bound, gam_bound, T_obs_eval, D_eval,
                             mesh, du_bin, eval_Z, n_eval,
                             lam_dispersion, dchis_gamma,
                             horizon, Q_comp, discrete) {
  if (discrete) {
    comp = compensator_discrete(lam_bound, gam_bound,
                                T_obs_eval, D_eval, mesh,
                                lam_dispersion, dchis_gamma, eval_Z)
  } else {
    comp = compensator_simpson(lam_bound, gam_bound,
                               T_obs_eval, D_eval, eval_Z, horizon, Q_comp,
                               lam_dispersion, dchis_gamma, du_bin)
  }

  M_train = length(mesh)
  dirac = rep(0, n_eval)
  events = which(D_eval == 1 & T_obs_eval <= horizon)
  if (length(events) > 0) {
    t_ev = T_obs_eval[events]
    bin_idx = pmin(ceiling(t_ev / du_bin), M_train)
    t_mesh = mesh[bin_idx]
    phi_ev = eval_phi_vec(gam_bound, t_mesh, events)
    Z_ev = cbind(t_mesh, eval_Z[events, , drop = FALSE])
    dirac[events] = dchis_gamma(Z_ev, phi_ev)
  }

  list(correction = dirac - comp$comp, noise_var = comp$noise_var,
       events = events)
}

# ============================================================
# Lambda and gamma fit functions (called per fold)
# ============================================================

fit_lambda_fold = function(stacked_lam, kern, eta_lam, lam_dispersion,
                            tgt_lam, K_lam, mesh, du_bin,
                            eval_Z, stacked_gam, mesh_u,
                            estimand, horizon, discrete, M_train,
                            alpha0 = NULL, beta0 = NULL) {
  lam_fit = kernel_bregman(stacked_lam$Z_pool, kern,
                           eta = eta_lam, dispersion = lam_dispersion,
                           target = tgt_lam$target,
                           target_null = tgt_lam$target_null,
                           K = K_lam,
                           alpha0 = alpha0, beta0 = beta0)

  predict_phi_fn = function(Z) predict_phi(lam_fit, Z)
  dr = compute_direct(predict_phi_fn, lam_dispersion, estimand, eval_Z,
                       M_train, du_bin, mesh, horizon, discrete)

  # Gamma target: mu_gam = sum r_{im} delta_{Z_{im}^dpsi}  (dpsi embedding)
  emb = gam_measure(estimand, dr$predict_haz, eval_Z, mesh_u, mesh, M_train, du_bin)
  c_gam = split_target(
    project_target(stacked_gam$Z_pool, kern, emb),
    nrow(stacked_gam$Z_pool))

  lam_bound = bind_subjects(lam_fit, eval_Z)

  list(lam_fit = lam_fit, direct = dr$direct, c_gam = c_gam,
       lam_bound = lam_bound)
}

fit_gamma_fold = function(lam_res, stacked_gam, kern, eta_gam, gam_disp_train,
                           eval_Z, T_obs_eval, D_eval, mesh, du_bin,
                           lam_dispersion, dchis_gamma,
                           horizon, Q_comp, n_eval, discrete,
                           alpha0 = NULL, beta0 = NULL) {
  gam_fit = kernel_bregman(stacked_gam$Z_pool, kern,
                           eta = eta_gam, dispersion = gam_disp_train,
                           target = lam_res$c_gam$target,
                           target_null = lam_res$c_gam$target_null,
                           alpha0 = alpha0, beta0 = beta0)

  gam_bound = bind_subjects(gam_fit, eval_Z)

  dr = correction_terms(lam_res$lam_bound, gam_bound,
                         T_obs_eval, D_eval, mesh, du_bin,
                         eval_Z, n_eval,
                         lam_dispersion, dchis_gamma,
                         horizon, Q_comp, discrete)

  list(direct = lam_res$direct, correction = dr$correction,
       noise_var = dr$noise_var,
       lam_fit = lam_res$lam_fit, gam_fit = gam_fit,
       lam_bound = lam_res$lam_bound, gam_bound = gam_bound,
       events = dr$events)
}

# ============================================================
# Init-based call (reuses fold setup from a previous call)
# ============================================================

bregbal_surv_init = function(init, eta_lam, eta_gam, bootstrap_reps, log) {
  n = init$n
  n_folds = init$n_folds

  ate_terms = rep(NA, n)
  ate_direct = rep(NA, n)
  ate_noise_var = rep(NA, n)

  if (bootstrap_reps > 0) {
    all_boot_ates = array(0, bootstrap_reps)
    all_boot_ses = array(0, bootstrap_reps)
  }

  for (ff in 1:n_folds) {
    fs = init$fold_state[[ff]]
    log("fold %d/%d", ff, n_folds)

    lam_res = fit_lambda_fold(fs$stacked_lam, init$kern, eta_lam, init$lam_dispersion,
                               fs$tgt_lam, fs$K_lam, fs$mesh, fs$du_bin,
                               fs$eval_Z, fs$stacked_gam, fs$mesh_u,
                               init$estimand, init$horizon,
                               init$discrete, init$M_train)

    res = fit_gamma_fold(lam_res, fs$stacked_gam, init$kern, eta_gam, fs$gam_disp_train,
                          fs$eval_Z, fs$T_obs_eval, fs$D_eval, fs$mesh, fs$du_bin,
                          init$lam_dispersion, fs$dchis_gamma,
                          init$horizon, init$Q_comp, fs$n_eval,
                          init$discrete)

    ate_terms[fs$eval_fold] = res$direct + res$correction
    ate_direct[fs$eval_fold] = res$direct
    ate_noise_var[fs$eval_fold] = res$noise_var

    if (bootstrap_reps > 0) {
      log("  bootstrap: %d reps", bootstrap_reps)
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
        fs$gam_disp_train, fs$dchis_gamma,
        fs$eval_Z,
        init$kern, init$lam_dispersion,
        init$estimand, init$horizon, init$Q_comp,
        init$discrete,
        eta_gam, init$bootstrap_gamma_onestep,
        init$w_boot)
      all_boot_ates = all_boot_ates + boot_res$ates
      all_boot_ses = all_boot_ses + boot_res$ses
    }
  }

  out = list(
    est = mean(ate_terms),
    se = sd(ate_terms) / sqrt(n),
    # noise_var can be negative when quadratic dispersion predicts negative hazard;
    # clamp per-subject noise_var at 0 since it estimates a conditional variance
    se_amle = sqrt((var(ate_direct) + mean(pmax(ate_noise_var, 0))) / n),
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
                          gam_disp_train, dchis_gamma,
                          eval_Z,
                          kern, lam_dispersion,
                          estimand, horizon, Q_comp,
                          discrete,
                          eta_gam, bootstrap_gamma_onestep,
                          w_boot) {
  boot_ates = array(NA_real_, n_boot)
  boot_ses = array(NA_real_, n_boot)
  M_train = length(mesh)

  # Precompute cross-kernel matrix for project_target (Z_dpsi is fixed across reps)
  # Use a dummy predict_haz to get Z layout, then compute K_cross once.
  dummy_haz = function(k, Z_ev) rep(0, nrow(atleast_2d(Z_ev)))
  emb_template = estimand$dpsi_grid(dummy_haz, eval_Z, mesh_u, mesh, M_train, du_bin)
  K_cross_gam = kernel_matrix(stacked_gam$Z_pool, emb_template$Z, kern)

  for (b in 1:n_boot) {
    w_subj = if (!is.null(w_boot)) w_boot[, b] else rpois(n, 1)
    w_pool_lam = w_subj[lam_fold][i_pool_lam]
    w_pool_gam = w_subj[gam_fold][i_pool_gam]
    w_eval = w_subj[eval_fold]

    # One-step Newton for lambda
    # Bootstrap target: mu_lam_b = sum w_i Y_i delta_{Z_i}
    c_lam_b = project_target(stacked_lam$Z_pool, kern,
                             list(Z = stacked_lam$Z_pool, r = w_pool_lam * stacked_lam$Y_pool),
                             K_cross = lam_fit$K)
    tgt_b = split_target(c_lam_b, nrow(stacked_lam$Z_pool))
    lam_boot = onestep_bregman(lam_fit,
                                tgt_b$target, tgt_b$target_null,
                                w_pool_lam)

    # Recompute direct term via compute_direct
    predict_phi_fn = function(Z) {
      predict_phi_alpha(lam_boot$alpha, lam_boot$beta, lam_fit, Z)
    }
    dr = compute_direct(predict_phi_fn, lam_dispersion, estimand, eval_Z,
                         M_train, du_bin, mesh, horizon, discrete)

    # Recompute gamma target from bootstrap lambda
    emb_b = gam_measure(estimand, dr$predict_haz, eval_Z, mesh_u, mesh,
                         M_train, du_bin, w_eval = w_eval)
    c_gam_b = split_target(
      project_target(stacked_gam$Z_pool, kern, emb_b, K_cross = K_cross_gam),
      nrow(stacked_gam$Z_pool))

    # Re-solve gamma
    if (bootstrap_gamma_onestep) {
      gam_boot = onestep_bregman(gam_fit,
                                  target_new = c_gam_b$target,
                                  target_null_new = c_gam_b$target_null,
                                  w_new = w_pool_gam)
    } else {
      gam_boot = kernel_bregman(stacked_gam$Z_pool, kern,
                                 eta = eta_gam,
                                 dispersion = gam_disp_train,
                                 target = c_gam_b$target,
                                 target_null = c_gam_b$target_null,
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

    # DR correction via correction_terms
    dr_b = correction_terms(lam_bound_b, gam_bound_b,
                             T_obs_eval, D_eval, mesh, du_bin,
                             eval_Z, n_eval,
                             lam_dispersion, dchis_gamma,
                             horizon, Q_comp, discrete)

    dr_terms_b = dr$direct + dr_b$correction

    sw = sum(w_eval)
    boot_ates[b] = sum(w_eval * dr_terms_b) / sw
    boot_ses[b] = sd(dr_terms_b) / sqrt(n_eval)
  }

  list(ates = boot_ates, ses = boot_ses)
}
