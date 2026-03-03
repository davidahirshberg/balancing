## Survival DR estimator with crossfitting and bootstrap.
##
## Per-fold primitives:
##   fit_lambda(mp, kern, eta, lambda_disp, warm, log)
##   fit_gamma(mp, lambda_fit, kern, eta, estimand, gamma_disp_fn, gamma_base, warm, log)
##   estimate(mp, lambda_fit, gamma_result, estimand, grid)
##
## Aggregation:
##   aggregate_dr(mps, fold_estimates, ...) -> dr_result
##
## Convenience wrapper:
##   survival_effect(observations, kern, estimand, ...) -> dr_result
##
## Bootstrap:
##   bootstrap(result, n_reps) -> boot_result

# ============================================================
# Internal helpers
# ============================================================

# Build the gamma target measure from the estimand's dpsi_grid.
gamma_measure = function(estimand, predict_haz, mp, w_eval = NULL) {
  emb = estimand$dpsi_grid(predict_haz, mp$eval_Z, mp$mesh_u, mp$train_grid)
  if (!is.null(w_eval)) {
    n_eval = nrow(atleast_2d(mp$eval_Z))
    n_per_subj = length(emb$r) %/% n_eval
    emb$r = emb$r * rep(w_eval, times = n_per_subj)
  }
  emb
}

# Compute plug-in (direct) estimand.
compute_direct = function(phi, lambda_disp, estimand, mp, grid) {
  dchis = \(p) .dchis(lambda_disp, p)

  # predict_haz on training mesh: model output = dLambda on training grid.
  # Used by gamma_measure / make_dchis_gamma_ts (always on training mesh).
  predict_haz = function(k, Z_ev) dchis(phi(cbind(mp$mesh[k], Z_ev)))

  # predict_dLambda on evaluation grid: scale by du_eval / du_bin.
  # For discrete, eval grid = train grid, so scale = 1.
  # For continuous, eval grid is finer than train grid.
  scale = grid$du / mp$du_bin
  predict_dLambda = function(k, Z_ev)
    pmax(dchis(phi(cbind(grid$points[k], Z_ev))) * scale, 0)

  direct = estimand$direct(predict_dLambda, mp$eval_Z, grid)
  list(direct = direct, predict_haz = predict_haz)
}

# DR correction terms: compensator integral + Dirac at event times.
correction_terms = function(mp, lambda_fit_, gamma_fit_, lambda_disp, dchis_gamma,
                            comp_grid) {
  comp = compensator(lambda_fit_, gamma_fit_,
                     mp$T_obs_eval, mp$eval_Z,
                     comp_grid, lambda_disp, dchis_gamma, mp$du_bin)

  M_train = length(mp$mesh)
  dirac = rep(0, mp$n_eval)
  events = which(mp$D_eval == 1 & mp$T_obs_eval <= mp$horizon)
  if (length(events) > 0) {
    t_ev = mp$T_obs_eval[events]
    bin_idx = pmin(ceiling(t_ev / mp$du_bin), M_train)
    t_mesh = mp$mesh[bin_idx]
    Z_ev = cbind(t_mesh, mp$eval_Z[events, , drop = FALSE])
    phi_events = \(Z) .phi(gamma_fit_[events], Z)
    phi_ev = phi_events(Z_ev)
    dirac[events] = dchis_gamma(Z_ev, phi_ev)
  }
  list(correction = dirac - comp$comp, noise_var = comp$noise_var,
       events = events)
}

#' Build dchis_gamma(Z, phi) for the gamma correction prediction.
#'
#' Looks up the per-observation offset r from the estimand's dpsi_grid,
#' constructs dispersion(gamma_base, offset = r, sigma = sign(r)),
#' and evaluates .dchis at phi.
#'
#' @param gamma_base  Base dispersion record. Default entropy_base (gives TSE).
make_dchis_gamma_ts = function(predict_haz, estimand_obj, eval_Z,
                                train_grid, gamma_base = entropy_base) {
  M_train = train_grid$M
  du_bin = train_grid$du
  n_eval = nrow(atleast_2d(eval_Z))
  eval_hash = apply(eval_Z, 1, paste0, collapse = "|")
  lookup = setNames(seq_len(n_eval), eval_hash)
  has_per_arm = !is.null(estimand_obj$dot_psi_arm)
  if (has_per_arm) {
    dp1 = estimand_obj$dot_psi_arm(predict_haz, eval_Z, train_grid, arm = 1)
    dp0 = estimand_obj$dot_psi_arm(predict_haz, eval_Z, train_grid, arm = 0)
    r1_mat = r0_mat = matrix(NA_real_, n_eval, M_train)
    for (k in 1:M_train) { r1_mat[, k] = dp1(k); r0_mat[, k] = -dp0(k) }
  } else {
    dpsi = estimand_obj$dot_psi_Z(predict_haz, eval_Z, train_grid)
    r_mat = matrix(NA_real_, n_eval, M_train)
    for (k in 1:M_train) r_mat[, k] = dpsi(k)
  }
  function(Z, phi) {
    u = Z[, 1]
    bin = pmin(pmax(ceiling(u / du_bin), 1L), M_train)
    idx = lookup[apply(Z[, -1, drop = FALSE], 1, paste0, collapse = "|")]
    if (has_per_arm) {
      A = Z[, 2]
      r = ifelse(A == 1, r1_mat[cbind(idx, bin)], r0_mat[cbind(idx, bin)])
    } else {
      r = r_mat[cbind(idx, bin)]
    }
    .dchis(dispersion(gamma_base, offset = r, sigma = sign(r)), phi)
  }
}

# ============================================================
# Mesh projection
# ============================================================

#' Project observations onto a mesh, one mesh_projection per crossfit fold.
#'
#' @param observations List with $T_obs, $D, $Z.
#' @param horizon Time horizon.
#' @param M_train Number of time bins.
#' @param n_folds Number of crossfit folds.
#' @param kern Kernel object.
#' @return List of mesh_projection objects, one per fold.
mesh_project = function(observations, horizon, M_train, n_folds, kern, discrete = FALSE) {
  T_obs = observations$T_obs; D = observations$D; Z = atleast_2d(observations$Z)
  n = length(T_obs)
  fold_id = (seq_len(n) %% n_folds) + 1L

  lapply(1:n_folds, function(ff) {
    if (n_folds == 2) {
      lambda_idx = which(fold_id != ff); gamma_idx = which(fold_id == ff); eval_idx = gamma_idx
    } else {
      lambda_idx = which(fold_id == ff)
      gamma_idx = which(fold_id == (ff %% n_folds) + 1)
      eval_idx = which(fold_id == ((ff + 1) %% n_folds) + 1)
    }

    lambda_train = build_training_mesh(T_obs[lambda_idx], D[lambda_idx],
                                       Z[lambda_idx, , drop = FALSE], horizon, M_train, discrete)
    stacked_lambda = at_risk_pairs(lambda_train)

    gamma_train = build_training_mesh(T_obs[gamma_idx], D[gamma_idx],
                                      Z[gamma_idx, , drop = FALSE], horizon, M_train, discrete)
    stacked_gamma = at_risk_pairs(gamma_train)

    K_lambda = kernel_matrix(stacked_lambda$Z_pool, stacked_lambda$Z_pool, kern)
    c_lambda = project_target(stacked_lambda$Z_pool, kern,
                              list(Z = stacked_lambda$Z_pool, r = stacked_lambda$Y_pool),
                              K_cross = K_lambda)
    tgt_lambda = split_target(c_lambda, nrow(stacked_lambda$Z_pool))

    eval_Z = Z[eval_idx, , drop = FALSE]

    # Training grid: always discrete (model output = dLambda on mesh points)
    train_grid = discrete_grid(lambda_train$mesh, du = lambda_train$du)

    list(lambda_idx = lambda_idx, gamma_idx = gamma_idx, eval_idx = eval_idx,
         stacked_lambda = stacked_lambda, stacked_gamma = stacked_gamma,
         mesh = lambda_train$mesh, du_bin = stacked_lambda$du_bin,
         mesh_u = unique(stacked_gamma$u_pool),
         K_lambda = K_lambda, tgt_lambda = tgt_lambda,
         train_grid = train_grid,
         eval_Z = eval_Z,
         T_obs_eval = T_obs[eval_idx], D_eval = D[eval_idx],
         n_eval = length(eval_idx),
         horizon = horizon,
         i_pool_lambda = stacked_lambda$i_pool, i_pool_gamma = stacked_gamma$i_pool)
  })
}

# ============================================================
# Per-fold primitives
# ============================================================

#' Fit the hazard model (lambda) on one mesh projection.
#' @param mp A mesh_projection (one element of mesh_project output).
#' @param warm Optional previous fit for warm-starting ($alpha, $beta).
#' @return A kernel_bregman fit object.
fit_lambda = function(mp, kern, eta, lambda_disp, warm = NULL, log = null_logger) {
  kernel_bregman(mp$stacked_lambda$Z_pool, kern,
                 eta = eta, dispersion = lambda_disp,
                 target = mp$tgt_lambda$target,
                 target_null = mp$tgt_lambda$target_null,
                 K = mp$K_lambda,
                 alpha0 = warm$alpha, beta0 = warm$beta,
                 log = log)
}

#' Fit the weight model (gamma) on one mesh projection, given a fitted lambda.
#'
#' @param mp A mesh_projection.
#' @param lambda_fit Output of fit_lambda.
#' @param gamma_disp_fn Factory: r -> dispersion object. Default target_scaled_entropy.
#' @param gamma_base Base dispersion record for correction weights. Default entropy_base (TSE).
#' @param warm Optional previous gamma fit for warm-starting.
#' @param gamma_target Pre-computed target (skips recomputation if provided).
#' @param gamma_target_null Pre-computed null-space target.
#' @param w Per-observation weights (for bootstrap).
#' @return List with $gamma_fit, $dchis_gamma, $predict_haz, $gamma_target, $gamma_target_null.
fit_gamma = function(mp, lambda_fit, kern, eta, estimand,
                     gamma_disp_fn = target_scaled_entropy,
                     gamma_base = entropy_base,
                     warm = NULL,
                     gamma_target = NULL, gamma_target_null = NULL,
                     w = NULL,
                     log = null_logger) {
  phi   = \(Z) .phi(lambda_fit, Z)
  dchis = \(p) .dchis(lambda_fit$dispersion, p)
  predict_haz = function(k, Z_ev) dchis(phi(cbind(mp$mesh[k], Z_ev)))

  dchis_gamma = make_dchis_gamma_ts(predict_haz, estimand, mp$eval_Z,
                                     mp$train_grid,
                                     gamma_base = gamma_base)

  if (is.null(gamma_target)) {
    emb = gamma_measure(estimand, predict_haz, mp)
    tgt = split_target(project_target(mp$stacked_gamma$Z_pool, kern, emb),
                       nrow(mp$stacked_gamma$Z_pool))
    gamma_target = tgt$target
    gamma_target_null = tgt$target_null
  }
  gamma_disp = gamma_disp_fn(gamma_target)

  gamma_fit = kernel_bregman(mp$stacked_gamma$Z_pool, kern,
                             eta = eta, dispersion = gamma_disp,
                             target = gamma_target, target_null = gamma_target_null,
                             w = w,
                             K      = warm$K,
                             alpha0 = warm$alpha, beta0 = warm$beta,
                             log = indent(log, "gamma"))

  list(gamma_fit = gamma_fit, dchis_gamma = dchis_gamma,
       predict_haz = predict_haz,
       gamma_target = gamma_target, gamma_target_null = gamma_target_null)
}

#' Compute per-subject DR terms on one mesh projection given fitted lambda and gamma.
#'
#' @param mp A mesh_projection.
#' @param lambda_fit Output of fit_lambda.
#' @param gamma_result Output of fit_gamma.
#' @param grid Evaluation grid (discrete_grid or continuous_grid).
#' @return List with $direct, $correction, $terms, $noise_var, $events.
estimate = function(mp, lambda_fit, gamma_result, estimand, grid) {
  lambda_disp = lambda_fit$dispersion
  lambda_fit_ = .bind(lambda_fit, mp$eval_Z)
  gamma_fit_ = .bind(gamma_result$gamma_fit, mp$eval_Z)
  phi = \(Z) .phi(lambda_fit, Z)

  dr = compute_direct(phi, lambda_disp, estimand, mp, grid)

  ct = correction_terms(mp, lambda_fit_, gamma_fit_,
                         lambda_disp, gamma_result$dchis_gamma,
                         grid)

  list(direct = dr$direct, correction = ct$correction,
       terms = dr$direct + ct$correction,
       noise_var = ct$noise_var, events = ct$events)
}

# ============================================================
# Aggregation
# ============================================================

#' Aggregate per-fold estimates into a dr_result.
#'
#' @param mps List of mesh_projections (from mesh_project).
#' @param lambda_fits List of per-fold lambda fits.
#' @param gamma_results List of per-fold gamma results (from fit_gamma).
#' @param fold_estimates List of per-fold estimates (from estimate).
#' @return A dr_result S3 object.
aggregate_dr = function(mps, lambda_fits, gamma_results, fold_estimates,
                        kern, eta_gam, lambda_disp, estimand,
                        horizon, grid,
                        gamma_disp_fn = target_scaled_entropy,
                        gamma_base = entropy_base) {
  n = sum(sapply(mps, `[[`, "n_eval"))
  terms = rep(NA_real_, n)
  direct = rep(NA_real_, n)
  noise_var = rep(NA_real_, n)

  for (ff in seq_along(mps)) {
    idx = mps[[ff]]$eval_idx
    terms[idx] = fold_estimates[[ff]]$terms
    direct[idx] = fold_estimates[[ff]]$direct
    noise_var[idx] = fold_estimates[[ff]]$noise_var
  }

  structure(
    list(est = mean(terms), se = sd(terms) / sqrt(n),
         se_amle = sqrt((var(direct) + mean(pmax(noise_var, 0))) / n),
         terms = terms, direct = direct, noise_var = noise_var,
         mps = mps,
         lambda_fits = lambda_fits,
         gamma_results = gamma_results,
         fold_estimates = fold_estimates,
         n = n,
         kern = kern, eta_gam = eta_gam, lambda_disp = lambda_disp,
         estimand = estimand,
         horizon = horizon, grid = grid,
         gamma_disp_fn = gamma_disp_fn, gamma_base = gamma_base),
    class = "dr_result")
}

# ============================================================
# Convenience wrapper
# ============================================================

#' Doubly-robust survival estimator with crossfitting.
#'
#' @param observations List with $T_obs, $D, $Z.
#' @param kern Kernel object.
#' @param estimand Estimand object (e.g. surv_ate()).
#' @param lambda_disp Dispersion for lambda (hazard model).
#' @param eta_lam Regularization for lambda.
#' @param eta_gam Regularization for gamma (weight model).
#' @param horizon Time horizon.
#' @param gamma_disp_fn Factory: r -> gamma dispersion. Default target_scaled_entropy.
#' @param gamma_base Base dispersion record for correction weights. Default entropy_base (TSE).
#' @param mps Optional pre-computed mesh projections (skips mesh_project).
#' @param lambda_fits Optional pre-fitted lambda results per fold (skips fit_lambda).
#' @param gamma_results Optional pre-fitted gamma results per fold (skips fit_gamma).
#' @return A dr_result S3 object.
survival_effect = function(observations, kern, estimand, lambda_disp,
                           eta_lam, eta_gam, horizon,
                           M_train = 15, n_folds = 2, Q_comp = 50,
                           discrete = FALSE, log = null_logger,
                           gamma_disp_fn = target_scaled_entropy,
                           gamma_base = entropy_base,
                           mps = NULL, lambda_fits = NULL, gamma_results = NULL) {
  if (is.null(mps))
    mps = mesh_project(observations, horizon, M_train, n_folds, kern, discrete)

  # Grid for plug-in estimand and compensator evaluation
  grid = if (discrete) discrete_grid(mps[[1]]$mesh, du = mps[[1]]$du_bin) else continuous_grid(horizon, Q_comp)

  if (is.null(lambda_fits))
    lambda_fits = lapply(seq_along(mps), function(ff)
      fit_lambda(mps[[ff]], kern, eta_lam, lambda_disp,
                 log = indent(log, sprintf("fold %d/%d lambda", ff, length(mps)))))

  if (is.null(gamma_results))
    gamma_results = lapply(seq_along(mps), function(ff)
      fit_gamma(mps[[ff]], lambda_fits[[ff]], kern, eta_gam, estimand,
                gamma_disp_fn = gamma_disp_fn, gamma_base = gamma_base,
                log = indent(log, sprintf("fold %d/%d", ff, length(mps)))))

  fold_estimates = lapply(seq_along(mps), function(ff)
    estimate(mps[[ff]], lambda_fits[[ff]], gamma_results[[ff]], estimand, grid))

  aggregate_dr(mps, lambda_fits, gamma_results, fold_estimates,
               kern = kern, eta_gam = eta_gam, lambda_disp = lambda_disp,
               estimand = estimand, horizon = horizon, grid = grid,
               gamma_disp_fn = gamma_disp_fn, gamma_base = gamma_base)
}

# ============================================================
# Bootstrap
# ============================================================

bootstrap = function(obj, ...) UseMethod("bootstrap")

bootstrap.dr_result = function(obj, n_reps, log = null_logger, ...) {
  n = obj$n
  kern = obj$kern; lambda_disp = obj$lambda_disp
  estimand = obj$estimand; horizon = obj$horizon
  grid = obj$grid; eta_gam = obj$eta_gam
  gamma_disp_fn = obj$gamma_disp_fn

  dchis_lam = \(p) .dchis(lambda_disp, p)
  boot_ates = rep(NA_real_, n_reps)
  boot_ses = rep(NA_real_, n_reps)
  gamma_iters = integer(0)
  gamma_status = character(0)

  # Pre-compute cross-kernel matrices for gamma target projection
  K_cross_gamma = lapply(seq_along(obj$mps), function(ff) {
    mp = obj$mps[[ff]]
    dummy_haz = function(k, Z_ev) rep(0, nrow(atleast_2d(Z_ev)))
    emb_template = estimand$dpsi_grid(dummy_haz, mp$eval_Z, mp$mesh_u, mp$train_grid)
    kernel_matrix(mp$stacked_gamma$Z_pool, emb_template$Z, kern)
  })

  for (b in 1:n_reps) {
    t_boot = proc.time()
    w_subj = rpois(n, 1)
    all_terms = c()
    blog = indent(log, sprintf("boot %d/%d", b, n_reps))

    for (ff in seq_along(obj$mps)) {
      mp = obj$mps[[ff]]
      lam_fit = obj$lambda_fits[[ff]]
      gam_res = obj$gamma_results[[ff]]
      flog = indent(blog, sprintf("fold %d/%d", ff, length(obj$mps)))

      w_pool_lambda = w_subj[mp$lambda_idx][mp$i_pool_lambda]
      w_pool_gamma = w_subj[mp$gamma_idx][mp$i_pool_gamma]
      w_eval = w_subj[mp$eval_idx]

      # One-step Newton for lambda
      c_lambda_b = project_target(mp$stacked_lambda$Z_pool, kern,
                                  list(Z = mp$stacked_lambda$Z_pool,
                                       r = w_pool_lambda * mp$stacked_lambda$Y_pool),
                                  K_cross = lam_fit$K)
      tgt_b = split_target(c_lambda_b, nrow(mp$stacked_lambda$Z_pool))
      lambda_boot = onestep_bregman(lam_fit, tgt_b$target, tgt_b$target_null,
                                    w_pool_lambda)

      # Gamma target from bootstrap lambda
      phi_boot = \(Z) .phi(lambda_boot, Z)
      predict_haz_b = function(k, Z_ev) dchis_lam(phi_boot(cbind(mp$mesh[k], Z_ev)))
      emb_b = gamma_measure(estimand, predict_haz_b, mp, w_eval = w_eval)
      c_gamma_b = split_target(
        project_target(mp$stacked_gamma$Z_pool, kern, emb_b,
                       K_cross = K_cross_gamma[[ff]]),
        nrow(mp$stacked_gamma$Z_pool))

      # Re-solve gamma with bootstrap weights, warm-started from original
      gam_b = fit_gamma(mp, lambda_boot, kern, eta_gam, estimand,
                        gamma_disp_fn = gamma_disp_fn,
                        gamma_base = obj$gamma_base,
                        warm = gam_res$gamma_fit,
                        gamma_target = c_gamma_b$target,
                        gamma_target_null = c_gamma_b$target_null,
                        w = w_pool_gamma,
                        log = flog)
      gamma_iters  = c(gamma_iters,  gam_b$gamma_fit$iter)
      gamma_status = c(gamma_status, gam_b$gamma_fit$status)

      # Estimate with bootstrap models but original dchis_gamma
      lambda_boot_ = .bind(lambda_boot, mp$eval_Z)
      gamma_boot_ = .bind(gam_b$gamma_fit, mp$eval_Z)

      dr_b = compute_direct(phi_boot, lambda_disp, estimand, mp, grid)

      ct_b = correction_terms(mp, lambda_boot_, gamma_boot_,
                               lambda_disp, gam_res$dchis_gamma,
                               grid)

      dr_terms_b = dr_b$direct + ct_b$correction
      all_terms = c(all_terms, w_eval * dr_terms_b)
    }

    sw = sum(w_subj)
    boot_ates[b] = sum(all_terms) / sw
    nonzero = all_terms != 0
    boot_ses[b] = if (sum(nonzero) > 1) sd(all_terms[nonzero]) / sqrt(sw) else NA
    elapsed_b = (proc.time() - t_boot)[3]
    blog$info("est=%.5f  %.1fs", boot_ates[b], elapsed_b)
  }

  log$info("gamma iters: min=%d, median=%.0f, max=%d",
           min(gamma_iters), median(gamma_iters), max(gamma_iters))
  log$data("gamma_iters", gamma_iters)
  log$data("gamma_status", gamma_status)

  structure(
    list(boot_ates = boot_ates, boot_ses = boot_ses,
         est = obj$est, se_amle = obj$se_amle,
         gamma_iters = gamma_iters, gamma_status = gamma_status),
    class = "boot_result")
}

# ============================================================
# boot_result methods
# ============================================================

se.boot_result = function(obj, ...) sd(obj$boot_ates, na.rm = TRUE)

ci_t = function(obj, ...) UseMethod("ci_t")
ci_t.boot_result = function(obj, level = 0.95, ...) {
  good = !is.na(obj$boot_ates) & !is.nan(obj$boot_ates) & abs(obj$boot_ates) < 100 &
         !is.na(obj$boot_ses) & obj$boot_ses > 0
  if (sum(good) < 10) return(c(NA, NA))
  alpha = 1 - level
  t_stats = (obj$boot_ates[good] - obj$est) / obj$boot_ses[good]
  q = quantile(t_stats, c(alpha / 2, 1 - alpha / 2))
  c(obj$est - q[2] * obj$se_amle, obj$est - q[1] * obj$se_amle)
}

ci_pct = function(obj, ...) UseMethod("ci_pct")
ci_pct.boot_result = function(obj, level = 0.95, ...) {
  good = !is.na(obj$boot_ates) & !is.nan(obj$boot_ates) & abs(obj$boot_ates) < 100
  if (sum(good) < 10) return(c(NA, NA))
  alpha = 1 - level
  unname(quantile(obj$boot_ates[good], c(alpha / 2, 1 - alpha / 2)))
}

n_bad = function(obj, ...) UseMethod("n_bad")
n_bad.boot_result = function(obj, ...) {
  sum(is.na(obj$boot_ates) | is.nan(obj$boot_ates) | abs(obj$boot_ates) > 100)
}
