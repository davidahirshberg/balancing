#' # Survival DR Estimator
#'
#' Doubly-robust estimation of survival causal effects with
#' crossfitting and Poisson multiplier bootstrap.
#'
#' **Paper** (spinoff3, eq:one-step-general): the DR estimator is
#' \deqn{\\hat\\psi = \\hat P\\bigl[\\psi_Z^1(\\hat\\lambda) - \\psi_Z^\{-1\}(\\hat\\lambda)\\bigr]
#'   + \\hat P\\,\\hat\\gamma(Z)\\\{Y - \\hat\\lambda(Z)\\\}}
#' The first term is the plug-in (direct) estimate using the
#' hazard model \eqn{\\hat\\lambda}; the second is the weighted-residual
#' correction using balancing weights \eqn{\\hat\\gamma}.
#'
#' With crossfitting (eq:crossfit):
#' \deqn{\\hat\\psi = \\frac\{1\}\{n\}\\sum_i \\hat V_i, \\quad
#'   \\hat V_i = \\psi_\{Z_i\}(\\hat\\lambda_\{-i\}) +
#'   E_U\\,\\hat\\gamma_i(U)\\\{Y_i^U - \\hat\\lambda_U(Z_i)\\\}}
#'
#' **Pipeline**: `mesh_project` \eqn{\\to} `fit_lambda` \eqn{\\to} `fit_gamma`
#' \eqn{\\to} `estimate` \eqn{\\to} `pool_folds`, wrapped by `survival_effect`.
#'
#' ## Per-fold primitives
#'
#' - `fit_lambda`: hazard model on fold \eqn{I_\\lambda}
#' - `fit_gamma`: weight model on fold \eqn{I_\\psi}, given \eqn{\\hat\\lambda}
#' - `estimate`: plug-in + correction on fold \eqn{I_\\psi}
#'
#' ## Aggregation
#'
#' - `pool_folds`: pool per-fold \eqn{\\hat V_i} into a `effect_estimate`
#'
#' ## Bootstrap
#'
#' - `bootstrap.effect_estimate`: Poisson multiplier bootstrap with
#'   one-step Newton refit

# ============================================================
# Internal helpers
# ============================================================

#' ### `gamma_measure`
#'
#' **Paper**: the embedding vector \eqn{c_\\gamma} for the weight
#' model encodes \eqn{\\nabla_\\phi\\,\\dot\\psi_Z}, the gradient of the
#' estimand functional w.r.t. the dual variable \eqn{\\phi_\\gamma}.
#' For survival ATE (eq:cgam):
#' \deqn{c_\{\\gamma,jm'\} = \\frac\{1\}\{|I_\\psi|\\,M\}\\sum_\{i \\in I_\\psi\}
#'   \\sum_w w\\,\\hat S_t(w, X_i) \\sum_m K\\bigl((u_\{m'\}, W_j, X_j),\\,
#'   (u_m, w, X_i)\\bigr)}
#'
#' **Code**: calls `estimand$dpsi_grid(predict_haz, ...)` to get
#' the signed measure \eqn{(Z^\{\\dot\\psi\}, r)}, then passes to
#' `project_target`.  Raw sum scale — no normalization by \eqn{|I_\\psi|},
#' matching the solver's sum-scale objective.
gamma_measure = function(estimand, predict_haz, mp, w_eval = NULL) {
  emb = estimand$dpsi_grid(predict_haz, mp$eval_Z, mp$mesh_u, mp$train_grid)
  if (!is.null(w_eval)) {
    n_eval = nrow(atleast_2d(mp$eval_Z))
    n_per_subj = length(emb$r) %/% n_eval
    emb$r = emb$r * rep(w_eval, times = n_per_subj)
  }
  emb
}

#' ### `compute_direct`
#'
#' **Paper**: the plug-in estimate
#' \eqn{\\hat P[\\psi_Z^w(\\hat\\lambda)]} evaluates the estimand
#' functional at the fitted hazard \eqn{\\hat\\lambda_u(w, X_i) =
#' \\dot\\chi^*_\\lambda(\\hat\\phi_\\lambda(u, w, X_i))}.
#'
#' **Code**: `predict_haz(k, Z_ev)` returns
#' \eqn{\\hat\\lambda_\{u_k\}(Z_\{\\mathrm\{ev\}\})} on the training mesh.
#' `predict_dLambda(k, Z_ev)` rescales to the evaluation grid
#' (for continuous: `scale = du_eval / du_bin`).
#' `estimand$direct(predict_dLambda, eval_Z, grid)` computes
#' the per-subject direct estimate.
compute_direct = function(phi, lambda_disp, estimand, mp, grid) {
  dchis = \(p) .dchis(lambda_disp, p)

  # predict_haz on training mesh: model output = dLambda on training grid.
  # Used by gamma_measure / build_gamma_correction (always on training mesh).
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

#' ### `correction_terms`
#'
#' **Paper**: the DR correction (second term of eq:one-step-general):
#' \deqn{E_U\\,\\hat\\gamma_i(U)\\\{Y_i^U - \\hat\\lambda_U(Z_i)\\\}
#'   = \\underbrace\{\\hat\\gamma_i(T_i)\\,\\mathbf\{1\}(D_i=1)\}_\{\\text\{Dirac\}\}
#'   - \\underbrace\{\\int_0^\{T_i^*\} \\hat\\gamma_u(W_i, X_i)\\,
#'     \\hat\\lambda_u(W_i, X_i)\\,du\}_\{\\text\{compensator\}\}}
#'
#' **Code**: `compensator(...)` computes the integral via
#' quadrature. The Dirac term evaluates \eqn{\\hat\\gamma} at each
#' subject's event time using `dchis_gamma(Z_ev, phi_ev)`.
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

#' ### `build_gamma_correction`
#'
#' **Paper**: the weight model's correction prediction is
#' \eqn{\\hat\\gamma_u(W_i, X_i) = \\dot\\chi^*_Z(\\hat\\phi_\\gamma(u, W_i, X_i))}
#' where the Z-dependent dispersion uses the Riesz representer
#' offset \eqn{r_\{im\} = \\dot\\psi_\{Z_i\}(\\hat\\lambda)(u_m)}.
#'
#' For the ATE with target-scaled entropy (TSE):
#' \deqn{\\hat\\gamma_\{im\} = r_\{im\} + \\mathrm\{sign\}(r_\{im\})\\,
#'   e^\{\\mathrm\{sign\}(r_\{im\})\\,\\hat\\phi_\{\\gamma,im\}\}}
#' Floors \eqn{|\\hat\\gamma|} at \eqn{|r|}, then exponential on the excess.
#'
#' **Code**: returns a closure `dchis_gamma(Z, phi)` that
#' (1) looks up \eqn{r_\{im\}} for each \eqn{(i, m)} pair from the
#' precomputed linearization grid, and
#' (2) evaluates `.dchis(dispersion(gamma_base, offset=r, sigma=sigma), phi)`.
#' For ATE estimands, \eqn{\\sigma = -W} (from `gamma_disp_sigma`);
#' for TSM, \eqn{\\sigma = \\mathrm\{sign\}(r)}.
build_gamma_correction = function(predict_haz, estimand_obj, eval_Z,
                                train_grid, gamma_base = entropy_base,
                                gamma_disp_sigma = estimand_obj$gamma_disp_sigma) {
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
    sigma = if (!is.null(gamma_disp_sigma)) gamma_disp_sigma(Z[, 2]) else sign(r)
    .dchis(dispersion(gamma_base, offset = r, sigma = sigma), phi)
  }
}

#' ## `mesh_project`
#'
#' **Paper**: discretize the time axis \eqn{[0, \\bar t]} into \eqn{M}
#' bins with mesh points \eqn{u_1, \\ldots, u_M}. For each subject
#' \eqn{i} with observed time \eqn{\\tilde T_i}, the at-risk set is
#' \eqn{\\mathcal\{D\}_i = \\\{m : \\tilde T_i \\ge u_m\\\}}. The stacked
#' data \eqn{\\mathcal\{D\} = \\\{(i, m) : i \\in I_\\lambda,\\, m \\in
#' \\mathcal\{D\}_i\\\}} has response \eqn{Y_\{im\} =
#' \\mathbf\{1\}(\\tilde T_i \\in [u_\{m-1\}, u_m),\\, D_i = 1)}.
#'
#' **Code**: one `mesh_projection` per crossfit fold, containing:
#' - `stacked_lambda`: at-risk pairs for hazard training
#' - `stacked_gamma`: at-risk pairs for weight estimation
#' - `K_lambda`, `K_gamma`: kernel matrices on stacked data
#' - `tgt_lambda`: projected target \eqn{c_\\lambda = K_\{\\mathcal\{D\}\}\\,Y}
#' - `eval_Z`, `T_obs_eval`, `D_eval`: evaluation fold data
mesh_project = function(observations, horizon, M_train, n_folds, kern,
                        kern_gamma = kern, discrete = FALSE) {
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

    K_gamma = kernel_matrix(stacked_gamma$Z_pool, stacked_gamma$Z_pool, kern_gamma)

    eval_Z = Z[eval_idx, , drop = FALSE]

    train_grid = if (discrete)
      discrete_grid(lambda_train$mesh, du = lambda_train$du)
    else
      continuous_training_grid(lambda_train$mesh, du = lambda_train$du)

    list(lambda_idx = lambda_idx, gamma_idx = gamma_idx, eval_idx = eval_idx,
         stacked_lambda = stacked_lambda, stacked_gamma = stacked_gamma,
         mesh = lambda_train$mesh, du_bin = stacked_lambda$du_bin,
         mesh_u = unique(stacked_gamma$u_pool),
         K_lambda = K_lambda, K_gamma = K_gamma, tgt_lambda = tgt_lambda,
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

#' ## `fit_lambda`
#'
#' **Paper** (eq:hazard-alpha): fit the hazard model on stacked
#' at-risk pairs \eqn{\\mathcal\{D\}}:
#' \deqn{\\min_\\alpha \\frac\{1\}\{|\\mathcal\{D\}|\}\\sum_\{(i,m) \\in \\mathcal\{D\}\}
#'   \\bigl[\\chi^*_\\lambda(\\phi_\{im\}) - Y_\{im\}\\,\\phi_\{im\}\\bigr]
#'   + \\eta_\\lambda\\,\\alpha^\\top K_\{\\mathcal\{D\}\}\\,\\alpha}
#' where \eqn{\\phi = K_\{\\mathcal\{D\}\}\\,\\alpha + B\\beta} and
#' \eqn{\\hat\\lambda_u(Z_i) = \\dot\\chi^*_\\lambda(\\phi_\{im\})}.
#'
#' **Code**: wraps `kernel_bregman` with `target` \eqn{= K_\{\\mathcal\{D\}\}\\,Y}
#' (precomputed in `mesh_project`).
fit_lambda = function(mp, kern, eta, lambda_disp, warm = NULL,
                      K_chol = NULL, maxiter = 100, log = null_logger) {
  kernel_bregman(mp$stacked_lambda$Z_pool, kern,
                 eta = eta, dispersion = lambda_disp,
                 target = mp$tgt_lambda$target,
                 target_null = mp$tgt_lambda$target_null,
                 K = mp$K_lambda, K_chol = K_chol,
                 alpha0 = warm$alpha, beta0 = warm$beta,
                 maxiter = maxiter,
                 log = log)
}

#' ## `fit_gamma`
#'
#' **Paper** (eq:weight-beta): fit the weight model on stacked
#' at-risk pairs \eqn{\\mathcal\{E\}}:
#' \deqn{\\min_\\beta \\frac\{1\}\{|\\mathcal\{E\}|\}\\sum_\{(i,m) \\in \\mathcal\{E\}\}
#'   \\chi^*_Z(\\phi_\{im\}) - c_\\gamma^\\top\\beta
#'   + \\eta_\\gamma\\,\\beta^\\top K_\{\\mathcal\{E\}\}\\,\\beta}
#' The target \eqn{c_\\gamma} embeds \eqn{\\nabla_\\phi\\,\\dot\\psi_Z}, which
#' evaluates at **counterfactual** arms \eqn{(u, w, X_i)} for both
#' \eqn{w \\in \\\{-1, +1\\\}}, weighted by \eqn{w\\,\\hat S_t(w, X_i)}.
#'
#' **Code**: calls `gamma_measure` to build the signed measure
#' \eqn{(Z^\{\\dot\\psi\}, r)}, projects onto the representer basis via
#' `project_target`, then solves with `kernel_bregman`.
#' The `gamma_disp_fn` maps the target \eqn{r} to a Z-dependent
#' dispersion (default: target-scaled entropy).
fit_gamma = function(mp, lambda_fit, kern, eta, estimand,
                     gamma_disp_fn = target_scaled_entropy_dispersion,
                     gamma_base = entropy_base,
                     warm = NULL,
                     gamma_target = NULL, gamma_target_null = NULL,
                     w = NULL, w_eval = NULL,
                     K_cross = NULL, dchis_gamma = NULL,
                     K_chol_gamma = NULL,
                     maxiter = 100,
                     log = null_logger) {
  lambda_fit = .bind(lambda_fit, mp$eval_Z)
  phi   = \(Z) .phi(lambda_fit, Z)
  dchis = \(p) .dchis(lambda_fit$dispersion, p)
  predict_haz = function(k, Z_ev) dchis(phi(cbind(mp$mesh[k], Z_ev)))

  if (is.null(dchis_gamma))
    dchis_gamma = build_gamma_correction(predict_haz, estimand, mp$eval_Z,
                                       mp$train_grid,
                                       gamma_base = gamma_base)

  if (is.null(gamma_target)) {
    emb = gamma_measure(estimand, predict_haz, mp, w_eval = w_eval)
    tgt = split_target(project_target(mp$stacked_gamma$Z_pool, kern, emb,
                                      K_cross = K_cross),
                       nrow(mp$stacked_gamma$Z_pool))
    gamma_target = tgt$target
    gamma_target_null = tgt$target_null
  }
  # TSE offset: per-observation scale (divide sum-scale target by n_eval).
  # The solver target is sum-scale to match sum(chi*), but the TSE offset
  # is the expected per-observation gamma and must be O(1).
  n_eval = nrow(atleast_2d(mp$eval_Z))
  gamma_offset = gamma_target / n_eval

  # For ATE estimands, use sigma = -W (from the Riesz representer
  # sign structure) instead of sign(gamma_target).  This makes the
  # TSE sign constraint robust to hazard estimates where
  # lambda_hat > 1 could flip sign(r_hat).  TSM estimands lack
  # gamma_disp_sigma and fall through to the default TSE.
  if (!is.null(estimand$gamma_disp_sigma)) {
    A_basis = mp$stacked_gamma$Z_pool[, 2]
    sigma = estimand$gamma_disp_sigma(A_basis)
    gamma_disp = dispersion(gamma_base, offset = gamma_offset, sigma = sigma)
  } else {
    gamma_disp = gamma_disp_fn(gamma_offset)
  }

  K_gam = if (!is.null(warm$K)) warm$K
          else if (!is.null(mp$K_gamma)) mp$K_gamma
          else NULL
  gamma_fit = kernel_bregman(mp$stacked_gamma$Z_pool, kern,
                             eta = eta, dispersion = gamma_disp,
                             target = gamma_target,
                             target_null = gamma_target_null,
                             w = w,
                             K = K_gam, K_chol = K_chol_gamma,
                             alpha0 = warm$alpha, beta0 = warm$beta,
                             maxiter = maxiter,
                             n = n_eval,
                             log = indent(log, "gamma"))

  list(gamma_fit = gamma_fit, dchis_gamma = dchis_gamma,
       predict_haz = predict_haz,
       gamma_target = gamma_target, gamma_target_null = gamma_target_null)
}

#' ## `estimate`
#'
#' **Paper**: compute the per-subject influence function terms
#' \deqn{\\hat V_i = \\psi_\{Z_i\}(\\hat\\lambda)
#'   + E_U\\,\\hat\\gamma_i(U)\\\{Y_i^U - \\hat\\lambda_U(Z_i)\\\}}
#' The first term is `direct` (plug-in at \eqn{\\hat\\lambda}); the
#' second is `correction` (Dirac minus compensator).
#' Returns `terms` \eqn{= \\hat V_i} for aggregation.
estimate = function(mp, lambda_fit, gamma_result, estimand, grid) {
  lambda_disp = lambda_fit$dispersion
  lambda_fit_ = .bind(lambda_fit, mp$eval_Z)
  gamma_fit_ = .bind(gamma_result$gamma_fit, mp$eval_Z)
  phi = \(Z) .phi(lambda_fit_, Z)

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

#' ## `pool_folds`
#'
#' **Paper** (eq:crossfit): pool crossfit folds:
#' \deqn{\\hat\\psi = \\frac\{1\}\{n\}\\sum_i \\hat V_i, \\quad
#'   \\hat\\sigma^2 = \\frac\{1\}\{n\}\\sum_i \\hat V_i^2}
#'
#' **Code**: collects per-fold `terms` \eqn{= \\hat V_i} into a
#' length-\eqn{n} vector. Returns `est` \eqn{= \\bar V}, `se` \eqn{=
#' \\mathrm\{sd\}(V)/\\sqrt\{n\}}, and `se_amle` \eqn{= \\sqrt\{(\\mathrm\{Var\}(
#' \\text\{direct\}) + \\bar\{\\text\{noise\\_var\}\})/n\}} (AMLE variance
#' estimate using plug-in variance + noise variance separately).
pool_folds = function(mps, lambda_fits, gamma_results, fold_estimates,
                        kern, eta_gam, lambda_disp, estimand,
                        horizon, grid,
                        kern_gamma = kern,
                        gamma_disp_fn = target_scaled_entropy_dispersion,
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
         kern = kern, kern_gamma = kern_gamma,
         eta_gam = eta_gam, lambda_disp = lambda_disp,
         estimand = estimand,
         horizon = horizon, grid = grid,
         gamma_disp_fn = gamma_disp_fn, gamma_base = gamma_base),
    class = "effect_estimate")
}

# ============================================================
# Convenience wrapper
# ============================================================

#' ## `survival_effect`
#'
#' Convenience wrapper: runs the full pipeline
#' `mesh_project` \eqn{\\to} `fit_lambda` \eqn{\\to} `fit_gamma` \eqn{\\to}
#' `estimate` \eqn{\\to} `pool_folds`.
#'
#' Accepts pre-computed intermediate results (`mps`,
#' `lambda_fits`, `gamma_results`) to allow reuse across
#' regularization grid sweeps.
survival_effect = function(observations, kern, estimand, lambda_disp,
                           eta_lam, eta_gam, horizon,
                           M_train = 15, n_folds = 2, Q_comp = 50,
                           discrete = FALSE, log = null_logger,
                           gamma_disp_fn = target_scaled_entropy_dispersion,
                           gamma_base = entropy_base,
                           kern_gamma = kern,
                           mps = NULL, lambda_fits = NULL, gamma_results = NULL) {
  if (is.null(mps))
    mps = mesh_project(observations, horizon, M_train, n_folds, kern,
                       kern_gamma = kern_gamma, discrete = discrete)

  # Grid for plug-in estimand and compensator evaluation
  grid = if (discrete) discrete_grid(mps[[1]]$mesh, du = mps[[1]]$du_bin) else continuous_grid(horizon, Q_comp)

  if (is.null(lambda_fits))
    lambda_fits = lapply(seq_along(mps), function(ff)
      fit_lambda(mps[[ff]], kern, eta_lam, lambda_disp,
                 log = indent(log, sprintf("fold %d/%d lambda", ff, length(mps)))))

  if (is.null(gamma_results))
    gamma_results = lapply(seq_along(mps), function(ff)
      fit_gamma(mps[[ff]], lambda_fits[[ff]], kern_gamma, eta_gam, estimand,
                gamma_disp_fn = gamma_disp_fn, gamma_base = gamma_base,
                log = indent(log, sprintf("fold %d/%d", ff, length(mps)))))

  fold_estimates = lapply(seq_along(mps), function(ff)
    estimate(mps[[ff]], lambda_fits[[ff]], gamma_results[[ff]], estimand, grid))

  pool_folds(mps, lambda_fits, gamma_results, fold_estimates,
               kern = kern, eta_gam = eta_gam, lambda_disp = lambda_disp,
               estimand = estimand, horizon = horizon, grid = grid,
               kern_gamma = kern_gamma,
               gamma_disp_fn = gamma_disp_fn, gamma_base = gamma_base)
}

#' ## Bootstrap
#'
#' ### `bootstrap.effect_estimate`
#'
#' **Paper** (spinoff3, sec:bootstrap): Poisson multiplier
#' bootstrap. For each replicate \eqn{b}:
#'
#' 1. Draw \eqn{\\xi_i \\sim \\mathrm\{Poisson\}(1)}, \eqn{i = 1,\\ldots,n}.
#' 2. One-step Newton refit of \eqn{\\hat\\lambda^*}: reweight the
#'    target \eqn{c_\\lambda^* = K_\{\\mathcal\{D\}\}(\\xi \\odot Y)} and
#'    observation weights \eqn{w^* = \\xi}, then
#'    \eqn{\\alpha_\\lambda^* = \\hat\\alpha_\\lambda - M_\\lambda^\{-1\}
#'    \\nabla L_\\lambda^*(\\hat\\alpha_\\lambda)} using cached \eqn{M}.
#' 3. Refit \eqn{\\hat\\gamma^*} with updated survival predictions
#'    from \eqn{\\hat\\lambda^*} and reweighted targets (`maxiter = 1`).
#' 4. Compute DR estimate \eqn{\\hat\\psi^* = n^\{-1\}\\sum_i \\xi_i\\hat V_i^*}.
#'
#' **Code**: pre-computes cross-kernel matrices `K_cross_gamma`
#' for target projection (avoids recomputing per replicate).
#' Lambda refit uses `onestep_bregman`; gamma refit uses
#' `fit_gamma` with `maxiter = 1` (warm-started from original).
bootstrap = function(obj, ...) UseMethod("bootstrap")

bootstrap.effect_estimate = function(obj, n_reps, log = null_logger, ...) {
  n = obj$n
  kern = obj$kern; lambda_disp = obj$lambda_disp
  kern_gamma = if (!is.null(obj$kern_gamma)) obj$kern_gamma else kern
  estimand = obj$estimand; horizon = obj$horizon
  grid = obj$grid; eta_gam = obj$eta_gam
  gamma_disp_fn = obj$gamma_disp_fn

  boot_ates = rep(NA_real_, n_reps)
  boot_ses = rep(NA_real_, n_reps)
  # Pre-compute cross-kernel matrices for gamma target projection
  K_cross_gamma = lapply(seq_along(obj$mps), function(ff) {
    mp = obj$mps[[ff]]
    dummy_haz = function(k, Z_ev) rep(0, nrow(atleast_2d(Z_ev)))
    emb_template = estimand$dpsi_grid(dummy_haz, mp$eval_Z, mp$mesh_u, mp$train_grid)
    kernel_matrix(mp$stacked_gamma$Z_pool, emb_template$Z, kern_gamma)
  })

  # Per-step timing accumulators
  t_acc = list(proj_lam = 0, onestep_lam = 0,
               fit_gamma = 0, estimate = 0)
  tick = function() proc.time()[3]

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
      t0 = tick()
      c_lambda_b = project_target(mp$stacked_lambda$Z_pool, kern,
                                  list(Z = mp$stacked_lambda$Z_pool,
                                       r = w_pool_lambda * mp$stacked_lambda$Y_pool),
                                  K_cross = lam_fit$K)
      tgt_b = split_target(c_lambda_b, nrow(mp$stacked_lambda$Z_pool))
      t1 = tick(); t_acc$proj_lam = t_acc$proj_lam + (t1 - t0)

      lambda_boot = onestep_bregman(lam_fit, tgt_b$target, tgt_b$target_null,
                                    w_pool_lambda)
      t2 = tick(); t_acc$onestep_lam = t_acc$onestep_lam + (t2 - t1)

      # One-step gamma via fit_gamma (reuses dchis_gamma, pre-computed K_cross)
      gam_res_b = fit_gamma(mp, lambda_boot, kern_gamma, eta_gam, estimand,
                            gamma_disp_fn = gamma_disp_fn,
                            warm = gam_res$gamma_fit,
                            w = w_pool_gamma, w_eval = w_eval,
                            K_cross = K_cross_gamma[[ff]],
                            dchis_gamma = gam_res$dchis_gamma,
                            K_chol_gamma = gam_res$gamma_fit$K_chol,
                            maxiter = 1, log = flog)
      t3 = tick(); t_acc$fit_gamma = t_acc$fit_gamma + (t3 - t2)

      # DR estimate via shared estimate() path
      est_b = estimate(mp, lambda_boot, gam_res_b, estimand, grid)
      t4 = tick(); t_acc$estimate = t_acc$estimate + (t4 - t3)

      dr_terms_b = est_b$terms
      all_terms = c(all_terms, w_eval * dr_terms_b)
    }

    sw = sum(w_subj)
    boot_ates[b] = sum(all_terms) / sw
    nonzero = all_terms != 0
    boot_ses[b] = if (sum(nonzero) > 1) sd(all_terms[nonzero]) / sqrt(sw) else NA
    elapsed_b = (proc.time() - t_boot)[3]
    blog$info("est=%.5f  %.1fs", boot_ates[b], elapsed_b)
  }

  # Report per-step timing
  total_acc = Reduce(`+`, t_acc)
  if (n_reps > 0 && total_acc > 0) {
    cat(sprintf("\n  Bootstrap step timing (%d reps × %d folds):\n", n_reps, length(obj$mps)))
    for (nm in names(t_acc))
      cat(sprintf("    %-16s %6.3fs  %5.1f%%  (%.4fs/rep)\n",
                  nm, t_acc[[nm]], 100 * t_acc[[nm]] / total_acc, t_acc[[nm]] / n_reps))
    cat(sprintf("    %-16s %6.3fs\n", "TOTAL", total_acc))
  }

  structure(
    list(boot_ates = boot_ates, boot_ses = boot_ses,
         est = obj$est, se_amle = obj$se_amle),
    class = "bootstrapped_estimate")
}

# ============================================================
# bootstrapped_estimate methods
# ============================================================

bootstrap_se = function(obj) sd(obj$boot_ates, na.rm = TRUE)

#' Bootstrap confidence interval (t-method).
#' @param obj A `bootstrapped_estimate` object.
#' @param level Confidence level (default 0.95).
#' @param ... Additional arguments (unused).
#' @return Length-2 vector (lower, upper).
bootstrap_t_interval = function(obj, ...) UseMethod("bootstrap_t_interval")
bootstrap_t_interval.bootstrapped_estimate = function(obj, level = 0.95, ...) {
  good = !is.na(obj$boot_ates) & !is.nan(obj$boot_ates) & abs(obj$boot_ates) < 100 &
         !is.na(obj$boot_ses) & obj$boot_ses > 0
  if (sum(good) < 10) return(c(NA, NA))
  alpha = 1 - level
  t_stats = (obj$boot_ates[good] - obj$est) / obj$boot_ses[good]
  q = quantile(t_stats, c(alpha / 2, 1 - alpha / 2))
  c(obj$est - q[2] * obj$se_amle, obj$est - q[1] * obj$se_amle)
}

#' Bootstrap confidence interval (percentile method).
#' @param obj A `bootstrapped_estimate` object.
#' @param level Confidence level (default 0.95).
#' @param ... Additional arguments (unused).
#' @return Length-2 vector (lower, upper).
percentile_bootstrap_interval = function(obj, ...) UseMethod("percentile_bootstrap_interval")
percentile_bootstrap_interval.bootstrapped_estimate = function(obj, level = 0.95, ...) {
  good = !is.na(obj$boot_ates) & !is.nan(obj$boot_ates) & abs(obj$boot_ates) < 100
  if (sum(good) < 10) return(c(NA, NA))
  alpha = 1 - level
  unname(quantile(obj$boot_ates[good], c(alpha / 2, 1 - alpha / 2)))
}

