## Cross-fitting scaffold.
##
## Splits data into folds, calls estimator on each train/eval split.
## The estimator logic doesn't know about folds; cross-fitting orchestrates.

#' Create fold assignments.
#' Returns integer vector of fold IDs (1:n_folds), balanced.
make_folds = function(n, n_folds = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  rep(1:n_folds, length.out = n)[sample(n)]
}

#' Three-fold cross-fitting scaffold for single-outcome estimation.
#'
#' Folds: I_fit (train outcome model), I_sel (Lepski selection), I_psi (DR evaluation).
#' Rotates so each unit serves as evaluation exactly once.
#'
#' @param Y Outcome vector.
#' @param W Treatment vector (binary 0/1).
#' @param X Covariate matrix.
#' @param kern Kernel object.
#' @param estimand Estimand object (single_outcome type).
#' @param dispersion Dispersion for outcome model.
#' @param eta_grid Regularization grid for outcome model.
#' @param tuning_outcome Tuning strategy for outcome eta.
#' @param eta_gam Regularization for weight model.
#' @param n_folds Number of folds (default 3).
#' @return List with $est, $se, $eif, $selected_eta.
crossfit_single = function(Y, W, X, kern, estimand, dispersion,
                           eta_grid, tuning_outcome,
                           eta_gam = 0.1, n_folds = 3, seed = NULL) {
  n = length(Y)
  Z = cbind(W, X)
  fold_id = make_folds(n, n_folds, seed)
  tuning_outcome = as_tuning(tuning_outcome)

  eif_all = rep(NA, n)
  est_all = numeric(n_folds)
  selected_eta = numeric(n_folds)

  for (ff in 1:n_folds) {
    # Three-fold rotation: fit on ff, select on (ff%%3)+1, evaluate on ((ff+1)%%3)+1
    fit_idx = fold_id == ff
    sel_idx = fold_id == (ff %% n_folds) + 1
    eval_idx = fold_id == ((ff + 1) %% n_folds) + 1

    Z_fit = Z[fit_idx, , drop = FALSE]
    Y_fit = Y[fit_idx]
    Z_sel = Z[sel_idx, , drop = FALSE]
    Z_eval = Z[eval_idx, , drop = FALSE]

    # Pre-compute kernel matrix and targets (shared across eta grid)
    K_fit = kernel_matrix(Z_fit, Z_fit, kern)
    B_fit = null_basis(Z_fit, kern)
    tgt_fit = dpsi_from_response(Y_fit, K_fit, B_fit)

    # Fit outcome model at each eta
    models = vector("list", length(eta_grid))
    mu1_sel = array(0, c(sum(sel_idx), length(eta_grid)))
    mu0_sel = array(0, c(sum(sel_idx), length(eta_grid)))

    Z1_sel = Z_sel; Z1_sel[, 1] = 1
    Z0_sel = Z_sel; Z0_sel[, 1] = 0

    for (j in seq_along(eta_grid)) {
      models[[j]] = kernel_bregman(Z_fit, kern, eta_grid[j], dispersion,
                                    target = tgt_fit$target,
                                    target_null = tgt_fit$target_null,
                                    K = K_fit)
      mu1_sel[, j] = predict(models[[j]], Z1_sel)
      mu0_sel[, j] = predict(models[[j]], Z0_sel)
    }

    # Compute dir_val / dir_se for tuning
    dir_val = numeric(length(eta_grid))
    dir_se = numeric(length(eta_grid))
    for (j in seq_along(eta_grid)) {
      dir_val[j] = estimand$theta(mu1_sel[, j], mu0_sel[, j])
      dir_se[j] = estimand$se(mu1_sel[, j], mu0_sel[, j], sum(sel_idx))
    }

    # Select eta via tuning context (lazy: each strategy pulls what it needs)
    ctx = tuning_context(eta_grid, n = sum(sel_idx),
      dir_val = function() dir_val,
      dir_se = function() dir_se,
      loo_scores = function() vapply(models, loo_loss, numeric(1)),
      models = function() models)
    sel = tuning_outcome$select(ctx)
    selected_eta[ff] = sel$eta
    model_sel = models[[sel$idx]]

    # Predict on evaluation fold
    Z1_eval = Z_eval; Z1_eval[, 1] = 1
    Z0_eval = Z_eval; Z0_eval[, 1] = 0
    mu1_eval = predict(model_sel, Z1_eval)
    mu0_eval = predict(model_sel, Z0_eval)

    # Fit weights on evaluation fold
    gam = balancing_weights(Z_eval, kern, eta_gam)

    # DR estimate
    dr = estimand$dr(mu1_eval, mu0_eval, Y[eval_idx], W[eval_idx],
                     gam$gamma1, gam$gamma0)
    eif_all[eval_idx] = dr$eif
    est_all[ff] = dr$est
  }

  est = mean(est_all)
  se = sd(eif_all) / sqrt(n)
  list(est = est, se = se, eif = eif_all, selected_eta = selected_eta)
}
