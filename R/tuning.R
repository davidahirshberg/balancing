## Tuning: pluggable regularization selection.
##
## A tuning strategy provides $select(fits, eta_grid, estimand, ...),
## which returns a list with $eta, $idx, $model, $mu1_hat, $mu0_hat.
##
## Or it's a bare number / function, coerced by as_tuning().

# ============================================================
# Coerce number/function to tuning strategy
# ============================================================

#' Coerce a tuning argument to a strategy object.
#' - numeric: fixed eta
#' - function: eta = fn(n, ...)
#' - list with $select: pass through
as_tuning = function(x) {
  if (is.list(x) && !is.null(x$select)) return(x)
  if (is.numeric(x)) return(fixed_tuning(x))
  if (is.function(x)) return(rule_tuning(x))
  stop("tuning must be a number, function, or tuning strategy (list with $select)")
}

#' Fixed eta (no selection).
fixed_tuning = function(eta) {
  list(
    name = "fixed",
    eta = eta,
    select = function(fits, eta_grid, estimand, ...) {
      # Find closest eta in grid, or just use the fixed value
      if (!missing(eta_grid) && !is.null(eta_grid)) {
        idx = which.min(abs(eta_grid - eta))
        list(eta = eta_grid[idx], idx = idx,
             model = fits$models[[idx]],
             mu1_hat = fits$mu1_hat[, idx],
             mu0_hat = fits$mu0_hat[, idx])
      } else {
        list(eta = eta, idx = NA, model = NULL, mu1_hat = NULL, mu0_hat = NULL)
      }
    }
  )
}

#' Rule-based eta: eta = fn(n, ...).
rule_tuning = function(fn) {
  list(
    name = "rule",
    fn = fn,
    select = function(fits, eta_grid, estimand, n = NULL, ...) {
      eta = fn(n, ...)
      if (!missing(eta_grid) && !is.null(eta_grid)) {
        idx = which.min(abs(eta_grid - eta))
        list(eta = eta_grid[idx], idx = idx,
             model = fits$models[[idx]],
             mu1_hat = fits$mu1_hat[, idx],
             mu0_hat = fits$mu0_hat[, idx])
      } else {
        list(eta = eta, idx = NA, model = NULL, mu1_hat = NULL, mu0_hat = NULL)
      }
    }
  )
}

# ============================================================
# Lepski: functional-stability selection
# ============================================================

#' Lepski tuning: select largest admissible eta where the functional is stable.
#'
#' Monitors estimand$theta(mu1, mu0) across the eta grid, selects the largest
#' eta where the functional hasn't drifted by more than C * estimand$se.
#'
#' Works for functionals whose plug-in "stabilizes then drifts" (TSM, ATE,
#' survival probability, risk ratio). Fails for variance-type functionals
#' (VCATE) where the plug-in monotonically shrinks with smoothing.
lepski_tuning = function(C = 3) {
  list(
    name = "lepski",
    C = C,
    select = function(fits, eta_grid, estimand, ...) {
      n = nrow(fits$mu1_hat)
      m = length(eta_grid)

      dir_val = numeric(m)
      dir_se = numeric(m)
      for (j in 1:m) {
        dir_val[j] = estimand$theta(fits$mu1_hat[, j], fits$mu0_hat[, j])
        dir_se[j] = estimand$se(fits$mu1_hat[, j], fits$mu0_hat[, j], n)
      }

      # Lepski: largest admissible eta
      selected = 1
      for (j in 2:m) {
        admissible = TRUE
        for (i in 1:(j - 1)) {
          if (abs(dir_val[i] - dir_val[j]) > C * dir_se[i]) {
            admissible = FALSE; break
          }
        }
        if (admissible) selected = j
      }

      list(eta = eta_grid[selected], idx = selected,
           model = fits$models[[selected]],
           mu1_hat = fits$mu1_hat[, selected],
           mu0_hat = fits$mu0_hat[, selected],
           all_dir_val = dir_val, all_dir_se = dir_se)
    }
  )
}

# ============================================================
# Oracle: pick best eta using known truth (simulation only)
# ============================================================

#' Oracle tuning: pick eta minimizing |estimate - truth|.
#' Only for simulations where truth is known.
oracle_tuning = function(truth) {
  list(
    name = "oracle",
    truth = truth,
    select = function(fits, eta_grid, estimand, ...) {
      m = length(eta_grid)
      dir_val = numeric(m)
      for (j in 1:m)
        dir_val[j] = estimand$theta(fits$mu1_hat[, j], fits$mu0_hat[, j])

      selected = which.min(abs(dir_val - truth))

      list(eta = eta_grid[selected], idx = selected,
           model = fits$models[[selected]],
           mu1_hat = fits$mu1_hat[, selected],
           mu0_hat = fits$mu0_hat[, selected],
           all_dir_val = dir_val)
    }
  )
}

# ============================================================
# CV: cross-validated loss
# ============================================================

#' CV tuning: select eta minimizing cross-validated prediction loss.
#' Uses the dispersion's loss (chis(phi) - Y*phi) as the CV criterion.
cv_tuning = function(folds = 5) {
  list(
    name = "cv",
    folds = folds,
    # CV selection requires refitting at each eta on each fold,
    # so it takes raw data rather than pre-fitted models.
    select_from_data = function(Y, Z, kern, eta_grid, dispersion,
                                intercept = TRUE, w = NULL) {
      n = length(Y)
      fold_id = sample(rep(1:folds, length.out = n))
      m = length(eta_grid)
      cv_loss = numeric(m)

      for (j in 1:m) {
        fold_losses = numeric(folds)
        for (ff in 1:folds) {
          train_idx = fold_id != ff
          test_idx = fold_id == ff
          fit = kernel_bregman(Y[train_idx], atleast_2d(Z)[train_idx, , drop = FALSE],
                               kern, eta_grid[j], dispersion,
                               intercept = intercept,
                               w = if (!is.null(w)) w[train_idx] else NULL)
          phi_test = predict_phi(fit, atleast_2d(Z)[test_idx, , drop = FALSE])
          # Loss: chis(phi) - Y*phi
          fold_losses[ff] = mean(dispersion$chis(phi_test) - Y[test_idx] * phi_test)
        }
        cv_loss[j] = mean(fold_losses)
      }

      selected = which.min(cv_loss)
      list(eta = eta_grid[selected], idx = selected, cv_loss = cv_loss)
    }
  )
}
