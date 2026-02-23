## Tuning: pluggable regularization selection.
##
## A tuning strategy provides $select(dir_val, dir_se, eta_grid),
## which returns a list with $idx and $eta (the selected grid point).
##
## Or it's a bare number / function, coerced by as_tuning().
##
## The pipeline computes dir_val[j] = plug-in functional at eta_grid[j]
## and dir_se[j] = its standard error, then hands these to the tuning
## strategy. The strategy knows nothing about survival vs single-outcome.

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
    needs_grid = FALSE,
    select = function(dir_val, dir_se, eta_grid, ...) {
      idx = which.min(abs(eta_grid - eta))
      list(eta = eta_grid[idx], idx = idx)
    }
  )
}

#' Rule-based eta: eta = fn(n, ...).
rule_tuning = function(fn) {
  list(
    name = "rule",
    fn = fn,
    needs_grid = FALSE,
    select = function(dir_val, dir_se, eta_grid, n = NULL, ...) {
      eta = fn(n, ...)
      idx = which.min(abs(eta_grid - eta))
      list(eta = eta_grid[idx], idx = idx)
    }
  )
}

# ============================================================
# Lepski: functional-stability selection
# ============================================================

#' Lepski tuning: select largest admissible eta where the functional is stable.
#'
#' Monitors dir_val across the eta grid, selects the largest eta where the
#' functional hasn't drifted by more than C * dir_se from any smaller eta.
#'
#' Works for functionals whose plug-in "stabilizes then drifts" (TSM, ATE,
#' survival probability, risk ratio). Fails for variance-type functionals
#' (VCATE) where the plug-in monotonically shrinks with smoothing.
lepski_tuning = function(C = 3) {
  list(
    name = "lepski",
    C = C,
    needs_grid = TRUE,
    select = function(dir_val, dir_se, eta_grid, ...) {
      m = length(eta_grid)

      # Lepski: largest admissible eta
      selected = 1
      for (j in 2:m) {
        if (is.na(dir_val[j])) next
        admissible = TRUE
        for (i in 1:(j - 1)) {
          if (is.na(dir_val[i])) next
          if (abs(dir_val[i] - dir_val[j]) > C * dir_se[i]) {
            admissible = FALSE; break
          }
        }
        if (admissible) selected = j
      }

      list(eta = eta_grid[selected], idx = selected)
    }
  )
}

# ============================================================
# Oracle: pick best eta using known truth (simulation only)
# ============================================================

#' Oracle tuning: pick eta minimizing |dir_val - truth|.
#' Only for simulations where truth is known.
oracle_tuning = function(truth) {
  list(
    name = "oracle",
    truth = truth,
    needs_grid = TRUE,
    select = function(dir_val, dir_se, eta_grid, ...) {
      ok = !is.na(dir_val)
      selected = which(ok)[which.min(abs(dir_val[ok] - truth))]
      list(eta = eta_grid[selected], idx = selected)
    }
  )
}

# ============================================================
# CV: cross-validated loss (works differently — needs raw data)
# ============================================================

#' CV tuning: select eta minimizing cross-validated prediction loss.
#' Uses the dispersion's loss (chis(phi) - Y*phi) as the CV criterion.
#'
#' Unlike other strategies, CV needs raw data to refit on sub-folds.
#' The pipeline calls $select_from_data when tuning$name == "cv".
cv_tuning = function(folds = 5) {
  list(
    name = "cv",
    folds = folds,
    needs_grid = TRUE,
    # Standard select interface — falls back to this if dir_val is provided
    select = function(dir_val, dir_se, eta_grid) {
      stop("CV tuning requires select_from_data, not select")
    },
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
