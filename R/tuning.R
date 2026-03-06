#' Tuning: pluggable regularization selection.
#'
#' A tuning strategy provides $select(ctx), where ctx is a tuning context
#' with precomputed vectors from the grid sweep.
#'
#' Context fields:
#'   ctx$eta_grid      — the regularization grid
#'   ctx$n             — sample size
#'   ctx$dir_val()     — plug-in functional at each grid point
#'   ctx$dir_se()      — SE of plug-in functional
#'   ctx$loo_scores()  — CV loss at each grid point
#'
#' Or the tuning argument is a bare number / function, coerced by as_tuning().

# ============================================================
# Tuning context
# ============================================================

#' Build a tuning context from precomputed grid sweep results.
#'
#' @param eta_grid The regularization grid.
#' @param n Sample size.
#' @param ... Named closures or values for the context.
#'   Common: dir_val, dir_se, loo_scores (closures returning vectors).
#' @return A tuning context (list).
tuning_context = function(eta_grid, n, ...) {
  c(list(eta_grid = eta_grid, n = n), list(...))
}

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
    select = function(ctx) {
      idx = which.min(abs(ctx$eta_grid - eta))
      list(eta = ctx$eta_grid[idx], idx = idx)
    }
  )
}

#' Rule-based eta: eta = fn(n, ...).
rule_tuning = function(fn) {
  list(
    name = "rule",
    fn = fn,
    needs_grid = FALSE,
    select = function(ctx) {
      eta = fn(ctx$n)
      idx = which.min(abs(ctx$eta_grid - eta))
      list(eta = ctx$eta_grid[idx], idx = idx)
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
    select = function(ctx) {
      dir_val = ctx$dir_val()
      dir_se = ctx$dir_se()
      m = length(ctx$eta_grid)
      if (m <= 1) return(list(eta = ctx$eta_grid[1], idx = 1L))

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

      list(eta = ctx$eta_grid[selected], idx = selected)
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
    select = function(ctx) {
      dir_val = ctx$dir_val()
      ok = !is.na(dir_val)
      selected = which(ok)[which.min(abs(dir_val[ok] - truth))]
      list(eta = ctx$eta_grid[selected], idx = selected)
    }
  )
}

# ============================================================
# CV: cross-validated loss on the dual objective
# ============================================================

#' CV tuning: select eta minimizing LOO cross-validated dual loss.
#'
#' Works for any model (lambda or gamma) because LOO is computed
#' from the solver's own objective via the FOC, not from a response Y.
cv_tuning = function() {
  list(
    name = "cv",
    needs_grid = TRUE,
    select = function(ctx) {
      loo = ctx$loo_scores()
      selected = which.min(loo)
      list(eta = ctx$eta_grid[selected], idx = selected, cv_loss = loo)
    }
  )
}
