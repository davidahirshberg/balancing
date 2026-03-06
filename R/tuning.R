#' @title Tuning strategies for regularization selection
#'
#' @description
#' Pluggable tuning strategies for selecting the regularization parameter
#' \eqn{\eta} in the kernel Bregman balancing solver. Each strategy is a
#' list with a \code{$select(ctx)} method that takes a tuning context and
#' returns the selected \eqn{\eta} and its grid index.
#'
#' Available strategies:
#' \itemize{
#'   \item \code{\link{lepski_tuning}}: functional-stability (Lepski) selection
#'   \item \code{\link{cv_tuning}}: leave-one-out cross-validated dual loss
#'   \item \code{\link{oracle_tuning}}: best \eqn{\eta} given known truth (simulation only)
#'   \item \code{\link{fixed_tuning}}: fixed \eqn{\eta}, no data-driven selection
#' }
#'
#' A bare numeric or function can also be passed as a tuning argument;
#' \code{as_tuning()} coerces it to a strategy object.

# ============================================================
# Tuning context
# ============================================================

#' Build a tuning context from precomputed grid sweep results
#'
#' Creates the context object passed to a tuning strategy's
#' \code{$select(ctx)} method. The context carries the regularization
#' grid and lazily-evaluated summaries from the solver sweep.
#'
#' @param eta_grid Numeric vector of candidate regularization values,
#'   ordered from smallest (least regularized) to largest.
#' @param n Integer; sample size.
#' @param ... Named closures or values to include in the context.
#'   Standard fields consumed by tuning strategies:
#'   \describe{
#'     \item{\code{dir_val}}{Closure returning a numeric vector: plug-in
#'       functional value at each grid point.}
#'     \item{\code{dir_se}}{Closure returning a numeric vector: standard
#'       error of the plug-in at each grid point.}
#'     \item{\code{loo_scores}}{Closure returning a numeric vector:
#'       leave-one-out cross-validated loss at each grid point.}
#'   }
#'
#' @return A list (tuning context) with components \code{eta_grid},
#'   \code{n}, and any additional named arguments from \code{...}.
#'
#' @examples
#' ctx <- tuning_context(
#'   eta_grid = 10^seq(-3, 0, length.out = 20),
#'   n = 200,
#'   dir_val = function() rnorm(20),
#'   dir_se  = function() rep(0.1, 20)
#' )
#'
#' @export
tuning_context = function(eta_grid, n, ...) {
  c(list(eta_grid = eta_grid, n = n), list(...))
}

# ============================================================
# Coerce number/function to tuning strategy
# ============================================================

#' Coerce a tuning argument to a strategy object
#'
#' Accepts a numeric value, a function, or a tuning strategy list, and
#' returns a valid tuning strategy. This allows callers to pass a bare
#' number (fixed \eqn{\eta}) or a function (\eqn{\eta = f(n)}) wherever
#' a tuning strategy is expected.
#'
#' @param x A tuning specification: a numeric scalar (coerced via
#'   \code{\link{fixed_tuning}}), a function of \code{n} (coerced via
#'   \code{rule_tuning}), or a list with a \code{$select} method
#'   (returned as-is).
#'
#' @return A tuning strategy list with at least a \code{$select(ctx)}
#'   method.
#'
#' @keywords internal
as_tuning = function(x) {
  if (is.list(x) && !is.null(x$select)) return(x)
  if (is.numeric(x)) return(fixed_tuning(x))
  if (is.function(x)) return(rule_tuning(x))
  stop("tuning must be a number, function, or tuning strategy (list with $select)")
}

#' Fixed regularization (no data-driven selection)
#'
#' Creates a tuning strategy that always selects the grid point closest
#' to a pre-specified \eqn{\eta}. Useful as a baseline or when
#' \eqn{\eta} has been determined externally.
#'
#' @param eta Numeric scalar; the fixed regularization parameter.
#'
#' @return A tuning strategy list with components:
#'   \describe{
#'     \item{\code{name}}{\code{"fixed"}.}
#'     \item{\code{eta}}{The specified \eqn{\eta}.}
#'     \item{\code{needs_grid}}{\code{FALSE}.}
#'     \item{\code{select(ctx)}}{Returns \code{list(eta, idx)} for the
#'       closest grid point.}
#'   }
#'
#' @examples
#' ft <- fixed_tuning(eta = 0.01)
#' ft$name
#' # "fixed"
#'
#' @seealso \code{\link{lepski_tuning}}, \code{\link{cv_tuning}},
#'   \code{\link{oracle_tuning}}
#' @export
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

#' Rule-based regularization
#'
#' Creates a tuning strategy where \eqn{\eta} is computed as a function
#' of the sample size: \eqn{\eta = f(n)}.
#'
#' @param fn A function taking \code{n} (sample size) and returning a
#'   numeric scalar \eqn{\eta}.
#'
#' @return A tuning strategy list with \code{$select(ctx)}.
#'
#' @keywords internal
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

#' Lepski tuning: functional-stability regularization selection
#'
#' Selects the largest admissible \eqn{\eta} where the plug-in functional
#' is stable. The criterion is Lepski's method applied to the regularization
#' path: \eqn{\eta_j} is admissible if for every smaller \eqn{\eta_i}
#' (\eqn{i < j}),
#' \deqn{|\hat\psi(\eta_i) - \hat\psi(\eta_j)| \le C \cdot
#'   \mathrm{SE}(\hat\psi(\eta_i))}
#' The selected \eqn{\eta} is the largest admissible value.
#'
#' @param C Numeric scalar; the tolerance multiplier for the stability
#'   criterion. Default \code{3}.
#'
#' @return A tuning strategy list with components:
#'   \describe{
#'     \item{\code{name}}{\code{"lepski"}.}
#'     \item{\code{C}}{The tolerance multiplier.}
#'     \item{\code{needs_grid}}{\code{TRUE}.}
#'     \item{\code{select(ctx)}}{Returns \code{list(eta, idx)}.}
#'   }
#'
#' @details
#' Works for functionals whose plug-in "stabilizes then drifts" (TSM, ATE,
#' survival probability, risk ratio). Fails for variance-type functionals
#' (VCATE) where the plug-in monotonically shrinks with smoothing.
#'
#' Requires the tuning context to provide \code{dir_val()} and
#' \code{dir_se()} closures.
#'
#' @examples
#' lt <- lepski_tuning(C = 3)
#' lt$name
#' # "lepski"
#'
#' @seealso \code{\link{cv_tuning}}, \code{\link{fixed_tuning}},
#'   \code{\link{oracle_tuning}}, \code{\link{tuning_context}}
#' @export
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

#' Oracle tuning: best regularization given known truth
#'
#' Selects the \eqn{\eta} that minimizes the absolute error of the
#' plug-in functional relative to the true parameter value. This is
#' only meaningful in simulation studies where the truth is known, and
#' serves as a benchmark for data-driven strategies.
#'
#' @param truth Numeric scalar; the true value of the estimand.
#'
#' @return A tuning strategy list with components:
#'   \describe{
#'     \item{\code{name}}{\code{"oracle"}.}
#'     \item{\code{truth}}{The true parameter value.}
#'     \item{\code{needs_grid}}{\code{TRUE}.}
#'     \item{\code{select(ctx)}}{Returns \code{list(eta, idx)} minimizing
#'       \eqn{|\hat\psi(\eta) - \psi_0|}.}
#'   }
#'
#' @examples
#' ot <- oracle_tuning(truth = 0.75)
#' ot$name
#' # "oracle"
#'
#' @seealso \code{\link{lepski_tuning}}, \code{\link{cv_tuning}},
#'   \code{\link{fixed_tuning}}
#' @export
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

#' Cross-validation tuning: LOO dual-loss selection
#'
#' Selects the \eqn{\eta} that minimizes the leave-one-out
#' cross-validated dual loss. The LOO scores are computed from the
#' solver's own first-order conditions, so this works for any model
#' (hazard or balancing weight) without requiring a held-out response.
#'
#' @return A tuning strategy list with components:
#'   \describe{
#'     \item{\code{name}}{\code{"cv"}.}
#'     \item{\code{needs_grid}}{\code{TRUE}.}
#'     \item{\code{select(ctx)}}{Returns \code{list(eta, idx, cv_loss)}
#'       where \code{cv_loss} is the full vector of LOO losses.}
#'   }
#'
#' @details
#' Requires the tuning context to provide a \code{loo_scores()} closure.
#'
#' @examples
#' ct <- cv_tuning()
#' ct$name
#' # "cv"
#'
#' @seealso \code{\link{lepski_tuning}}, \code{\link{fixed_tuning}},
#'   \code{\link{oracle_tuning}}, \code{\link{tuning_context}}
#' @export
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
