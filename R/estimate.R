# ============================================================
# Top-level estimator: dispatches on estimand class
# ============================================================

#' Estimate a causal estimand using kernel Bregman balancing
#'
#' Top-level entry point for the balancing package. Dispatches on the
#' class of \code{estimand}: survival estimands (class
#' \code{"survival_estimand"}) are handled by \code{\link{survival_effect}},
#' which runs the full crossfit pipeline (hazard model, balancing weights,
#' DR correction, bootstrap). Single-outcome estimands (class
#' \code{"effect_estimand"}) are not yet wired into the pipeline; use the
#' estimand's \code{$dr} method directly for now.
#'
#' @param estimand An estimand object created by one of the estimand
#'   constructors: \code{\link{treatment_specific_mean}},
#'   \code{\link{vcate}}, \code{\link{risk_ratio}},
#'   \code{\link{survival_probability_ate}},
#'   \code{\link{survival_probability_tsm}}, or \code{\link{rmst_ate}}.
#' @param ... Arguments passed to the dispatched method. For survival
#'   estimands these are forwarded to \code{\link{survival_effect}} and
#'   include:
#'   \describe{
#'     \item{\code{observations}}{A list with components \code{U} (event/censoring
#'       time), \code{Delta} (event indicator), \code{A} (treatment arm),
#'       \code{X} (covariate matrix).}
#'     \item{\code{kern}}{A kernel object, e.g. from \code{\link{matern_kernel}}
#'       or \code{\link{direct_product_kernel}}.}
#'     \item{\code{lambda_disp}}{A dispersion object for the hazard model,
#'       e.g. \code{\link{entropy_dispersion}}.}
#'     \item{\code{eta_lam}}{Tuning strategy or numeric regularization
#'       parameter for the hazard model.}
#'     \item{\code{eta_gam}}{Tuning strategy or numeric regularization
#'       parameter for the balancing weights.}
#'     \item{\code{horizon}}{Numeric; the time horizon \eqn{\bar t} for
#'       survival or RMST estimands.}
#'     \item{\code{discrete}}{Logical; if \code{TRUE}, use a discrete-time
#'       grid. Default \code{FALSE}.}
#'   }
#' @return For survival estimands, the return value of
#'   \code{\link{survival_effect}}: an \code{effect_estimate} object
#'   containing the DR point estimate, estimated standard error,
#'   per-subject influence function values, and fitted nuisance models
#'   (for bootstrap). For single-outcome estimands, not yet implemented.
#'
#' @examples
#' \dontrun{
#' # --- Survival estimand (ATE of survival probability) ---
#' est <- balancing_estimate(
#'   estimand    = survival_probability_ate(),
#'   observations = list(U = U, Delta = Delta, A = A, X = X),
#'   kern        = matern_kernel(nu = 2.5),
#'   lambda_disp = entropy_dispersion(),
#'   eta_lam     = lepski_tuning(),
#'   eta_gam     = cv_tuning(),
#'   horizon     = 5
#' )
#'
#' # --- Single-outcome estimand (use $dr directly) ---
#' tsm <- treatment_specific_mean(arm = 1)
#' result <- tsm$dr(mu1, mu0, Y, W, gamma1, gamma0)
#' }
#'
#' @seealso \code{\link{survival_effect}} for the survival pipeline,
#'   \code{\link{treatment_specific_mean}}, \code{\link{survival_probability_ate}}
#'   for estimand constructors, \code{\link{lepski_tuning}} and
#'   \code{\link{cv_tuning}} for tuning strategies.
#' @export
balancing_estimate = function(estimand, ...) UseMethod("balancing_estimate")

#' @describeIn balancing_estimate Method for survival estimands; delegates
#'   to \code{\link{survival_effect}}.
#' @export
balancing_estimate.survival_estimand = function(estimand, ...) {
  survival_effect(estimand = estimand, ...)
}

#' @describeIn balancing_estimate Method for single-outcome effect estimands
#'   (not yet implemented).
#' @export
balancing_estimate.effect_estimand = function(estimand, ...) {
  stop("single-outcome pipeline not yet implemented; use the estimand's $dr method directly")
}
