# ============================================================
# Top-level estimator: dispatches on estimand class
# ============================================================

#' Estimate a causal estimand using kernel Bregman balancing.
#'
#' @param estimand  An estimand object (from estimands.R constructors).
#' @param ...       Arguments passed to the method (observations, kern, etc.).
#' @return Depends on the estimand type.
balancing_estimate = function(estimand, ...) UseMethod("balancing_estimate")

#' Survival estimands: delegates to survival_effect.
balancing_estimate.survival_estimand = function(estimand, ...) {
  survival_effect(estimand = estimand, ...)
}

#' Effect estimands: not yet implemented.
balancing_estimate.effect_estimand = function(estimand, ...) {
  stop("single-outcome pipeline not yet implemented; use the estimand's $dr method directly")
}
