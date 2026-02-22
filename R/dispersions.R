## Dispersions: the pluggable chi / chi* choice.
##
## Each dispersion carries:
##   chi(g)       — primal penalty (for reference/documentation)
##   chis(phi)    — conjugate chi*
##   dchis(phi)   — dot chi* (maps dual -> primal: gamma = dchis(phi))
##   ddchis(phi)  — ddot chi* (curvature, IRLS weight)
##   resid(g, Y, w) — weighted gradient of loss in primal terms
##   linear       — if TRUE, solver uses direct solve instead of IRLS
##
## Names match the paper macros: \chis, \dchis, \ddchis.

# ============================================================
# Base dispersions
# ============================================================

#' Entropy dispersion: chi(g) = g log g - g on R+.
#' Conjugate chi*(phi) = exp(phi). Maps dual -> primal via exp.
#' Natural for: ATT weights (positive), continuous hazard (Poisson loss).
#' The current code's exp_link() is this dispersion.
entropy_dispersion = function() {
  list(
    name = "entropy",
    linear = FALSE,
    chi = function(g) ifelse(g > 0, g * log(g) - g, Inf),
    chis = function(phi) exp(phi),
    dchis = function(phi) exp(phi),
    ddchis = function(phi) exp(phi),
    # Loss gradient: w * (g - Y) where g = dchis(phi) = exp(phi)
    resid = function(g, Y, w) w * (g - Y),
    # IRLS RHS for direct solver: w*g*(phi - 1) + w*Y
    rhs = function(phi, g, Y, w) w * g * (phi - 1) + w * Y,
    predict = function(phi) exp(phi)
  )
}

#' Softplus dispersion: chi(g) = g log g + (1-g) log(1-g) (binary entropy) on (0,1).
#' Conjugate chi*(phi) = log(1 + exp(phi)). Maps dual -> primal via expit.
#' Natural for: discrete hazard probabilities, binary outcomes.
#' The current code's logistic_link() is this dispersion.
softplus_dispersion = function() {
  list(
    name = "softplus",
    linear = FALSE,
    chi = function(g) {
      g = pmin(pmax(g, 1e-15), 1 - 1e-15)
      g * log(g) + (1 - g) * log(1 - g)
    },
    chis = function(phi) log(1 + exp(phi)),
    dchis = function(phi) plogis(phi),
    ddchis = function(phi) { p = plogis(phi); p * (1 - p) },
    # Loss gradient: w * (g - Y) where g = expit(phi)
    resid = function(g, Y, w) w * (g - Y),
    # IRLS RHS
    rhs = function(phi, g, Y, w) {
      wc = w * g * (1 - g)
      wc * phi + w * (Y - g)
    },
    predict = function(phi) plogis(phi)
  )
}

#' Quadratic dispersion: chi(g) = g^2/2 on R.
#' Conjugate chi*(phi) = phi^2/2. Linear: dchis(phi) = phi, ddchis = 1.
#' The solver reduces to a single linear system (no IRLS).
#' Natural for: ridge regression, linear balancing weights.
quadratic_dispersion = function() {
  list(
    name = "quadratic",
    linear = TRUE,
    chi = function(g) g^2 / 2,
    chis = function(phi) phi^2 / 2,
    dchis = function(phi) phi,
    ddchis = function(phi) rep(1, length(phi)),
    resid = function(g, Y, w) w * (g - Y),
    rhs = function(phi, g, Y, w) w * Y,
    predict = function(phi) phi
  )
}

#' Positivity-constrained quadratic: chi(g) = (g_+)^2/2 on R+.
#' Conjugate chi*(phi) = (phi_+)^2/2. dchis(phi) = phi_+, ddchis = 1{phi>0}.
#' Like quadratic but forces gamma >= 0. Not linear (active set changes).
pos_quadratic_dispersion = function() {
  list(
    name = "pos_quadratic",
    linear = FALSE,
    chi = function(g) ifelse(g > 0, g^2 / 2, Inf),
    chis = function(phi) pmax(phi, 0)^2 / 2,
    dchis = function(phi) pmax(phi, 0),
    ddchis = function(phi) as.numeric(phi > 0),
    resid = function(g, Y, w) w * (g - Y),
    rhs = function(phi, g, Y, w) {
      active = phi > 0
      ifelse(active, w * phi + w * (Y - g), 0)
    },
    predict = function(phi) pmax(phi, 0)
  )
}

# ============================================================
# Sign-flip modifier
# ============================================================

#' Sign-flip a base dispersion: chi_Z(gamma) = chi(W * gamma).
#' W in {-1, +1} forces gamma > 0 for treated, gamma < 0 for control.
#'
#' Conjugate: chis_Z(phi) = chis(W * phi).
#' Derivatives: dchis_Z(phi) = W * dchis(W * phi), ddchis_Z(phi) = ddchis(W * phi).
#' The second equality uses W^2 = 1: curvature is preserved.
#'
#' This is the paper's chi_Z(gamma) = chi(W gamma) construction (spinoff1 eq:sign-flip-chi).
#' Applied to any R+-domain base dispersion to get sign-correct ATE weights.
signflip = function(base, W) {
  list(
    name = paste0("signflip_", base$name),
    linear = base$linear,
    chi = function(g) base$chi(W * g),
    chis = function(phi) base$chis(W * phi),
    dchis = function(phi) W * base$dchis(W * phi),
    ddchis = function(phi) base$ddchis(W * phi),
    resid = function(g, Y, w) {
      # g = W * dchis_base(W*phi), Y is the target
      # gradient in phi: w * (chis_Z'(phi) - Y) but we work through the base
      w * (g - Y)
    },
    rhs = function(phi, g, Y, w) {
      # IRLS working response: delegate to base with W*phi
      base$rhs(W * phi, base$dchis(W * phi), W * Y, w) * W
    },
    predict = function(phi) W * base$dchis(W * phi),
    base = base,
    W = W
  )
}
