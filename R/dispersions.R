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
# Shifted dispersions (domain [1, infty))
# ============================================================

#' Shifted entropy: chi(gamma) = (gamma-1) log(gamma-1) - (gamma-1) + 1 on [1, infty).
#' Conjugate chi*(phi) = phi + exp(phi). dchis = 1 + exp(phi), ddchis = exp(phi).
#' Natural for: ATE/TSM inverse propensity weights (gamma = 1/pi >= 1).
#' With sign-flip W = 2A-1: weights are 1+exp(phi) for treated, -(1+exp(-phi)) for control.
shifted_entropy_dispersion = function() {
  list(
    name = "shifted_entropy",
    linear = FALSE,
    chi = function(g) ifelse(g >= 1, (g - 1) * log(g - 1) - (g - 1) + 1, Inf),
    chis = function(phi) phi + exp(phi),
    dchis = function(phi) 1 + exp(phi),
    ddchis = function(phi) exp(phi),
    resid = function(g, Y, w) w * (g - Y),
    rhs = function(phi, g, Y, w) {
      # g = 1 + exp(phi), ddchis = exp(phi)
      ep = exp(phi)
      w * ep * (phi - 1) + w * (Y - 1)
    },
    predict = function(phi) 1 + exp(phi)
  )
}

#' Shifted positive-part quadratic: chi(gamma) = (gamma-1)^2/2 on [1, infty).
#' Conjugate chi*(phi) = phi + max(phi, 0)^2/2. dchis = 1 + max(phi, 0), ddchis = I(phi > 0).
#' Like shifted entropy but with quadratic tail behavior.
shifted_pos_quadratic_dispersion = function() {
  list(
    name = "shifted_pos_quadratic",
    linear = FALSE,
    chi = function(g) ifelse(g >= 1, (g - 1)^2 / 2, Inf),
    chis = function(phi) phi + pmax(phi, 0)^2 / 2,
    dchis = function(phi) 1 + pmax(phi, 0),
    ddchis = function(phi) as.numeric(phi > 0),
    resid = function(g, Y, w) w * (g - Y),
    rhs = function(phi, g, Y, w) {
      active = phi > 0
      ifelse(active, w * phi + w * (Y - 1), 0)
    },
    predict = function(phi) 1 + pmax(phi, 0)
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

# ============================================================
# Target-sign dispersions
# ============================================================
#
# Standard sign-flip uses sigma = W = 2A-1 (sign from treatment).
# These use sigma = sign(r) (sign from the target derivative).
#
# Why: entropy + sign-flip + per-level null space is infeasible when the
# per-arm target has the wrong sign for the dispersion's range.  With
# sigma = sign(r), the weight sign matches the target sign by construction,
# so per-level null-space constraints are always feasible.
#
# The ATE subtraction must live in r (negate arm 0's derivative), not in W.

#' Sign-from-target entropy: sigma_Z = sign(r_Z), shift s = 0.
#'
#' chi*_Z(phi) = exp(sigma phi) - 1
#' gamma  = sigma exp(sigma phi):  sign(gamma) = sign(r),  |gamma| > 0
#'
#' Equivalent to signflip(entropy, sign(r)).  Per-level feasibility
#' when sign(r) is constant within each level (always true for survival ATE).
#'
#' @param r  Per-observation target derivative (n-vector, nonzero).
#'   For ATE survival: r < 0 at arm-1 points, r > 0 at arm-0 points
#'   (ATE subtraction already in r).
signflip_r = function(r) {
  signflip(entropy_dispersion(), sign(r))
}

#' Target-scaled sign-flip entropy: sigma_Z = sign(r_Z), shift s_Z = |r_Z|.
#'
#' chi*_Z(phi)  = r phi + exp(sigma phi) - 1
#' gamma        = r + sigma exp(sigma phi)
#' |gamma|     >= |r|,  sign(gamma) = sign(r)
#'
#' The shift floors |gamma| at |r|: each weight is at least as large as the
#' target it needs to recover.  The excess |gamma| - |r| = exp(sigma phi)
#' is the adaptive entropy part.
#'
#' Per-level feasibility: guaranteed when sign(r) is constant within each
#' level, because each arm's null-space target c_a = sum(r_i) has the same
#' sign as the weights.
#'
#' @param r  Per-observation target derivative (n-vector, nonzero).
target_scaled_entropy = function(r) {
  sigma = sign(r)
  abs_r = abs(r)
  list(
    name = "target_scaled_entropy",
    linear = FALSE,
    chi = function(g) {
      # chi(gamma) = chi_{|r|}(sigma * gamma)  where chi_s is shifted entropy
      # Domain: sigma*gamma > |r|, i.e., sign(gamma) = sigma and |gamma| > |r|
      u = sigma * g - abs_r
      ifelse(u > 0, u * log(u) - u + 1, Inf)
    },
    # chi*(phi) = r*phi + exp(sigma*phi) - 1
    chis = function(phi) r * phi + exp(sigma * phi) - 1,
    # gamma = r + sigma * exp(sigma*phi)
    dchis = function(phi) r + sigma * exp(sigma * phi),
    # curvature = exp(sigma*phi), same as sign-flipped entropy
    ddchis = function(phi) exp(sigma * phi),
    resid = function(g, Y, w) w * (g - Y),
    rhs = function(phi, g, Y, w) {
      # IRLS: ddchis * phi + (Y - g)
      ep = exp(sigma * phi)
      w * ep * phi + w * (Y - g)
    },
    predict = function(phi) r + sigma * exp(sigma * phi),
    sigma = sigma,
    r = r
  )
}

#' Variance-weighted quadratic with sign-flip and target shift.
#'
#' chi_Z(gamma) = (v_Z / 2)(gamma - sigma_Z |r_Z|)^2
#'   on sign(gamma) = sigma_Z, |gamma| >= |r_Z|
#'
#' where sigma_Z = sign(r_Z), v_Z = lambda_hat(Z)(1 - lambda_hat(Z)) 1(T >= U)
#' is the conditional Bernoulli variance of Y_U.
#'
#' Conjugate:
#'   chi*_Z(phi) = r_Z phi + (sigma_Z phi)_+^2 / (2 v_Z)
#'   gamma       = r_Z + sigma_Z (sigma_Z phi)_+ / v_Z
#'   ddchis      = 1{sigma_Z phi > 0} / v_Z
#'
#' The primal penalty P_hat[chi_Z(gamma)] = P_hat[v_Z (gamma - r)^2 / 2]
#' equals the Bernoulli-variance term in the estimator's MSE, so eta
#' directly controls the bias-variance exchange rate in the right units.
#'
#' Excess weight |gamma| - |r| = (sigma phi)_+ / v_Z: large when the RKHS
#' says to upweight and the conditional variance is small (precise obs).
#'
#' @param r  Per-observation target derivative (n-vector, nonzero).
#' @param v  Per-observation conditional variance (n-vector, positive).
#'   Typically v_Z = lambda_hat(1 - lambda_hat) * at_risk.
variance_weighted_quadratic = function(r, v) {
  sigma = sign(r)
  abs_r = abs(r)
  list(
    name = "variance_weighted_quadratic",
    linear = FALSE,
    chi = function(g) {
      u = sigma * g - abs_r
      ifelse(u >= 0, v / 2 * u^2, Inf)
    },
    # chi*(phi) = r phi + (sigma phi)_+^2 / (2v)
    chis = function(phi) {
      sp = pmax(sigma * phi, 0)
      r * phi + sp^2 / (2 * v)
    },
    # gamma = r + sigma (sigma phi)_+ / v
    dchis = function(phi) r + sigma * pmax(sigma * phi, 0) / v,
    # curvature = 1{sigma phi > 0} / v
    ddchis = function(phi) as.numeric(sigma * phi > 0) / v,
    resid = function(g, Y, w) w * (g - Y),
    rhs = function(phi, g, Y, w) {
      # On active set (sigma phi > 0): ddchis*phi + (Y - g)
      #   = phi/v + Y - r - phi/v = Y - r.  Clean cancellation.
      active = sigma * phi > 0
      ifelse(active, w * (Y - r), 0)
    },
    predict = function(phi) r + sigma * pmax(sigma * phi, 0) / v,
    sigma = sigma,
    r = r,
    v = v
  )
}
