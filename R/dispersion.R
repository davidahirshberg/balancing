# ============================================================
# Dispersions: S3 class with composable base records
# ============================================================
#
#' **Primal.** The Z-dependent dispersion shifts and sign-flips a base chi:
#'   chi_Z(gamma) = v * chi(sigma * (gamma - o))
#' where chi is the base, o is the offset (balancing target),
#' sigma = sign(o) is the sign-flip, v >= 0 is the penalty scale.
#'
#' **Conjugate.** With psi = sigma * phi / v:
#'   chi*_Z(phi)      = o * phi + v * chi*(psi)
#'   dot chi*_Z(phi)   = o + sigma * dot chi*(psi)        [v cancels in chain rule]
#'   ddot chi*_Z(phi)  = ddot chi*(psi) / v
#'
#' Base records carry the irreducible math. The S3 class handles composition.

# ---- S3 generics (dot-prefix to avoid clash with local lambdas) ----

.dchis  = function(disp, phi) UseMethod(".dchis")
.ddchis = function(disp, phi) UseMethod(".ddchis")
.chis   = function(disp, phi) UseMethod(".chis")

# ---- Base dispersion records ----
#' Each base record defines chi* and its first two derivatives on the
#' "canonical" domain (no shift, no sign-flip, no scale).
#' A custom base only needs these three functions.

#' Entropy: chi(g) = g log g - g on R+.
#' chi*(phi) = exp(phi). Bregman divergence is KL / Poisson.
entropy_base = list(
  dchis  = function(phi) exp(phi),
  ddchis = function(phi) exp(phi),
  chis   = function(phi) exp(phi)
)

#' Quadratic: chi(g) = g^2/2 on R.
#' chi*(phi) = phi^2/2. Linear: one Newton step exact.
quadratic_base = list(
  dchis    = function(phi) phi,
  ddchis   = function(phi) rep(1, length(phi)),
  chis     = function(phi) phi^2 / 2,
  quadratic = TRUE
)

#' Softplus: chi(g) = g log g + (1-g) log(1-g) on (0,1).
#' chi*(phi) = log(1 + exp(phi)). Maps dual -> primal via expit.
softplus_base = list(
  dchis  = function(phi) plogis(phi),
  ddchis = function(phi) { p = plogis(phi); p * (1 - p) },
  chis   = function(phi) log(1 + exp(phi))
)

#' Shifted entropy: chi(g) = (g-1) log(g-1) - (g-1) + 1 on [1, inf).
#' chi*(phi) = phi + exp(phi). dchis = 1 + exp(phi).
shifted_entropy_base = list(
  dchis  = function(phi) 1 + exp(phi),
  ddchis = function(phi) exp(phi),
  chis   = function(phi) phi + exp(phi)
)

#' Positivity-constrained quadratic: chi(g) = (g_+)^2/2 on R+.
#' chi*(phi) = (phi_+)^2/2. Active set changes with sign of phi.
pos_quadratic_base = list(
  dchis  = function(phi) pmax(phi, 0),
  ddchis = function(phi) as.numeric(phi > 0),
  chis   = function(phi) pmax(phi, 0)^2 / 2
)

# ---- Constructor ----

#' Build a dispersion from a base record with optional composition.
#'
#' @param base  Named list with $dchis, $ddchis, $chis.
#' @param offset  Per-observation offset o (balancing target). Default 0.
#' @param sigma  Per-observation sign-flip in {-1, +1}. Default 1 (no flip).
#' @param v  Per-observation penalty scale (v >= 0). Default 1.
dispersion = function(base, offset = 0, sigma = 1, v = 1)
  structure(list(base = base, offset = offset, sigma = sigma, v = v),
            class = "dispersion")

# ---- S3 methods ----

#' dot chi*_Z(phi) = o + sigma * dot chi*(psi),  psi = sigma * phi / v
.dchis.dispersion = function(disp, phi) {
  psi = disp$sigma * phi / disp$v
  disp$offset + disp$sigma * disp$base$dchis(psi)
}

#' ddot chi*_Z(phi) = ddot chi*(psi) / v
.ddchis.dispersion = function(disp, phi) {
  psi = disp$sigma * phi / disp$v
  disp$base$ddchis(psi) / disp$v
}

#' chi*_Z(phi) = o * phi + v * chi*(psi)
.chis.dispersion = function(disp, phi) {
  psi = disp$sigma * phi / disp$v
  disp$offset * phi + disp$v * disp$base$chis(psi)
}

# ---- Convenience aliases ----
#' These preserve existing call signatures. Each returns a dispersion S3 object.

entropy_dispersion   = function() dispersion(entropy_base)
softplus_dispersion  = function() dispersion(softplus_base)
quadratic_dispersion = function() dispersion(quadratic_base)
shifted_entropy_dispersion = function() dispersion(shifted_entropy_base)

#' Sign-flip: chi_Z(gamma) = chi(W * gamma). W in {-1, +1}.
signflip = function(base, W) dispersion(base, sigma = W)

#' Target-scaled entropy: gamma = r + sigma * exp(sigma * phi), sigma = sign(r).
#' Floors |gamma| at |r|, then adaptive entropy on the excess.
target_scaled_entropy = function(r)
  dispersion(entropy_base, offset = r, sigma = sign(r))

#' Variance-weighted quadratic with sign-flip and target shift.
#' chi_Z(gamma) = (v/2)(gamma - r)^2 on sign(gamma) = sign(r), |gamma| >= |r|.
variance_weighted_quadratic = function(r, v)
  dispersion(pos_quadratic_base, offset = r, sigma = sign(r), v = v)
