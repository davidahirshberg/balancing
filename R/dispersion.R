#' # Dispersions
#'
#' S3 class with composable base records for the Bregman
#' divergence penalty.
#'
#' ## Z-Dependent Dispersion
#'
#' **Paper**: the primal dispersion \eqn{\chi_Z} shifts and
#' sign-flips a base \eqn{\chi}:
#' \deqn{\chi_Z(\gamma) = v\,\chi\bigl(\sigma(\gamma - o)\bigr)}
#' where \eqn{o} is the offset (balancing target), \eqn{\sigma =
#' \mathrm{sign}(o)} is the sign-flip, and \eqn{v \ge 0} is the
#' penalty scale.
#'
#' **Code** (conjugate, used by the dual solver):
#' with \eqn{\psi = \sigma\,\phi / v},
#' \deqn{\chi^*_Z(\phi) = o\,\phi + v\,\chi^*(\psi)}
#' \deqn{\dot\chi^*_Z(\phi) = o + \sigma\,\dot\chi^*(\psi)}
#' \deqn{\ddot\chi^*_Z(\phi) = \ddot\chi^*(\psi) / v}
#' The \eqn{v} cancels in \eqn{\dot\chi^*_Z} because \eqn{\partial\psi /
#' \partial\phi = \sigma/v} and \eqn{\dot\chi^*(\psi)\cdot\sigma/v
#' \cdot v = \sigma\,\dot\chi^*(\psi)}.
#'
#' ## Dispersion Table
#'
#' | Name | \eqn{\chi^*(\phi)} | \eqn{\dot\chi^*(\phi)} | \eqn{\ddot\chi^*(\phi)} |
#' |------|----------------|--------------------|--------------------|
#' | Quadratic | \eqn{\phi^2/2} | \eqn{\phi} | \eqn{1} |
#' | Entropy | \eqn{e^\phi} | \eqn{e^\phi} | \eqn{e^\phi} |
#' | Softplus | \eqn{\log(1+e^\phi)} | \eqn{\sigma(\phi)} | \eqn{\sigma(\phi)(1-\sigma(\phi))} |
#' | Shifted entropy | \eqn{\phi + e^\phi} | \eqn{1 + e^\phi} | \eqn{e^\phi} |
#' | Pos. quadratic | \eqn{(\phi_+)^2/2} | \eqn{\phi_+} | \eqn{\mathbf{1}(\phi>0)} |

# ---- S3 generics (dot-prefix to avoid clash with local lambdas) ----

.dchis  = function(disp, phi) UseMethod(".dchis")
.ddchis = function(disp, phi) UseMethod(".ddchis")
.chis   = function(disp, phi) UseMethod(".chis")

# ---- Base dispersion records ----
#' Each base record defines chi* and its first two derivatives on the
#' "canonical" domain (no shift, no sign-flip, no scale).
#' A custom base only needs these three functions.

#' Entropy: \eqn{\chi(g) = g \log g - g} on \eqn{R_+}.
#' \eqn{\chi^*(\phi) = \exp(\phi)}. Bregman divergence is KL / Poisson.
entropy_base = list(
  dchis  = function(phi) exp(phi),
  ddchis = function(phi) exp(phi),
  chis   = function(phi) exp(phi)
)

#' Quadratic: \eqn{\chi(g) = g^2/2} on \eqn{R}.
#' \eqn{\chi^*(\phi) = \phi^2/2}. Linear: one Newton step exact.
quadratic_base = list(
  dchis    = function(phi) phi,
  ddchis   = function(phi) rep(1, length(phi)),
  chis     = function(phi) phi^2 / 2,
  quadratic = TRUE
)

#' Softplus: \eqn{\chi(g) = g \log g + (1-g) \log(1-g)} on \eqn{(0,1)}.
#' \eqn{\chi^*(\phi) = \log(1 + \exp(\phi))}. Maps dual -> primal via expit.
softplus_base = list(
  dchis  = function(phi) plogis(phi),
  ddchis = function(phi) { p = plogis(phi); p * (1 - p) },
  chis   = function(phi) log(1 + exp(phi))
)

#' Shifted entropy: \eqn{\chi(g) = (g-1) \log(g-1) - (g-1) + 1} on \eqn{[1, \infty)}.
#' \eqn{\chi^*(\phi) = \phi + \exp(\phi)}. \eqn{\dot\chi^* = 1 + \exp(\phi)}.
shifted_entropy_base = list(
  dchis  = function(phi) 1 + exp(phi),
  ddchis = function(phi) exp(phi),
  chis   = function(phi) phi + exp(phi)
)

#' Positivity-constrained quadratic: \eqn{\chi(g) = (g_+)^2/2} on \eqn{R_+}.
#' \eqn{\chi^*(\phi) = (\phi_+)^2/2}. Active set changes with sign of \eqn{\phi}.
pos_quadratic_base = list(
  dchis  = function(phi) pmax(phi, 0),
  ddchis = function(phi) as.numeric(phi > 0),
  chis   = function(phi) pmax(phi, 0)^2 / 2
)

#' ## Constructor
#'
#' `dispersion(base, offset, sigma, v)` builds a Z-dependent
#' dispersion from a base record. The S3 methods `.dchis`,
#' `.ddchis`, `.chis` apply the composition formulas above.
dispersion = function(base, offset = 0, sigma = 1, v = 1)
  structure(list(base = base, offset = offset, sigma = sigma, v = v),
            class = "dispersion")

# ---- S3 methods ----

#' \eqn{\dot\chi^*_Z(\phi) = o + \sigma \cdot \dot\chi^*(\psi)}, \eqn{\psi = \sigma \cdot \phi / v}
.dchis.dispersion = function(disp, phi) {
  psi = disp$sigma * phi / disp$v
  disp$offset + disp$sigma * disp$base$dchis(psi)
}

#' \eqn{\ddot\chi^*_Z(\phi) = \ddot\chi^*(\psi) / v}
.ddchis.dispersion = function(disp, phi) {
  psi = disp$sigma * phi / disp$v
  disp$base$ddchis(psi) / disp$v
}

#' \eqn{\chi^*_Z(\phi) = o \cdot \phi + v \cdot \chi^*(\psi)}
.chis.dispersion = function(disp, phi) {
  psi = disp$sigma * phi / disp$v
  disp$offset * phi + disp$v * disp$base$chis(psi)
}

#' ## Convenience Constructors
#'
#' Wrappers that build common dispersions with readable names.

entropy_dispersion   = function() dispersion(entropy_base)
#' @rdname entropy_dispersion
softplus_dispersion  = function() dispersion(softplus_base)
#' @rdname entropy_dispersion
quadratic_dispersion = function() dispersion(quadratic_base)
shifted_entropy_dispersion = function() dispersion(shifted_entropy_base)

#' ### `sign_flip_dispersion`
#'
#' **Paper**: for ATE weights, \eqn{\chi_Z(\gamma) = \chi(W\gamma)}
#' with \eqn{W \in \{-1, +1\}}. Treated units (\eqn{W = +1}) get positive
#' weights; control units (\eqn{W = -1}) get negative weights.
#'
#' **Code**: `dispersion(base, sigma = W)` — the sign-flip goes
#' into the \eqn{\sigma} slot, which maps \eqn{\phi \mapsto W\phi}
#' before applying the base \eqn{\dot\chi^*}.
sign_flip_dispersion = function(base, W) dispersion(base, sigma = W)

#' ### `target_scaled_entropy_dispersion`
#'
#' **Paper**: \eqn{\chi^*_Z(\phi) = r\,\phi + e^{\mathrm{sign}(r)\,\phi}}
#' giving \eqn{\hat\gamma = r + \mathrm{sign}(r)\,e^{\mathrm{sign}(r)\,\phi}}.
#' Floors \eqn{|\hat\gamma|} at \eqn{|r|} (the Riesz representer
#' approximation), then exponential on the excess.
#'
#' **Code**: `dispersion(entropy_base, offset = r, sigma = sign(r))`.
#' This is the default weight dispersion in the survival pipeline.
target_scaled_entropy_dispersion = function(r)
  dispersion(entropy_base, offset = r, sigma = sign(r))

#' ### `variance_weighted_quadratic_dispersion`
#'
#' **Paper**: \eqn{\chi_Z(\gamma) = \frac{v}{2}(\gamma - r)^2} on
#' \eqn{\mathrm{sign}(\gamma) = \mathrm{sign}(r)}, \eqn{|\gamma| \ge |r|}.
#' The scale \eqn{v} allows heteroskedastic weighting.
#'
#' **Code**: `dispersion(pos_quadratic_base, offset = r,
#' sigma = sign(r), v = v)`.
variance_weighted_quadratic_dispersion = function(r, v)
  dispersion(pos_quadratic_base, offset = r, sigma = sign(r), v = v)
