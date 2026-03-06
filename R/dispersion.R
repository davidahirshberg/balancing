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

#' Entropy dispersion
#'
#' Construct a dispersion based on the entropy (exponential) conjugate.
#'
#' @details
#' A dispersion encodes the penalty structure for the balancing weights
#' problem.
#' The solver works with the conjugate \eqn{\chi^*} (the Legendre
#' transform of the primal penalty \eqn{\chi}), so each dispersion is
#' defined by \eqn{\chi^*} and its derivatives.
#'
#' The entropy dispersion corresponds to the primal
#' \deqn{\chi(\gamma) = \gamma \log \gamma - \gamma}
#' on \eqn{\gamma > 0}, with conjugate
#' \deqn{\chi^*(\phi) = \exp(\phi).}
#' All three derivatives coincide:
#' \eqn{\dot\chi^*(\phi) = \ddot\chi^*(\phi) = \exp(\phi)}.
#' The primal weights are therefore
#' \eqn{\hat\gamma = \exp(\hat\phi)}, which are strictly positive.
#'
#' The associated Bregman divergence is the KL divergence (Poisson
#' deviance).  This is the default choice when weights must be
#' nonnegative and the scale of \eqn{\gamma} is unconstrained.
#'
#' @return A \code{dispersion} S3 object (a list with slots \code{base},
#'   \code{offset}, \code{sigma}, \code{v}) whose S3 methods provide
#'   \eqn{\chi^*_Z}, \eqn{\dot\chi^*_Z}, and \eqn{\ddot\chi^*_Z} for
#'   the dual solver.
#'
#' @examples
#' d <- entropy_dispersion()
#'
#' @export
entropy_dispersion   = function() dispersion(entropy_base)

#' Softplus dispersion
#'
#' Construct a dispersion based on the softplus (logistic) conjugate.
#'
#' @details
#' The softplus dispersion corresponds to the primal
#' \deqn{\chi(\gamma) = \gamma \log \gamma + (1 - \gamma) \log(1 - \gamma)}
#' on \eqn{\gamma \in (0, 1)}, with conjugate
#' \deqn{\chi^*(\phi) = \log(1 + \exp(\phi)).}
#' The first derivative is the logistic (expit) function
#' \eqn{\dot\chi^*(\phi) = \sigma(\phi) = 1/(1 + \exp(-\phi))}, and the
#' second derivative is
#' \eqn{\ddot\chi^*(\phi) = \sigma(\phi)(1 - \sigma(\phi))}.
#'
#' The primal weights are therefore
#' \eqn{\hat\gamma = \sigma(\hat\phi) \in (0, 1)}, making this
#' dispersion appropriate when weights must lie in the unit interval.
#' The associated Bregman divergence is the binary KL divergence.
#'
#' @return A \code{dispersion} S3 object.  See \code{\link{entropy_dispersion}}
#'   for the structure.
#'
#' @examples
#' d <- softplus_dispersion()
#'
#' @export
softplus_dispersion  = function() dispersion(softplus_base)

#' Quadratic dispersion
#'
#' Construct a dispersion based on the quadratic conjugate.
#'
#' @details
#' The quadratic dispersion corresponds to the primal
#' \deqn{\chi(\gamma) = \gamma^2 / 2}
#' on \eqn{\gamma \in \mathbf{R}}, with conjugate
#' \deqn{\chi^*(\phi) = \phi^2 / 2.}
#' The first derivative is \eqn{\dot\chi^*(\phi) = \phi} and the
#' second derivative is \eqn{\ddot\chi^*(\phi) = 1} (constant).
#'
#' Because the conjugate is itself quadratic, Newton's method converges
#' in a single step.  The primal weights are
#' \eqn{\hat\gamma = \hat\phi}, which are unconstrained in sign and
#' magnitude.  This is the natural choice for ATE-type estimands where
#' negative weights are acceptable.
#'
#' @return A \code{dispersion} S3 object.  See \code{\link{entropy_dispersion}}
#'   for the structure.
#'
#' @examples
#' d <- quadratic_dispersion()
#'
#' @export
quadratic_dispersion = function() dispersion(quadratic_base)

shifted_entropy_dispersion = function() dispersion(shifted_entropy_base)

#' Sign-flipped dispersion
#'
#' Construct a dispersion whose primal penalty is sign-flipped by a
#' treatment indicator, giving opposite-sign weights for treated and
#' control units.
#'
#' @param base A base dispersion record (e.g. \code{entropy_base},
#'   \code{quadratic_base}).  Each base is a list with elements
#'   \code{dchis}, \code{ddchis}, and \code{chis} defining
#'   \eqn{\chi^*} and its first two derivatives on the canonical
#'   (unshifted, unflipped) domain.
#' @param W A numeric vector of signs, typically a treatment indicator
#'   taking values in \eqn{\{-1, +1\}}.  Treated units (\eqn{W = +1})
#'   receive positive weights; control units (\eqn{W = -1}) receive
#'   negative weights.
#'
#' @details
#' For ATE-type estimands the primal penalty is
#' \deqn{\chi_Z(\gamma) = \chi(W \gamma)}
#' where \eqn{W \in \{-1, +1\}} is the treatment indicator.
#' The sign-flip enters the conjugate via
#' \deqn{\dot\chi^*_Z(\phi) = W \, \dot\chi^*(W \phi),}
#' so the dual-to-primal map
#' \eqn{\hat\gamma = \dot\chi^*_Z(\hat\phi)} automatically produces
#' weights with the correct sign for each arm.
#'
#' Internally this sets the \eqn{\sigma} slot of the dispersion to
#' \code{W}, which maps \eqn{\phi \mapsto W\phi} before applying the
#' base \eqn{\dot\chi^*}.
#'
#' @return A \code{dispersion} S3 object with \code{sigma = W}.
#'   See \code{\link{entropy_dispersion}} for the structure.
#'
#' @examples
#' # ATE with entropy penalty: positive weights for treated,
#' # negative for control
#' W <- c(1, 1, -1, -1, 1)
#' d <- sign_flip_dispersion(entropy_base, W)
#'
#' @export
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
