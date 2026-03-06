#' # Dispersions
#'
#' S3 class with composable base records for the Bregman
#' divergence penalty.
#'
#' ## Z-Dependent Dispersion
#'
#' **Paper**: the primal dispersion $\chi_Z$ shifts and
#' sign-flips a base $\chi$:
#' $$\chi_Z(\gamma) = v\,\chi\bigl(\sigma(\gamma - o)\bigr)$$
#' where $o$ is the offset (balancing target), $\sigma =
#' \mathrm{sign}(o)$ is the sign-flip, and $v \ge 0$ is the
#' penalty scale.
#'
#' **Code** (conjugate, used by the dual solver):
#' with $\psi = \sigma\,\phi / v$,
#' $$\chi^*_Z(\phi) = o\,\phi + v\,\chi^*(\psi)$$
#' $$\dot\chi^*_Z(\phi) = o + \sigma\,\dot\chi^*(\psi)$$
#' $$\ddot\chi^*_Z(\phi) = \ddot\chi^*(\psi) / v$$
#' The $v$ cancels in $\dot\chi^*_Z$ because $\partial\psi /
#' \partial\phi = \sigma/v$ and $\dot\chi^*(\psi)\cdot\sigma/v
#' \cdot v = \sigma\,\dot\chi^*(\psi)$.
#'
#' ## Dispersion Table
#'
#' | Name | $\chi^*(\phi)$ | $\dot\chi^*(\phi)$ | $\ddot\chi^*(\phi)$ |
#' |------|----------------|--------------------|--------------------|
#' | Quadratic | $\phi^2/2$ | $\phi$ | $1$ |
#' | Entropy | $e^\phi$ | $e^\phi$ | $e^\phi$ |
#' | Softplus | $\log(1+e^\phi)$ | $\sigma(\phi)$ | $\sigma(\phi)(1-\sigma(\phi))$ |
#' | Shifted entropy | $\phi + e^\phi$ | $1 + e^\phi$ | $e^\phi$ |
#' | Pos. quadratic | $(\phi_+)^2/2$ | $\phi_+$ | $\mathbf{1}(\phi>0)$ |

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

#' ## Constructor
#'
#' `dispersion(base, offset, sigma, v)` builds a Z-dependent
#' dispersion from a base record. The S3 methods `.dchis`,
#' `.ddchis`, `.chis` apply the composition formulas above.
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
#' **Paper**: for ATE weights, $\chi_Z(\gamma) = \chi(W\gamma)$
#' with $W \in \{-1, +1\}$. Treated units ($W = +1$) get positive
#' weights; control units ($W = -1$) get negative weights.
#'
#' **Code**: `dispersion(base, sigma = W)` — the sign-flip goes
#' into the $\sigma$ slot, which maps $\phi \mapsto W\phi$
#' before applying the base $\dot\chi^*$.
sign_flip_dispersion = function(base, W) dispersion(base, sigma = W)

#' ### `target_scaled_entropy_dispersion`
#'
#' **Paper**: $\chi^*_Z(\phi) = r\,\phi + e^{\mathrm{sign}(r)\,\phi}$
#' giving $\hat\gamma = r + \mathrm{sign}(r)\,e^{\mathrm{sign}(r)\,\phi}$.
#' Floors $|\hat\gamma|$ at $|r|$ (the Riesz representer
#' approximation), then exponential on the excess.
#'
#' **Code**: `dispersion(entropy_base, offset = r, sigma = sign(r))`.
#' This is the default weight dispersion in the survival pipeline.
target_scaled_entropy_dispersion = function(r)
  dispersion(entropy_base, offset = r, sigma = sign(r))

#' ### `variance_weighted_quadratic_dispersion`
#'
#' **Paper**: $\chi_Z(\gamma) = \frac{v}{2}(\gamma - r)^2$ on
#' $\mathrm{sign}(\gamma) = \mathrm{sign}(r)$, $|\gamma| \ge |r|$.
#' The scale $v$ allows heteroskedastic weighting.
#'
#' **Code**: `dispersion(pos_quadratic_base, offset = r,
#' sigma = sign(r), v = v)`.
variance_weighted_quadratic_dispersion = function(r, v)
  dispersion(pos_quadratic_base, offset = r, sigma = sign(r), v = v)
