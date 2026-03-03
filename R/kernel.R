## Kernel Bregman infrastructure.
##
## General penalized dual solver, kernel constructors (isotropic + product),
## dispersions, and S3 interface for time-varying evaluation.
##
## Replaces: balancing/R/{kernels,bregman,dispersions}.R

source("R/utils.R")


# ============================================================
# Utilities
# ============================================================

atleast_2d = function(Z) {
  if (is.null(dim(Z))) Z = matrix(Z, ncol = 1)
  if (length(dim(Z)) == 1) dim(Z) = c(dim(Z), 1)
  Z
}

with.treatment = function(z, a) { z2 = z; z2[, 1] = a; z2 }

.chol_solve = function(R, b) {
  backsolve(R, forwardsolve(R, b, upper.tri = TRUE, transpose = TRUE))
}

.cg_solve = function(matvec, b, x0 = NULL, precond = NULL, tol = 1e-8, maxiter = 200) {
  x = if (!is.null(x0)) x0 else numeric(length(b))
  r = b - matvec(x)
  z = if (!is.null(precond)) precond(r) else r
  p = z
  rz = sum(r * z)
  bnorm = sqrt(sum(b^2))
  if (bnorm == 0) return(list(x = x, iter = 0, resid = 0))
  for (k in 1:maxiter) {
    Ap = matvec(p)
    pAp = sum(p * Ap)
    if (pAp <= 0) break
    a = rz / pAp
    x = x + a * p
    r = r - a * Ap
    if (sqrt(sum(r^2)) < tol * bnorm) break
    z = if (!is.null(precond)) precond(r) else r
    rz_new = sum(r * z)
    p = z + (rz_new / rz) * p
    rz = rz_new
  }
  list(x = x, iter = k, resid = sqrt(sum(r^2)))
}


# ============================================================
# Kernel constructors
# ============================================================

#' Matern kernel. null_basis: Z -> matrix(n, d) for ker(rho).
#' Default: constant functions (intercept). No null space:
#'   function(Z) matrix(nrow = nrow(atleast_2d(Z)), ncol = 0)
matern_kernel = function(sigma = 1, nu = 3/2,
                         null_basis = function(Z) matrix(1, nrow(atleast_2d(Z)), 1)) {
  pointwise = function(z, zp) {
    z = atleast_2d(z)
    d2 = colSums((t(z) - c(zp))^2)
    r = (sqrt(2 * nu) / sigma) * sqrt(d2)
    if (nu == 1/2) return(exp(-r))
    if (nu == 3/2) return((1 + r) * exp(-r))
    if (nu == 5/2) return((1 + r + r^2 / 3) * exp(-r))
    ifelse(r == 0, 1, (2^(1 - nu) / gamma(nu)) * r^nu * besselK(r, nu))
  }
  structure(pointwise, class = "kernel", type = "matern", sigma = sigma, nu = nu,
            null_basis = null_basis)
}

gaussian_kernel = function(sigma = 1,
                           null_basis = function(Z) matrix(1, nrow(atleast_2d(Z)), 1)) {
  pointwise = function(z, zp) {
    z = atleast_2d(z)
    d2 = colSums((t(z) - c(zp))^2)
    exp(-d2 / (2 * sigma^2))
  }
  structure(pointwise, class = "kernel", type = "gaussian", sigma = sigma,
            null_basis = null_basis)
}

#' Direct product: k((w,x),(w',x')) = k_x(x,x') * 1{w=w'}.
#' One copy of ker(rho_x) per grouping level.
#' iw: column indices in Z that hold grouping variable(s).
#' levels (required): for single iw, a vector of values. For multiple iw,
#'   a list of per-column level vectors (Cartesian product taken automatically).
#'   Defines the direct product decomposition — data must only contain
#'   declared levels.
direct_product = function(k_x, iw = 1, levels) {
  # Expand list of per-column levels to pasted Cartesian product
  if (is.list(levels)) {
    grid = do.call(expand.grid, levels)
    levels = apply(grid, 1, paste, collapse = "_")
  }
  nb_base = attr(k_x, "null_basis")
  nb = function(Z) {
    Z = atleast_2d(Z)
    niw = setdiff(seq_len(ncol(Z)), iw)
    B_base = if (!is.null(nb_base)) nb_base(Z[, niw, drop = FALSE])
             else matrix(1, nrow(Z), 1)
    d_base = ncol(B_base)
    gv = if (length(iw) == 1) Z[, iw]
         else apply(Z[, iw, drop = FALSE], 1, paste, collapse = "_")
    bad = setdiff(unique(gv), levels)
    if (length(bad) > 0)
      stop("direct_product: data contains undeclared levels: ",
           paste(head(bad, 5), collapse = ", "))
    n = nrow(Z)
    B = matrix(0, n, length(levels) * d_base)
    if (d_base > 0) {
      for (l in seq_along(levels)) {
        idx = which(gv == levels[l])
        cols = ((l - 1) * d_base + 1):(l * d_base)
        B[idx, cols] = B_base[idx, ]
      }
    }
    B
  }
  structure(list(k_x = k_x, iw = iw, levels = levels),
            class = c("product_kernel", "kernel"),
            null_basis = nb)
}

#' Evaluate null-space basis at points Z.
null_basis = function(Z, kern) {
  Z = atleast_2d(Z)
  nb = attr(kern, "null_basis")
  if (is.null(nb)) return(matrix(1, nrow(Z), 1))
  nb(Z)
}


# ============================================================
# Kernel matrix computation
# ============================================================

kernel_matrix = function(A, B, kern) UseMethod("kernel_matrix", kern)

kernel_matrix.product_kernel = function(A, B, kern) {
  A = atleast_2d(A); B = atleast_2d(B)
  iw = kern$iw
  niw = setdiff(1:ncol(A), iw)
  Kx = kernel_matrix(A[, niw, drop = FALSE], B[, niw, drop = FALSE], kern$k_x)
  mask = matrix(TRUE, nrow(A), nrow(B))
  for (j in iw) mask = mask & outer(A[, j], B[, j], "==")
  Kx * mask
}

kernel_matrix.kernel = function(A, B, kern) {
  A = atleast_2d(A); B = atleast_2d(B)
  ktype = attr(kern, "type")
  if (!is.null(ktype) && ktype %in% c("matern", "gaussian")) {
    D2 = .fast_sqdist(A, B)
    return(.kernel_from_D2(D2, kern))
  }
  outerproduct(A, B, kern)
}

kernel_matrix.default = kernel_matrix.kernel

outerproduct = function(A, B, K) {
  A = atleast_2d(A); B = atleast_2d(B)
  O = matrix(0, nrow(A), nrow(B))
  for (ib in 1:nrow(B)) O[, ib] = K(A, B[ib, , drop = FALSE])
  O
}

.kernel_from_D2 = function(D2, kern) {
  ktype = attr(kern, "type")
  sig = attr(kern, "sigma")
  if (ktype == "matern") {
    nu = attr(kern, "nu")
    R = (sqrt(2 * nu) / sig) * sqrt(D2)
    if (nu == 1/2) return(exp(-R))
    if (nu == 3/2) return((1 + R) * exp(-R))
    if (nu == 5/2) return((1 + R + R^2 / 3) * exp(-R))
    return(ifelse(R == 0, 1, (2^(1 - nu) / gamma(nu)) * R^nu * besselK(R, nu)))
  }
  if (ktype == "gaussian") return(exp(-D2 / (2 * sig^2)))
  stop("Unknown kernel type: ", ktype)
}

.fast_sqdist = function(A, B) {
  ssA = rowSums(A^2); ssB = rowSums(B^2)
  D2 = outer(ssA, ssB, "+") - 2 * tcrossprod(A, B)
  D2[D2 < 0] = 0
  D2
}


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

#' Target-shifted quadratic: gamma = r + phi. No sign constraint.
#' Works with mixed-sign r (ATE) where TSE's sign constraint breaks.
target_shifted_quadratic = function(r)
  dispersion(quadratic_base, offset = r)


# ============================================================
# Link functions (survival-specific, used by cts-gridfree.R)
# ============================================================

exp_link = function() {
  list(
    name = "exp",
    Q_surv = function(horizon) 100,
    mu = function(eta) exp(eta),
    curv = function(eta, mu) mu,
    resid = function(mu, Y, w) w * (mu - Y),
    predict = function(phi) exp(phi),
    cum_surv = function(haz_fn, Q, du) {
      running = rep(0, length(haz_fn(1)))
      h_prev = haz_fn(1)
      for (k in 2:(Q + 1)) {
        h_curr = haz_fn(k)
        running = running + (h_prev + h_curr) / 2 * du
        h_prev = h_curr
      }
      exp(-running)
    },
    cum_rmst = function(haz_fn, Q, du) {
      n = length(haz_fn(1))
      running_cumhaz = rep(0, n); running_rmst = rep(0, n)
      S_prev = rep(1, n); h_prev = haz_fn(1)
      for (k in 2:(Q + 1)) {
        h_curr = haz_fn(k)
        running_cumhaz = running_cumhaz + (h_prev + h_curr) / 2 * du
        S_curr = exp(-running_cumhaz)
        running_rmst = running_rmst + (S_prev + S_curr) / 2 * du
        S_prev = S_curr; h_prev = h_curr
      }
      running_rmst
    },
    surv_curve = function(haz_fn, Q, du) {
      n = length(haz_fn(1))
      S = matrix(1, n, Q + 1); running = rep(0, n); h_prev = haz_fn(1)
      for (k in 2:(Q + 1)) {
        h_curr = haz_fn(k)
        running = running + (h_prev + h_curr) / 2 * du
        S[, k] = exp(-running); h_prev = h_curr
      }
      S
    }
  )
}

logistic_link = function() {
  list(
    name = "logistic",
    Q_surv = function(horizon) as.integer(horizon),
    mu = function(eta) plogis(eta),
    curv = function(eta, mu) mu * (1 - mu),
    resid = function(mu, Y, w) w * (mu - Y),
    predict = function(phi) plogis(phi),
    cum_surv = function(haz_fn, Q, du) {
      log_surv = rep(0, length(haz_fn(2)))
      for (k in 2:(Q + 1)) {
        log_surv = log_surv + log(pmax(1 - haz_fn(k), 1e-15))
      }
      exp(log_surv)
    },
    cum_rmst = function(haz_fn, Q, du) {
      n = length(haz_fn(2))
      log_surv = rep(0, n); running_rmst = rep(0, n); S_prev = rep(1, n)
      for (k in 2:(Q + 1)) {
        log_surv = log_surv + log(pmax(1 - haz_fn(k), 1e-15))
        S_curr = exp(log_surv)
        running_rmst = running_rmst + (S_prev + S_curr) / 2 * du
        S_prev = S_curr
      }
      running_rmst
    },
    surv_curve = function(haz_fn, Q, du) {
      n = length(haz_fn(2))
      S = matrix(1, n, Q + 1); log_surv = rep(0, n)
      for (k in 2:(Q + 1)) {
        log_surv = log_surv + log(pmax(1 - haz_fn(k), 1e-15))
        S[, k] = exp(log_surv)
      }
      S
    }
  )
}


# ============================================================
# General kernel Bregman solver
# ============================================================
#
# Solves:  min_f  P_n[chi*_Z(f(Z))] + (eta/2) ||f||_rho^2 - dpsiZ(f)
#
# phi = K alpha + B beta, where B = null_basis(Z, kern) spans ker(rho).
# Two-step Newton: H = blkdiag(K, I_d) M where
#   M = [diag(wc)K + n*eta*I, diag(wc)B; B'diag(wc)K, B'diag(wc)B]
# Stores lazy M_inv for one-step bootstrap.

kernel_bregman = function(Z, kern, eta, dispersion,
                          target, target_null = NULL,
                          w = NULL, K = NULL,
                          alpha0 = NULL, beta0 = NULL,
                          maxiter = 25, tol = 1e-3,
                          log = null_logger) {
  Z = atleast_2d(Z)
  n = nrow(Z)
  if (is.null(K)) K = kernel_matrix(Z, Z, kern)
  K = K + (1e-8 * sum(diag(K)) / n) * diag(n)
  if (is.null(w)) w = rep(1, n)

  B = null_basis(Z, kern)
  d = ncol(B)
  if (is.null(target_null)) target_null = rep(0, d)

  alpha = if (!is.null(alpha0)) alpha0 else rep(0, n)
  beta = if (d > 0) {
    if (!is.null(beta0)) rep_len(as.numeric(beta0), d) else rep(0, d)
  } else numeric(0)

  K_chol = chol(K)

  #' Lambda bindings: dchis(phi) reads like dot chi*(phi).
  disp = dispersion
  dchis  = function(phi) .dchis(disp, phi)
  ddchis = function(phi) .ddchis(disp, phi)
  chis   = function(phi) .chis(disp, phi)

  fg_kernel = function(theta) {
    a = theta[1:n]
    b = if (d > 0) theta[(n+1):(n+d)] else numeric(0)
    phi = as.vector(K %*% a)
    if (d > 0) phi = phi + as.vector(B %*% b)
    mu = dchis(phi)
    if (!all(is.finite(mu))) return(list(val = Inf, grad = theta * 0))
    Ka = as.vector(K %*% a)
    val = sum(w * chis(phi)) + (n * eta / 2) * sum(a * Ka) -
          sum(target * a) - if (d > 0) sum(target_null * b) else 0
    g_a = as.vector(K %*% (w * mu)) + n * eta * Ka - target
    g_b = if (d > 0) as.vector(crossprod(B, w * mu)) - target_null else numeric(0)
    list(val = val, grad = c(g_a, g_b))
  }

  hv_kernel = function(theta, v) {
    a = theta[1:n]
    b = if (d > 0) theta[(n+1):(n+d)] else numeric(0)
    phi = as.vector(K %*% a)
    if (d > 0) phi = phi + as.vector(B %*% b)
    wc = w * ddchis(phi)
    wc = pmax(wc, 1e-10 * max(abs(wc), 1))
    v_a = v[1:n]
    v_b = if (d > 0) v[(n+1):(n+d)] else numeric(0)
    Kv_a = as.vector(K %*% v_a)
    inner = wc * (Kv_a + if (d > 0) as.vector(B %*% v_b) else 0)
    Hv_a = as.vector(K %*% inner) + n * eta * Kv_a
    Hv_b = if (d > 0) as.vector(crossprod(B, inner)) else numeric(0)
    c(Hv_a, Hv_b)
  }

  #' Bregman divergence D_chi*(phi_new, phi_old) = sum w_i [chi*(phi_new_i)
  #' - chi*(phi_old_i) - dchi*(phi_old_i) * (phi_new_i - phi_old_i)]
  bregman_conv = function(theta_old, theta_new) {
    phi_fn = function(th) {
      a = th[1:n]; b = if (d > 0) th[(n+1):(n+d)] else numeric(0)
      phi = as.vector(K %*% a)
      if (d > 0) phi = phi + as.vector(B %*% b)
      phi
    }
    phi_old = phi_fn(theta_old)
    phi_new = phi_fn(theta_new)
    gam_old = dchis(phi_old)
    sum(w * (chis(phi_new) - chis(phi_old) -
             gam_old * (phi_new - phi_old)))
  }

  phi_from_theta = function(th) {
    a = th[1:n]; b = if (d > 0) th[(n+1):(n+d)] else numeric(0)
    phi = as.vector(K %*% a)
    if (d > 0) phi = phi + as.vector(B %*% b)
    phi
  }
  weight_fn = function(th) dchis(phi_from_theta(th))

  delta0 = if (isTRUE(disp$base$quadratic)) Inf else 1.0
  tr_res = .tr_newton(fg_kernel, hv_kernel, c(alpha, beta),
                      maxiter = maxiter, tol = tol, delta0 = delta0,
                      conv_fn = bregman_conv, weight_fn = weight_fn)

  alpha = tr_res$theta[1:n]
  beta  = if (d > 0) tr_res$theta[(n+1):(n+d)] else numeric(0)
  iter  = tr_res$iter
  status = if (tr_res$iter < maxiter) "converged"
           else sprintf("maxiter (%d)", maxiter)

  log$info("iter=%d (%s)", iter, status)
  log$data("trace", tr_res$trace)

  # Lazy M_inv for one-step bootstrap
  M_lu_cache = NULL
  get_M_lu = function() {
    if (!is.null(M_lu_cache)) return(M_lu_cache)
    phi_f = as.vector(K %*% alpha)
    if (d > 0) phi_f = phi_f + as.vector(B %*% beta)
    wc_f = w * ddchis(phi_f)
    wc_f = pmax(wc_f, 1e-10 * max(abs(wc_f), 1))
    wcK_f = sweep(K, 1, wc_f, "*")
    M_f = if (d > 0) {
      wcB_f = wc_f * B
      rbind(cbind(wcK_f + n * eta * diag(n), wcB_f),
            cbind(crossprod(B, wcK_f), crossprod(B, wcB_f)))
    } else {
      wcK_f + n * eta * diag(n)
    }
    m = nrow(M_f); M_f = M_f + (1e-8 * sum(diag(M_f)) / m) * diag(m)
    M_lu_cache <<- solve(M_f)
    M_lu_cache
  }

  structure(list(
    alpha = alpha, beta = beta,
    kern = kern, Z = Z, K = K, B = B, eta = eta,
    K_chol = K_chol, get_M_lu = get_M_lu,
    dispersion = dispersion, iter = iter, status = status,
    tr_trace = tr_res$trace,
    w = w, target = target, target_null = target_null
  ), class = "kernel_bregman")
}


# ============================================================
# Prediction
# ============================================================

predict.kernel_bregman = function(obj, newZ, ...) {
  phi   = \(Z) .phi(obj, Z)
  dchis = \(p) .dchis(obj$dispersion, p)
  dchis(phi(newZ))
}



# ============================================================
# Target projection
# ============================================================

#' Project signed measure mu = sum r_i delta_{z_i} onto representer basis.
#' Returns c = [target; target_null] of length n + d.
project_target = function(basis_Z, kern, measure, K_cross = NULL) {
  basis_Z = atleast_2d(basis_Z)
  measure$Z = atleast_2d(measure$Z)
  if (is.null(K_cross))
    K_cross = kernel_matrix(basis_Z, measure$Z, kern)
  B = null_basis(measure$Z, kern)
  c(as.vector(K_cross %*% measure$r),
    if (ncol(B) > 0) as.vector(crossprod(B, measure$r)) else numeric(0))
}

#' Split project_target output into target (n) and target_null (d).
split_target = function(c_vec, n) {
  list(target = c_vec[1:n],
       target_null = if (length(c_vec) > n) c_vec[(n + 1):length(c_vec)] else numeric(0))
}


# ============================================================
# Cross-validation
# ============================================================

#' Out-of-sample dual loss: mean[ chi*(phi_test) - target_test * phi_test ].
cv_dual_loss = function(fit, Z_test, target_test, dispersion_test = NULL) {
  phi  = \(Z) .phi(fit, Z)
  disp = if (!is.null(dispersion_test)) dispersion_test else fit$dispersion
  chis = \(p) .chis(disp, p)
  phi_test = phi(Z_test)
  mean(chis(phi_test) - target_test * phi_test)
}

#' LOO CV loss (linearized at convergence).
loo_loss = function(fit) {
  K = fit$K; n = nrow(K)
  phi = as.vector(K %*% fit$alpha)
  d = length(fit$beta)
  if (d > 0) phi = phi + as.vector(fit$B %*% fit$beta)
  w = if (!is.null(fit$w)) fit$w else rep(1, n)
  ddchis = \(p) .ddchis(fit$dispersion, p)
  wc = w * ddchis(phi)
  wc = pmax(wc, 1e-10 * max(abs(wc), 1))
  raw_resid = n * fit$eta * fit$alpha / wc
  sqw = sqrt(wc)
  Kt = sweep(sweep(K, 1, sqw, "*"), 2, sqw, "*")
  evd = eigen(Kt, symmetric = TRUE)
  shrink = pmax(evd$values, 0) / (pmax(evd$values, 0) + n * fit$eta)
  H_diag = rowSums(evd$vectors^2 * rep(shrink, each = n))
  loo_resid = raw_resid / pmax(1 - H_diag, 1e-10)
  mean(wc * loo_resid^2)
}


# ============================================================
# One-step Newton (bootstrap)
# ============================================================

#' One-step Newton: gradient at current (alpha, beta) with new targets/weights,
#' solved using stored Hessian factorization (M_inv from get_M_lu).
onestep_bregman = function(model, target_new, target_null_new = NULL, w_new = NULL) {
  K = model$K; n = nrow(K)
  B = model$B; d = ncol(B)
  alpha = model$alpha
  beta = if (d > 0) model$beta else numeric(0)
  eta = model$eta; disp = model$dispersion

  if (is.null(w_new)) w_new = model$w
  if (is.null(target_null_new)) target_null_new = model$target_null

  dchis = \(p) .dchis(disp, p)
  phi = as.vector(K %*% alpha)
  if (d > 0) phi = phi + as.vector(B %*% beta)
  mu = dchis(phi)
  Ka = as.vector(K %*% alpha)
  g_a = as.vector(K %*% (w_new * mu)) + n * eta * Ka - target_new
  v_a = .chol_solve(model$K_chol, -g_a)

  if (d > 0) {
    g_b = as.vector(crossprod(B, w_new * mu)) - target_null_new
    rhs = c(v_a, -g_b)
  } else {
    rhs = v_a
  }
  delta = as.vector(model$get_M_lu() %*% rhs)
  d_a = delta[1:n]
  d_b = if (d > 0) delta[(n + 1):(n + d)] else numeric(0)
  model$alpha = alpha + d_a
  model$beta  = if (d > 0) beta + d_b else model$beta
  model
}


# ============================================================
# Z accessors: Z = (u, W, X) with time in col 1
# ============================================================

u = function(Z) Z[, 1]
w = function(Z) Z[, 2]
x = function(Z) Z[, -c(1, 2), drop = FALSE]


# ============================================================
# S3 generics
# ============================================================

.bind    = function(model, Z, ...) UseMethod(".bind")
.lambda  = function(model, Z, ...) UseMethod(".lambda")
.gamma   = function(model, Z, ...) UseMethod(".gamma")
.phi     = function(model, Z, ...) UseMethod(".phi")
one_step = function(model, ...) UseMethod("one_step")
.gamma_vec = function(model, ...) UseMethod(".gamma_vec")


# ============================================================
# bind: cache covariate distances. Returns same type.
# ============================================================

.bind.kernel_bregman = function(model, Z, ...) {
  Z = atleast_2d(Z)
  Z_centers = model$Z
  u_centers = Z_centers[, 1]
  cov_centers = Z_centers[, -1, drop = FALSE]
  kern = model$kern

  if (inherits(kern, "product_kernel")) {
    iw_cov = kern$iw[kern$iw > 1] - 1
    niw_cov = setdiff(1:ncol(cov_centers), iw_cov)
    model$.K_x_cache = kernel_matrix(Z[, niw_cov, drop = FALSE],
                                      cov_centers[, niw_cov, drop = FALSE], kern$k_x)
    model$.cov_centers = cov_centers
    model$.iw_cov = iw_cov
    model$.bind_method = "product"
  } else {
    model$.D2_cov_cache = .fast_sqdist(Z, cov_centers)
    model$.bind_method = "isotropic"
  }
  model$.u_centers = u_centers
  model$.n_bound = nrow(Z)
  model
}


# ============================================================
# [ : subset a bound model
# ============================================================

`[.kernel_bregman` = function(x, i, ...) {
  if (!is.null(x$.D2_cov_cache))
    x$.D2_cov_cache = x$.D2_cov_cache[i, , drop = FALSE]
  if (!is.null(x$.K_x_cache))
    x$.K_x_cache = x$.K_x_cache[i, , drop = FALSE]
  if (!is.null(x$.n_bound))
    x$.n_bound = sum(rep(TRUE, x$.n_bound)[i])
  x
}


# ============================================================
# phi: linear predictor. Z = (u, covariates), time in col 1.
# ============================================================

.phi.kernel_bregman = function(model, Z, ...) {
  Z = atleast_2d(Z)
  kern = model$kern

  if (!is.null(model$.bind_method) && model$.bind_method == "isotropic") {
    delta_u2 = outer(u(Z), model$.u_centers, function(a, b) (a - b)^2)
    D2 = model$.D2_cov_cache + delta_u2
    D2[D2 < 0] = 0
    K_eval = .kernel_from_D2(D2, kern)
  } else if (!is.null(model$.bind_method) && model$.bind_method == "product") {
    cov = Z[, -1, drop = FALSE]
    K_x = model$.K_x_cache
    u_mask = outer(u(Z), model$.u_centers, "==")
    mask = u_mask
    for (j in model$.iw_cov)
      mask = mask & outer(cov[, j], model$.cov_centers[, j], "==")
    K_eval = K_x * mask
  } else {
    K_eval = kernel_matrix(Z, model$Z, kern)
  }

  p = as.vector(K_eval %*% model$alpha)
  d = length(model$beta)
  if (d > 0) {
    if (d == 1) {
      p = p + model$beta[1]
    } else {
      p = p + as.vector(null_basis(Z, kern) %*% model$beta)
    }
  }
  p
}


# ============================================================
# lambda: hazard rate = dchis(phi) / du_bin
# Model must have $du_bin (set by fit_lambda).
# ============================================================

.lambda.kernel_bregman = function(model, Z, ...) {
  pmax(.dchis(model$dispersion, .phi(model, Z)) / model$du_bin, 0)
}


# ============================================================
# gamma: weight = dchis_gamma(phi, Z)
# Model must have $predict_gamma (set by fit_gamma).
# ============================================================

.gamma.kernel_bregman = function(model, Z, ...) {
  model$predict_gamma(.phi(model, Z), Z)
}

.gamma_vec.kernel_bregman = function(model, Z, ...) {
  model$predict_gamma(.phi(model, Z), Z)
}


# ============================================================
# one_step: S3 wrapper for onestep_bregman
# ============================================================

one_step.kernel_bregman = function(model, target_new, target_null_new = NULL, w_new = NULL, ...) {
  onestep_bregman(model, target_new, target_null_new, w_new)
}

