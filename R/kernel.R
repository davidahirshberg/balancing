## Kernel Bregman infrastructure.
##
## General penalized dual solver, kernel constructors (isotropic + product),
## dispersions, and S3 interface for time-varying evaluation.
##
## Replaces: balancing/R/{kernels,bregman,dispersions}.R


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
direct_product = function(k_x, iw = 1, levels = NULL) {
  nb_base = attr(k_x, "null_basis")
  nb = function(Z) {
    Z = atleast_2d(Z)
    niw = setdiff(seq_len(ncol(Z)), iw)
    B_base = if (!is.null(nb_base)) nb_base(Z[, niw, drop = FALSE])
             else matrix(1, nrow(Z), 1)
    d_base = ncol(B_base)
    gv = if (length(iw) == 1) Z[, iw]
         else apply(Z[, iw, drop = FALSE], 1, paste, collapse = "_")
    lvls = if (!is.null(levels)) levels else sort(unique(gv))
    n = nrow(Z)
    B = matrix(0, n, length(lvls) * d_base)
    for (l in seq_along(lvls)) {
      idx = which(gv == lvls[l])
      cols = ((l - 1) * d_base + 1):(l * d_base)
      B[idx, cols] = B_base[idx, ]
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
# Dispersions
# ============================================================

entropy_dispersion = function() {
  list(name = "entropy", linear = FALSE,
       chi = function(g) ifelse(g > 0, g * log(g) - g, Inf),
       chis = function(phi) exp(phi),
       dchis = function(phi) exp(phi),
       ddchis = function(phi) exp(phi),
       resid = function(g, Y, w) w * (g - Y),
       predict = function(phi) exp(phi))
}

softplus_dispersion = function() {
  list(name = "softplus", linear = FALSE,
       chi = function(g) { g = pmin(pmax(g, 1e-15), 1 - 1e-15); g * log(g) + (1 - g) * log(1 - g) },
       chis = function(phi) log(1 + exp(phi)),
       dchis = function(phi) plogis(phi),
       ddchis = function(phi) { p = plogis(phi); p * (1 - p) },
       resid = function(g, Y, w) w * (g - Y),
       predict = function(phi) plogis(phi))
}

quadratic_dispersion = function() {
  list(name = "quadratic", linear = TRUE,
       chi = function(g) g^2 / 2,
       chis = function(phi) phi^2 / 2,
       dchis = function(phi) phi,
       ddchis = function(phi) rep(1, length(phi)),
       resid = function(g, Y, w) w * (g - Y),
       predict = function(phi) phi)
}

shifted_entropy_dispersion = function() {
  list(name = "shifted_entropy", linear = FALSE,
       chi = function(g) ifelse(g >= 1, (g - 1) * log(g - 1) - (g - 1) + 1, Inf),
       chis = function(phi) phi + exp(phi),
       dchis = function(phi) 1 + exp(phi),
       ddchis = function(phi) exp(phi),
       resid = function(g, Y, w) w * (g - Y),
       predict = function(phi) 1 + exp(phi))
}

#' Sign-flip: chi_Z(gamma) = chi(W * gamma). W in {-1,+1}.
signflip = function(base, W) {
  list(name = paste0("signflip_", base$name), linear = base$linear,
       chi = function(g) base$chi(W * g),
       chis = function(phi) base$chis(W * phi),
       dchis = function(phi) W * base$dchis(W * phi),
       ddchis = function(phi) base$ddchis(W * phi),
       resid = function(g, Y, w) w * (g - Y),
       predict = function(phi) W * base$dchis(W * phi),
       base = base, W = W)
}

#' Target-scaled entropy: gamma = r + sigma*exp(sigma*phi), sigma = sign(r).
target_scaled_entropy = function(r) {
  sigma = sign(r)
  list(name = "target_scaled_entropy", linear = FALSE,
       chi = function(g) { u = sigma * g - abs(r); ifelse(u > 0, u * log(u) - u + 1, Inf) },
       chis = function(phi) r * phi + exp(sigma * phi) - 1,
       dchis = function(phi) r + sigma * exp(sigma * phi),
       ddchis = function(phi) exp(sigma * phi),
       resid = function(g, Y, w) w * (g - Y),
       predict = function(phi) r + sigma * exp(sigma * phi),
       sigma = sigma, r = r)
}

#' Variance-weighted quadratic with sign-flip and target shift.
variance_weighted_quadratic = function(r, v) {
  sigma = sign(r)
  list(name = "variance_weighted_quadratic", linear = FALSE,
       chi = function(g) { u = sigma * g - abs(r); ifelse(u >= 0, v / 2 * u^2, Inf) },
       chis = function(phi) { sp = pmax(sigma * phi, 0); r * phi + sp^2 / (2 * v) },
       dchis = function(phi) r + sigma * pmax(sigma * phi, 0) / v,
       ddchis = function(phi) as.numeric(sigma * phi > 0) / v,
       resid = function(g, Y, w) w * (g - Y),
       predict = function(phi) r + sigma * pmax(sigma * phi, 0) / v,
       sigma = sigma, r = r, v = v)
}


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
                          maxiter = 25, tol = 1e-6) {
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
  if (isTRUE(dispersion$linear)) maxiter = 1

  grad_norm_0 = NULL
  status = sprintf("maxiter (%d)", maxiter)
  for (iter in 1:maxiter) {
    phi = as.vector(K %*% alpha)
    if (d > 0) phi = phi + as.vector(B %*% beta)

    mu = dispersion$dchis(phi)
    wc = w * dispersion$ddchis(phi)
    if (any(!is.finite(mu)) || any(!is.finite(wc))) {
      status = sprintf("non-finite at iter %d", iter); break
    }
    wc = pmax(wc, 1e-10 * max(abs(wc), 1))

    Ka = as.vector(K %*% alpha)
    g_a = as.vector(K %*% (w * mu)) + n * eta * Ka - target
    v_a = .chol_solve(K_chol, -g_a)

    wcK = sweep(K, 1, wc, "*")
    M_aa = wcK + n * eta * diag(n)

    if (d > 0) {
      g_b = as.vector(crossprod(B, w * mu)) - target_null
      wcB = wc * B
      M = rbind(cbind(M_aa, wcB), cbind(crossprod(B, wcK), crossprod(B, wcB)))
      m = nrow(M); M = M + (1e-8 * sum(diag(M)) / m) * diag(m)
      delta = solve(M, c(v_a, -g_b))
      d_a = delta[1:n]; d_b = delta[(n + 1):(n + d)]
    } else {
      d_a = as.vector(solve(M_aa, v_a)); d_b = numeric(0)
    }

    step_size = max(abs(d_a), if (d > 0) max(abs(d_b)) else 0)
    grad_norm = sqrt(sum(g_a^2) + if (d > 0) sum(g_b^2) else 0)
    if (!is.finite(step_size) || !is.finite(grad_norm)) {
      status = sprintf("non-finite step at iter %d", iter); break
    }
    alpha = alpha + d_a; beta = beta + d_b
    if (is.null(grad_norm_0)) grad_norm_0 = grad_norm
    if (step_size < tol || grad_norm < tol * grad_norm_0) {
      status = "converged"; break
    }
  }

  # Lazy M_inv for one-step bootstrap
  M_lu_cache = NULL
  get_M_lu = function() {
    if (!is.null(M_lu_cache)) return(M_lu_cache)
    phi_f = as.vector(K %*% alpha)
    if (d > 0) phi_f = phi_f + as.vector(B %*% beta)
    wc_f = w * dispersion$ddchis(phi_f)
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
    w = w, target = target, target_null = target_null
  ), class = "kernel_bregman")
}


# ============================================================
# Prediction
# ============================================================

predict.kernel_bregman = function(obj, newZ, ...) {
  obj$dispersion$dchis(predict_phi(obj, newZ))
}

predict_phi = function(obj, newZ) {
  newZ = atleast_2d(newZ)
  Knew = kernel_matrix(newZ, obj$Z, obj$kern)
  phi = as.vector(Knew %*% obj$alpha)
  d = length(obj$beta)
  if (d > 0) {
    B_new = null_basis(newZ, obj$kern)
    phi = phi + as.vector(B_new %*% obj$beta)
  }
  phi
}

predict_phi_alpha = function(alpha, beta, model, newZ) {
  newZ = atleast_2d(newZ)
  Knew = kernel_matrix(newZ, model$Z, model$kern)
  phi = as.vector(Knew %*% alpha)
  d = length(beta)
  if (d > 0) {
    B_new = null_basis(newZ, model$kern)
    phi = phi + as.vector(B_new %*% beta)
  }
  phi
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
  phi_test = predict_phi(fit, Z_test)
  disp = if (!is.null(dispersion_test)) dispersion_test else fit$dispersion
  mean(disp$chis(phi_test) - target_test * phi_test)
}

#' LOO CV loss (linearized at convergence).
loo_loss = function(fit) {
  K = fit$K; n = nrow(K)
  phi = as.vector(K %*% fit$alpha)
  d = length(fit$beta)
  if (d > 0) phi = phi + as.vector(fit$B %*% fit$beta)
  w = if (!is.null(fit$w)) fit$w else rep(1, n)
  wc = w * fit$dispersion$ddchis(phi)
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

  phi = as.vector(K %*% alpha)
  if (d > 0) phi = phi + as.vector(B %*% beta)
  mu = disp$dchis(phi)
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
  list(alpha = alpha + d_a, beta = beta + d_b)
}


# ============================================================
# S3 generics
# ============================================================

bind     = function(model, Z, ...) UseMethod("bind")
lambda   = function(model, u, Z, ...) UseMethod("lambda")
gamma    = function(model, u, Z, ...) UseMethod("gamma")
phi      = function(model, ...) UseMethod("phi")
one_step = function(model, ...) UseMethod("one_step")
phi_vec  = function(model, ...) UseMethod("phi_vec")
gamma_vec = function(model, ...) UseMethod("gamma_vec")


# ============================================================
# bind: cache covariate distances. Returns same type.
# ============================================================

bind.kernel_bregman = function(model, Z, ...) {
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
# phi: linear predictor at scalar time u
# ============================================================

phi.kernel_bregman = function(model, u, Z, ...) {
  Z = atleast_2d(Z)
  kern = model$kern

  if (!is.null(model$.bind_method) && model$.bind_method == "isotropic") {
    delta_u2 = (u - model$.u_centers)^2
    D2 = sweep(model$.D2_cov_cache, 2, delta_u2, "+")
    D2[D2 < 0] = 0
    K_eval = .kernel_from_D2(D2, kern)
  } else if (!is.null(model$.bind_method) && model$.bind_method == "product") {
    K_x = model$.K_x_cache
    u_mask = outer(rep(u, nrow(K_x)), model$.u_centers, "==")
    mask = u_mask
    for (j in model$.iw_cov)
      mask = mask & outer(Z[, j], model$.cov_centers[, j], "==")
    K_eval = K_x * mask
  } else {
    K_eval = kernel_matrix(cbind(u, Z), model$Z, kern)
  }

  p = as.vector(K_eval %*% model$alpha)
  d = length(model$beta)
  if (d > 0) {
    if (d == 1) {
      p = p + model$beta[1]
    } else {
      p = p + as.vector(null_basis(cbind(u, Z), kern) %*% model$beta)
    }
  }
  p
}

#' phi_vec: per-subject times. u_vec[j] is the time for subject j.
phi_vec.kernel_bregman = function(model, u_vec, Z, ...) {
  Z = atleast_2d(Z)
  kern = model$kern

  if (!is.null(model$.bind_method) && model$.bind_method == "isotropic") {
    delta_u2 = outer(u_vec, model$.u_centers, function(a, b) (a - b)^2)
    D2 = model$.D2_cov_cache + delta_u2
    D2[D2 < 0] = 0
    K_eval = .kernel_from_D2(D2, kern)
  } else if (!is.null(model$.bind_method) && model$.bind_method == "product") {
    K_x = model$.K_x_cache
    u_mask = outer(u_vec, model$.u_centers, "==")
    mask = u_mask
    for (j in model$.iw_cov)
      mask = mask & outer(Z[, j], model$.cov_centers[, j], "==")
    K_eval = K_x * mask
  } else {
    K_eval = kernel_matrix(cbind(u_vec, Z), model$Z, kern)
  }

  p = as.vector(K_eval %*% model$alpha)
  d = length(model$beta)
  if (d > 0) {
    if (d == 1) {
      p = p + model$beta[1]
    } else {
      p = p + as.vector(null_basis(cbind(u_vec, Z), kern) %*% model$beta)
    }
  }
  p
}


# ============================================================
# lambda: hazard rate = dispersion$dchis(phi) / du_bin
# Model must have $du_bin (set by fit_lambda).
# ============================================================

lambda.kernel_bregman = function(model, u, Z, ...) {
  pmax(model$dispersion$dchis(phi(model, u, Z)) / model$du_bin, 0)
}


# ============================================================
# gamma: weight = dchis_gamma(phi, Z)
# Model must have $predict_gamma (set by fit_gamma).
# ============================================================

gamma.kernel_bregman = function(model, u, Z, ...) {
  model$predict_gamma(phi(model, u, Z), Z)
}

gamma_vec.kernel_bregman = function(model, u_vec, Z, ...) {
  model$predict_gamma(phi_vec(model, u_vec, Z), Z)
}


# ============================================================
# one_step: S3 wrapper for onestep_bregman
# ============================================================

one_step.kernel_bregman = function(model, target_new, target_null_new = NULL, w_new = NULL, ...) {
  onestep_bregman(model, target_new, target_null_new, w_new)
}


# ============================================================
# Legacy compatibility: bind_subjects, eval_phi, eval_phi_vec
# These support estimate.R / survival.R from the balancing package.
# ============================================================

bind_subjects = function(model, Z_eval) {
  Z_eval = atleast_2d(Z_eval)
  Z_centers = model$Z
  u_centers = Z_centers[, 1]
  cov_centers = Z_centers[, -1, drop = FALSE]
  kern = model$kern

  bound = list(u_centers = u_centers, alpha = model$alpha,
               beta = model$beta, kern = kern, Z_eval = Z_eval)

  if (inherits(kern, "product_kernel")) {
    iw_cov = kern$iw[kern$iw > 1] - 1
    niw_cov = setdiff(1:ncol(cov_centers), iw_cov)
    bound$K_x = kernel_matrix(Z_eval[, niw_cov, drop = FALSE],
                              cov_centers[, niw_cov, drop = FALSE], kern$k_x)
    bound$cov_centers = cov_centers
    bound$iw_cov = iw_cov
    bound$method = "product"
  } else {
    ktype = attr(kern, "type")
    if (!is.null(ktype) && ktype %in% c("matern", "gaussian")) {
      bound$D2_cov = .fast_sqdist(Z_eval, cov_centers)
      bound$method = "isotropic"
    } else {
      bound$Z_centers = Z_centers
      bound$method = "generic"
    }
  }
  bound
}

.null_phi = function(bound, Z_pred) {
  d = length(bound$beta)
  if (d == 0) return(0)
  if (d == 1) return(bound$beta[1])
  as.vector(null_basis(Z_pred, bound$kern) %*% bound$beta)
}

eval_phi = function(bound, u, idx = NULL) {
  if (bound$method == "isotropic") {
    D2_cov = if (is.null(idx)) bound$D2_cov else bound$D2_cov[idx, , drop = FALSE]
    delta_u2 = (u - bound$u_centers)^2
    D2 = sweep(D2_cov, 2, delta_u2, "+"); D2[D2 < 0] = 0
    K = .kernel_from_D2(D2, bound$kern)
  } else if (bound$method == "product") {
    K_x = if (is.null(idx)) bound$K_x else bound$K_x[idx, , drop = FALSE]
    n_eval = nrow(K_x)
    u_mask = outer(rep(u, n_eval), bound$u_centers, "==")
    mask = u_mask
    Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
    for (j in bound$iw_cov)
      mask = mask & outer(Z_ev[, j], bound$cov_centers[, j], "==")
    K = K_x * mask
  } else {
    Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
    K = kernel_matrix(cbind(u, Z_ev), bound$Z_centers, bound$kern)
  }
  phi = as.vector(K %*% bound$alpha)
  Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
  phi + .null_phi(bound, cbind(u, Z_ev))
}

eval_phi_vec = function(bound, u_vec, idx = NULL) {
  if (bound$method == "isotropic") {
    D2_cov = if (is.null(idx)) bound$D2_cov else bound$D2_cov[idx, , drop = FALSE]
    delta_u2 = outer(u_vec, bound$u_centers, function(a, b) (a - b)^2)
    D2 = D2_cov + delta_u2; D2[D2 < 0] = 0
    K = .kernel_from_D2(D2, bound$kern)
  } else if (bound$method == "product") {
    K_x = if (is.null(idx)) bound$K_x else bound$K_x[idx, , drop = FALSE]
    u_mask = outer(u_vec, bound$u_centers, "==")
    mask = u_mask
    Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
    for (j in bound$iw_cov)
      mask = mask & outer(Z_ev[, j], bound$cov_centers[, j], "==")
    K = K_x * mask
  } else {
    Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
    K = kernel_matrix(cbind(u_vec, Z_ev), bound$Z_centers, bound$kern)
  }
  phi = as.vector(K %*% bound$alpha)
  Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
  phi + .null_phi(bound, cbind(u_vec, Z_ev))
}
