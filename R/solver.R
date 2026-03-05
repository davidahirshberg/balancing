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
                          w = NULL, K = NULL, K_chol = NULL,
                          alpha0 = NULL, beta0 = NULL,
                          maxiter = 25, tol = 1e-3,
                          log = null_logger) {
  Z = atleast_2d(Z)
  n = nrow(Z)
  if (is.null(K)) K = kernel_matrix(Z, Z, kern)
  if (is.null(K_chol)) {
    K_chol = pchol(K)
  } else if (!inherits(K_chol, "pchol")) {
    # Legacy: caller passed a plain Cholesky factor. Wrap it.
    K_chol = structure(list(R = K_chol, pivot = seq_len(n),
                            rank = n, n = n), class = "pchol")
  }
  if (is.null(w)) w = rep(1, n)

  B = null_basis(Z, kern)
  d = ncol(B)
  if (is.null(target_null)) target_null = rep(0, d)

  alpha = if (!is.null(alpha0)) alpha0 else rep(0, n)
  beta = if (d > 0) {
    if (!is.null(beta0)) rep_len(as.numeric(beta0), d) else rep(0, d)
  } else numeric(0)

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


# ============================================================
# S3 generics
# ============================================================

.bind    = function(model, Z, ...) UseMethod(".bind")
.phi     = function(model, Z, ...) UseMethod(".phi")


# ============================================================
# bind: cache covariate distances. Returns same type.
# ============================================================

.bind.kernel_bregman = function(model, Z, ...) {
  Z = atleast_2d(Z)
  # Already bound to same data? Return as-is.
  if (!is.null(model$.bind_method) && !is.null(model$.n_bound) &&
      model$.n_bound == nrow(Z) && !is.null(model$.bind_cov)) {
    kern = model$kern
    if (inherits(kern, "product_kernel")) {
      iw_cov = kern$iw[kern$iw > 1] - 1
      niw_cov = setdiff(1:ncol(Z), iw_cov)
      Z_cov = Z[, niw_cov, drop = FALSE]
    } else {
      Z_cov = Z
    }
    if (ncol(Z_cov) == ncol(model$.bind_cov) && max(abs(Z_cov - model$.bind_cov)) < 1e-12)
      return(model)
  }
  Z_centers = model$Z
  u_centers = Z_centers[, 1]
  cov_centers = Z_centers[, -1, drop = FALSE]
  kern = model$kern

  if (inherits(kern, "product_kernel")) {
    iw_cov = kern$iw[kern$iw > 1] - 1
    niw_cov = setdiff(1:ncol(cov_centers), iw_cov)
    Z_cov = Z[, niw_cov, drop = FALSE]
    K_x = kernel_matrix(Z_cov, cov_centers[, niw_cov, drop = FALSE], kern$k_x)
    model$.K_x_cache = K_x
    model$.cov_centers = cov_centers
    model$.iw_cov = iw_cov
    model$.bind_cov = Z_cov
    model$.bind_method = "product"
    # Precompute center groups keyed by (time, iw) for block-diagonal .phi.
    # Each group's K_x columns are extracted once here; .phi only does
    # cheap row matching at call time.
    key_ctr = as.character(u_centers)
    for (j in iw_cov)
      key_ctr = paste0(key_ctr, "|", cov_centers[, j])
    ctr_idx = split(seq_along(key_ctr), key_ctr)
    model$.ctr_blocks = lapply(ctr_idx, function(cols)
      list(cols = cols, K_x = K_x[, cols, drop = FALSE]))
  } else {
    model$.D2_cov_cache = .fast_sqdist(Z, cov_centers)
    model$.bind_cov = Z
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
  if (!is.null(x$.bind_cov))
    x$.bind_cov = x$.bind_cov[i, , drop = FALSE]
  if (!is.null(x$.n_bound))
    x$.n_bound = sum(rep(TRUE, x$.n_bound)[i])
  if (!is.null(x$.ctr_blocks))
    x$.ctr_blocks = lapply(x$.ctr_blocks, function(b)
      list(cols = b$cols, K_x = b$K_x[i, , drop = FALSE]))
  x
}


# ============================================================
# phi: linear predictor. Z = (u, covariates), time in col 1.
# ============================================================

.phi.kernel_bregman = function(model, Z, ...) {
  Z = atleast_2d(Z)
  kern = model$kern

  if (!is.null(model$.bind_method) && model$.bind_method == "isotropic") {
    if (nrow(Z) != model$.n_bound || ncol(Z) != ncol(model$.bind_cov) ||
        max(abs(Z - model$.bind_cov)) > 1e-12)
      stop(".phi: Z covariates do not match bound covariates")
    delta_u2 = outer(u(Z), model$.u_centers, function(a, b) (a - b)^2)
    D2 = model$.D2_cov_cache + delta_u2
    D2[D2 < 0] = 0
    K_eval = .kernel_from_D2(D2, kern)
  } else if (!is.null(model$.bind_method) && model$.bind_method == "product") {
    cov = Z[, -1, drop = FALSE]
    Z_cov = cov[, setdiff(1:ncol(cov), model$.iw_cov), drop = FALSE]
    if (nrow(Z_cov) != model$.n_bound || ncol(Z_cov) != ncol(model$.bind_cov) ||
        max(abs(Z_cov - model$.bind_cov)) > 1e-12)
      stop(".phi: Z covariates do not match bound covariates")
    # Block-diagonal multiply: match new-point (time, iw) to precomputed
    # center blocks. Each block is a dense K_x submatrix; rest is zero.
    key_new = as.character(u(Z))
    for (j in model$.iw_cov)
      key_new = paste0(key_new, "|", cov[, j])
    n_new = nrow(Z)
    blocks = list()
    for (g in names(model$.ctr_blocks)) {
      rows = which(key_new == g)
      if (length(rows) > 0) {
        b = model$.ctr_blocks[[g]]
        blocks[[length(blocks) + 1L]] =
          list(rows = rows, cols = b$cols, K = b$K_x[rows, , drop = FALSE])
      }
    }
    K_bd = block_diag(blocks, n_new)
    p_kern = K_bd %*% model$alpha
  } else {
    K_eval = kernel_matrix(Z, model$Z, kern)
  }

  p = if (exists("p_kern")) p_kern else as.vector(K_eval %*% model$alpha)
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
