## Kernel Bregman solver.
##
## Solves the penalized dual problem (eq:dual-problem from the paper):
##
##   min_f  P_n[chi*_Z(f(Z))] + (eta/2) ||f||_rho^2 - dpsiZ(f)
##
## Lambda (hazard) and gamma (weights) are both instances — they differ
## only in dpsiZ. The solver doesn't know which.
##
## In the representer basis phi = K alpha + B beta, where B = null_basis(Z, kern)
## spans ker(rho), the sum-scaled objective is:
##
##   min  sum w_i chi*(phi_i) + (n*eta/2) alpha'K alpha - [target'alpha + target_null'beta]
##
## where target and target_null are the representer-coefficient projections
## of the linear functional dpsiZ:
##
##   target (n-vector):     d/d(alpha) dpsiZ(phi)
##   target_null (d-vector): d/d(beta) dpsiZ(phi)
##
## Newton iteration uses the two-step factorization H = blkdiag(K, I_d) M
## where M = [diag(wc)K + n*eta*I, diag(wc)B; B'diag(wc)K, B'diag(wc)B].
## M is not symmetric but has condition number kappa(K), not kappa(K)^2,
## so the solve stays well-conditioned. M_inv is stored for O(n^2) one-step
## bootstrap solves.

# ============================================================
# Cholesky back-solve helper
# ============================================================

.chol_solve = function(R, b) {
  backsolve(R, forwardsolve(R, b, upper.tri = TRUE, transpose = TRUE))
}

# ============================================================
# dpsiZ helpers
# ============================================================

#' Compute representer-coefficient projections of dpsiZ from an
#' observation-space response Y.
#'
#' When dpsiZ(f) = P_n[Y f(Z)] (e.g. the hazard model's counting-process
#' increment), the representer projections are:
#'   target      = K(w * Y)     (n-vector)
#'   target_null = B'(w * Y)    (d-vector)
#'
#' @param Y Response vector (n).
#' @param K Kernel matrix (n x n), with nugget already added.
#' @param B Null-space basis at training points (n x d).
#' @param w Observation weights (n). Default: all ones.
#' @return List with $target (n-vector) and $target_null (d-vector).
dpsi_from_response = function(Y, K, B, w = NULL) {
  if (is.null(w)) w = rep(1, length(Y))
  wY = w * Y
  list(target = as.vector(K %*% wY),
       target_null = if (ncol(B) > 0) as.vector(crossprod(B, wY)) else numeric(0))
}

# ============================================================
# Main solver
# ============================================================

#' Solve the penalized dual problem.
#'
#' @param Z Covariates (n x p).
#' @param kern Kernel object (carries null_basis for ker(rho)).
#' @param eta Regularization strength (paper scale). Penalty = (eta/2)||f||^2.
#'   The solver handles the n-scaling internally: ridge term = n*eta*I.
#' @param dispersion Dispersion object (chis/dchis/ddchis).
#' @param target n-vector: representer projection of dpsiZ onto the RKHS
#'   component. For the hazard model with response Y: target = K(w*Y)
#'   (use dpsi_from_response). For the weight model: target = kernel
#'   mean embedding of dpsiZ (from gam_embedding).
#' @param target_null d-vector: projection of dpsiZ onto ker(rho).
#'   Default: zero (no null-space component in the linear functional).
#' @param w Observation weights (n). Default: all ones.
#' @param K Precomputed kernel matrix (optional; nugget already added).
#' @param alpha0 Warm start for alpha (n).
#' @param beta0 Warm start for beta (d-vector or scalar recycled to d).
#' @param maxiter Max Newton iterations.
#' @param tol Convergence tolerance on max |step|.
kernel_bregman = function(Z, kern, eta, dispersion,
                          target,
                          target_null = NULL,
                          w = NULL, K = NULL,
                          alpha0 = NULL, beta0 = NULL,
                          maxiter = 25, tol = 1e-6) {
  Z = atleast_2d(Z)
  n = nrow(Z)

  if (is.null(K)) K = kernel_matrix(Z, Z, kern)
  # Numerical regularization for chol(K) in the two-step factorization.
  # g_a is always in im(K), so this doesn't affect the solution.
  K = K + (1e-8 * sum(diag(K)) / n) * diag(n)
  if (is.null(w)) w = rep(1, n)

  # Null-space basis from the kernel's seminorm
  B = null_basis(Z, kern)
  d = ncol(B)

  if (is.null(target_null)) target_null = rep(0, d)

  alpha = if (!is.null(alpha0)) alpha0 else rep(0, n)
  beta = if (d > 0) {
    if (!is.null(beta0)) rep_len(as.numeric(beta0), d) else rep(0, d)
  } else numeric(0)

  K_chol = chol(K)

  # Linear dispersion (e.g. quadratic): one Newton step is exact.
  if (isTRUE(dispersion$linear)) maxiter = 1

  for (iter in 1:maxiter) {
    phi = as.vector(K %*% alpha)
    if (d > 0) phi = phi + as.vector(B %*% beta)

    mu = dispersion$dchis(phi)
    wc = w * dispersion$ddchis(phi)
    wc = pmax(wc, 1e-10 * max(abs(wc), 1))

    # --- Gradient ---
    Ka = as.vector(K %*% alpha)
    g_a = as.vector(K %*% (w * mu)) + n * eta * Ka - target

    # --- Two-step Newton solve ---
    #
    # The full Hessian factors as H = blkdiag(K, I_d) M where
    #   M = [diag(wc) K + n*eta*I,   diag(wc) B]
    #       [B' diag(wc) K,          B' diag(wc) B]
    #
    # M is NOT symmetric but has condition number kappa(K), not kappa(K)^2.
    # Solve: delta = M^{-1} blkdiag(K^{-1}, I_d) (-g).
    v_a = .chol_solve(K_chol, -g_a)

    wcK = sweep(K, 1, wc, "*")            # diag(wc) K
    M_aa = wcK + n * eta * diag(n)

    if (d > 0) {
      g_b = as.vector(crossprod(B, w * mu)) - target_null

      wcB = wc * B                         # diag(wc) B, n x d
      M_ba = crossprod(B, wcK)             # d x n
      M_bb = crossprod(B, wcB)             # d x d

      M = rbind(cbind(M_aa, wcB), cbind(M_ba, M_bb))
      delta = solve(M, c(v_a, -g_b))
      d_a = delta[1:n]
      d_b = delta[(n + 1):(n + d)]
    } else {
      d_a = as.vector(solve(M_aa, v_a))
      d_b = numeric(0)
    }

    step_size = max(abs(d_a), if (d > 0) max(abs(d_b)) else 0)
    alpha = alpha + d_a
    beta = beta + d_b
    if (step_size < tol) break
  }

  # --- Pre-invert M at convergence (for one-step bootstrap) ---
  #
  # The one-step bootstrap solves H delta = -g_new using the original Hessian.
  # Store K_chol and M_inv so each bootstrap rep is O(n^2) multiply, not O(n^3) solve.
  phi_f = as.vector(K %*% alpha)
  if (d > 0) phi_f = phi_f + as.vector(B %*% beta)
  wc_f = w * dispersion$ddchis(phi_f)
  wc_f = pmax(wc_f, 1e-10 * max(abs(wc_f), 1))

  wcK_f = sweep(K, 1, wc_f, "*")
  M_aa_f = wcK_f + n * eta * diag(n)
  if (d > 0) {
    wcB_f = wc_f * B
    M_f = rbind(cbind(M_aa_f, wcB_f),
                cbind(crossprod(B, wcK_f), crossprod(B, wcB_f)))
  } else {
    M_f = M_aa_f
  }
  M_inv = solve(M_f)

  structure(list(
    alpha = alpha, beta = beta,
    kern = kern, Z = Z, K = K, B = B, eta = eta,
    K_chol = K_chol, M_inv = M_inv,
    dispersion = dispersion, iter = iter,
    w = w, target = target, target_null = target_null
  ), class = "kernel_bregman")
}

# ============================================================
# One-step Newton update (bootstrap)
# ============================================================

#' One-step Newton: compute gradient at current (alpha, beta) with new
#' targets/weights, then solve using the stored Hessian factorization.
#'
#' @param model A fitted kernel_bregman object.
#' @param target_new New target vector (n). Required.
#' @param target_null_new New target_null vector (d). Default: model's.
#' @param w_new New observation weights (n). Default: model's.
#' @return List with $alpha (n), $beta (d).
onestep_bregman = function(model, target_new, target_null_new = NULL, w_new = NULL) {
  K = model$K; n = nrow(K)
  B = model$B; d = ncol(B)
  alpha = model$alpha
  beta = if (d > 0) model$beta else numeric(0)
  eta = model$eta
  disp = model$dispersion

  if (is.null(w_new)) w_new = model$w
  if (is.null(target_null_new)) target_null_new = model$target_null

  phi = as.vector(K %*% alpha)
  if (d > 0) phi = phi + as.vector(B %*% beta)
  mu = disp$dchis(phi)
  Ka = as.vector(K %*% alpha)

  # Gradient with new targets/weights
  g_a = as.vector(K %*% (w_new * mu)) + n * eta * Ka - target_new

  # Two-step: delta = M_inv %*% blkdiag(K^{-1}, I_d) (-g)
  v_a = .chol_solve(model$K_chol, -g_a)

  if (d > 0) {
    g_b = as.vector(crossprod(B, w_new * mu)) - target_null_new
    rhs = c(v_a, -g_b)
  } else {
    rhs = v_a
  }

  delta = as.vector(model$M_inv %*% rhs)
  d_a = delta[1:n]
  d_b = if (d > 0) delta[(n + 1):(n + d)] else numeric(0)

  list(alpha = alpha + d_a, beta = beta + d_b)
}

# ============================================================
# Predict
# ============================================================

#' Predict from kernel_bregman: returns dchis(phi) (primal prediction).
predict.kernel_bregman = function(obj, newZ, ...) {
  obj$dispersion$dchis(predict_phi(obj, newZ))
}

#' Predict linear predictor phi (dual), not transformed.
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

#' Predict phi with new alpha/beta, reusing the model's kernel centers.
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
# LOO CV loss (linearized)
# ============================================================

#' Leave-one-out CV loss for a fitted kernel_bregman model.
#'
#' Works from the fit object alone — no observation-space Y needed.
#' At convergence the FOC gives K(w*mu + n*eta*alpha) = target, so the
#' IRLS working residual is n*eta*alpha / wc. This applies to both
#' lambda (hazard) and gamma (weight) models identically.
#'
#' For quadratic dispersion: exact LOO.
#' For other dispersions: approximate LOO via linearization at the
#' converged solution (same approximation as glm).
#'
#' @param fit A kernel_bregman object.
#' @return Scalar LOO CV loss.
loo_loss = function(fit) {
  K = fit$K
  n = nrow(K)
  phi = as.vector(K %*% fit$alpha)
  d = length(fit$beta)
  if (d > 0) phi = phi + as.vector(fit$B %*% fit$beta)
  w = if (!is.null(fit$w)) fit$w else rep(1, n)

  # Curvature weights at convergence
  wc = w * fit$dispersion$ddchis(phi)
  wc = pmax(wc, 1e-10 * max(abs(wc), 1))

  # Working residual from the FOC. At convergence,
  #   K(w*mu + n*eta*alpha) = target
  # so the IRLS working residual is n*eta*alpha / wc.
  # For lambda with target = K*(w*Y), this equals (Y - mu) / ddchis(phi).
  raw_resid = n * fit$eta * fit$alpha / wc

  # Hat matrix diagonal: H = K (Wc K + n*eta*I)^{-1} Wc
  sqw = sqrt(wc)
  Kt = sweep(sweep(K, 1, sqw, "*"), 2, sqw, "*")
  evd = eigen(Kt, symmetric = TRUE)
  V = evd$vectors
  l = pmax(evd$values, 0)

  shrink = l / (l + n * fit$eta)
  H_diag = rowSums(V^2 * rep(shrink, each = n))

  # LOO residual
  loo_resid = raw_resid / pmax(1 - H_diag, 1e-10)

  # Weighted LOO MSE
  mean(wc * loo_resid^2)
}

#' Factory for LOO scoring function. Returns a closure
#' score_fn(eta_grid) -> numeric(m) of LOO losses.
#'
#' Used to decouple CV tuning from the solver: the pipeline constructs
#' the closure (which knows about K, targets, etc.), and tuning just picks argmin.
#' Works for both lambda (hazard) and gamma (weight) models — the LOO loss
#' is computed from the solver's FOC, not from a response Y.
#'
#' @param Z Covariates (n x p).
#' @param kern Kernel object.
#' @param dispersion Dispersion object.
#' @param target Representer-space target (n).
#' @param target_null Null-space target (d). Default: zero.
#' @param w Observation weights. Default: all ones.
#' @return Function(eta_grid) -> numeric(length(eta_grid)).
loo_score_fn = function(Z, kern, dispersion, target,
                        target_null = NULL, w = NULL) {
  Z = atleast_2d(Z)
  n = nrow(Z)
  K = kernel_matrix(Z, Z, kern)
  B = null_basis(Z, kern)
  if (is.null(w)) w = rep(1, n)
  if (is.null(target_null)) target_null = rep(0, ncol(B))

  function(eta_grid) {
    vapply(eta_grid, function(eta) {
      fit = tryCatch(
        kernel_bregman(Z, kern, eta, dispersion,
                       target = target, target_null = target_null,
                       w = w, K = K),
        error = function(e) NULL)
      if (!is.null(fit)) loo_loss(fit) else Inf
    }, numeric(1))
  }
}

# ============================================================
# Subject binding: precompute covariate distances for fast
# time-varying evaluation of kernel expansions
# ============================================================

#' Precompute covariate distances between evaluation subjects and model centers.
#' Model centers have form (u, z_1, ..., z_p). Covariate columns z are fixed
#' across time evaluations; only the time column u changes.
#'
#' For isotropic kernels (matern, gaussian): precomputes D2_cov for fast
#' eval_phi via D2 = D2_cov + delta_u^2.
#' For product kernels: precomputes K_x (covariate kernel) and stores
#' center metadata for mask construction at eval time.
#' For other kernels: stores Z_eval and Z_centers for full kernel_matrix fallback.
bind_subjects = function(model, Z_eval) {
  Z_eval = atleast_2d(Z_eval)
  Z_centers = model$Z
  u_centers = Z_centers[, 1]
  cov_centers = Z_centers[, -1, drop = FALSE]
  kern = model$kern

  bound = list(
    u_centers = u_centers,
    alpha = model$alpha,
    beta = model$beta,
    kern = kern,
    dispersion = model$dispersion,
    Z_eval = Z_eval
  )

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

# Null-space contribution to phi for bound objects.
# For d=1 (intercept): beta[1] (scalar broadcast, avoids constructing Z_pred).
# For d>1: evaluates null_basis at prediction points.
.null_phi = function(bound, Z_pred) {
  d = length(bound$beta)
  if (d == 0) return(0)
  if (d == 1) return(bound$beta[1])
  as.vector(null_basis(Z_pred, bound$kern) %*% bound$beta)
}

#' Evaluate phi at scalar time u for all bound subjects (or a subset).
eval_phi = function(bound, u, idx = NULL) {
  if (bound$method == "isotropic") {
    D2_cov = if (is.null(idx)) bound$D2_cov else bound$D2_cov[idx, , drop = FALSE]
    delta_u2 = (u - bound$u_centers)^2
    D2 = sweep(D2_cov, 2, delta_u2, "+")
    D2[D2 < 0] = 0
    K = .kernel_from_D2(D2, bound$kern)
  } else if (bound$method == "product") {
    K_x = if (is.null(idx)) bound$K_x else bound$K_x[idx, , drop = FALSE]
    u_mask = outer(rep(u, if (is.null(idx)) nrow(K_x) else length(idx)),
                   bound$u_centers, "==")
    mask = u_mask
    Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
    for (j in bound$iw_cov)
      mask = mask & outer(Z_ev[, j], bound$cov_centers[, j], "==")
    K = K_x * mask
  } else {
    Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
    Z_aug = cbind(u, Z_ev)
    K = kernel_matrix(Z_aug, cbind(bound$u_centers, bound$Z_eval), bound$kern)
  }
  phi = as.vector(K %*% bound$alpha)
  Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
  phi + .null_phi(bound, cbind(u, Z_ev))
}

#' Evaluate phi at per-subject times (u_vec[j] for subject idx[j]).
eval_phi_vec = function(bound, u_vec, idx = NULL) {
  if (bound$method == "isotropic") {
    D2_cov = if (is.null(idx)) bound$D2_cov else bound$D2_cov[idx, , drop = FALSE]
    delta_u2 = outer(u_vec, bound$u_centers, function(a, b) (a - b)^2)
    D2 = D2_cov + delta_u2
    D2[D2 < 0] = 0
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
    Z_aug = cbind(u_vec, Z_ev)
    Z_cen = cbind(bound$u_centers, bound$Z_eval)
    K = kernel_matrix(Z_aug, Z_cen, bound$kern)
  }
  phi = as.vector(K %*% bound$alpha)
  Z_ev = if (is.null(idx)) bound$Z_eval else bound$Z_eval[idx, , drop = FALSE]
  phi + .null_phi(bound, cbind(u_vec, Z_ev))
}
