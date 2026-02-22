## Unified Bregman regression solver.
##
## Minimizes: sum w_i * chis(phi_i) - sum w_i * Y_i * phi_i + eta * alpha' K alpha
## where phi = K alpha + beta0 (when intercept = TRUE).
##
## The dispersion provides chis/dchis/ddchis. The solver calls dchis (for gradient)
## and ddchis (for Hessian/IRLS curvature). For linear dispersions (quadratic),
## it reduces to a single solve.
##
## This unifies the current kernel_bregman (outcome), kernel_bregman_gamma (sign-flip
## weights), and balancing_weights (quadratic/seminorm).

# ============================================================
# Main solver
# ============================================================

#' Fit kernel Bregman regression with given dispersion.
#'
#' Solves: min_phi sum w * chis(phi) - sum w * Y * phi + eta * ||f||_K^2
#' where phi_i = beta0 + sum_j alpha_j K(z_i, z_j).
#'
#' For linear dispersions (quadratic): single direct solve.
#' For nonlinear dispersions: IRLS (Newton on the augmented system).
#' First iteration uses direct Cholesky; subsequent use CG with that
#' Cholesky as preconditioner.
#'
#' @param Y Target vector (n). For outcome regression: observed Y.
#'   For weight regression: dot_psi targets.
#' @param Z Covariate matrix (n x p).
#' @param kern Kernel object.
#' @param eta Regularization parameter (scalar). Penalty is eta * alpha'K*alpha.
#' @param dispersion Dispersion object (from dispersions.R).
#' @param intercept Include unpenalized intercept beta0?
#' @param w Observation weights (n). Default: all ones.
#' @param K Precomputed kernel matrix (optional).
#' @param alpha0 Initial alpha (optional, for warm start).
#' @param beta0_init Initial beta0 (optional).
#' @param maxiter Max IRLS iterations (ignored for linear dispersions).
#' @param tol Convergence tolerance on max |delta alpha|.
#' @param solver "auto", "direct", or "cg". Auto picks CG for n > 500.
#' @return A kernel_bregman object with alpha, beta0, K, R_chol, dispersion, etc.
kernel_bregman = function(Y, Z, kern, eta, dispersion, intercept = TRUE,
                          w = NULL, K = NULL, alpha0 = NULL, beta0_init = NULL,
                          maxiter = 25, tol = 1e-6, solver = "auto") {
  Z = atleast_2d(Z)
  if (is.null(K)) K = kernel_matrix(Z, Z, kern)
  n = nrow(K)
  if (is.null(w)) w = rep(1, n)
  alpha = if (!is.null(alpha0)) alpha0 else rep(0, n)
  beta0 = if (!is.null(beta0_init)) beta0_init else 0

  if (solver == "auto") solver = if (n > 500) "cg" else "direct"

  # Linear dispersion: single direct solve
  if (isTRUE(dispersion$linear)) {
    return(.bregman_linear(Y, Z, kern, eta, dispersion, intercept, w, K, alpha, beta0))
  }

  # Nonlinear: IRLS
  R_chol = NULL
  for (iter in 1:maxiter) {
    phi = as.vector(K %*% alpha) + beta0
    g = dispersion$dchis(phi)        # primal prediction
    wc = w * dispersion$ddchis(phi)  # curvature weight
    # Floor curvature to avoid singular Hessian (e.g. pos_quadratic at phi=0)
    wc = pmax(wc, 1e-10 * max(abs(wc), 1))
    wr = w * (g - Y)                 # gradient of data term w.r.t. phi

    if (solver == "direct") {
      # IRLS working response: linearize around current phi
      rhs_a = wc * phi - wr  # = wc * phi + w * (Y - g)
      if (intercept) {
        rhs_b = sum(rhs_a)
        H_aa = sweep(K, 1, wc, "*") + 2 * eta * diag(n)
        M_aug = rbind(cbind(H_aa, wc), c(wc, sum(wc)))
        sol = solve(M_aug, c(rhs_a, rhs_b))
        alpha_new = sol[1:n]; beta0_new = sol[n + 1]
      } else {
        H = sweep(K, 1, wc, "*") + 2 * eta * diag(n)
        alpha_new = as.vector(solve(H, rhs_a)); beta0_new = 0
      }
    } else {
      # Newton with CG on symmetric Hessian
      Kwc = as.vector(K %*% wc)
      s = sum(wc)
      g_a = as.vector(K %*% wr) + 2 * eta * as.vector(K %*% alpha)

      if (intercept) {
        g_b = sum(wr)
        neg_grad = -c(g_a, g_b)

        matvec = function(v) {
          va = v[1:n]; vb = v[n + 1]
          Kva = as.vector(K %*% va)
          c(as.vector(K %*% (wc * Kva)) + 2 * eta * Kva + vb * Kwc,
            sum(Kwc * va) + s * vb)
        }

        if (is.null(R_chol)) {
          B = sweep(K, 1, sqrt(pmax(wc, 0)), "*")
          H_aa = crossprod(B) + 2 * eta * K
          H_aug = rbind(cbind(H_aa, Kwc), c(Kwc, s))
          R_chol = chol(H_aug)
          step = backsolve(R_chol,
                           forwardsolve(R_chol, neg_grad,
                                        upper.tri = TRUE, transpose = TRUE))
        } else {
          precond_fn = function(v)
            backsolve(R_chol,
                      forwardsolve(R_chol, v,
                                   upper.tri = TRUE, transpose = TRUE))
          cg_res = .cg_solve(matvec, neg_grad, precond = precond_fn,
                             tol = tol * 0.1, maxiter = 50)
          step = cg_res$x
        }
        alpha_new = alpha + step[1:n]
        beta0_new = beta0 + step[n + 1]
      } else {
        neg_grad = -g_a

        matvec = function(v) {
          Kv = as.vector(K %*% v)
          as.vector(K %*% (wc * Kv)) + 2 * eta * Kv
        }

        if (is.null(R_chol)) {
          B = sweep(K, 1, sqrt(pmax(wc, 0)), "*")
          H = crossprod(B) + 2 * eta * K
          R_chol = chol(H)
          step = backsolve(R_chol,
                           forwardsolve(R_chol, neg_grad,
                                        upper.tri = TRUE, transpose = TRUE))
        } else {
          precond_fn = function(v)
            backsolve(R_chol,
                      forwardsolve(R_chol, v,
                                   upper.tri = TRUE, transpose = TRUE))
          cg_res = .cg_solve(matvec, neg_grad, precond = precond_fn,
                             tol = tol * 0.1, maxiter = 50)
          step = cg_res$x
        }
        alpha_new = alpha + step; beta0_new = 0
      }
    }

    delta = max(abs(alpha_new - alpha), abs(beta0_new - beta0))
    alpha = alpha_new; beta0 = beta0_new
    if (delta < tol) break
  }

  # Hessian Cholesky at convergence (for one-step bootstrap)
  phi = as.vector(K %*% alpha) + beta0
  wc_final = dispersion$ddchis(phi)
  B = sweep(K, 1, sqrt(pmax(wc_final, 0)), "*")
  if (intercept) {
    Kwc = as.vector(K %*% wc_final)
    H_aug = rbind(cbind(crossprod(B) + 2 * eta * K, Kwc), c(Kwc, sum(wc_final)))
    R_chol = chol(H_aug)
  } else {
    R_chol = chol(crossprod(B) + 2 * eta * K)
  }

  structure(list(alpha = alpha, beta0 = beta0, intercept = intercept,
                 kern = kern, Z = Z, K = K, eta = eta, R_chol = R_chol,
                 dispersion = dispersion, iter = iter, w = w),
            class = "kernel_bregman")
}

# ============================================================
# Linear dispersion: single direct solve
# ============================================================

# For quadratic dispersion: dchis(phi) = phi, ddchis = 1.
# Loss = sum w * phi^2/2 - sum w * Y * phi + eta * alpha'K*alpha.
# FOC: (diag(w) K + 2 eta I) alpha + w * beta0 = w * Y.
.bregman_linear = function(Y, Z, kern, eta, dispersion, intercept, w, K, alpha, beta0) {
  n = nrow(K)

  if (intercept) {
    H = sweep(K, 1, w, "*") + 2 * eta * diag(n)
    rhs = w * Y
    M_aug = rbind(cbind(H, w), c(w, sum(w)))
    sol = solve(M_aug, c(rhs, sum(rhs)))
    alpha = sol[1:n]; beta0 = sol[n + 1]
    R_chol = chol(M_aug)
  } else {
    H = sweep(K, 1, w, "*") + 2 * eta * diag(n)
    alpha = as.vector(solve(H, w * Y))
    beta0 = 0
    R_chol = chol(H)
  }

  structure(list(alpha = alpha, beta0 = beta0, intercept = intercept,
                 kern = kern, Z = Z, K = K, eta = eta, R_chol = R_chol,
                 dispersion = dispersion, iter = 1, w = w),
            class = "kernel_bregman")
}

# ============================================================
# One-step Newton update (bootstrap)
# ============================================================

#' One-step Newton update: new weights and/or targets.
#' Uses precomputed Cholesky of Hessian from original fit.
onestep_bregman = function(model, Y_new, w_new) {
  K = model$K; n = nrow(K)
  alpha = model$alpha; beta0 = model$beta0
  eta_reg = model$eta; R = model$R_chol
  disp = model$dispersion

  phi = as.vector(K %*% alpha) + beta0
  g = disp$dchis(phi)
  resid = w_new * (g - Y_new)
  g_a = as.vector(K %*% resid) + 2 * eta_reg * as.vector(K %*% alpha)

  if (model$intercept) {
    grad = c(g_a, sum(resid))
    step = backsolve(R, forwardsolve(R, -grad, upper.tri = TRUE, transpose = TRUE))
    list(alpha = alpha + step[1:n], beta0 = beta0 + step[n + 1])
  } else {
    step = backsolve(R, forwardsolve(R, -g_a, upper.tri = TRUE, transpose = TRUE))
    list(alpha = alpha + step, beta0 = 0)
  }
}

# ============================================================
# Predict
# ============================================================

#' Predict from kernel_bregman: returns dchis(phi) (primal prediction).
predict.kernel_bregman = function(obj, newZ, ...) {
  newZ = atleast_2d(newZ)
  Knew = kernel_matrix(newZ, obj$Z, obj$kern)
  phi = as.vector(Knew %*% obj$alpha) + obj$beta0
  obj$dispersion$dchis(phi)
}

#' Predict linear predictor phi (dual), not transformed.
predict_phi = function(obj, newZ) {
  newZ = atleast_2d(newZ)
  Knew = kernel_matrix(newZ, obj$Z, obj$kern)
  as.vector(Knew %*% obj$alpha) + obj$beta0
}

# ============================================================
# Subject binding: precompute covariate distances for fast
# time-varying evaluation of kernel expansions
# ============================================================

#' Precompute covariate distances between evaluation subjects and model centers.
#' Model centers have form (u, z_1, ..., z_p). Covariate columns z are fixed
#' across time evaluations; only the time column u changes.
bind_subjects = function(model, Z_eval) {
  Z_eval = atleast_2d(Z_eval)
  Z_centers = model$Z
  u_centers = Z_centers[, 1]
  cov_centers = Z_centers[, -1, drop = FALSE]

  D2_cov = .fast_sqdist(Z_eval, cov_centers)

  list(
    D2_cov = D2_cov,
    u_centers = u_centers,
    alpha = model$alpha,
    beta0 = if (!is.null(model$beta0)) model$beta0 else 0,
    kern = model$kern,
    dispersion = model$dispersion
  )
}

#' Evaluate phi at scalar time u for all bound subjects (or a subset).
eval_phi = function(bound, u, idx = NULL) {
  D2_cov = if (is.null(idx)) bound$D2_cov else bound$D2_cov[idx, , drop = FALSE]
  delta_u2 = (u - bound$u_centers)^2
  D2 = sweep(D2_cov, 2, delta_u2, "+")
  D2[D2 < 0] = 0
  K = .kernel_from_D2(D2, bound$kern)
  as.vector(K %*% bound$alpha) + bound$beta0
}

#' Evaluate phi at per-subject times (u_vec[j] for subject idx[j]).
eval_phi_vec = function(bound, u_vec, idx = NULL) {
  D2_cov = if (is.null(idx)) bound$D2_cov else bound$D2_cov[idx, , drop = FALSE]
  delta_u2 = outer(u_vec, bound$u_centers, function(a, b) (a - b)^2)
  D2 = D2_cov + delta_u2
  D2[D2 < 0] = 0
  K = .kernel_from_D2(D2, bound$kern)
  as.vector(K %*% bound$alpha) + bound$beta0
}
