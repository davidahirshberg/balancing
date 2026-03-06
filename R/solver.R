#' # Kernel Bregman Solver
#'
#' ## `kernel_bregman`
#'
#' **Paper**: solve the regularized dual problem (eq:dual in spinoff3)
#' $$\hat\phi = \arg\min_{\rho(\phi)<\infty} \hat P\,\chi^*_Z(\phi)
#'   + \eta\,\rho(\phi)^2 - \dot\psi_Z(\phi)$$
#' The weights are $\hat\gamma = \dot\chi^*_Z(\hat\phi)$.
#' For a seminorm $\rho$ with null space $\ker(\rho) = \mathrm{span}(B)$,
#' the representer theorem gives $\hat\phi = K\alpha + B\beta$.
#'
#' **Code**: with $\phi = K\alpha + B\beta$, the finite-dimensional
#' objective over $\theta = (\alpha, \beta)$ is
#' $$L(\theta) = \sum_i w_i\,\chi^*(\phi_i)
#'   + \tfrac{n\eta}{2}\,\alpha^\top K\alpha
#'   - c^\top\alpha - c_0^\top\beta$$
#' where $n$ is the sample size (passed by the caller; defaults to
#' `nrow(Z)` which is correct when basis = sample), $c = $ `target`
#' $= K r$ (from `project_target`) and $c_0 = $ `target_null`
#' $= B^\top r$.
#'
#' The gradient is
#' $$\nabla_\alpha L = K(w \odot \dot\chi^*(\phi))
#'   + n\eta\,K\alpha - c, \quad
#'   \nabla_\beta L = B^\top(w \odot \dot\chi^*(\phi)) - c_0$$
#'
#' The Hessian-vector product $H v$ uses
#' $w_c = w \odot \ddot\chi^*(\phi)$ (curvature weights):
#' $$H \begin{pmatrix} v_\alpha \\ v_\beta \end{pmatrix} =
#'   \begin{pmatrix} K(w_c \odot (K v_\alpha + B v_\beta)) + n\eta\,K v_\alpha \\
#'   B^\top(w_c \odot (K v_\alpha + B v_\beta)) \end{pmatrix}$$
#'
#' Solved by trust-region Newton (Steihaug-Toint CG). Convergence
#' criterion: Bregman divergence $D_{\chi^*}(\phi^{(t)}, \phi^{(t-1)}) < $ `tol`.
#' For quadratic $\chi^*$ the Hessian is constant and the trust
#' region is infinite (`delta0 = Inf`), so one Newton step is exact.
#'
#' Stores a lazy Hessian inverse `get_M_lu` for `onestep_bregman`.
kernel_bregman = function(Z, kern, eta, dispersion,
                          target, target_null = NULL,
                          w = NULL, K = NULL, K_chol = NULL,
                          alpha0 = NULL, beta0 = NULL,
                          maxiter = 25, tol = 1e-3,
                          n = NULL,
                          log = null_logger) {
  Z = atleast_2d(Z)
  n_row = nrow(Z)
  if (is.null(n)) n = n_row
  if (is.null(K)) K = kernel_matrix(Z, Z, kern)
  if (is.null(K_chol)) {
    K_chol = pchol(K)
  } else if (!inherits(K_chol, "pchol")) {
    # Legacy: caller passed a plain Cholesky factor. Wrap it.
    K_chol = structure(list(R = K_chol, pivot = seq_len(n_row),
                            rank = n_row, n = n_row), class = "pchol")
  }
  if (is.null(w)) w = rep(1, n_row)

  B = null_basis(Z, kern)
  d = ncol(B)
  if (is.null(target_null)) target_null = rep(0, d)

  alpha = if (!is.null(alpha0)) alpha0 else rep(0, n_row)
  beta = if (d > 0) {
    if (!is.null(beta0)) rep_len(as.numeric(beta0), d) else rep(0, d)
  } else numeric(0)

  #' Lambda bindings: dchis(phi) reads like dot chi*(phi).
  disp = dispersion
  dchis  = function(phi) .dchis(disp, phi)
  ddchis = function(phi) .ddchis(disp, phi)
  chis   = function(phi) .chis(disp, phi)

  fg_kernel = function(theta) {
    a = theta[1:n_row]
    b = if (d > 0) theta[(n_row+1):(n_row+d)] else numeric(0)
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
    a = theta[1:n_row]
    b = if (d > 0) theta[(n_row+1):(n_row+d)] else numeric(0)
    phi = as.vector(K %*% a)
    if (d > 0) phi = phi + as.vector(B %*% b)
    wc = w * ddchis(phi)
    wc = pmax(wc, 1e-10 * max(abs(wc), 1))
    v_a = v[1:n_row]
    v_b = if (d > 0) v[(n_row+1):(n_row+d)] else numeric(0)
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
      a = th[1:n_row]; b = if (d > 0) th[(n_row+1):(n_row+d)] else numeric(0)
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
    a = th[1:n_row]; b = if (d > 0) th[(n_row+1):(n_row+d)] else numeric(0)
    phi = as.vector(K %*% a)
    if (d > 0) phi = phi + as.vector(B %*% b)
    phi
  }
  weight_fn = function(th) dchis(phi_from_theta(th))

  # Trust region in dual geometry: ||ddchi*(phi) · (K·da + B·db)||.
  # The Hessian of chi* is the metric tensor of the dual; its action
  # on phi-steps gives the induced primal (gamma) change.  Permissive
  # where gamma is insensitive (tails), tight where it's sensitive.
  # For quadratic chi* (ddchi* = const), this is proportional to the
  # phi-space norm, and delta0 = Inf makes it a no-op anyway.
  dual_map = function(theta, p) {
    phi = phi_from_theta(theta)
    da = p[1:n_row]
    db = if (d > 0) p[(n_row+1):(n_row+d)] else numeric(0)
    phi_step = as.vector(K %*% da)
    if (d > 0) phi_step = phi_step + as.vector(B %*% db)
    ddchis(phi) * phi_step
  }

  delta0 = if (isTRUE(disp$base$quadratic)) Inf else 1.0
  tr_res = .tr_newton(fg_kernel, hv_kernel, c(alpha, beta),
                      maxiter = maxiter, tol = tol, delta0 = delta0,
                      conv_fn = bregman_conv, weight_fn = weight_fn,
                      dual_map = dual_map)

  alpha = tr_res$theta[1:n_row]
  beta  = if (d > 0) tr_res$theta[(n_row+1):(n_row+d)] else numeric(0)
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
      rbind(cbind(wcK_f + n * eta * diag(n_row), wcB_f),
            cbind(crossprod(B, wcK_f), crossprod(B, wcB_f)))
    } else {
      wcK_f + n * eta * diag(n_row)
    }
    m = nrow(M_f); M_f = M_f + (1e-8 * sum(diag(M_f)) / m) * diag(m)
    M_lu_cache <<- solve(M_f)
    M_lu_cache
  }

  structure(list(
    alpha = alpha, beta = beta,
    kern = kern, Z = Z, K = K, B = B, n = n, eta = eta,
    K_chol = K_chol, get_M_lu = get_M_lu,
    dispersion = dispersion, iter = iter, status = status,
    tr_trace = tr_res$trace,
    w = w, target = target, target_null = target_null
  ), class = "kernel_bregman")
}


#' ## `project_target`
#'
#' **Paper**: the linear functional $\dot\psi_Z$ acts on the
#' representer $\phi = K\alpha + B\beta$ as
#' $$\dot\psi_Z(\phi) = \sum_i r_i\,\phi(Z_i^{\dot\psi})
#'   = c^\top\alpha + c_0^\top\beta$$
#' where $r_i$ are the signed measure weights and
#' $Z_i^{\dot\psi}$ are the evaluation points (possibly
#' counterfactual: $(u, w, X_i)$ for each arm $w$).
#'
#' **Code**: $c_j = \sum_i r_i\,K(Z_j^{\mathrm{basis}}, Z_i^{\dot\psi})
#'   = (K_{\mathrm{cross}}\, r)_j$ and
#' $c_0 = B^\top r$. Returns `c(c, c_0)`.
project_target = function(basis_Z, kern, measure, K_cross = NULL) {
  basis_Z = atleast_2d(basis_Z)
  measure$Z = atleast_2d(measure$Z)
  if (is.null(K_cross))
    K_cross = kernel_matrix(basis_Z, measure$Z, kern)
  B = null_basis(measure$Z, kern)
  c(as.vector(K_cross %*% measure$r),
    if (ncol(B) > 0) as.vector(crossprod(B, measure$r)) else numeric(0))
}

#' ## `split_target`
#'
#' Splits the concatenated output of `project_target` into
#' `target` (first $n$ entries) and `target_null` (remaining $d$ entries).
split_target = function(c_vec, n) {
  list(target = c_vec[1:n],
       target_null = if (length(c_vec) > n) c_vec[(n + 1):length(c_vec)] else numeric(0))
}


#' ## Cross-Validation
#'
#' ### `cv_dual_loss`
#'
#' **Paper**: the dual objective evaluated on held-out data:
#' $$\hat L_{\mathrm{cv}}(\eta) = \frac{1}{|I_{\mathrm{test}}|}
#'   \sum_{i \in I_{\mathrm{test}}} \bigl[\chi^*(\hat\phi_\eta(Z_i))
#'   - c_i\,\hat\phi_\eta(Z_i)\bigr]$$
#' where $c_i$ is the per-observation target from `project_target`
#' and $\hat\phi_\eta$ is the model fitted on the training fold.
#'
#' **Code**: `phi_test` $= \hat\phi_\eta(Z_{\mathrm{test}})$ via `.phi`,
#' then `mean(chis(phi_test) - target_test * phi_test)`.
cv_dual_loss = function(fit, Z_test, target_test, dispersion_test = NULL) {
  phi  = \(Z) .phi(fit, Z)
  disp = if (!is.null(dispersion_test)) dispersion_test else fit$dispersion
  chis = \(p) .chis(disp, p)
  phi_test = phi(Z_test)
  mean(chis(phi_test) - target_test * phi_test)
}

#' ### `loo_loss`
#'
#' **Paper**: leave-one-out cross-validation loss, linearized at
#' convergence. At the optimum $\hat\phi$, the LOO residual for
#' observation $i$ is
#' $$e_i^{\mathrm{loo}} = \frac{e_i^{\mathrm{raw}}}{1 - H_{ii}}$$
#' where $e_i^{\mathrm{raw}} = n\eta\,\alpha_i / w_{c,i}$ is the
#' raw residual ($w_{c,i} = w_i\,\ddot\chi^*(\phi_i)$ is the
#' curvature weight) and $H_{ii}$ is the $i$-th diagonal of the
#' hat matrix $H = W_c^{1/2} K (W_c^{1/2} K W_c^{1/2} + n\eta I)^{-1} W_c^{1/2}$.
#'
#' **Code**: $\tilde K = W_c^{1/2} K W_c^{1/2}$ has eigendecomposition
#' $V \Lambda V^\top$. Then $H_{ii} = \sum_j V_{ij}^2\,\lambda_j /
#' (\lambda_j + n\eta)$ and
#' `loo_loss` $= n^{-1}\sum_i w_{c,i}\,(e_i^{\mathrm{loo}})^2$.
loo_loss = function(fit) {
  K = fit$K; n_row = nrow(K); n = if (!is.null(fit$n)) fit$n else n_row
  phi = as.vector(K %*% fit$alpha)
  d = length(fit$beta)
  if (d > 0) phi = phi + as.vector(fit$B %*% fit$beta)
  w = if (!is.null(fit$w)) fit$w else rep(1, n_row)
  ddchis = \(p) .ddchis(fit$dispersion, p)
  wc = w * ddchis(phi)
  wc = pmax(wc, 1e-10 * max(abs(wc), 1))
  raw_resid = n * fit$eta * fit$alpha / wc
  sqw = sqrt(wc)
  Kt = sweep(sweep(K, 1, sqw, "*"), 2, sqw, "*")
  evd = eigen(Kt, symmetric = TRUE)
  shrink = pmax(evd$values, 0) / (pmax(evd$values, 0) + n * fit$eta)
  H_diag = rowSums(evd$vectors^2 * rep(shrink, each = n_row))
  loo_resid = raw_resid / pmax(1 - H_diag, 1e-10)
  mean(wc * loo_resid^2)
}


#' ## `onestep_bregman`
#'
#' **Paper**: Poisson multiplier bootstrap. Draw $\xi_i \sim
#' \mathrm{Poisson}(1)$, reweight, and take one Newton step from
#' the original solution $(\hat\alpha, \hat\beta)$:
#' $$\theta^* = \hat\theta - M^{-1}\,\nabla L^*(\hat\theta)$$
#' where $M = H(\hat\theta)$ is the Hessian at the original
#' solution (cached via `get_M_lu`) and $\nabla L^*$ is the
#' gradient with bootstrap weights $w^* = \xi$ and possibly
#' updated targets $c^*$.
#'
#' **Code**: the gradient $\nabla_\alpha L^* = K(w^* \odot
#' \dot\chi^*(\hat\phi)) + n\eta\,K\hat\alpha - c^*$ is
#' pre-multiplied by $K^{-1}$ (via `pchol_solve`), then
#' $$\delta = M^{-1}\begin{pmatrix} K^{-1}(-\nabla_\alpha L^*) \\
#'   -\nabla_\beta L^* \end{pmatrix}$$
#' gives $\alpha^* = \hat\alpha + \delta_\alpha$,
#' $\beta^* = \hat\beta + \delta_\beta$.
#'
#' One step is correct because the multiplier bootstrap targets
#' the linear approximation; iterating would undo the perturbation.
onestep_bregman = function(model, target_new, target_null_new = NULL, w_new = NULL) {
  K = model$K; n_row = nrow(K); n = if (!is.null(model$n)) model$n else n_row
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
  d_a = delta[1:n_row]
  d_b = if (d > 0) delta[(n_row + 1):(n_row + d)] else numeric(0)
  model$alpha = alpha + d_a
  model$beta  = if (d > 0) beta + d_b else model$beta
  model
}


# ============================================================
# S3 generics
# ============================================================

.bind    = function(model, Z, ...) UseMethod(".bind")
.phi     = function(model, Z, ...) UseMethod(".phi")


#' ## `.bind.kernel_bregman`
#'
#' Cache for product kernels that block on time (column 1 in `iw`).
#' When `k_x` only sees `X`, `K_x(X_eval, X_centers)` is
#' time-independent and can be precomputed once. `.phi` then
#' assembles the full kernel matrix as `K_x * mask`.
#'
#' For kernels where time is in `k_x` (anisotropic), no cache — `.phi`
#' falls through to `kernel_matrix`.
#'
#' Callers pass `Z = eval_Z = (A, X)` (no time column).
#' Centers are `model$Z = (u, A, X)` (with time).
#' `k_x` sees the non-indicator columns: for `iw = c(1, 2)`, that's `X`.
.bind.kernel_bregman = function(model, Z, ...) {
  Z = atleast_2d(Z)
  kern = model$kern
  if (!inherits(kern, "product_kernel") || !(1 %in% kern$iw))
    return(model)

  # Product kernel with time blocked: k_x sees only X.
  # niw of centers (u, A, X) with iw = c(1, 2) gives X columns.
  Z_centers = model$Z
  niw = setdiff(1:ncol(Z_centers), kern$iw)
  X_centers = Z_centers[, niw, drop = FALSE]

  # eval_Z = (A, X). Strip iw columns that are > 1 (arm), shifted by -1
  # since eval_Z lacks the time column.
  iw_eval = kern$iw[kern$iw > 1] - 1
  niw_eval = setdiff(1:ncol(Z), iw_eval)
  X_eval = Z[, niw_eval, drop = FALSE]

  model$.Kx_cache = kernel_matrix(X_eval, X_centers, kern$k_x)
  model$.bind_n = nrow(Z)
  model
}


`[.kernel_bregman` = function(x, i, ...) x  # no-op: .phi computes K fresh


#' ## `.phi.kernel_bregman`
#'
#' Evaluate the linear predictor
#' $\hat\phi(Z) = \sum_j \alpha_j\,k(Z, Z_j) + B(Z)^\top\beta$
#' at new points $Z$.
#'
#' One path: `kernel_matrix(Z, model$Z, kern)`. If `.bind` cached
#' `K_x` for a time-blocked product kernel, assembles the full
#' kernel matrix as `K_x * mask` instead. Same result, avoids
#' recomputing `K_x` each call.
.phi.kernel_bregman = function(model, Z, ...) {
  Z = atleast_2d(Z)
  kern = model$kern

  if (!is.null(model$.Kx_cache) && nrow(Z) == model$.bind_n) {
    # Cached: K = K_x * mask. K_x is precomputed, mask is cheap.
    iw = kern$iw
    mask = matrix(TRUE, nrow(Z), nrow(model$Z))
    for (j in iw) mask = mask & outer(Z[, j], model$Z[, j], "==")
    K_eval = model$.Kx_cache * mask
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
