## Utilities

#' Ensure Z is a matrix (at least 2d).
atleast_2d = function(Z) {
  if (is.null(dim(Z))) Z = matrix(Z, ncol = 1)
  if (length(dim(Z)) == 1) dim(Z) = c(dim(Z), 1)
  Z
}

#' Set treatment column to a given value.
with.treatment = function(z, a) { z2 = z; z2[, 1] = a; z2 }

#' Legacy outerproduct: loop over rows of B. Slow but general.
outerproduct = function(A, B, K) {
  A = atleast_2d(A)
  B = atleast_2d(B)
  O = matrix(0, nrow(A), nrow(B))
  for (ib in 1:nrow(B)) O[, ib] = K(A, B[ib, , drop = FALSE])
  O
}

#' CG solver for symmetric PD systems.
#' matvec: function v -> A*v. precond: optional preconditioner (function).
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
    if (pAp <= 0) break  # not PD in this direction; bail
    a = rz / pAp
    x = x + a * p
    r = r - a * Ap
    rnorm = sqrt(sum(r^2))
    if (rnorm < tol * bnorm) break
    z = if (!is.null(precond)) precond(r) else r
    rz_new = sum(r * z)
    p = z + (rz_new / rz) * p
    rz = rz_new
  }
  list(x = x, iter = k, resid = sqrt(sum(r^2)))
}

#' Trust-region truncated Newton (Steihaug-Toint CG).
#'
#' Solves: min f(θ) using exact Hessian-vector products via CG.
#' Each outer iteration approximately solves H·p = -g by CG, truncating
#' when ||p|| hits the trust radius (Steihaug) or CG detects negative
#' curvature. Trust region handles overflow (val=Inf → shrink radius).
#'
#' @param fg Function theta -> list(val, grad). May return val=Inf.
#' @param hv Function (theta, v) -> H(theta)·v. Hessian-vector product.
#' @param theta0 Initial parameter vector.
#' @param precond Optional function v -> M⁻¹v (preconditioner).
#' @param maxiter Maximum outer (Newton) iterations.
#' @param cg_maxiter Maximum CG iterations per Newton step.
#' @param tol Convergence tolerance on gradient norm.
#' @param delta0 Initial trust-region radius.
#' @return list(theta, val, grad, iter)
.tr_newton = function(fg, hv, theta0, precond = NULL,
                      maxiter = 50, cg_maxiter = 50, tol = 1e-6,
                      delta0 = 1.0) {
  theta = theta0
  res = fg(theta)
  val = res$val; grad = res$grad
  delta = delta0
  trace = data.frame(iter = integer(0), val = numeric(0),
                     gnorm = numeric(0), delta = numeric(0),
                     accepted = logical(0))

  for (iter in 1:maxiter) {
    gnorm = sqrt(sum(grad^2))
    if (gnorm < tol) break

    # Steihaug-Toint CG: approximately solve H·p = -g within trust region.
    # If we hit the boundary or negative curvature, return the boundary point.
    p = .steihaug_cg(function(v) hv(theta, v), -grad, delta,
                     precond = precond, maxiter = cg_maxiter,
                     tol = min(0.5, sqrt(gnorm)) * gnorm)

    # Evaluate candidate
    theta_new = theta + p
    res_new = fg(theta_new)
    val_new = res_new$val

    # Predicted reduction: -g'p - 0.5 p'Hp (quadratic model)
    Hp = hv(theta, p)
    predicted = -sum(grad * p) - 0.5 * sum(p * Hp)
    actual = val - val_new

    if (predicted > 0 && is.finite(val_new) && actual / predicted > 0.1) {
      # Accept step
      theta = theta_new; val = val_new; grad = res_new$grad
      ratio = actual / predicted
      pnorm = sqrt(sum(p^2))
      trace[nrow(trace) + 1, ] = list(iter, val, gnorm, delta, TRUE)
      if (ratio > 0.75 && pnorm > 0.8 * delta)
        delta = min(delta * 2, 1e8)
      else if (ratio < 0.25)
        delta = delta * 0.25
    } else {
      # Reject step, shrink
      trace[nrow(trace) + 1, ] = list(iter, val, gnorm, delta, FALSE)
      delta = delta * 0.25
      if (delta < 1e-15 * (1 + sqrt(sum(theta^2)))) break
    }
  }
  list(theta = theta, val = val, grad = grad, iter = iter, trace = trace)
}

#' Steihaug-Toint truncated CG for trust-region subproblem.
#'
#' Approximately solves H·p = rhs subject to ||p|| <= delta.
#' Returns early if: (a) negative curvature detected, (b) step hits boundary,
#' (c) residual small enough.
#'
#' @param hv Function v -> H·v.
#' @param rhs Right-hand side (typically -grad).
#' @param delta Trust-region radius.
#' @param precond Optional preconditioner M⁻¹.
#' @param maxiter Max CG iterations.
#' @param tol Residual tolerance.
#' @return Step vector p with ||p|| <= delta.
.steihaug_cg = function(hv, rhs, delta, precond = NULL, maxiter = 50, tol = 1e-8) {
  p = numeric(length(rhs))
  r = rhs  # residual = rhs - H·0 = rhs
  z = if (!is.null(precond)) precond(r) else r
  d = z
  rz = sum(r * z)
  rhs_norm = sqrt(sum(rhs^2))

  for (k in 1:maxiter) {
    Hd = hv(d)
    dHd = sum(d * Hd)

    if (dHd <= 0) {
      # Negative curvature: go to trust-region boundary along d
      return(.to_boundary(p, d, delta))
    }

    alpha = rz / dHd
    p_new = p + alpha * d

    if (sqrt(sum(p_new^2)) >= delta) {
      # Step exceeds trust region: go to boundary along d
      return(.to_boundary(p, d, delta))
    }

    p = p_new
    r = r - alpha * Hd
    if (sqrt(sum(r^2)) < tol) break

    z = if (!is.null(precond)) precond(r) else r
    rz_new = sum(r * z)
    d = z + (rz_new / rz) * d
    rz = rz_new
  }
  p
}

#' Move from p along direction d to the trust-region boundary ||p + t*d|| = delta.
.to_boundary = function(p, d, delta) {
  # Solve ||p + t*d||² = delta² for t > 0
  pp = sum(p^2); pd = sum(p * d); dd = sum(d^2)
  disc = pd^2 - dd * (pp - delta^2)
  if (disc < 0) return(p)  # shouldn't happen
  t = (-pd + sqrt(disc)) / dd
  p + t * d
}
