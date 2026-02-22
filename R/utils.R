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
