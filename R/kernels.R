## Kernel infrastructure
## Constructors return closure-lists with type/sigma/nu attributes.
## kernel_matrix dispatches on type for fast BLAS paths.

# ============================================================
# Kernel constructors
# ============================================================

#' Matern kernel with lengthscale sigma and smoothness nu.
#' Closed forms for nu = 1/2, 3/2, 5/2; general Bessel fallback.
matern_kernel = function(sigma = 1, nu = 3/2) {
  pointwise = function(z, zp) {
    z = atleast_2d(z)
    d2 = colSums((t(z) - c(zp))^2)
    d = sqrt(d2)
    r = (sqrt(2 * nu) / sigma) * d
    if (nu == 1/2) return(exp(-r))
    if (nu == 3/2) return((1 + r) * exp(-r))
    if (nu == 5/2) return((1 + r + r^2 / 3) * exp(-r))
    ifelse(r == 0, 1, (2^(1 - nu) / gamma(nu)) * r^nu * besselK(r, nu))
  }
  structure(pointwise, class = "kernel", type = "matern", sigma = sigma, nu = nu)
}

#' Gaussian (RBF) kernel with lengthscale sigma.
gaussian_kernel = function(sigma = 1) {
  pointwise = function(z, zp) {
    z = atleast_2d(z)
    d2 = colSums((t(z) - c(zp))^2)
    exp(-d2 / (2 * sigma^2))
  }
  structure(pointwise, class = "kernel", type = "gaussian", sigma = sigma)
}

#' Direct product kernel: K_x(x,x') * 1{w=w'}.
#' For binary treatment: zero similarity between treatment groups.
#' Restricting to one arm gives K_x.
product_kernel = function(k_x, iw = 1) {
  structure(list(k_x = k_x, iw = iw), class = c("product_kernel", "kernel"))
}

# ============================================================
# Kernel matrix computation
# ============================================================

#' Compute kernel matrix K[i,j] = k(A[i,], B[j,]).
#' Dispatches to BLAS-based code for matern/gaussian,
#' product kernel logic for product_kernel,
#' pointwise loop fallback for custom kernels.
kernel_matrix = function(A, B, kern) {
  UseMethod("kernel_matrix", kern)
}

kernel_matrix.product_kernel = function(A, B, kern) {
  A = atleast_2d(A); B = atleast_2d(B)
  iw = kern$iw
  niw = setdiff(1:ncol(A), iw)
  Kx = kernel_matrix(A[, niw, drop = FALSE], B[, niw, drop = FALSE], kern$k_x)
  if (length(iw) == 1) {
    mask = outer(A[, iw], B[, iw], "==")
  } else {
    # Row-wise equality across multiple grouping columns
    mask = matrix(TRUE, nrow(A), nrow(B))
    for (j in iw) mask = mask & outer(A[, j], B[, j], "==")
  }
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

kernel_matrix.default = function(A, B, kern) {
  A = atleast_2d(A); B = atleast_2d(B)
  outerproduct(A, B, kern)
}

# ============================================================
# Fast internals
# ============================================================

#' Apply kernel function to a squared distance matrix.
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

  if (ktype == "gaussian") {
    return(exp(-D2 / (2 * sig^2)))
  }

  stop("Unknown kernel type: ", ktype)
}

#' Squared distance matrix via BLAS.
#' ||a - b||^2 = ||a||^2 + ||b||^2 - 2 a'b
.fast_sqdist = function(A, B) {
  ssA = rowSums(A^2)
  ssB = rowSums(B^2)
  D2 = outer(ssA, ssB, "+") - 2 * tcrossprod(A, B)
  D2[D2 < 0] = 0  # numerical floor
  D2
}
