## Kernel infrastructure.
##
## A kernel defines a seminorm rho on an RKHS: the reproducing kernel k(z,z')
## and a basis for ker(rho) (the unpenalized null space). Constructors return
## closure-lists with type/sigma/nu attributes. kernel_matrix dispatches on
## type for fast BLAS paths.

# ============================================================
# Kernel constructors
# ============================================================

#' Matern kernel with lengthscale sigma and smoothness nu.
#' Closed forms for nu = 1/2, 3/2, 5/2; general Bessel fallback.
#'
#' @param null_basis Function Z -> matrix(n, d) evaluating the null-space basis
#'   of the seminorm at arbitrary points. Default: constant functions (intercept).
#'   For no null space: function(Z) matrix(nrow = nrow(Z), ncol = 0).
matern_kernel = function(sigma = 1, nu = 3/2,
                         null_basis = function(Z) matrix(1, nrow(atleast_2d(Z)), 1)) {
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
  structure(pointwise, class = "kernel", type = "matern", sigma = sigma, nu = nu,
            null_basis = null_basis)
}

#' Gaussian (RBF) kernel with lengthscale sigma.
#' @param null_basis See matern_kernel.
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

#' Direct product of RKHS copies, one per grouping level.
#'
#' Takes a kernel k_x on X and returns the direct product space on W x X:
#'   Kernel:   k((w,x),(w',x')) = k_x(x,x') * 1{w=w'}
#'   Null space: one copy of ker(rho_x) per level — if ker(rho_x) = span{1},
#'               then ker(rho) = span{1_{w=0}, 1_{w=1}}
#'   Seminorm: rho(f)^2 = sum_g rho_x(f_g)^2
#'
#' @param k_x Base kernel on X.
#' @param iw Column index(es) in Z that hold the grouping variable(s).
#' @param levels The grouping levels (e.g. c(0, 1) for binary treatment).
#'   Part of the space definition — not discovered from data.
direct_product = function(k_x, iw = 1, levels = NULL) {
  nb_base = attr(k_x, "null_basis")
  nb = function(Z) {
    Z = atleast_2d(Z)
    niw = setdiff(seq_len(ncol(Z)), iw)
    B_base = if (!is.null(nb_base)) {
      nb_base(Z[, niw, drop = FALSE])
    } else {
      matrix(1, nrow(Z), 1)
    }
    d_base = ncol(B_base)
    if (length(iw) == 1) {
      gv = Z[, iw]
    } else {
      gv = apply(Z[, iw, drop = FALSE], 1, paste, collapse = "_")
    }
    # Use stored levels if provided at construction, otherwise auto-detect.
    # Stored levels eliminate threading; auto-detect is fallback for
    # cases where levels depend on the data (e.g. discrete time mesh).
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

# ============================================================
# Null-space basis
# ============================================================

#' Evaluate the null-space basis of a kernel's seminorm at points Z.
#' Returns matrix(n, d) where d = dim(ker rho). If the kernel has no
#' null_basis attribute, defaults to constant functions (intercept).
null_basis = function(Z, kern) {
  Z = atleast_2d(Z)
  nb = attr(kern, "null_basis")
  if (is.null(nb)) return(matrix(1, nrow(Z), 1))
  nb(Z)
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
