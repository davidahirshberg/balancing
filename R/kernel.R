#' # Kernel Infrastructure
#'
#' Kernel constructors (isotropic + product), kernel matrix
#' computation, pivoted Cholesky, block-diagonal matrix, and
#' link functions.

# When sourced directly (not via load_all), pull in utils.R.
if (!exists("atleast_2d", mode = "function")) source("R/utils.R")


#' ## Pivoted Cholesky
#'
#' `pchol(K)` computes \eqn{K[\\pi, \\pi] = R^\\top R} where \eqn{R} is
#' upper triangular and \eqn{\\pi} is a permutation that improves
#' numerical stability. `pchol_solve(pc, b)` solves \eqn{Kx = b}
#' via forward/back substitution on the permuted system.
pchol = function(K, eps = 1e-8) {
  n = ncol(K)
  K = K + (eps * sum(diag(K)) / n) * diag(n)
  R = suppressWarnings(chol(K, pivot = TRUE))
  structure(list(R = R, pivot = attr(R, "pivot"),
                 rank = attr(R, "rank"), n = n),
            class = "pchol")
}

#' Solve K x = b using a pivoted Cholesky factorization.
pchol_solve = function(pc, b) {
  bp = b[pc$pivot]
  z = forwardsolve(t(pc$R), bp, upper.tri = FALSE)
  x_piv = backsolve(pc$R, z)
  x = numeric(pc$n)
  x[pc$pivot] = x_piv
  x
}

# Legacy wrapper for backward compatibility
.chol_solve = function(R, b) {
  if (inherits(R, "pchol")) return(pchol_solve(R, b))
  backsolve(R, forwardsolve(R, b, upper.tri = TRUE, transpose = TRUE))
}


#' ## Block-Diagonal Matrix
#'
#' **Paper**: for a product kernel \eqn{k((w,x),(w',x')) =
#' k_X(x,x')\\,\\mathbf\{1\}\\\{w=w'\\\}}, the kernel matrix is
#' block-diagonal with one dense block per treatment level.
#'
#' **Code**: `block_diag` stores the dense sub-blocks without
#' materializing the full \eqn{n \\times m} matrix. The `%*%` method
#' performs block-diagonal matrix-vector multiply:
#' \eqn{y_\{\\mathrm\{rows\}_b\} \\mathrel\{+\}= K_b\\, x_\{\\mathrm\{cols\}_b\}}.
block_diag = function(blocks, n_row) {
  structure(list(blocks = blocks, n_row = n_row), class = "block_diag")
}

#' Matrix-vector multiply for block_diag: y = sum_b K_b %*% x[cols_b],
#' placed into result[rows_b].
`%*%.block_diag` = function(x, y) {
  result = numeric(x$n_row)
  for (b in x$blocks) {
    if (length(b$rows) > 0)
      result[b$rows] = result[b$rows] + b$K %*% y[b$cols]
  }
  result
}


#' ## Kernel Constructors
#'
#' ### `matern_kernel`
#'
#' **Paper**: the Matern kernel with smoothness \eqn{\\nu} and
#' length scale \eqn{\\sigma}:
#' \deqn{k(z, z') = \\frac\{2^\{1-\\nu\}\}\{\\Gamma(\\nu)\}
#'   \\Bigl(\\frac\{\\sqrt\{2\\nu\}\}\{\\sigma\}\\|z-z'\\|\\Bigr)^\\nu
#'   K_\\nu\\!\\Bigl(\\frac\{\\sqrt\{2\\nu\}\}\{\\sigma\}\\|z-z'\\|\\Bigr)}
#' Special cases: \eqn{\\nu = 1/2} (exponential), \eqn{\\nu = 3/2}
#' (default), \eqn{\\nu = 5/2}.
#'
#' The `null_basis` attribute defines \eqn{\\ker(\\rho)}. Default:
#' constant functions (intercept, \eqn{B = \\mathbf\{1\}}). For norm
#' penalty (no null space): `function(Z) matrix(nrow=nrow(Z), ncol=0)`.
matern_kernel = function(sigma = 1, nu = 3/2, scale = NULL,
                         null_basis = function(Z) matrix(1, nrow(atleast_2d(Z)), 1)) {
  pointwise = function(z, zp) {
    z = atleast_2d(z)
    diff = t(z) - c(zp)
    if (!is.null(scale)) diff = diff / scale  # per-dimension length scales
    d2 = colSums(diff^2)
    r = sqrt(2 * nu) * sqrt(d2)
    if (nu == 1/2) return(exp(-r))
    if (nu == 3/2) return((1 + r) * exp(-r))
    if (nu == 5/2) return((1 + r + r^2 / 3) * exp(-r))
    ifelse(r == 0, 1, (2^(1 - nu) / gamma(nu)) * r^nu * besselK(r, nu))
  }
  if (is.null(scale) && sigma != 1) scale = NULL  # isotropic: handled via sigma in fast path
  structure(pointwise, class = "kernel", type = "matern", sigma = sigma, nu = nu,
            scale = scale, null_basis = null_basis)
}


#' ### `direct_product_kernel`
#'
#' **Paper**: the product kernel
#' \deqn{k\\bigl((w,x),(w',x')\\bigr) = k_X(x,x')\\,\\mathbf\{1\}\\\{w=w'\\\}}
#' decomposes the RKHS into one copy of \eqn{\\mathcal\{H\}_\{k_X\}} per
#' treatment arm. The null space of the product seminorm has one
#' copy of \eqn{\\ker(\\rho_X)} per level: \eqn{B = \\mathrm\{blkdiag\}(B_X,
#' \\ldots, B_X)}.
#'
#' **Code**: `iw` gives the column indices in \eqn{Z} that hold the
#' grouping variable(s). `levels` declares the set of
#' treatment values (data must only contain declared levels).
#' The kernel matrix is computed as \eqn{K_X \\odot \\mathbf\{1\}\\\{W_i = W_j\\\}}.
direct_product_kernel = function(k_x, iw = 1, levels) {
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
      stop("direct_product_kernel: data contains undeclared levels: ",
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


#' ## Kernel Matrix Computation
#'
#' `kernel_matrix(A, B, kern)` computes \eqn{K_\{ij\} = k(A_i, B_j)}.
#' S3 dispatch on `kern`:
#' - `product_kernel`: \eqn{K_X \\odot \\mathrm\{mask\}}, mask from grouping cols
#' - `kernel` (isotropic): fast squared-distance + `.kernel_from_D2`
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
    sc = attr(kern, "scale")
    if (!is.null(sc)) {
      A = t(t(A) / sc); B = t(t(B) / sc)
    }
    D2 = .fast_sqdist(A, B)
    return(.kernel_from_D2(D2, kern))
  }
  outerproduct(A, B, kern)
}

.kernel_from_D2 = function(D2, kern) {
  ktype = attr(kern, "type")
  sig = attr(kern, "sigma")
  sc = attr(kern, "scale")
  if (ktype == "matern") {
    nu = attr(kern, "nu")
    # When scale is set, data is already rescaled; otherwise use isotropic sigma
    R = if (!is.null(sc)) sqrt(2 * nu) * sqrt(D2) else (sqrt(2 * nu) / sig) * sqrt(D2)
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


#' ## Link Functions (Survival)
#'
#' **Paper**: the hazard model uses a link function \eqn{\\sigma}
#' mapping the dual variable \eqn{\\phi_\\lambda} to the hazard rate:
#' \eqn{\\hat\\lambda_u = \\sigma(\\hat\\phi_\{\\lambda,u\})}.
#'
#' - **Exponential**: \eqn{\\sigma(\\phi) = e^\\phi}. Continuous-time
#'   survival via \eqn{S_t = \\exp(-\\int_0^t \\lambda_u\\,du)}
#'   (trapezoidal quadrature).
#' - **Logistic**: \eqn{\\sigma(\\phi) = e^\\phi/(1+e^\\phi)}.
#'   Discrete-time survival via \eqn{S_t = \\prod_u (1 - \\lambda_u)}.
#'
#' Each link provides `cum_surv(haz_fn, Q, du)` for
#' \eqn{\\hat S_t}, `cum_rmst(...)` for \eqn{\\int_0^t \\hat S_u\\,du},
#' and `surv_curve(...)` for the full \eqn{[\\hat S_\{u_1\}, \\ldots,
#' \\hat S_\{u_Q\}]} matrix.
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
