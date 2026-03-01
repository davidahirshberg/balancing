# Test the unified Bregman solver.

test_that("signflip(entropy) produces correctly-signed weights", {
  pkg_dir = file.path(Sys.getenv("BALANCING_DIR", "/Users/skip/work/balancing"), "R")
  for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R"))
    source(file.path(pkg_dir, f))

  set.seed(42)
  n = 100; p = 3
  X = matrix(rnorm(n * p), n, p)
  A = rbinom(n, 1, plogis(0.5 * X[, 1]))
  W = 2 * A - 1
  Z = cbind(A, X)
  kern = matern_kernel(sigma = 2, nu = 3/2)
  K = kernel_matrix(Z, Z, kern)


  # Target: dot_psi for ATE (sum of kernel columns with counterfactual treatment)
  Z1 = Z; Z1[, 1] = 1
  c_gam = colSums(kernel_matrix(Z1, Z, kern)) / n

  disp = signflip(entropy_dispersion(), W)
  fit = kernel_bregman(Z, kern, eta = 1,
                       dispersion = disp, target = c_gam, K = K)

  gamma = predict(fit, Z)

  # Treated (W=+1) should have positive weights
  expect_true(all(gamma[W == 1] > 0))
  # Control (W=-1) should have negative weights
  expect_true(all(gamma[W == -1] < 0))
})

test_that("signflip(pos_quadratic) produces correctly-signed weights", {
  pkg_dir = file.path(Sys.getenv("BALANCING_DIR", "/Users/skip/work/balancing"), "R")
  for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R"))
    source(file.path(pkg_dir, f))

  set.seed(42)
  n = 100; p = 3
  X = matrix(rnorm(n * p), n, p)
  A = rbinom(n, 1, plogis(0.5 * X[, 1]))
  W = 2 * A - 1
  Z = cbind(A, X)
  kern = matern_kernel(sigma = 2, nu = 3/2)
  K = kernel_matrix(Z, Z, kern)


  Z1 = Z; Z1[, 1] = 1
  c_gam = colSums(kernel_matrix(Z1, Z, kern)) / n

  disp = signflip(pos_quadratic_dispersion(), W)
  fit = kernel_bregman(Z, kern, eta = 1,
                       dispersion = disp, target = c_gam, K = K)

  gamma = predict(fit, Z)
  expect_true(all(gamma[W == 1] >= -1e-10))
  expect_true(all(gamma[W == -1] <= 1e-10))
})

test_that("quadratic_dispersion gives correct phi", {
  pkg_dir = file.path(Sys.getenv("BALANCING_DIR", "/Users/skip/work/balancing"), "R")
  for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R"))
    source(file.path(pkg_dir, f))

  set.seed(42)
  n = 50; p = 3
  Z = matrix(rnorm(n * p), n, p)
  Y = rnorm(n)
  kern = matern_kernel(sigma = 2, nu = 3/2)
  K = kernel_matrix(Z, Z, kern)

  eta_val = 0.5

  B = null_basis(Z, kern)
  tgt = dpsi_from_response(Y, K, B)
  fit = kernel_bregman(Z, kern, eta = eta_val,
                       dispersion = quadratic_dispersion(),
                       target = tgt$target, target_null = tgt$target_null,
                       K = K)

  # For quadratic dispersion: the optimal phi satisfies
  #   phi_i + 2eta (K^{-1} (phi - B beta))_i = Y_i
  # and B'(phi - Y) = 0 (null space condition)
  #
  # Equivalently: (K + 2eta I)^{-1} (Y - B beta) gives alpha,
  # where 1'K alpha + n beta = sum(Y).
  #
  # Check: predictions at training data match an independent solve.
  # Use fit$K (which includes the solver's nugget) for the reference computation.
  Ks = fit$K
  phi = as.vector(Ks %*% fit$alpha) + as.vector(B %*% fit$beta)

  # Independent computation: solve the full Hessian system
  H = rbind(cbind(crossprod(sweep(Ks, 1, 1, "*")) + n * eta_val * Ks, Ks %*% B),
            cbind(t(B) %*% Ks, crossprod(B)))
  rhs = c(Ks %*% Y, crossprod(B, Y))
  sol = solve(H, rhs)
  phi_check = as.vector(Ks %*% sol[1:n]) + as.vector(B %*% sol[(n+1):(n+ncol(B))])

  expect_equal(phi, phi_check, tolerance = 1e-7)

  # Also: predictions at new points should match
  Z_test = matrix(rnorm(10 * p), 10, p)
  pred = predict_phi(fit, Z_test)
  K_test = kernel_matrix(Z_test, Z, kern)
  B_test = null_basis(Z_test, kern)
  pred_check = as.vector(K_test %*% sol[1:n]) + as.vector(B_test %*% sol[(n+1):(n+ncol(B))])
  expect_equal(pred, pred_check, tolerance = 1e-7)
})

test_that("no-null-space kernel matches pure ridge", {
  pkg_dir = file.path(Sys.getenv("BALANCING_DIR", "/Users/skip/work/balancing"), "R")
  for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R"))
    source(file.path(pkg_dir, f))

  set.seed(42)
  n = 50; p = 3
  Z = matrix(rnorm(n * p), n, p)
  Y = rnorm(n)
  no_null = function(Z) matrix(nrow = nrow(atleast_2d(Z)), ncol = 0)
  kern = matern_kernel(sigma = 2, nu = 3/2, null_basis = no_null)
  K = kernel_matrix(Z, Z, kern)
  eta_val = 0.5

  B = null_basis(Z, kern)
  tgt = dpsi_from_response(Y, K, B)
  fit = kernel_bregman(Z, kern, eta = eta_val,
                       dispersion = quadratic_dispersion(),
                       target = tgt$target, K = K)

  # FOC: Ks(Ks + n*eta*I) alpha = target, where Ks is the solver's nuggeted K
  Ks = fit$K
  alpha_ref = solve(Ks + n * eta_val * diag(n), solve(Ks, tgt$target))
  expect_equal(fit$alpha, alpha_ref, tolerance = 1e-7)
  expect_equal(length(fit$beta), 0)
})
