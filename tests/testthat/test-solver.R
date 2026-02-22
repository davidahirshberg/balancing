# Test the unified Bregman solver.

test_that("entropy_dispersion matches old exp_link kernel_bregman", {
  old_env = new.env()
  source(file.path(Sys.getenv("SPINOFFS_DIR", "/Users/skip/work/spinoffs"),
                   "code/spinoff3/kernel.R"), local = old_env)

  pkg_dir = file.path(Sys.getenv("BALANCING_DIR", "/Users/skip/work/balancing"), "R")
  for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R"))
    source(file.path(pkg_dir, f))

  set.seed(42)
  n = 50; p = 3
  Z = matrix(rnorm(n * p), n, p)
  Y = rpois(n, lambda = exp(0.5 * Z[, 1]))
  kern = matern_kernel(sigma = 2, nu = 3/2)
  K = kernel_matrix(Z, Z, kern)
  sigma = 1

  old_fit = old_env$kernel_bregman(Y = Y, Z = Z, kern = kern, sigma = sigma,
                                    intercept = TRUE, link = old_env$exp_link(), K = K)

  new_fit = kernel_bregman(Y = Y, Z = Z, kern = kern, eta = sigma^2,
                           dispersion = entropy_dispersion(), intercept = TRUE, K = K)

  expect_equal(new_fit$alpha, old_fit$alpha, tolerance = 1e-6)
  expect_equal(new_fit$beta0, old_fit$beta0, tolerance = 1e-6)
})

test_that("softplus_dispersion matches old logistic_link kernel_bregman", {
  old_env = new.env()
  source(file.path(Sys.getenv("SPINOFFS_DIR", "/Users/skip/work/spinoffs"),
                   "code/spinoff3/kernel.R"), local = old_env)

  pkg_dir = file.path(Sys.getenv("BALANCING_DIR", "/Users/skip/work/balancing"), "R")
  for (f in c("utils.R", "kernels.R", "dispersions.R", "bregman.R"))
    source(file.path(pkg_dir, f))

  set.seed(42)
  n = 50; p = 3
  Z = matrix(rnorm(n * p), n, p)
  Y = rbinom(n, 1, plogis(Z[, 1]))
  kern = matern_kernel(sigma = 2, nu = 3/2)
  K = kernel_matrix(Z, Z, kern)
  sigma = 1

  old_fit = old_env$kernel_bregman(Y = Y, Z = Z, kern = kern, sigma = sigma,
                                    intercept = TRUE, link = old_env$logistic_link(), K = K)

  new_fit = kernel_bregman(Y = Y, Z = Z, kern = kern, eta = sigma^2,
                           dispersion = softplus_dispersion(), intercept = TRUE, K = K)

  expect_equal(new_fit$alpha, old_fit$alpha, tolerance = 1e-6)
  expect_equal(new_fit$beta0, old_fit$beta0, tolerance = 1e-6)
})

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
  Y = colSums(kernel_matrix(Z1, Z, kern)) / n

  disp = signflip(entropy_dispersion(), W)
  fit = kernel_bregman(Y = Y, Z = Z, kern = kern, eta = 1,
                       dispersion = disp, intercept = TRUE, K = K)

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
  Y = colSums(kernel_matrix(Z1, Z, kern)) / n

  disp = signflip(pos_quadratic_dispersion(), W)
  fit = kernel_bregman(Y = Y, Z = Z, kern = kern, eta = 1,
                       dispersion = disp, intercept = TRUE, K = K)

  gamma = predict(fit, Z)
  expect_true(all(gamma[W == 1] >= -1e-10))
  expect_true(all(gamma[W == -1] <= 1e-10))
})

test_that("quadratic_dispersion matches direct solve", {
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

  fit = kernel_bregman(Y = Y, Z = Z, kern = kern, eta = eta_val,
                       dispersion = quadratic_dispersion(), intercept = TRUE, K = K)

  # Manual solve: (K + 2*eta*I) alpha + beta0 * 1 = Y, 1'alpha + n*beta0 = sum(Y)
  # With intercept via augmented system
  H = K + 2 * eta_val * diag(n)
  w = rep(1, n)
  M = rbind(cbind(H, w), c(w, n))
  sol = solve(M, c(Y, sum(Y)))
  expect_equal(fit$alpha, sol[1:n], tolerance = 1e-10)
  expect_equal(fit$beta0, sol[n + 1], tolerance = 1e-10)
})
