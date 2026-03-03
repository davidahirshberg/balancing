## Smoke test: verify S3 dispersion refactor + surv_estimate plumbing.
## Tiny dataset, no bootstrap. Just checks that everything sources
## and produces a finite DR estimate.

setwd("/Users/skip/work/balancing")
source("R/kernel.R")
source("R/dispersions.R")
source("R/estimands.R")
source("R/survival.R")
source("R/surv_estimate.R")
source("R/dgp.R")

cat("=== S3 dispersion smoke test ===\n\n")

# --- 1. Verify S3 dispatch on base dispersions ---
cat("1. S3 dispatch on base dispersions\n")
d = entropy_dispersion()
stopifnot(inherits(d, "dispersion"))
stopifnot(abs(.dchis(d, 0) - 1) < 1e-10)       # exp(0) = 1
stopifnot(abs(.dchis(d, 1) - exp(1)) < 1e-10)   # exp(1)
cat("   entropy: OK\n")

d = softplus_dispersion()
stopifnot(abs(.dchis(d, 0) - 0.5) < 1e-10)      # expit(0) = 0.5
cat("   softplus: OK\n")

d = quadratic_dispersion()
stopifnot(abs(.dchis(d, 3) - 3) < 1e-10)         # identity
stopifnot(isTRUE(d$base$quadratic))
cat("   quadratic: OK\n")

# --- 2. Verify composition (target_scaled_entropy) ---
cat("2. Composition: target_scaled_entropy\n")
r = c(0.5, -0.3, 0.8)
d = target_scaled_entropy(r)
test_phi = c(0, 0, 0)
gamma = .dchis(d, test_phi)
# gamma = r + sign(r) * exp(sign(r) * phi) = r + sign(r) * 1
expected = r + sign(r)
stopifnot(max(abs(gamma - expected)) < 1e-10)
cat("   dchis at phi=0: OK\n")

# Check sign constraint: |gamma| >= |r|
stopifnot(all(abs(gamma) >= abs(r) - 1e-10))
cat("   sign constraint: OK\n")

# --- 3. Verify variance_weighted_quadratic ---
cat("3. Composition: variance_weighted_quadratic\n")
r = c(1, -2); v = c(0.5, 0.5)
d = variance_weighted_quadratic(r, v)
test_phi = c(0.1, -0.1)
gamma = .dchis(d, test_phi)
# psi = sign(r) * phi / v, dchis = r + sign(r) * pmax(psi, 0)
psi = sign(r) * test_phi / v  # c(0.2, 0.2)
expected = r + sign(r) * pmax(psi, 0)  # c(1 + 0.2, -2 - 0.2) = c(1.2, -2.2)
stopifnot(max(abs(gamma - expected)) < 1e-10)
cat("   dchis: OK\n")

# --- 4. Verify vectorized dispersion in make_dchis_gamma_ts context ---
cat("4. Vectorized dispersion (dispersion with vector offset)\n")
r_vec = c(0.5, -0.3, 0.1, -0.8)
d = dispersion(entropy_base, offset = r_vec, sigma = sign(r_vec))
test_phi = c(0.1, -0.2, 0.3, -0.1)
gamma = .dchis(d, test_phi)
expected = r_vec + sign(r_vec) * exp(sign(r_vec) * test_phi)
stopifnot(max(abs(gamma - expected)) < 1e-10)
cat("   OK\n")

# --- 5. End-to-end: survival_effect ---
cat("5. End-to-end: survival_effect on paper_cts DGP (n=50)\n")
set.seed(42)
dgp = paper_cts_dgp(horizon = 1)
dat = dgp$generate(50, p = 2)
kern = direct_product(matern_kernel(sigma = 2, nu = 3/2), iw = 2, levels = c(0, 1))

obs = list(T_obs = dat$T_obs, D = dat$D, Z = cbind(dat$A, dat$X))
estimand = surv_prob_estimand()

dr = survival_effect(obs, kern, estimand,
                     lambda_disp = entropy_dispersion(),
                     eta_lam = 1/50, eta_gam = 1/50,
                     horizon = 1, M_train = 10, n_folds = 2, Q_comp = 30)

cat("   est=", sprintf("%.4f", dr$est), "\n")
cat("   se=", sprintf("%.4f", dr$se), "\n")
cat("   se_amle=", sprintf("%.4f", dr$se_amle), "\n")
stopifnot(is.finite(dr$est))
stopifnot(is.finite(dr$se))
stopifnot(dr$se > 0)
cat("   finite estimate + positive SE: OK\n")

truth = dat$psi1_true - dat$psi0_true
cat("   truth=", sprintf("%.4f", truth), "\n")
cat("   error=", sprintf("%+.4f", dr$est - truth), "\n")

# --- 6. gamma_base = quadratic_base (mixed-sign ATE) ---
cat("6. gamma_base = quadratic_base\n")
dr_q = survival_effect(obs, kern, estimand,
                       lambda_disp = entropy_dispersion(),
                       eta_lam = 1/50, eta_gam = 1/50,
                       horizon = 1, M_train = 10, n_folds = 2, Q_comp = 30,
                       gamma_base = quadratic_base)

cat("   est=", sprintf("%.4f", dr_q$est), "\n")
stopifnot(is.finite(dr_q$est))
cat("   finite estimate: OK\n")

# --- 7. Discrete mesh uses grid points, not midpoints ---
cat("7. Discrete mesh\n")
tm_cts = build_training_mesh(1:8, rep(1,8), matrix(0,8,2), 8, 8)
tm_disc = build_training_mesh(1:8, rep(1,8), matrix(0,8,2), 8, 8, discrete = TRUE)
stopifnot(all(tm_cts$mesh == (1:8 - 0.5)))   # midpoints
stopifnot(all(tm_disc$mesh == 1:8))           # grid points
cat("   continuous midpoints:", paste(tm_cts$mesh, collapse=", "), "\n")
cat("   discrete grid points:", paste(tm_disc$mesh, collapse=", "), "\n")
cat("   OK\n")

cat("\n=== All tests passed ===\n")
