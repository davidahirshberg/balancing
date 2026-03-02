#!/usr/bin/env Rscript
## Coverage test: full simulation with method grid, resolution sweep,
## bootstrap coverage for balancing methods, point estimates for GRF/1S.
##
## Cell = DGP × resolution. 15 cells total.
##   1-5   paper_cts × {2, 4, 8, 16, cts}
##   6-10  grf_type1 × {2, 4, 8, 16, cts}
##   11-15 grf_type2 × {2, 4, 8, 16, cts}
##
## Lambda (hazard model) is shared across estimands within a method.
##
## Usage:
##   Rscript test-coverage.R --cell 4 --nreps 200 --boot 200 --workers 50

library(future)
library(furrr)

# ============================================================
# Parse arguments
# ============================================================
args = commandArgs(trailingOnly = TRUE)
n_reps = 50; n_obs = 200; boot_reps = 200; n_workers = 1; start_rep = 1; cell_id = 1

i = 1
while (i <= length(args)) {
  if (args[i] == "--nreps")      { n_reps = as.integer(args[i + 1]); i = i + 2 }
  else if (args[i] == "--n")     { n_obs = as.integer(args[i + 1]); i = i + 2 }
  else if (args[i] == "--boot")  { boot_reps = as.integer(args[i + 1]); i = i + 2 }
  else if (args[i] == "--workers") { n_workers = as.integer(args[i + 1]); i = i + 2 }
  else if (args[i] == "--start-rep") { start_rep = as.integer(args[i + 1]); i = i + 2 }
  else if (args[i] == "--cell")  { cell_id = as.integer(args[i + 1]); i = i + 2 }
  else { i = i + 1 }
}

# ============================================================
# Cell grid: DGP × resolution
# ============================================================
dgp_defs = list(
  paper_cts = list(dgp = "paper_cts", horizon = 1,   p = 2),
  grf_type1 = list(dgp = "grf_type1", horizon = 0.8, p = 5),
  grf_type2 = list(dgp = "grf_type2", horizon = 1.2, p = 5)
)
resolutions = c(2, 4, 8, 16, Inf)  # Inf = continuous

cells = expand.grid(dgp_name = names(dgp_defs), res = resolutions,
                    stringsAsFactors = FALSE)
cells = cells[order(match(cells$dgp_name, names(dgp_defs)), cells$res), ]
rownames(cells) = NULL

if (cell_id < 1 || cell_id > nrow(cells))
  stop(sprintf("Invalid cell %d (max %d)", cell_id, nrow(cells)))
cell = cells[cell_id, ]
dgp_def = dgp_defs[[cell$dgp_name]]
is_cts = is.infinite(cell$res)
cell_name = sprintf("%s_%s", cell$dgp_name,
                    if (is_cts) "cts" else as.integer(cell$res))

# ============================================================
# Method grid
# ============================================================
# Balancing methods: get bootstrap + coverage
# kern: "pham" = product_kernel (block-diagonal), "smooth" = matern (shared time)
# lam_disp: dispersion for lambda (hazard model)
# time: "discrete" or "continuous"
bal_methods = list(
  Pham_quad     = list(kern = "pham",        lam_disp = "quadratic", time = "discrete", n_folds = 2, res = c(2, 4, 8, 16)),
  Pham_normed   = list(kern = "pham_normed", lam_disp = "quadratic", time = "discrete", n_folds = 2, res = c(2, 4, 8, 16)),
  Pham_logistic = list(kern = "pham",        lam_disp = "softplus",  time = "discrete", n_folds = 2, res = c(2, 4, 8, 16)),
  KBW_disc      = list(kern = "smooth", lam_disp = "softplus",  time = "discrete", n_folds = 2, res = c(2, 4, 8, 16)),
  KBW_cts       = list(kern = "smooth", lam_disp = "entropy",   time = "continuous", n_folds = 2, res = c(Inf)),
  PI            = list(kern = "pham",   lam_disp = "quadratic", time = "discrete", n_folds = 3, res = c(2, 4, 8, 16),
                       skip_estimands = "rmst")
)

# External methods: point estimate + SE only, no bootstrap
ext_methods = list(
  GRF = list(res = c(Inf)),
  `1S` = list(res = c(2, 4, 8, 16))
)

active_bal = names(bal_methods)[sapply(bal_methods, function(m) cell$res %in% m$res)]
active_ext = names(ext_methods)[sapply(ext_methods, function(m) cell$res %in% m$res)]

# Estimands: surv and rmst always; tsm1 for discrete only
estimand_names = c("surv_prob", "rmst")
if (!is_cts) estimand_names = c(estimand_names, "tsm1")

# ============================================================
# Setup
# ============================================================
source("R/kernel.R")
for (f in c("R/estimands.R", "R/survival.R", "R/tuning.R", "R/dgp.R"))
  source(f)

# ============================================================
# Helpers (from estimate.R — the pieces that do real work)
# ============================================================

# Build the gamma target measure from the estimand's dpsi_grid.
# Returns a measure list(Z, r) suitable for project_target.
gam_measure = function(estimand, predict_haz, eval_Z, mesh_u, mesh,
                       M_train, du_bin, w_eval = NULL) {
  emb = estimand$dpsi_grid(predict_haz, eval_Z, mesh_u, mesh, M_train, du_bin)
  if (!is.null(w_eval)) {
    n_eval = nrow(atleast_2d(eval_Z))
    n_per_subj = length(emb$r) %/% n_eval
    emb$r = emb$r * rep(w_eval, times = n_per_subj)
  }
  emb
}

# Compute direct estimand from a predict_phi function.
# Returns list(direct, predict_haz).
compute_direct = function(predict_phi_fn, lam_disp, estimand, eval_Z,
                          M_train, du_bin, mesh, horizon, discrete) {
  predict_haz = function(k, Z_ev)
    lam_disp$dchis(predict_phi_fn(cbind(mesh[k], Z_ev)))

  if (discrete) {
    direct = estimand$direct_discrete(predict_haz, eval_Z, M_train, du_bin)
  } else {
    Q_surv = 100
    predict_haz_cts = function(k, Z_ev) {
      phi = predict_phi_fn(cbind((k - 1) * horizon / Q_surv, Z_ev))
      pmax(lam_disp$dchis(phi) / du_bin, 0)
    }
    direct = estimand$direct_cts(predict_haz_cts, eval_Z, horizon, Q_surv)
  }
  list(direct = direct, predict_haz = predict_haz)
}

# Compute DR correction terms (compensator + dirac) and noise variance.
correction_terms = function(lam_bound, gam_bound, T_obs_eval, D_eval,
                            mesh, du_bin, eval_Z, n_eval,
                            lam_disp, dchis_gamma,
                            horizon, Q_comp, discrete) {
  if (discrete) {
    comp = compensator_discrete(lam_bound, gam_bound,
                                T_obs_eval, D_eval, mesh,
                                lam_disp, dchis_gamma, eval_Z)
  } else {
    comp = compensator_simpson(lam_bound, gam_bound,
                               T_obs_eval, D_eval, eval_Z, horizon, Q_comp,
                               lam_disp, dchis_gamma, du_bin)
  }

  M_train = length(mesh)
  dirac = rep(0, n_eval)
  events = which(D_eval == 1 & T_obs_eval <= horizon)
  if (length(events) > 0) {
    t_ev = T_obs_eval[events]
    bin_idx = pmin(ceiling(t_ev / du_bin), M_train)
    t_mesh = mesh[bin_idx]
    phi_ev = eval_phi_vec(gam_bound, t_mesh, events)
    Z_ev = cbind(t_mesh, eval_Z[events, , drop = FALSE])
    dirac[events] = dchis_gamma(Z_ev, phi_ev)
  }

  list(correction = dirac - comp$comp, noise_var = comp$noise_var,
       events = events)
}

# ============================================================
# DR estimator + crossfit + bootstrap (S3)
# ============================================================

# Compute DR estimate on one fold's data.
# Fits gamma, computes direct + correction, caches state for bootstrap.
dr_fold = function(fd, lam_fit, kern, eta_gam, gam_disp, dchis_gamma,
                   lam_disp, estimand, horizon, Q_comp, discrete,
                   gam_target = NULL, gam_target_null = NULL) {
  n_eval = fd$n_eval
  lam_bound = bind_subjects(lam_fit, fd$eval_Z)
  M_train = length(fd$mesh)

  # Direct estimate
  predict_phi_fn = function(Z) predict_phi(lam_fit, Z)
  dr = compute_direct(predict_phi_fn, lam_disp, estimand, fd$eval_Z,
                       M_train, fd$du_bin, fd$mesh, horizon, discrete)

  # Gamma target and dispersion
  if (is.null(gam_target)) {
    emb = gam_measure(estimand, dr$predict_haz, fd$eval_Z,
                      fd$mesh_u, fd$mesh, M_train, fd$du_bin)
    tgt = split_target(project_target(fd$stacked_gam$Z_pool, kern, emb),
                       nrow(fd$stacked_gam$Z_pool))
    gam_target = tgt$target
    gam_target_null = tgt$target_null
    gam_disp = target_scaled_entropy(gam_target)
  }

  # Fit gamma
  gam_fit = kernel_bregman(fd$stacked_gam$Z_pool, kern,
                           eta = eta_gam, dispersion = gam_disp,
                           target = gam_target, target_null = gam_target_null)
  gam_bound = bind_subjects(gam_fit, fd$eval_Z)

  # DR correction
  ct = correction_terms(lam_bound, gam_bound, fd$T_obs_eval, fd$D_eval,
                         fd$mesh, fd$du_bin, fd$eval_Z, n_eval,
                         lam_disp, dchis_gamma, horizon, Q_comp, discrete)

  list(direct = dr$direct, correction = ct$correction,
       noise_var = ct$noise_var, events = ct$events,
       lam_fit = lam_fit, gam_fit = gam_fit,
       lam_bound = lam_bound, gam_bound = gam_bound,
       predict_haz = dr$predict_haz, dchis_gamma = dchis_gamma)
}

# Cross-fit a DR estimator over fold partitions.
# partition_fn(fd) -> per-fold result with $direct, $correction, $noise_var + model state.
# Stores problem params so bootstrap() needs only n_reps.
crossfit = function(n, fold_data, partition_fn,
                    kern, eta_gam, lam_disp, gam_disp,
                    estimand, horizon, Q_comp, discrete) {
  n_folds = length(fold_data)
  fold_results = lapply(fold_data, partition_fn)

  terms = rep(NA_real_, n)
  direct = rep(NA_real_, n)
  noise_var = rep(NA_real_, n)
  for (ff in seq_along(fold_data)) {
    idx = fold_data[[ff]]$eval_idx
    terms[idx] = fold_results[[ff]]$direct + fold_results[[ff]]$correction
    direct[idx] = fold_results[[ff]]$direct
    noise_var[idx] = fold_results[[ff]]$noise_var
  }

  structure(
    list(est = mean(terms), se = sd(terms) / sqrt(n),
         se_amle = sqrt((var(direct) + mean(pmax(noise_var, 0))) / n),
         terms = terms, direct = direct, noise_var = noise_var,
         folds = fold_data, fold_results = fold_results, n = n,
         kern = kern, eta_gam = eta_gam, lam_disp = lam_disp,
         gam_disp = gam_disp, estimand = estimand,
         horizon = horizon, Q_comp = Q_comp, discrete = discrete),
    class = "dr_result")
}

# S3 generic
bootstrap = function(obj, ...) UseMethod("bootstrap")

# Multiplier bootstrap for dr_result.
# Returns a boot_result with $boot_ates, $boot_ses, and CI methods.
bootstrap.dr_result = function(obj, n_reps, ...) {
  n = obj$n
  kern = obj$kern; lam_disp = obj$lam_disp; gam_disp = obj$gam_disp
  estimand = obj$estimand; horizon = obj$horizon
  Q_comp = obj$Q_comp; discrete = obj$discrete; eta_gam = obj$eta_gam

  boot_ates = rep(NA_real_, n_reps)
  boot_ses = rep(NA_real_, n_reps)

  # Pre-compute cross-kernel matrices for gamma target projection
  K_cross_gam = lapply(seq_along(obj$folds), function(ff) {
    fd = obj$folds[[ff]]
    M_train = length(fd$mesh)
    dummy_haz = function(k, Z_ev) rep(0, nrow(atleast_2d(Z_ev)))
    emb_template = estimand$dpsi_grid(dummy_haz, fd$eval_Z, fd$mesh_u,
                                       fd$mesh, M_train, fd$du_bin)
    kernel_matrix(fd$stacked_gam$Z_pool, emb_template$Z, kern)
  })

  for (b in 1:n_reps) {
    w_subj = rpois(n, 1)
    all_terms = c()

    for (ff in seq_along(obj$folds)) {
      fd = obj$folds[[ff]]
      fr = obj$fold_results[[ff]]
      M_train = length(fd$mesh)

      w_pool_lam = w_subj[fd$lam_idx][fd$i_pool_lam]
      w_pool_gam = w_subj[fd$gam_idx][fd$i_pool_gam]
      w_eval = w_subj[fd$eval_idx]

      # One-step Newton for lambda
      c_lam_b = project_target(fd$stacked_lam$Z_pool, kern,
                               list(Z = fd$stacked_lam$Z_pool,
                                    r = w_pool_lam * fd$stacked_lam$Y_pool),
                               K_cross = fr$lam_fit$K)
      tgt_b = split_target(c_lam_b, nrow(fd$stacked_lam$Z_pool))
      lam_boot = onestep_bregman(fr$lam_fit, tgt_b$target, tgt_b$target_null,
                                  w_pool_lam)

      # Direct from bootstrap lambda
      predict_phi_fn = function(Z)
        predict_phi_alpha(lam_boot$alpha, lam_boot$beta, fr$lam_fit, Z)
      dr_b = compute_direct(predict_phi_fn, lam_disp, estimand, fd$eval_Z,
                             M_train, fd$du_bin, fd$mesh, horizon, discrete)

      # Gamma target from bootstrap lambda
      emb_b = gam_measure(estimand, dr_b$predict_haz, fd$eval_Z,
                           fd$mesh_u, fd$mesh, M_train, fd$du_bin, w_eval = w_eval)
      c_gam_b = split_target(
        project_target(fd$stacked_gam$Z_pool, kern, emb_b,
                       K_cross = K_cross_gam[[ff]]),
        nrow(fd$stacked_gam$Z_pool))

      # Re-solve gamma (warm-start from original fit)
      gam_disp_b = target_scaled_entropy(c_gam_b$target)
      gam_boot = kernel_bregman(fd$stacked_gam$Z_pool, kern,
                                 eta = eta_gam, dispersion = gam_disp_b,
                                 target = c_gam_b$target,
                                 target_null = c_gam_b$target_null,
                                 w = w_pool_gam,
                                 K = fr$gam_fit$K,
                                 alpha0 = fr$gam_fit$alpha,
                                 beta0 = fr$gam_fit$beta)

      # Patch bound objects with bootstrap alphas
      lam_bound_b = fr$lam_bound
      lam_bound_b$alpha = lam_boot$alpha
      lam_bound_b$beta = lam_boot$beta

      gam_bound_b = fr$gam_bound
      gam_bound_b$alpha = gam_boot$alpha
      gam_bound_b$beta = gam_boot$beta

      ct_b = correction_terms(lam_bound_b, gam_bound_b,
                               fd$T_obs_eval, fd$D_eval,
                               fd$mesh, fd$du_bin, fd$eval_Z, fd$n_eval,
                               lam_disp, fr$dchis_gamma,
                               horizon, Q_comp, discrete)

      dr_terms_b = dr_b$direct + ct_b$correction
      all_terms = c(all_terms, w_eval * dr_terms_b)
    }

    sw = sum(w_subj)
    boot_ates[b] = sum(all_terms) / sw
    boot_ses[b] = sd(all_terms[all_terms != 0]) / sqrt(n)
  }

  structure(
    list(boot_ates = boot_ates, boot_ses = boot_ses,
         est = obj$est, se_amle = obj$se_amle),
    class = "boot_result")
}

# boot_result methods
se.boot_result = function(obj, ...) sd(obj$boot_ates)

ci_t = function(obj, ...) UseMethod("ci_t")
ci_t.boot_result = function(obj, level = 0.95, ...) {
  good = !is.na(obj$boot_ates) & !is.nan(obj$boot_ates) & abs(obj$boot_ates) < 100
  if (sum(good) < 10) return(c(NA, NA))
  alpha = 1 - level
  t_stats = (obj$boot_ates[good] - obj$est) / obj$boot_ses[good]
  q = quantile(t_stats, c(alpha / 2, 1 - alpha / 2))
  c(obj$est - q[2] * obj$se_amle, obj$est - q[1] * obj$se_amle)
}

ci_pct = function(obj, ...) UseMethod("ci_pct")
ci_pct.boot_result = function(obj, level = 0.95, ...) {
  good = !is.na(obj$boot_ates) & !is.nan(obj$boot_ates) & abs(obj$boot_ates) < 100
  if (sum(good) < 10) return(c(NA, NA))
  alpha = 1 - level
  unname(quantile(obj$boot_ates[good], c(alpha / 2, 1 - alpha / 2)))
}

n_bad.boot_result = function(obj) {
  sum(is.na(obj$boot_ates) | is.nan(obj$boot_ates) | abs(obj$boot_ates) > 100)
}

get_dgp_fn = function(dgp_name) {
  switch(dgp_name,
    paper_cts = paper_cts_dgp,
    grf_type1 = grf_type1_dgp,
    grf_type2 = grf_type2_dgp,
    stop("Unknown DGP: ", dgp_name))
}
get_estimand = function(est_name) {
  switch(est_name,
    surv_prob = surv_prob_estimand(),
    rmst = rmst_estimand(),
    tsm1 = surv_tsm(1),
    stop("Unknown estimand: ", est_name))
}
make_kern = function(kern_name) {
  no_null = function(Z) matrix(nrow = nrow(atleast_2d(Z)), ncol = 0)
  base = matern_kernel(sigma = 2, nu = 3/2)
  base_normed = matern_kernel(sigma = 2, nu = 3/2, null_basis = no_null)
  # levels = c(0, 1) ensures null_basis has consistent columns for prediction
  # (even when evaluating at a single treatment arm)
  if (kern_name == "pham") direct_product(base, iw = c(1, 2))              # block time + treatment, seminorm
  else if (kern_name == "pham_normed") direct_product(base_normed, iw = c(1, 2))  # block time + treatment, norm
  else if (kern_name == "smooth") direct_product(base, iw = 2, levels = c(0, 1)) # block treatment, smooth over time
  else stop("Unknown kernel: ", kern_name)
}
make_lam_disp = function(name) {
  switch(name, quadratic = quadratic_dispersion(),
         softplus = softplus_dispersion(), entropy = entropy_dispersion())
}

# Build dchis_gamma(Z, phi) for target_scaled_entropy prediction.
# Z = (u, treatment, X...), phi = dual variable.
# Returns r(Z) + sign(r(Z)) * exp(sign(r(Z)) * phi) where r comes from
# the estimand's dot_psi_Z evaluated at the fitted lambda.
#
# dpsi_grid always uses discrete building blocks (M_train mesh bins), even
# for continuous time — the hazard is piecewise-constant on the training mesh.
# So r is always looked up by mapping u to a mesh bin.
#
# For ATE estimands, r is the PER-ARM derivative (dp1 at arm=1, dp0 at arm=0),
# not the combined ATE derivative. This matches the dpsi_grid training stacking.
make_dchis_gamma_ts = function(predict_haz, estimand_obj, eval_Z,
                                M_train, du_bin) {
  n_eval = nrow(atleast_2d(eval_Z))
  eval_hash = apply(eval_Z, 1, paste0, collapse = "|")
  lookup = setNames(seq_len(n_eval), eval_hash)
  has_per_arm = !is.null(estimand_obj$dot_psi_arm_discrete)

  if (has_per_arm) {
    # ATE: precompute per-arm r matrices
    dp1 = estimand_obj$dot_psi_arm_discrete(predict_haz, eval_Z, M_train, du_bin, arm = 1)
    dp0 = estimand_obj$dot_psi_arm_discrete(predict_haz, eval_Z, M_train, du_bin, arm = 0)
    r1_mat = r0_mat = matrix(NA_real_, n_eval, M_train)
    for (k in 1:M_train) { r1_mat[, k] = dp1(k); r0_mat[, k] = dp0(k) }
  } else {
    # TSM: single arm
    dpsi = estimand_obj$dot_psi_Z_discrete(predict_haz, eval_Z, M_train, du_bin)
    r_mat = matrix(NA_real_, n_eval, M_train)
    for (k in 1:M_train) r_mat[, k] = dpsi(k)
  }

  function(Z, phi) {
    u = Z[, 1]
    bin = pmin(pmax(ceiling(u / du_bin), 1L), M_train)
    idx = lookup[apply(Z[, -1, drop = FALSE], 1, paste0, collapse = "|")]
    if (has_per_arm) {
      A = Z[, 2]
      r = ifelse(A == 1, r1_mat[cbind(idx, bin)], r0_mat[cbind(idx, bin)])
    } else {
      r = r_mat[cbind(idx, bin)]
    }
    sigma = sign(r)
    r + sigma * exp(sigma * phi)
  }
}

# Population truth
cts_dgp = get_dgp_fn(cell$dgp_name)(horizon = dgp_def$horizon)
set.seed(42)
big = cts_dgp$generate(10000, p = dgp_def$p)
truths = list(surv_prob = big$psi1_true - big$psi0_true,
              rmst = big$rmst1_true - big$rmst0_true,
              tsm1 = big$psi1_true)

cat(sprintf("=== Cell %d: %s ===\n", cell_id, cell_name))
cat(sprintf("DGP=%s, res=%s, n=%d, p=%d\n",
            cell$dgp_name, if (is_cts) "cts" else as.integer(cell$res),
            n_obs, dgp_def$p))
cat(sprintf("Balancing methods: %s\n", paste(active_bal, collapse = ", ")))
cat(sprintf("External methods: %s\n", paste(active_ext, collapse = ", ")))
cat(sprintf("Estimands: %s\n", paste(estimand_names, collapse = ", ")))
for (en in estimand_names) cat(sprintf("  truth[%s] = %.5f\n", en, truths[[en]]))
cat(sprintf("n_reps=%d, n=%d, boot=%d, workers=%d\n\n", n_reps, n_obs, boot_reps, n_workers))

if (n_workers > 1) plan(multisession, workers = n_workers)

log_dir = sprintf("code/spinoff3/coverage-logs-%s-n%d-%s",
                  cell_name, n_obs, format(Sys.time(), "%Y%m%d-%H%M"))
dir.create(log_dir, showWarnings = FALSE)
cat(sprintf("Logs: %s/\n", log_dir))

make_logger = function(rep_id) {
  f = file.path(log_dir, sprintf("rep-%03d.log", rep_id))
  function(...) cat(sprintf(...), "\n", file = f, append = TRUE)
}

# ============================================================
# One rep
# ============================================================
one_rep = function(rep_id) {
  # Re-source for worker processes (multisession)
  source("R/kernel.R")
  for (f in c("R/estimands.R", "R/survival.R", "R/tuning.R", "R/dgp.R"))
    source(f)

  log = make_logger(rep_id)
  log("rep %d started (cell %s)", rep_id, cell_name)

  # Generate continuous data (same seed regardless of resolution)
  cts_dgp = get_dgp_fn(cell$dgp_name)(horizon = dgp_def$horizon)
  set.seed(rep_id)
  cts_dat = cts_dgp$generate(n_obs, p = dgp_def$p)
  n = n_obs

  # Discretize if needed
  disc_dat = NULL
  if (!is_cts) {
    disc_dgp_obj = discretize(cts_dgp, as.integer(cell$res))
    disc_dat = disc_dgp_obj$generate_from(cts_dat)
  }

  sigma2_grid = 2^(-6:6)
  m_grid = length(sigma2_grid)

  na_row = function(mname, est_name) data.frame(
    method = mname, estimand = est_name, rep = rep_id,
    resolution = if (is_cts) "cts" else as.character(as.integer(cell$res)),
    est = NA, bias = NA,
    selected_sigma2_lam = NA, selected_sigma2_gam = NA,
    tuning = NA,
    covers_eif = NA, covers_amle = NA,
    covers_boot_t = NA, covers_pct = NA,
    se_eif = NA, se_amle = NA, se_boot = NA,
    ci_width_eif = NA, ci_width_amle = NA,
    ci_width_boot_t = NA, ci_width_pct = NA,
    n_boot_bad = NA)

  all_coverage = list()
  all_paths = list()

  # ==========================================================
  # External methods (no grid sweep, no bootstrap)
  # ==========================================================

  if ("GRF" %in% active_ext) {
    for (en in setdiff(estimand_names, "tsm1")) {
      grf_res = tryCatch({
        if (!requireNamespace("grf", quietly = TRUE)) stop("grf not installed")
        Y_grf = pmin(cts_dat$T_obs, dgp_def$horizon)
        D_grf = cts_dat$D; D_grf[cts_dat$T_obs > dgp_def$horizon] = 0
        h = dgp_def$horizon - 0.01
        grf_target = if (en == "rmst") "RMST" else "survival.probability"
        csf = grf::causal_survival_forest(
          X = cts_dat$X, Y = Y_grf, W = cts_dat$A, D = D_grf,
          target = grf_target, horizon = h,
          W.hat = cts_dat$pi_x, num.trees = 2000)
        ate = grf::average_treatment_effect(csf)
        list(est = ate[1], se = ate[2])
      }, error = function(e) { log("GRF %s error: %s", en, e$message); NULL })

      if (!is.null(grf_res)) {
        truth = truths[[en]]
        ci = grf_res$est + c(-1, 1) * 1.96 * grf_res$se
        row = na_row("GRF", en)
        row$est = grf_res$est; row$bias = grf_res$est - truth
        row$tuning = "none"
        row$se_eif = grf_res$se
        row$covers_eif = truth >= ci[1] & truth <= ci[2]
        row$ci_width_eif = diff(ci)
        all_coverage[[length(all_coverage) + 1]] = row
      }
    }
    log("GRF done")
  }

  if ("1S" %in% active_ext) {
    os_dir = if (dir.exists("/Users/skip/work/spinoffs/code/spinoff3"))
      "/Users/skip/work/spinoffs/code/spinoff3" else
      path.expand("~/work/spinoffs/code/spinoff3")
    source(file.path(os_dir, "one-step-estimator.R"))

    for (en in setdiff(estimand_names, "rmst")) {
      os_estimand = if (en == "tsm1") "tsm1" else "surv"
      os_res = tryCatch(
        one_step_surv(disc_dat$T_obs, disc_dat$D, disc_dat$A, disc_dat$X,
                      n_steps = as.integer(cell$res),
                      horizon = dgp_def$horizon,
                      estimand_name = os_estimand),
        error = function(e) { log("1S %s error: %s", en, e$message); NULL })

      if (!is.null(os_res) && !is.na(os_res$est)) {
        truth = truths[[en]]
        se = if (!is.null(os_res$se)) os_res$se else NA
        row = na_row("1S", en)
        row$est = os_res$est; row$bias = os_res$est - truth
        row$tuning = "none"
        if (!is.na(se)) {
          ci = os_res$est + c(-1, 1) * 1.96 * se
          row$se_eif = se
          row$covers_eif = truth >= ci[1] & truth <= ci[2]
          row$ci_width_eif = diff(ci)
        }
        all_coverage[[length(all_coverage) + 1]] = row
      }
    }
    log("1S done")
  }

  # ==========================================================
  # Balancing methods (grid sweep + bootstrap)
  # ==========================================================

  for (mname in active_bal) {
    t_method = proc.time()
    mdef = bal_methods[[mname]]
    kern = make_kern(mdef$kern)
    lam_disp = make_lam_disp(mdef$lam_disp)
    gam_disp = entropy_dispersion()
    time_type = mdef$time

    # Data for this method
    if (time_type == "discrete") {
      T_obs = disc_dat$T_obs; D = disc_dat$D; A = disc_dat$A; X = disc_dat$X
      M_tr = as.integer(cell$res)
    } else {
      T_obs = cts_dat$T_obs; D = cts_dat$D; A = cts_dat$A; X = cts_dat$X
      M_tr = 15
    }

    log("--- method: %s (kern=%s, disp=%s, time=%s) ---",
        mname, mdef$kern, mdef$lam_disp, time_type)

    discrete = (time_type == "discrete")
    Z = cbind(A, X)
    n_folds = mdef$n_folds
    horizon = dgp_def$horizon
    Q_comp = 50

    # Fold setup: build training meshes, kernel matrices, targets
    fold_id = (seq_len(n) %% n_folds) + 1L
    fold_data = lapply(1:n_folds, function(ff) {
      if (n_folds == 2) {
        lam_idx = which(fold_id != ff); gam_idx = which(fold_id == ff); eval_idx = gam_idx
      } else {
        lam_idx = which(fold_id == ff)
        gam_idx = which(fold_id == (ff %% n_folds) + 1)
        eval_idx = which(fold_id == ((ff + 1) %% n_folds) + 1)
      }

      lam_train = build_training_mesh(T_obs[lam_idx], D[lam_idx],
                                       Z[lam_idx, , drop = FALSE], horizon, M_tr)
      stacked_lam = at_risk_pairs(lam_train)

      gam_train = build_training_mesh(T_obs[gam_idx], D[gam_idx],
                                       Z[gam_idx, , drop = FALSE], horizon, M_tr)
      stacked_gam = at_risk_pairs(gam_train)

      K_lam = kernel_matrix(stacked_lam$Z_pool, stacked_lam$Z_pool, kern)
      c_lam = project_target(stacked_lam$Z_pool, kern,
                             list(Z = stacked_lam$Z_pool, r = stacked_lam$Y_pool),
                             K_cross = K_lam)
      tgt_lam = split_target(c_lam, nrow(stacked_lam$Z_pool))

      eval_Z = Z[eval_idx, , drop = FALSE]

      list(lam_idx = lam_idx, gam_idx = gam_idx, eval_idx = eval_idx,
           stacked_lam = stacked_lam, stacked_gam = stacked_gam,
           mesh = lam_train$mesh, du_bin = stacked_lam$du_bin,
           mesh_u = unique(stacked_gam$u_pool),
           K_lam = K_lam, tgt_lam = tgt_lam,
           eval_Z = eval_Z,
           T_obs_eval = T_obs[eval_idx], D_eval = D[eval_idx],
           n_eval = length(eval_idx),
           i_pool_lam = stacked_lam$i_pool, i_pool_gam = stacked_gam$i_pool)
    })

    # --- Lambda grid sweep (estimand-independent) ---
    t_lam = proc.time()
    lam_cache = vector("list", m_grid)
    cv_lam = rep(NA_real_, m_grid)
    lam_ok_vec = logical(m_grid)

    prev_lam_alpha = vector("list", n_folds)
    prev_lam_beta = vector("list", n_folds)

    sweep_order = order(sigma2_grid, decreasing = TRUE)
    for (jl in sweep_order) {
      eta_lam = sigma2_grid[jl] / n
      fold_fits = vector("list", n_folds)

      for (ff in 1:n_folds) {
        fd = fold_data[[ff]]
        lam_fit = kernel_bregman(fd$stacked_lam$Z_pool, kern,
                         eta = eta_lam, dispersion = lam_disp,
                         target = fd$tgt_lam$target,
                         target_null = fd$tgt_lam$target_null,
                         K = fd$K_lam,
                         alpha0 = prev_lam_alpha[[ff]],
                         beta0 = prev_lam_beta[[ff]])

        lam_bound = bind_subjects(lam_fit, fd$eval_Z)
        fold_fits[[ff]] = list(lam_fit = lam_fit, lam_bound = lam_bound)
        prev_lam_alpha[[ff]] = lam_fit$alpha
        prev_lam_beta[[ff]] = lam_fit$beta
      }

      lam_cache[[jl]] = fold_fits
      lam_ok_vec[jl] = TRUE

      cv_lam[jl] = mean(sapply(1:n_folds, function(ff) {
        fd = fold_data[[ff]]
        cv_dual_loss(fold_fits[[ff]]$lam_fit,
                     fd$stacked_gam$Z_pool,
                     fd$stacked_gam$Y_pool)
      }))
    }
    log("%s lambda sweep: %.1fs", mname, (proc.time() - t_lam)[3])

    # --- Per-estimand: gamma sweep, selection, bootstrap ---
    method_ests = setdiff(estimand_names, mdef$skip_estimands)
    for (en in method_ests) {
      t_est = proc.time()
      estimand_obj = get_estimand(en)
      truth = truths[[en]]
      grid_est = array(NA_real_, c(m_grid, m_grid))
      grid_dir_val = array(NA_real_, c(m_grid, m_grid))
      grid_dir_se = array(NA_real_, c(m_grid, m_grid))
      grid_dr_se = array(NA_real_, c(m_grid, m_grid))
      grid_cv_lam = array(NA_real_, c(m_grid, m_grid))
      grid_cv_gam = array(NA_real_, c(m_grid, m_grid))

      for (jl in 1:m_grid) {
        if (!lam_ok_vec[jl]) next
        fold_fits = lam_cache[[jl]]

        # Direct estimate per fold (via compute_direct)
        dir_per_fold = lapply(1:n_folds, function(ff) {
          fd = fold_data[[ff]]
          lf = fold_fits[[ff]]$lam_fit
          predict_phi_fn = function(Z) predict_phi(lf, Z)
          compute_direct(predict_phi_fn, lam_disp, estimand_obj, fd$eval_Z,
                         length(fd$mesh), fd$du_bin, fd$mesh, horizon, discrete)
        })
        dir_vals = unlist(lapply(dir_per_fold, `[[`, "direct"))
        dir_mean = mean(dir_vals)
        dir_se = sd(dir_vals) / sqrt(n)

        # Gamma measures and targets for training and test folds
        predict_haz_folds = lapply(1:n_folds, function(ff) {
          fd = fold_data[[ff]]
          lf = fold_fits[[ff]]$lam_fit
          function(k, Z_ev) lam_disp$dchis(predict_phi(lf, cbind(fd$mesh[k], Z_ev)))
        })

        emb_gam_train = lapply(1:n_folds, function(ff) {
          fd = fold_data[[ff]]
          gam_measure(estimand_obj, predict_haz_folds[[ff]], fd$eval_Z,
                      fd$mesh_u, fd$mesh, length(fd$mesh), fd$du_bin)
        })
        c_gam_train = lapply(1:n_folds, function(ff) {
          fd = fold_data[[ff]]
          split_target(project_target(fd$stacked_gam$Z_pool, kern,
                                      emb_gam_train[[ff]]),
                       nrow(fd$stacked_gam$Z_pool))
        })

        emb_gam_test = lapply(1:n_folds, function(ff) {
          fd = fold_data[[ff]]
          gam_measure(estimand_obj, predict_haz_folds[[ff]], fd$eval_Z,
                      fd$mesh_u, fd$mesh, length(fd$mesh), fd$du_bin)
        })
        c_gam_test = lapply(1:n_folds, function(ff) {
          fd = fold_data[[ff]]
          split_target(project_target(fd$stacked_lam$Z_pool, kern,
                                      emb_gam_test[[ff]]),
                       nrow(fd$stacked_lam$Z_pool))
        })

        gam_disp_test = lapply(1:n_folds, function(ff)
          target_scaled_entropy(c_gam_test[[ff]]$target))
        gam_disp_train = lapply(1:n_folds, function(ff)
          target_scaled_entropy(c_gam_train[[ff]]$target))

        dchis_gamma_folds = lapply(1:n_folds, function(ff) {
          fd = fold_data[[ff]]
          make_dchis_gamma_ts(predict_haz_folds[[ff]], estimand_obj, fd$eval_Z,
                              length(fd$mesh), fd$du_bin)
        })

        prev_gam_alpha = vector("list", n_folds)
        prev_gam_beta = vector("list", n_folds)

        for (jg in sweep_order) {
          eta_gam = sigma2_grid[jg] / n
          all_terms = c(); all_direct = c(); all_noise_var = c()
          gam_fits = vector("list", n_folds)

          for (ff in 1:n_folds) {
            fd = fold_data[[ff]]
            # Fit gamma directly
            gam_fit = kernel_bregman(fd$stacked_gam$Z_pool, kern,
                             eta = eta_gam, dispersion = gam_disp_train[[ff]],
                             target = c_gam_train[[ff]]$target,
                             target_null = c_gam_train[[ff]]$target_null,
                             alpha0 = prev_gam_alpha[[ff]],
                             beta0 = prev_gam_beta[[ff]])
            gam_bound = bind_subjects(gam_fit, fd$eval_Z)

            # DR correction directly
            ct = correction_terms(fold_fits[[ff]]$lam_bound, gam_bound,
                                   fd$T_obs_eval, fd$D_eval,
                                   fd$mesh, fd$du_bin, fd$eval_Z, fd$n_eval,
                                   lam_disp, dchis_gamma_folds[[ff]],
                                   horizon, Q_comp, discrete)

            all_terms = c(all_terms, dir_per_fold[[ff]]$direct + ct$correction)
            all_direct = c(all_direct, dir_per_fold[[ff]]$direct)
            all_noise_var = c(all_noise_var, ct$noise_var)
            gam_fits[[ff]] = gam_fit
            prev_gam_alpha[[ff]] = gam_fit$alpha
            prev_gam_beta[[ff]] = gam_fit$beta
          }
          # Sanity check: if gamma phi is extreme, the solve is degenerate
          gam_degenerate = any(sapply(1:n_folds, function(ff) {
            gf = gam_fits[[ff]]
            phi = as.vector(gf$K %*% gf$alpha)
            if (ncol(gf$B) > 0) phi = phi + as.vector(gf$B %*% gf$beta)
            max(abs(phi)) > 50
          }))
          if (gam_degenerate) break

          cv_gam_vals = sapply(1:n_folds, function(ff)
            cv_dual_loss(gam_fits[[ff]],
                         fold_data[[ff]]$stacked_lam$Z_pool,
                         c_gam_test[[ff]]$target,
                         dispersion_test = gam_disp_test[[ff]]))

          grid_est[jl, jg] = mean(all_terms)
          grid_dir_val[jl, jg] = dir_mean
          grid_dir_se[jl, jg] = dir_se
          grid_dr_se[jl, jg] = sd(all_terms) / sqrt(n)
          grid_cv_lam[jl, jg] = cv_lam[jl]
          grid_cv_gam[jl, jg] = mean(cv_gam_vals)
        }
      }

      # Build path
      path = data.frame(
        method = mname, estimand = en,
        resolution = if (is_cts) "cts" else as.character(as.integer(cell$res)),
        sigma2_lam = rep(sigma2_grid, each = m_grid),
        sigma2_gam = rep(sigma2_grid, times = m_grid),
        dir_val = as.vector(t(grid_dir_val)),
        dir_se = as.vector(t(grid_dir_se)),
        dr_est = as.vector(t(grid_est)),
        dr_se = as.vector(t(grid_dr_se)),
        loo_lam = as.vector(t(grid_cv_lam)),
        loo_gam = as.vector(t(grid_cv_gam)),
        rep = rep_id)
      all_paths[[length(all_paths) + 1]] = path

      # --- Selection strategies ---
      jg_mid = ceiling(m_grid / 2)
      lam_slice = path[path$sigma2_gam == sigma2_grid[jg_mid], ]
      lam_slice = lam_slice[order(lam_slice$sigma2_lam), ]

      select_lam = function(tuning) {
        ctx = tuning_context(sigma2_grid / n, n = n,
          dir_val = function() lam_slice$dir_val,
          dir_se = function() lam_slice$dir_se,
          loo_scores = function() lam_slice$loo_lam)
        tuning$select(ctx)
      }

      select_gam = function(sel_lam_val) {
        gam_slice = path[path$sigma2_lam == sel_lam_val, ]
        gam_slice = gam_slice[order(gam_slice$sigma2_gam), ]
        idx = which.min(gam_slice$loo_gam)
        if (length(idx) == 0) return(NA_real_)
        sigma2_grid[idx]
      }

      strategies = list(
        list(label = "fixed",  tuning = fixed_tuning(1 / n)),
        list(label = "lepski", tuning = lepski_tuning()),
        list(label = "cv",     tuning = cv_tuning()))

      # Estimate + bootstrap at a selected (s2_lam, s2_gam).
      # Uses crossfit to get a dr_result, then bootstrap for CIs.
      estimate_at = function(s2_lam, s2_gam) {
        eta_l = s2_lam / n; eta_g = s2_gam / n
        crossfit(n, fold_data, function(fd) {
          lam_fit = kernel_bregman(fd$stacked_lam$Z_pool, kern,
                           eta = eta_l, dispersion = lam_disp,
                           target = fd$tgt_lam$target,
                           target_null = fd$tgt_lam$target_null,
                           K = fd$K_lam)
          # Build dchis_gamma for this fold
          predict_haz_fd = function(k, Z_ev)
            lam_disp$dchis(predict_phi(lam_fit, cbind(fd$mesh[k], Z_ev)))
          dchis_gam = make_dchis_gamma_ts(predict_haz_fd, estimand_obj,
                                           fd$eval_Z, length(fd$mesh), fd$du_bin)
          dr_fold(fd, lam_fit, kern, eta_g, NULL, dchis_gam,
                  lam_disp, estimand_obj, horizon, Q_comp, discrete)
        }, kern = kern, eta_gam = eta_g, lam_disp = lam_disp,
           gam_disp = gam_disp, estimand = estimand_obj,
           horizon = horizon, Q_comp = Q_comp, discrete = discrete)
      }

      # Build a coverage row from a dr_result + optional boot_result.
      coverage_row = function(result, label, s2_lam, s2_gam, boot = NULL) {
        row = na_row(mname, en)
        row$est = result$est; row$bias = result$est - truth
        row$selected_sigma2_lam = s2_lam; row$selected_sigma2_gam = s2_gam
        row$tuning = label
        row$se_eif = result$se; row$se_amle = result$se_amle

        ci_eif = result$est + c(-1, 1) * 1.96 * result$se
        ci_amle = result$est + c(-1, 1) * 1.96 * result$se_amle
        row$covers_eif = truth >= ci_eif[1] & truth <= ci_eif[2]
        row$covers_amle = truth >= ci_amle[1] & truth <= ci_amle[2]
        row$ci_width_eif = diff(ci_eif); row$ci_width_amle = diff(ci_amle)

        if (!is.null(boot)) {
          nb = n_bad.boot_result(boot)
          if (nb > 0)
            log("%s %s %s: %d/%d bootstrap reps degenerate",
                mname, en, label, nb, length(boot$boot_ates))
          row$se_boot = se.boot_result(boot)
          row$n_boot_bad = nb

          bt = ci_t(boot)
          if (!any(is.na(bt))) {
            row$covers_boot_t = truth >= bt[1] & truth <= bt[2]
            row$ci_width_boot_t = diff(bt)
          }
          bp = ci_pct(boot)
          if (!any(is.na(bp))) {
            row$covers_pct = truth >= bp[1] & truth <= bp[2]
            row$ci_width_pct = diff(bp)
          }
        }
        row
      }

      for (strat in strategies) {
        label = strat$label
        sel = select_lam(strat$tuning)
        sel_sigma2_lam = sel$eta * n
        sel_sigma2_gam = select_gam(sel_sigma2_lam)
        if (is.na(sel_sigma2_gam)) {
          log("%s %s %s: no valid gamma at sigma2_lam=%.3g", mname, en, label, sel_sigma2_lam)
          all_coverage[[length(all_coverage) + 1]] = na_row(mname, en)
          next
        }
        log("%s %s %s selected: sigma2_lam=%.3g, sigma2_gam=%.3g",
            mname, en, label, sel_sigma2_lam, sel_sigma2_gam)

        result = estimate_at(sel_sigma2_lam, sel_sigma2_gam)
        boot = if (boot_reps > 0) bootstrap(result, boot_reps) else NULL
        all_coverage[[length(all_coverage) + 1]] = coverage_row(result, label,
                                                                 sel_sigma2_lam, sel_sigma2_gam, boot)

        if (label == "lepski") {
          sel_idx = which(sigma2_grid == sel_sigma2_lam)
          us_idx = sel_idx - 2
          if (us_idx >= 1) {
            us_sigma2_gam = select_gam(sigma2_grid[us_idx])
            result2 = estimate_at(sigma2_grid[us_idx], us_sigma2_gam)
            boot2 = if (boot_reps > 0) bootstrap(result2, boot_reps) else NULL
            all_coverage[[length(all_coverage) + 1]] = coverage_row(result2, "lepski-2",
                                                                     sigma2_grid[us_idx], us_sigma2_gam, boot2)
          }
        }
      }
      log("%s %s: %.1fs", mname, en, (proc.time() - t_est)[3])
    }
    log("%s total: %.1fs", mname, (proc.time() - t_method)[3])
  }

  list(coverage = do.call(rbind, all_coverage),
       paths = if (length(all_paths) > 0) do.call(rbind, all_paths) else NULL)
}

# ============================================================
# Run
# ============================================================
cat("Running...\n")
t0 = proc.time()

rep_ids = start_rep:(start_rep + n_reps - 1)

if (n_workers > 1) {
  all_results = future_map(rep_ids, one_rep,
                           .options = furrr_options(seed = TRUE),
                           .progress = TRUE)
} else {
  all_results = lapply(rep_ids, function(r) {
    set.seed(r)
    cat(sprintf("  rep %d/%d\n", r - start_rep + 1, n_reps))
    one_rep(r)
  })
}

elapsed = (proc.time() - t0)[3]
df = do.call(rbind, lapply(all_results, `[[`, "coverage"))
all_paths = do.call(rbind, Filter(Negate(is.null), lapply(all_results, `[[`, "paths")))

cat(sprintf("\nDone in %.1f seconds (%.1f s/rep)\n\n", elapsed, elapsed / n_reps))

# ============================================================
# Summary
# ============================================================
for (en in estimand_names) {
  cat(sprintf("\n========== %s (truth=%.5f) ==========\n", en, truths[[en]]))
  sub_en = df[df$estimand == en, ]
  for (m in unique(sub_en$method)) {
    sub_m = sub_en[sub_en$method == m, ]
    for (tn in unique(sub_m$tuning)) {
      sub = sub_m[sub_m$tuning == tn & !is.na(sub_m$est), ]
      if (nrow(sub) == 0) next
      cat(sprintf("  [%s/%s] (%d reps)\n", m, tn, nrow(sub)))
      cat(sprintf("    bias: %.4f (RMSE: %.4f)\n", mean(sub$bias), sqrt(mean(sub$bias^2))))
      if (!all(is.na(sub$covers_eif)))
        cat(sprintf("    cov EIF: %.2f  (width: %.4f)\n",
                    mean(sub$covers_eif, na.rm = TRUE), mean(sub$ci_width_eif, na.rm = TRUE)))
      if (!all(is.na(sub$covers_boot_t)))
        cat(sprintf("    cov boot-t: %.2f  (width: %.4f)\n",
                    mean(sub$covers_boot_t, na.rm = TRUE), mean(sub$ci_width_boot_t, na.rm = TRUE)))
    }
  }
}

# Save
out_file = sprintf("code/spinoff3/coverage-results-%s-n%d-%s.rds",
                   cell_name, n_obs, format(Sys.time(), "%Y%m%d-%H%M"))
saveRDS(list(results = df, paths = all_paths,
             params = list(n_obs = n_obs, n_reps = n_reps,
                           boot_reps = boot_reps, truths = truths,
                           cell = cell, dgp_def = dgp_def)), out_file)
cat(sprintf("Saved: %s\n", out_file))
