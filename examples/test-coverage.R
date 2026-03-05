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
  paper_cts     = list(dgp = "paper_cts",     horizon = 1,   p = 2, valid_res = c(2, 4, 8, 16, Inf)),
  grf_type1     = list(dgp = "grf_type1",     horizon = 0.8, p = 5, valid_res = c(2, 4, 8, 16, Inf)),
  grf_type2     = list(dgp = "grf_type2",     horizon = 1.2, p = 5, valid_res = c(2, 4, 8, 16, Inf)),
  pham_discrete = list(dgp = "pham_discrete", horizon = 8,   p = 2, valid_res = c(2, 4, 8, 16),
                       native_discrete = TRUE)
)
resolutions = c(2, 4, 8, 16, Inf)  # Inf = continuous

cells = expand.grid(dgp_name = names(dgp_defs), res = resolutions,
                    stringsAsFactors = FALSE)
cells = cells[mapply(function(d, r) r %in% dgp_defs[[d]]$valid_res,
                     cells$dgp_name, cells$res), ]
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
# lambda_disp: dispersion for lambda (hazard model)
# time: "discrete" or "continuous"
bal_methods = list(
  Pham_quad     = list(kern = "pham",        lambda_disp = "quadratic", time = "discrete", n_folds = 2, res = c(2, 4, 8, 16)),
  Pham_normed   = list(kern = "pham_normed", lambda_disp = "quadratic", time = "discrete", n_folds = 2, res = c(2, 4, 8, 16)),
  Pham_logistic = list(kern = "pham",        lambda_disp = "softplus",  time = "discrete", n_folds = 2, res = c(2, 4, 8, 16)),
  KBW_disc      = list(kern = "smooth", lambda_disp = "softplus",  time = "discrete", n_folds = 2, res = c(2, 4, 8, 16)),
  KBW_cts       = list(kern = "smooth", lambda_disp = "entropy",   time = "continuous", n_folds = 2, res = c(Inf)),
  PI            = list(kern = "pham",   lambda_disp = "quadratic", time = "discrete", n_folds = 3, res = c(2, 4, 8, 16),
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
source("R/all.R")
for (f in c("R/tuning.R", "R/dgp.R"))
  source(f)

# Helpers: gamma_measure, compute_direct, correction_terms, build_gamma_correction,
# dr_fold/crossfit/bootstrap now come from R/surv_estimate.R.


get_dgp_fn = function(dgp_name) {
  switch(dgp_name,
    paper_cts     = paper_cts_dgp,
    grf_type1     = grf_type1_dgp,
    grf_type2     = grf_type2_dgp,
    pham_discrete = function(horizon, ...) pham_discrete_dgp(t_max = 8),
    stop("Unknown DGP: ", dgp_name))
}
get_estimand = function(est_name) {
  switch(est_name,
    surv_prob = survival_probability_ate(),
    rmst = rmst_ate(),
    tsm1 = survival_probability_tsm(1),
    stop("Unknown estimand: ", est_name))
}
make_kern = function(kern_name, mesh = NULL) {
  no_null = function(Z) matrix(nrow = nrow(atleast_2d(Z)), ncol = 0)
  base = matern_kernel(sigma = 2, nu = 3/2)
  base_normed = matern_kernel(sigma = 2, nu = 3/2, null_basis = no_null)
  # levels defines the direct product decomposition of the RKHS.
  # Must be explicit for seminorm kernels — otherwise null_basis dimensions
  # depend on which (u, W) combos appear in the data, breaking prediction.
  if (kern_name == "pham") direct_product(base, iw = c(1, 2), levels = list(mesh, c(0, 1)))
  else if (kern_name == "pham_normed") direct_product(base_normed, iw = c(1, 2), levels = list(mesh, c(0, 1)))
  else if (kern_name == "smooth") direct_product(base, iw = 2, levels = c(0, 1))
  else stop("Unknown kernel: ", kern_name)
}
make_lambda_disp = function(name) {
  switch(name, quadratic = quadratic_dispersion(),
         softplus = softplus_dispersion(), entropy = entropy_dispersion())
}

# Population truth — discretized at the estimator's resolution
cts_dgp = get_dgp_fn(cell$dgp_name)(horizon = dgp_def$horizon)
set.seed(42)
big = cts_dgp$generate(10000, p = dgp_def$p)
if (is_cts) {
  truths = list(surv_prob = big$psi1_true - big$psi0_true,
                rmst = big$rmst1_true - big$rmst0_true,
                tsm1 = big$psi1_true)
} else {
  disc_dgp = discretize(cts_dgp, n_steps = as.integer(cell$res))
  disc_truth = compute_truth(disc_dgp$hazard, big$Z, horizon, "discrete", disc_dgp$grid)
  truths = list(surv_prob = disc_truth$psi1 - disc_truth$psi0,
                rmst = disc_truth$rmst1 - disc_truth$rmst0,
                tsm1 = disc_truth$psi1)
}

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

log_dir = sprintf("examples/coverage-logs-%s-n%d-%s",
                  cell_name, n_obs, format(Sys.time(), "%Y%m%d-%H%M"))
dir.create(log_dir, showWarnings = FALSE)
cat(sprintf("Logs: %s/\n", log_dir))


# ============================================================
# One rep
# ============================================================
one_rep = function(rep_id) {
  # Define before re-sourcing: source() writes to globalenv and would overwrite
  # make_logger with utils.R's list-returning version, breaking findFun.
  log_file = file.path(log_dir, sprintf("rep-%03d.log", rep_id))
  logger = function(...) cat(sprintf(...), "\n", file = log_file, append = TRUE)

  # Re-source for worker processes (multisession)
  source("R/kernel.R")
  for (f in c("R/estimands.R", "R/survival.R", "R/tuning.R", "R/dgp.R"))
    source(f)

  logger("rep %d started (cell %s)", rep_id, cell_name)

  dgp_obj = get_dgp_fn(cell$dgp_name)(horizon = dgp_def$horizon)
  set.seed(rep_id)
  n = n_obs

  if (isTRUE(dgp_def$native_discrete)) {
    # Native discrete DGP: generate directly, no continuous→discrete step
    disc_dat = dgp_obj$generate(n_obs, p = dgp_def$p)
    cts_dat = disc_dat
  } else {
    # Continuous DGP: generate continuous, discretize if needed
    cts_dat = dgp_obj$generate(n_obs, p = dgp_def$p)
    disc_dat = NULL
    if (!is_cts) {
      disc_dgp_obj = discretize(dgp_obj, as.integer(cell$res))
      disc_dat = disc_dgp_obj$generate_from(cts_dat)
    }
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
  all_boot = list()

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
      }, error = function(e) { logger("GRF %s error: %s", en, e$message); NULL })

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
    logger("GRF done")
  }

  if ("1S" %in% active_ext) {
    os_dir = if (dir.exists("/Users/skip/work/spinoffs/code/spinoff3"))
      "/Users/skip/work/spinoffs/code/spinoff3" else
      path.expand("~/work/spinoffs/code/spinoff3")
    source(file.path(os_dir, "one-step-estimator.R"))

    for (en in setdiff(estimand_names, "rmst")) {
      os_estimand = if (en == "tsm1") "tsm1" else "surv_prob"
      os_res = tryCatch(
        one_step_surv(disc_dat$T_obs, disc_dat$D, disc_dat$A, disc_dat$X,
                      n_steps = as.integer(cell$res),
                      horizon = dgp_def$horizon,
                      estimand_name = os_estimand),
        error = function(e) { logger("1S %s error: %s", en, e$message); NULL })

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
    logger("1S done")
  }

  # ==========================================================
  # Balancing methods (grid sweep + bootstrap)
  # ==========================================================

  for (mname in active_bal) {
    t_method = proc.time()
    mdef = bal_methods[[mname]]
    lambda_disp = make_lambda_disp(mdef$lambda_disp)
    gamma_disp_fn = if (!is.null(mdef$gamma_disp_fn)) mdef$gamma_disp_fn else target_scaled_entropy
    time_type = mdef$time

    # Data for this method
    if (time_type == "discrete") {
      T_obs = disc_dat$T_obs; D = disc_dat$D; A = disc_dat$A; X = disc_dat$X
      M_tr = as.integer(cell$res)
    } else {
      T_obs = cts_dat$T_obs; D = cts_dat$D; A = cts_dat$A; X = cts_dat$X
      M_tr = 15
    }

    discrete = (time_type == "discrete")
    horizon = dgp_def$horizon
    mesh = ((1:M_tr) - 0.5) * (horizon / M_tr)
    kern = make_kern(mdef$kern, mesh = mesh)

    logger("--- method: %s (kern=%s, disp=%s, time=%s) ---",
        mname, mdef$kern, mdef$lambda_disp, time_type)

    Z = cbind(A, X)
    n_folds = mdef$n_folds
    Q_comp = 50
    observations = list(T_obs = T_obs, D = D, Z = Z)
    mps = mesh_project(observations, horizon, M_tr, n_folds, kern)

    # --- Lambda grid sweep (estimand-independent) ---
    t_lam = proc.time()
    lambda_cache = vector("list", m_grid)
    cv_lambda = rep(NA_real_, m_grid)
    lambda_ok_vec = logical(m_grid)

    prev_lambda_alpha = vector("list", n_folds)
    prev_lambda_beta = vector("list", n_folds)

    sweep_order = order(sigma2_grid, decreasing = TRUE)
    for (jl in sweep_order) {
      eta_lam = sigma2_grid[jl] / n
      fold_fits = vector("list", n_folds)

      for (ff in 1:n_folds) {
        mp = mps[[ff]]
        warm = if (!is.null(prev_lambda_alpha[[ff]]))
          list(alpha = prev_lambda_alpha[[ff]], beta = prev_lambda_beta[[ff]]) else NULL
        lambda_fit = fit_lambda(mp, kern, eta_lam, lambda_disp, warm = warm)
        lambda_fit_ = .bind(lambda_fit, mp$eval_Z)
        fold_fits[[ff]] = list(lambda_fit = lambda_fit, lambda_fit_ = lambda_fit_)
        prev_lambda_alpha[[ff]] = lambda_fit$alpha
        prev_lambda_beta[[ff]] = lambda_fit$beta
      }

      lambda_cache[[jl]] = fold_fits
      lambda_ok_vec[jl] = TRUE

      cv_lambda[jl] = mean(sapply(1:n_folds, function(ff) {
        mp = mps[[ff]]
        cv_dual_loss(fold_fits[[ff]]$lambda_fit,
                     mp$stacked_gamma$Z_pool,
                     mp$stacked_gamma$Y_pool)
      }))
    }
    logger("%s lambda sweep: %.1fs", mname, (proc.time() - t_lam)[3])

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
      grid_cv_lambda = array(NA_real_, c(m_grid, m_grid))
      grid_cv_gamma = array(NA_real_, c(m_grid, m_grid))

      for (jl in 1:m_grid) {
        if (!lambda_ok_vec[jl]) next
        fold_fits = lambda_cache[[jl]]

        # Direct estimate per fold (via compute_direct)
        dir_per_fold = lapply(1:n_folds, function(ff) {
          mp = mps[[ff]]
          lf = fold_fits[[ff]]$lambda_fit
          phi = function(Z) .phi(lf, Z)
          compute_direct(phi, lambda_disp, estimand_obj, mp, discrete)
        })
        dir_vals = unlist(lapply(dir_per_fold, `[[`, "direct"))
        dir_mean = mean(dir_vals)
        dir_se = sd(dir_vals) / sqrt(n)

        # Gamma measures and targets for training and test folds
        predict_haz_folds = lapply(1:n_folds, function(ff) {
          mp = mps[[ff]]
          lf = fold_fits[[ff]]$lambda_fit
          phi   = \(Z) .phi(lf, Z)
          dchis = \(p) .dchis(lambda_disp, p)
          function(k, Z_ev) dchis(phi(cbind(mp$mesh[k], Z_ev)))
        })

        emb_gamma_train = lapply(1:n_folds, function(ff)
          gamma_measure(estimand_obj, predict_haz_folds[[ff]], mps[[ff]]))
        c_gamma_train = lapply(1:n_folds, function(ff) {
          mp = mps[[ff]]
          split_target(project_target(mp$stacked_gamma$Z_pool, kern,
                                      emb_gamma_train[[ff]]),
                       nrow(mp$stacked_gamma$Z_pool))
        })

        emb_gamma_test = lapply(1:n_folds, function(ff)
          gamma_measure(estimand_obj, predict_haz_folds[[ff]], mps[[ff]]))
        c_gamma_test = lapply(1:n_folds, function(ff) {
          mp = mps[[ff]]
          split_target(project_target(mp$stacked_lambda$Z_pool, kern,
                                      emb_gamma_test[[ff]]),
                       nrow(mp$stacked_lambda$Z_pool))
        })

        gamma_disp_test = lapply(1:n_folds, function(ff)
          gamma_disp_fn(c_gamma_test[[ff]]$target))
        gamma_disp_train = lapply(1:n_folds, function(ff)
          gamma_disp_fn(c_gamma_train[[ff]]$target))

        dchis_gamma_folds = lapply(1:n_folds, function(ff) {
          mp = mps[[ff]]
          build_gamma_correction(predict_haz_folds[[ff]], estimand_obj, mp$eval_Z,
                              length(mp$mesh), mp$du_bin)
        })

        prev_gamma_alpha = vector("list", n_folds)
        prev_gamma_beta = vector("list", n_folds)

        for (jg in sweep_order) {
          eta_gam = sigma2_grid[jg] / n
          all_terms = c(); all_direct = c(); all_noise_var = c()
          gamma_fits = vector("list", n_folds)

          for (ff in 1:n_folds) {
            mp = mps[[ff]]
            gamma_fit = kernel_bregman(mp$stacked_gamma$Z_pool, kern,
                             eta = eta_gam, dispersion = gamma_disp_train[[ff]],
                             target = c_gamma_train[[ff]]$target,
                             target_null = c_gamma_train[[ff]]$target_null,
                             alpha0 = prev_gamma_alpha[[ff]],
                             beta0 = prev_gamma_beta[[ff]])
            gamma_fit_ = .bind(gamma_fit, mp$eval_Z)

            ct = correction_terms(mp, fold_fits[[ff]]$lambda_fit_, gamma_fit_,
                                   lambda_disp, dchis_gamma_folds[[ff]],
                                   Q_comp, discrete)

            all_terms = c(all_terms, dir_per_fold[[ff]]$direct + ct$correction)
            all_direct = c(all_direct, dir_per_fold[[ff]]$direct)
            all_noise_var = c(all_noise_var, ct$noise_var)
            gamma_fits[[ff]] = gamma_fit
            prev_gamma_alpha[[ff]] = gamma_fit$alpha
            prev_gamma_beta[[ff]] = gamma_fit$beta
          }
          # Sanity check: gamma = dchis(phi) has non-finite values.
          # phi is an implementation detail (not identifiable); only gamma matters.
          gamma_degenerate = any(sapply(1:n_folds, function(ff) {
            gf = gamma_fits[[ff]]
            phi = as.vector(gf$K %*% gf$alpha)
            if (ncol(gf$B) > 0) phi = phi + as.vector(gf$B %*% gf$beta)
            gamma_vals = .dchis(gf$dispersion, phi)
            any(!is.finite(gamma_vals))
          }))
          if (gamma_degenerate) break

          cv_gamma_vals = sapply(1:n_folds, function(ff)
            cv_dual_loss(gamma_fits[[ff]],
                         mps[[ff]]$stacked_lambda$Z_pool,
                         c_gamma_test[[ff]]$target,
                         dispersion_test = gamma_disp_test[[ff]]))

          grid_est[jl, jg] = mean(all_terms)
          grid_dir_val[jl, jg] = dir_mean
          grid_dir_se[jl, jg] = dir_se
          grid_dr_se[jl, jg] = sd(all_terms) / sqrt(n)
          grid_cv_lambda[jl, jg] = cv_lambda[jl]
          grid_cv_gamma[jl, jg] = mean(cv_gamma_vals)
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
        loo_lam = as.vector(t(grid_cv_lambda)),
        loo_gam = as.vector(t(grid_cv_gamma)),
        rep = rep_id)
      all_paths[[length(all_paths) + 1]] = path

      # --- Selection strategies ---
      jg_mid = ceiling(m_grid / 2)
      lambda_slice = path[path$sigma2_gam == sigma2_grid[jg_mid], ]
      lambda_slice = lambda_slice[order(lambda_slice$sigma2_lam), ]

      select_lam = function(tuning) {
        ctx = tuning_context(sigma2_grid / n, n = n,
          dir_val = function() lambda_slice$dir_val,
          dir_se = function() lambda_slice$dir_se,
          loo_scores = function() lambda_slice$loo_lam)
        tuning$select(ctx)
      }

      select_gamma = function(sel_lam_val) {
        gamma_slice = path[path$sigma2_lam == sel_lam_val, ]
        gamma_slice = gamma_slice[order(gamma_slice$sigma2_gam), ]
        idx = which.min(gamma_slice$loo_gam)
        if (length(idx) == 0) return(NA_real_)
        sigma2_grid[idx]
      }

      strategies = list(
        list(label = "fixed",  tuning = fixed_tuning(1 / n)),
        list(label = "lepski", tuning = lepski_tuning()),
        list(label = "cv",     tuning = cv_tuning()))

      # Estimate + bootstrap at a selected (s2_lam, s2_gam).
      # Uses survival_effect from the package to get a dr_result.
      estimate_at = function(s2_lam, s2_gam) {
        survival_effect(observations, kern, estimand_obj, lambda_disp,
                        eta_lam = s2_lam / n, eta_gam = s2_gam / n,
                        horizon = horizon, M_train = M_tr, n_folds = n_folds,
                        Q_comp = Q_comp, discrete = discrete,
                        gamma_disp_fn = gamma_disp_fn)
      }

      # Build a coverage row from a dr_result + optional boot_result.
      coverage_row = function(result, label, s2_lam, s2_gam, boot = NULL) {
        all_boot[[length(all_boot) + 1]] <<- list(
          rep = rep_id, method = mname, estimand = en, tuning = label,
          result = result, boot = boot)
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
            logger("%s %s %s: %d/%d bootstrap reps degenerate",
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
        sel_sigma2_gam = select_gamma(sel_sigma2_lam)
        if (is.na(sel_sigma2_gam)) {
          logger("%s %s %s: no valid gamma at sigma2_lam=%.3g", mname, en, label, sel_sigma2_lam)
          all_coverage[[length(all_coverage) + 1]] = na_row(mname, en)
          next
        }
        logger("%s %s %s selected: sigma2_lam=%.3g, sigma2_gam=%.3g",
            mname, en, label, sel_sigma2_lam, sel_sigma2_gam)

        result = estimate_at(sel_sigma2_lam, sel_sigma2_gam)
        boot = if (boot_reps > 0) bootstrap(result, boot_reps) else NULL
        all_coverage[[length(all_coverage) + 1]] = coverage_row(result, label,
                                                                 sel_sigma2_lam, sel_sigma2_gam, boot)

        if (label == "lepski") {
          sel_idx = which(sigma2_grid == sel_sigma2_lam)
          us_idx = sel_idx - 2
          if (us_idx >= 1) {
            us_sigma2_gam = select_gamma(sigma2_grid[us_idx])
            result2 = estimate_at(sigma2_grid[us_idx], us_sigma2_gam)
            boot2 = if (boot_reps > 0) bootstrap(result2, boot_reps) else NULL
            all_coverage[[length(all_coverage) + 1]] = coverage_row(result2, "lepski-2",
                                                                     sigma2_grid[us_idx], us_sigma2_gam, boot2)
          }
        }
      }
      logger("%s %s: %.1fs", mname, en, (proc.time() - t_est)[3])
    }
    logger("%s total: %.1fs", mname, (proc.time() - t_method)[3])
  }

  list(coverage = do.call(rbind, all_coverage),
       paths = if (length(all_paths) > 0) do.call(rbind, all_paths) else NULL,
       boot = all_boot)
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
all_boot = do.call(c, lapply(all_results, `[[`, "boot"))

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
out_file = sprintf("examples/coverage-results-%s-n%d-%s-r%d-%d.rds",
                   cell_name, n_obs, format(Sys.time(), "%Y%m%d-%H%M"),
                   start_rep, start_rep + n_reps - 1)
saveRDS(list(results = df, paths = all_paths, boot = all_boot,
             params = list(n_obs = n_obs, n_reps = n_reps,
                           boot_reps = boot_reps, truths = truths,
                           cell = cell, dgp_def = dgp_def)), out_file)
cat(sprintf("Saved: %s\n", out_file))
