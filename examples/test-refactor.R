## Quick test: verify the refactored test-coverage.R plumbing works.
## Runs one rep of KBW_cts on paper_cts with no bootstrap.

# Source everything the same way test-coverage.R does
source("R/kernel.R")
for (f in c("R/estimands.R", "R/survival.R", "R/tuning.R", "R/dgp.R"))
  source(f)

# Source just the helpers from test-coverage.R (lines 96-310ish)
# For now, inline the key ones:

# --- gam_measure ---
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

# --- compute_direct ---
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

# --- correction_terms ---
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
  list(correction = dirac - comp$comp, noise_var = comp$noise_var, events = events)
}

# --- make_dchis_gamma_ts ---
make_dchis_gamma_ts = function(predict_haz, estimand_obj, eval_Z,
                                M_train, du_bin) {
  n_eval = nrow(atleast_2d(eval_Z))
  eval_hash = apply(eval_Z, 1, paste0, collapse = "|")
  lookup = setNames(seq_len(n_eval), eval_hash)
  has_per_arm = !is.null(estimand_obj$dot_psi_arm_discrete)
  if (has_per_arm) {
    dp1 = estimand_obj$dot_psi_arm_discrete(predict_haz, eval_Z, M_train, du_bin, arm = 1)
    dp0 = estimand_obj$dot_psi_arm_discrete(predict_haz, eval_Z, M_train, du_bin, arm = 0)
    r1_mat = r0_mat = matrix(NA_real_, n_eval, M_train)
    for (k in 1:M_train) { r1_mat[, k] = dp1(k); r0_mat[, k] = dp0(k) }
  } else {
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

# --- dr_fold ---
dr_fold = function(fd, lam_fit, kern, eta_gam, gam_disp, dchis_gamma,
                   lam_disp, estimand, horizon, Q_comp, discrete,
                   gam_target = NULL, gam_target_null = NULL) {
  n_eval = fd$n_eval
  lam_bound = bind_subjects(lam_fit, fd$eval_Z)
  M_train = length(fd$mesh)
  predict_phi_fn = function(Z) predict_phi(lam_fit, Z)
  dr = compute_direct(predict_phi_fn, lam_disp, estimand, fd$eval_Z,
                       M_train, fd$du_bin, fd$mesh, horizon, discrete)
  if (is.null(gam_target)) {
    emb = gam_measure(estimand, dr$predict_haz, fd$eval_Z,
                      fd$mesh_u, fd$mesh, M_train, fd$du_bin)
    tgt = split_target(project_target(fd$stacked_gam$Z_pool, kern, emb),
                       nrow(fd$stacked_gam$Z_pool))
    gam_target = tgt$target; gam_target_null = tgt$target_null
    gam_disp = target_scaled_entropy(gam_target)
  }
  gam_fit = kernel_bregman(fd$stacked_gam$Z_pool, kern,
                           eta = eta_gam, dispersion = gam_disp,
                           target = gam_target, target_null = gam_target_null)
  gam_bound = bind_subjects(gam_fit, fd$eval_Z)
  ct = correction_terms(lam_bound, gam_bound, fd$T_obs_eval, fd$D_eval,
                         fd$mesh, fd$du_bin, fd$eval_Z, n_eval,
                         lam_disp, dchis_gamma, horizon, Q_comp, discrete)
  list(direct = dr$direct, correction = ct$correction,
       noise_var = ct$noise_var, events = ct$events,
       lam_fit = lam_fit, gam_fit = gam_fit,
       lam_bound = lam_bound, gam_bound = gam_bound,
       predict_haz = dr$predict_haz, dchis_gamma = dchis_gamma)
}

# --- crossfit ---
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

# ===========================================================
# Test: paper_cts DGP, continuous time, 2-fold, surv_prob ATE
# ===========================================================

cat("=== Refactor plumbing test ===\n")

set.seed(1)
n = 100; p = 2; M_tr = 15; horizon = 1
cts_dgp = paper_cts_dgp(horizon = horizon)
dat = cts_dgp$generate(n, p = p)
n_folds = 2
discrete = FALSE
Q_comp = 50

base = matern_kernel(sigma = 2, nu = 3/2)
kern = direct_product(base, iw = 2, levels = c(0, 1))
lam_disp = entropy_dispersion()
gam_disp = entropy_dispersion()
estimand_obj = surv_prob_estimand()

T_obs = dat$T_obs; D = dat$D; A = dat$A; X = dat$X
Z = cbind(A, X)

cat("Data: n=", n, "events=", sum(D), "ATE_true=", sprintf("%.4f", dat$psi1_true - dat$psi0_true), "\n")

# Build fold data
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

cat("Fold data built: ", length(fold_data), " folds\n")
cat("  fold 1: lam=", length(fold_data[[1]]$lam_idx),
    "gam=", length(fold_data[[1]]$gam_idx),
    "eval=", length(fold_data[[1]]$eval_idx), "\n")

# Test dr_fold on fold 1
cat("\n--- Testing dr_fold on fold 1 ---\n")
fd = fold_data[[1]]
eta_lam = 1 / n; eta_gam = 5 / n

lam_fit = kernel_bregman(fd$stacked_lam$Z_pool, kern,
                         eta = eta_lam, dispersion = lam_disp,
                         target = fd$tgt_lam$target,
                         target_null = fd$tgt_lam$target_null,
                         K = fd$K_lam)
cat("Lambda fit: iter=", lam_fit$iter, "beta0=", sprintf("%.4f", lam_fit$beta), "\n")

# Build dchis_gamma
predict_haz_fn = function(k, Z_ev)
  lam_disp$dchis(predict_phi(lam_fit, cbind(fd$mesh[k], Z_ev)))
dchis_gam = make_dchis_gamma_ts(predict_haz_fn, estimand_obj, fd$eval_Z,
                                 length(fd$mesh), fd$du_bin)
res = dr_fold(fd, lam_fit, kern, eta_gam, NULL, dchis_gam,
              lam_disp, estimand_obj, horizon, Q_comp, discrete)
cat("dr_fold result:\n")
cat("  direct: mean=", sprintf("%.4f", mean(res$direct)), "\n")
cat("  correction: mean=", sprintf("%.4f", mean(res$correction)), "\n")
cat("  fold estimate: ", sprintf("%.4f", mean(res$direct + res$correction)), "\n")

# Test crossfit
cat("\n--- Testing crossfit ---\n")
result = crossfit(n, fold_data, function(fd) {
  lam_fit = kernel_bregman(fd$stacked_lam$Z_pool, kern,
                           eta = eta_lam, dispersion = lam_disp,
                           target = fd$tgt_lam$target,
                           target_null = fd$tgt_lam$target_null,
                           K = fd$K_lam)
  predict_haz_fn = function(k, Z_ev)
    lam_disp$dchis(predict_phi(lam_fit, cbind(fd$mesh[k], Z_ev)))
  dchis_gam = make_dchis_gamma_ts(predict_haz_fn, estimand_obj, fd$eval_Z,
                                   length(fd$mesh), fd$du_bin)
  dr_fold(fd, lam_fit, kern, eta_gam, NULL, dchis_gam,
          lam_disp, estimand_obj, horizon, Q_comp, discrete)
}, kern = kern, eta_gam = eta_gam, lam_disp = lam_disp,
   gam_disp = gam_disp, estimand = estimand_obj,
   horizon = horizon, Q_comp = Q_comp, discrete = discrete)

truth = dat$psi1_true - dat$psi0_true
cat("crossfit result:\n")
cat("  est=", sprintf("%.4f", result$est), "\n")
cat("  se=", sprintf("%.4f", result$se), "\n")
cat("  se_amle=", sprintf("%.4f", result$se_amle), "\n")
cat("  truth=", sprintf("%.4f", truth), "\n")
cat("  error=", sprintf("%+.4f", result$est - truth), "\n")
cat("  class=", class(result), "\n")
cat("  n_fold_results=", length(result$fold_results), "\n")

# Test bootstrap
cat("\n--- Testing bootstrap ---\n")

# S3 generic + method (inline from test-coverage.R)
bootstrap = function(obj, ...) UseMethod("bootstrap")

bootstrap.dr_result = function(obj, n_reps, ...) {
  n = obj$n
  kern = obj$kern; lam_disp = obj$lam_disp
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

ci_pct = function(obj, ...) UseMethod("ci_pct")
ci_pct.boot_result = function(obj, level = 0.95, ...) {
  good = !is.na(obj$boot_ates) & !is.nan(obj$boot_ates) & abs(obj$boot_ates) < 100
  if (sum(good) < 10) return(c(NA, NA))
  alpha = 1 - level
  unname(quantile(obj$boot_ates[good], c(alpha / 2, 1 - alpha / 2)))
}

set.seed(42)
boot = bootstrap(result, n_reps = 20)
cat("bootstrap result:\n")
cat("  boot_se=", sprintf("%.4f", sd(boot$boot_ates)), "\n")
cat("  n_bad=", sum(is.na(boot$boot_ates) | abs(boot$boot_ates) > 100), "\n")
ci = ci_pct(boot)
cat("  ci_pct=", sprintf("[%.4f, %.4f]", ci[1], ci[2]), "\n")
cat("  truth in CI?", truth >= ci[1] & truth <= ci[2], "\n")

cat("\n=== Test passed ===\n")
