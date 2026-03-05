## Estimands: the (psi, dot_psi) pair that defines what we estimate.
##
## Single-outcome estimands provide:
##   theta(mu1, mu0)                   — plug-in functional
##   se(mu1, mu0, n)                   — SE of direct term (for Lepski)
##   dr(mu1, mu0, Y, W, gamma1, gamma0) — DR-corrected estimate
##
## Survival estimands provide:
##
##   direct(predict_dLambda, eval_Z, grid) — plug-in estimand per subject
##   dot_psi_Z(predict_dLambda, eval_Z, grid)  — function(k) giving n-vector
##   dpsi_grid(predict_dLambda, eval_Z, mesh_u, train_grid)
##                                              — dpsi evaluation grid (Z, r)
##
##   predict_dLambda(k, Z) takes grid index k and covariate matrix Z (with
##   treatment in col 1), returns n-vector of hazard increments dLambda_k.
##   The estimand evaluates at whatever counterfactual Z it needs.
##
##   All survival operations use the product integral and grid operations
##   from grid.R, unifying discrete and continuous time.
##
##   The estimand provides W_fn(A) returning the per-subject sign vector.
##   The pipeline constructs signflip(gamma_dispersion, W_fn(A)) automatically.
##
## Both carry a $type field ("single_outcome" or "survival") for dispatch.

# ============================================================
# Single-outcome estimands (spinoff4)
# ============================================================

#' Treatment-specific mean: E[mu_arm(X)].
#' DR: direct + mean(gamma_arm * (Y - mu_arm)).
tsm_estimand = function(arm = 1) {
  structure(list(
    name = paste0("tsm", arm),
    type = "single_outcome",

    theta = function(mu1, mu0) {
      mean(if (arm == 1) mu1 else mu0)
    },

    se = function(mu1, mu0, n) {
      vals = if (arm == 1) mu1 else mu0
      sd(vals) / sqrt(n)
    },

    dr = function(mu1, mu0, Y, W, gamma1, gamma0) {
      n = length(Y)
      if (arm == 1) {
        mu_arm = mu1; gam = gamma1
      } else {
        mu_arm = mu0; gam = gamma0
      }
      resid_arm = Y - mu_arm

      direct = mean(mu_arm)
      correction = mean(gam * resid_arm)
      est = direct + correction

      eif = mu_arm - est + gam * resid_arm
      se = sd(eif) / sqrt(n)
      list(est = est, se = se, eif = eif)
    }
  ), class = c("single_outcome_estimand", "estimand"))
}

#' Variance of CATE: Var(mu1(X) - mu0(X)).
vcate_estimand = function() {
  structure(list(
    name = "vcate",
    type = "single_outcome",

    theta = function(mu1, mu0) {
      tau = mu1 - mu0
      mean(tau^2) - mean(tau)^2
    },

    se = function(mu1, mu0, n) {
      tau = mu1 - mu0
      vcate = mean(tau^2) - mean(tau)^2
      phi = (tau - mean(tau))^2 - vcate
      sd(phi) / sqrt(n)
    },

    dr = function(mu1, mu0, Y, W, gamma1, gamma0) {
      n = length(Y)
      tau_hat = mu1 - mu0
      resid1 = Y - mu1; resid0 = Y - mu0
      gamma_resid = gamma1 * resid1 - gamma0 * resid0

      ate_dr = mean(tau_hat) + mean(gamma_resid)
      e_tau2_dr = mean(tau_hat^2) + 2 * mean(tau_hat * gamma_resid)
      vcate = e_tau2_dr - ate_dr^2

      eif = (tau_hat - ate_dr)^2 + 2 * (tau_hat - ate_dr) * gamma_resid - vcate
      se = sd(eif) / sqrt(n)
      list(est = vcate, se = se, eif = eif, ate_dr = ate_dr)
    }
  ), class = c("single_outcome_estimand", "estimand"))
}

#' Risk ratio: E[mu1(X)] / E[mu0(X)].
rr_estimand = function() {
  structure(list(
    name = "rr",
    type = "single_outcome",

    theta = function(mu1, mu0) {
      mean(mu1) / mean(mu0)
    },

    se = function(mu1, mu0, n) {
      a = mean(mu1); b = mean(mu0)
      va = var(mu1) / n; vb = var(mu0) / n; cab = cov(mu1, mu0) / n
      sqrt(va / b^2 + a^2 * vb / b^4 - 2 * a * cab / b^3)
    },

    dr = function(mu1, mu0, Y, W, gamma1, gamma0) {
      n = length(Y)
      resid1 = Y - mu1; resid0 = Y - mu0

      tsm1 = mean(mu1) + mean(gamma1 * resid1)
      tsm0 = mean(mu0) + mean(gamma0 * resid0)
      rr = tsm1 / tsm0

      eif1 = mu1 - tsm1 + gamma1 * resid1
      eif0 = mu0 - tsm0 + gamma0 * resid0
      eif = eif1 / tsm0 - tsm1 / tsm0^2 * eif0

      se = sd(eif) / sqrt(n)
      list(est = rr, se = se, eif = eif, tsm1 = tsm1, tsm0 = tsm0)
    }
  ), class = c("single_outcome_estimand", "estimand"))
}

# ============================================================
# Survival estimands (spinoff3)
# ============================================================

#' Survival probability TSM: E[S(tau | a, X)].
#'
#' direct: S(tau | a, X) via product integral.
#' dot_psi_Z(k): -S(tau | a, X) / (1 - dLambda_k(a, X)).
survival_probability_tsm = function(arm) {
  structure(list(
    name = paste0("survival.probability.tsm", arm),
    type = "survival",
    arm = arm,
    W_fn = function(A) as.numeric(if (arm == 1) A else 1 - A),

    direct = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      Za = cbind(arm, eval_X)
      dL = materialize_dLambda(\(k) predict_dLambda(k, Za), grid$M)
      surv_prob(dL, grid)
    },

    dot_psi_Z = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      Za = cbind(arm, eval_X)
      dL = materialize_dLambda(\(k) predict_dLambda(k, Za), grid$M)
      surv_prob_dot(dL, grid)
    },

    dpsi_grid = function(predict_dLambda, eval_Z, mesh_u, train_grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      Za = cbind(arm, eval_X)
      dL = materialize_dLambda(\(k) predict_dLambda(k, Za), train_grid$M)
      dp = surv_prob_dot(dL, train_grid)
      mesh_k = match(mesh_u, train_grid$points)
      n_eval = nrow(eval_Z)
      n_mu = length(mesh_u)
      r_subj = array(NA_real_, c(n_eval, n_mu))
      for (j in seq_along(mesh_u)) r_subj[, j] = dp(mesh_k[j])
      Z_ctf = do.call(rbind, lapply(seq_along(mesh_u), function(j) {
        cbind(mesh_u[j], arm, eval_X)
      }))
      list(Z = Z_ctf, r = as.vector(r_subj))
    },

    grf_target = "survival.probability"
  ), class = c("survival_estimand", "estimand"))
}

#' Survival probability ATE: E[S(tau | 1, X)] - E[S(tau | 0, X)].
#'
#' direct: S(tau|1,X) - S(tau|0,X) via product integral.
#' dot_psi_Z(k): -S(tau|1,x)/(1-dL_k(1,x)) + S(tau|0,x)/(1-dL_k(0,x)).
survival_probability_ate = function() {
  structure(list(
    name = "survival.probability",
    type = "survival",
    W_fn = function(A) 2 * A - 1,

    direct = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), grid$M)
      surv_prob(dL1, grid) - surv_prob(dL0, grid)
    },

    dot_psi_Z = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), grid$M)
      d1 = surv_prob_dot(dL1, grid)
      d0 = surv_prob_dot(dL0, grid)
      function(k) d1(k) - d0(k)
    },

    # Per-arm dot_psi (for gamma embedding with ATE sign-flip).
    dot_psi_arm = function(predict_dLambda, eval_Z, grid, arm) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL = materialize_dLambda(\(k) predict_dLambda(k, cbind(arm, eval_X)), grid$M)
      surv_prob_dot(dL, grid)
    },

    dpsi_grid = function(predict_dLambda, eval_Z, mesh_u, train_grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), train_grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), train_grid$M)
      dp1 = surv_prob_dot(dL1, train_grid)
      dp0 = surv_prob_dot(dL0, train_grid)
      mesh_k = match(mesh_u, train_grid$points)
      n_eval = nrow(eval_Z)
      n_mu = length(mesh_u)
      r1_subj = r0_subj = array(NA_real_, c(n_eval, n_mu))
      for (j in seq_along(mesh_u)) {
        r1_subj[, j] = dp1(mesh_k[j])
        r0_subj[, j] = -dp0(mesh_k[j])  # negate: ATE = psi1 - psi0
      }
      Z_ctf1 = do.call(rbind, lapply(seq_along(mesh_u), function(j) cbind(mesh_u[j], 1, eval_X)))
      Z_ctf0 = do.call(rbind, lapply(seq_along(mesh_u), function(j) cbind(mesh_u[j], 0, eval_X)))
      list(Z = rbind(Z_ctf1, Z_ctf0),
           r = c(as.vector(r1_subj), as.vector(r0_subj)))
    },

    grf_target = "survival.probability"
  ), class = c("survival_estimand", "estimand"))
}

# ---- RMST estimands ----

#' RMST TSM: E[int_0^tau S(u | a, X) du].
rmst_tsm = function(arm) {
  structure(list(
    name = paste0("rmst.tsm", arm),
    type = "survival",
    arm = arm,
    W_fn = function(A) as.numeric(if (arm == 1) A else 1 - A),

    direct = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      Za = cbind(arm, eval_X)
      dL = materialize_dLambda(\(k) predict_dLambda(k, Za), grid$M)
      rmst(dL, grid)
    },

    dot_psi_Z = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      Za = cbind(arm, eval_X)
      dL = materialize_dLambda(\(k) predict_dLambda(k, Za), grid$M)
      rmst_dot(dL, grid)
    },

    dpsi_grid = function(predict_dLambda, eval_Z, mesh_u, train_grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      Za = cbind(arm, eval_X)
      dL = materialize_dLambda(\(k) predict_dLambda(k, Za), train_grid$M)
      dp = rmst_dot(dL, train_grid)
      mesh_k = match(mesh_u, train_grid$points)
      n_eval = nrow(eval_Z)
      n_mu = length(mesh_u)
      r_subj = array(NA_real_, c(n_eval, n_mu))
      for (j in seq_along(mesh_u)) r_subj[, j] = dp(mesh_k[j])
      Z_ctf = do.call(rbind, lapply(seq_along(mesh_u), function(j) {
        cbind(mesh_u[j], arm, eval_X)
      }))
      list(Z = Z_ctf, r = as.vector(r_subj))
    },

    grf_target = "RMST"
  ), class = c("survival_estimand", "estimand"))
}

#' RMST ATE: E[int_0^tau S(u|1,X) du] - E[int_0^tau S(u|0,X) du].
rmst_ate = function() {
  structure(list(
    name = "rmst",
    type = "survival",
    W_fn = function(A) 2 * A - 1,

    direct = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), grid$M)
      rmst(dL1, grid) - rmst(dL0, grid)
    },

    dot_psi_Z = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), grid$M)
      d1 = rmst_dot(dL1, grid)
      d0 = rmst_dot(dL0, grid)
      function(k) d1(k) - d0(k)
    },

    dot_psi_arm = function(predict_dLambda, eval_Z, grid, arm) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL = materialize_dLambda(\(k) predict_dLambda(k, cbind(arm, eval_X)), grid$M)
      rmst_dot(dL, grid)
    },

    dpsi_grid = function(predict_dLambda, eval_Z, mesh_u, train_grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), train_grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), train_grid$M)
      dp1 = rmst_dot(dL1, train_grid)
      dp0 = rmst_dot(dL0, train_grid)
      mesh_k = match(mesh_u, train_grid$points)
      n_eval = nrow(eval_Z)
      n_mu = length(mesh_u)
      r1_subj = r0_subj = array(NA_real_, c(n_eval, n_mu))
      for (j in seq_along(mesh_u)) {
        r1_subj[, j] = dp1(mesh_k[j])
        r0_subj[, j] = -dp0(mesh_k[j])  # negate: ATE = psi1 - psi0
      }
      Z_ctf1 = do.call(rbind, lapply(seq_along(mesh_u), function(j) cbind(mesh_u[j], 1, eval_X)))
      Z_ctf0 = do.call(rbind, lapply(seq_along(mesh_u), function(j) cbind(mesh_u[j], 0, eval_X)))
      list(Z = rbind(Z_ctf1, Z_ctf0),
           r = c(as.vector(r1_subj), as.vector(r0_subj)))
    },

    grf_target = "RMST"
  ), class = c("survival_estimand", "estimand"))
}
