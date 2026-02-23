## Estimands: the (psi, dot_psi) pair that defines what we estimate.
##
## Single-outcome estimands provide:
##   theta(mu1, mu0)                   — plug-in functional
##   se(mu1, mu0, n)                   — SE of direct term (for Lepski)
##   dr(mu1, mu0, Y, W, gamma1, gamma0) — DR-corrected estimate
##
## Survival estimands provide two levels of interface:
##
##   LOW-LEVEL (building blocks, don't know about treatment arms):
##     psi(haz_fn, Q, du)               — outcome functional given a haz_fn
##     dot_psi(haz_fn, Q, du)            — Gateaux derivative given a haz_fn
##     psi_discrete(h_fn, M)             — discrete analog
##     dot_psi_discrete(h_fn, M)         — discrete analog
##
##   HIGH-LEVEL (what bregbal_surv calls; knows about treatment arms):
##     direct_discrete(predict_haz, X, M)       — plug-in estimand per subject
##     dot_psi_Z_discrete(predict_haz, X, M)    — function(k) giving n-vector
##     direct_cts(predict_haz, X, horizon, Q, du_bin)
##     dot_psi_Z_cts(predict_haz, X, horizon, Q, du_bin)
##
##   predict_haz(k, Z) takes mesh index k and covariate matrix Z (with
##   treatment in col 1), returns n-vector of hazard values.
##   The estimand evaluates at whatever counterfactual Z it needs.
##
##   The user passes an appropriate dispersion for the gamma model.
##   (Entropy for TSM, sign-flip entropy for ATE, quadratic for Pham-style, etc.)
##
## Both carry a $type field ("single_outcome" or "survival") for dispatch.

# ============================================================
# Single-outcome estimands (spinoff4)
# ============================================================

#' Treatment-specific mean: E[mu_arm(X)].
#' DR: direct + mean(gamma_arm * (Y - mu_arm)).
tsm_estimand = function(arm = 1) {
  list(
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
  )
}

#' Variance of CATE: Var(mu1(X) - mu0(X)).
vcate_estimand = function() {
  list(
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
  )
}

#' Risk ratio: E[mu1(X)] / E[mu0(X)].
rr_estimand = function() {
  list(
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
  )
}

# ============================================================
# Survival estimands (spinoff3)
# ============================================================

# ---- Low-level survival building blocks ----
# These compute psi and dot_psi given an h_fn. They don't know about
# treatment arms — the caller provides h_fn evaluated at whatever Z it wants.

.surv_prob_psi = function(haz_fn, Q, du) {
  running = rep(0, length(haz_fn(1)))
  h_prev = haz_fn(1)
  for (k in 2:(Q + 1)) {
    h_curr = haz_fn(k)
    running = running + (h_prev + h_curr) / 2 * du
    h_prev = h_curr
  }
  exp(-running)
}

.surv_prob_psi_discrete = function(h_fn, M) {
  S = rep(1, length(h_fn(1)))
  for (k in 1:M) S = S * (1 - h_fn(k))
  S
}

.surv_prob_dot_psi_discrete = function(h_fn, M) {
  S = .surv_prob_psi_discrete(h_fn, M)
  h_vals = matrix(0, length(S), M)
  for (k in 1:M) h_vals[, k] = h_fn(k)
  function(k) -S / (1 - h_vals[, k])
}

.surv_prob_dot_psi = function(haz_fn, Q, du) {
  S = .surv_prob_psi(haz_fn, Q, du)
  function(u) -S
}

.surv_prob_surv_curve = function(haz_fn, Q, du) {
  n = length(haz_fn(1))
  S = matrix(1, n, Q + 1)
  running = rep(0, n)
  h_prev = haz_fn(1)
  for (k in 2:(Q + 1)) {
    h_curr = haz_fn(k)
    running = running + (h_prev + h_curr) / 2 * du
    S[, k] = exp(-running)
    h_prev = h_curr
  }
  S
}

# ---- High-level survival estimand constructors ----
# These know about treatment arms and provide the interface bregbal_surv calls.
#
# predict_haz(k, Z): takes mesh index k and covariate matrix Z (treatment in
#   col 1), returns n-vector of hazard values. The estimand evaluates at
#   whatever counterfactual Z it needs.

#' Survival probability TSM: E[S(tau | a, X)].
#'
#' dot_psi_Z(u, w, x) = -S(tau | a, x) / (1 - h_u(a, x)).
#' Does not depend on observed treatment w.
surv_prob_tsm = function(arm) {
  list(
    name = paste0("survival.probability.tsm", arm),
    type = "survival",
    arm = arm,
    contrast = FALSE,

    # Low-level (for direct use when h_fn already evaluated at the right Z)
    psi = .surv_prob_psi,
    psi_discrete = .surv_prob_psi_discrete,
    dot_psi = .surv_prob_dot_psi,
    dot_psi_discrete = .surv_prob_dot_psi_discrete,
    surv_curve = .surv_prob_surv_curve,

    # High-level: plug-in E[S(tau | a, X)]
    direct_discrete = function(predict_haz, X, M) {
      Za = cbind(arm, X)
      h_fn = function(k) predict_haz(k, Za)
      .surv_prob_psi_discrete(h_fn, M)
    },

    # High-level: dot_psi_Z evaluated at (a, X) for all subjects
    dot_psi_Z_discrete = function(predict_haz, X, M) {
      Za = cbind(arm, X)
      h_fn = function(k) predict_haz(k, Za)
      .surv_prob_dot_psi_discrete(h_fn, M)
    },

    # Continuous versions
    direct_cts = function(predict_haz, X, horizon, Q) {
      Za = cbind(arm, X)
      du = horizon / Q
      h_fn = function(k) predict_haz(k, Za)
      .surv_prob_psi(h_fn, Q, du)
    },

    dot_psi_Z_cts = function(predict_haz, X, horizon, Q) {
      Za = cbind(arm, X)
      du = horizon / Q
      h_fn = function(k) predict_haz(k, Za)
      .surv_prob_dot_psi(h_fn, Q, du)
    },

    grf_target = "survival.probability"
  )
}

#' Survival probability ATE: E[S(tau | 1, X)] - E[S(tau | 0, X)].
#'
#' dot_psi_Z(u, w, x) = -S(tau|1,x)/(1-h_u(1,x)) + S(tau|0,x)/(1-h_u(0,x)).
#' Does not depend on observed treatment w.
surv_prob_ate = function() {
  list(
    name = "survival.probability",
    type = "survival",
    contrast = TRUE,

    # Low-level (kept for backward compat and building blocks)
    psi = .surv_prob_psi,
    psi_discrete = .surv_prob_psi_discrete,
    dot_psi = .surv_prob_dot_psi,
    dot_psi_discrete = .surv_prob_dot_psi_discrete,
    surv_curve = .surv_prob_surv_curve,

    # High-level: E[S(tau|1,X)] - E[S(tau|0,X)]
    direct_discrete = function(predict_haz, X, M) {
      h_fn_1 = function(k) predict_haz(k, cbind(1, X))
      h_fn_0 = function(k) predict_haz(k, cbind(0, X))
      .surv_prob_psi_discrete(h_fn_1, M) - .surv_prob_psi_discrete(h_fn_0, M)
    },

    # dot_psi_Z = deriv of S(tau|1,X) minus deriv of S(tau|0,X)
    dot_psi_Z_discrete = function(predict_haz, X, M) {
      h_fn_1 = function(k) predict_haz(k, cbind(1, X))
      h_fn_0 = function(k) predict_haz(k, cbind(0, X))
      d1 = .surv_prob_dot_psi_discrete(h_fn_1, M)
      d0 = .surv_prob_dot_psi_discrete(h_fn_0, M)
      function(k) d1(k) - d0(k)
    },

    # Per-arm dot_psi (for gamma embedding with ATE).
    dot_psi_arm_discrete = function(predict_haz, X, M, arm) {
      h_fn = function(k) predict_haz(k, cbind(arm, X))
      .surv_prob_dot_psi_discrete(h_fn, M)
    },

    direct_cts = function(predict_haz, X, horizon, Q) {
      du = horizon / Q
      h_fn_1 = function(k) predict_haz(k, cbind(1, X))
      h_fn_0 = function(k) predict_haz(k, cbind(0, X))
      .surv_prob_psi(h_fn_1, Q, du) - .surv_prob_psi(h_fn_0, Q, du)
    },

    dot_psi_Z_cts = function(predict_haz, X, horizon, Q) {
      du = horizon / Q
      h_fn_1 = function(k) predict_haz(k, cbind(1, X))
      h_fn_0 = function(k) predict_haz(k, cbind(0, X))
      d1 = .surv_prob_dot_psi(h_fn_1, Q, du)
      d0 = .surv_prob_dot_psi(h_fn_0, Q, du)
      function(u) d1(u) - d0(u)
    },

    grf_target = "survival.probability"
  )
}

#' Backward-compatible wrapper. Returns surv_prob_ate().
surv_prob_estimand = function() surv_prob_ate()

# ---- Low-level RMST building blocks ----

.rmst_psi = function(haz_fn, Q, du) {
  n = length(haz_fn(1))
  running_cumhaz = rep(0, n)
  running_rmst = rep(0, n)
  S_prev = rep(1, n)
  h_prev = haz_fn(1)
  for (k in 2:(Q + 1)) {
    h_curr = haz_fn(k)
    running_cumhaz = running_cumhaz + (h_prev + h_curr) / 2 * du
    S_curr = exp(-running_cumhaz)
    running_rmst = running_rmst + (S_prev + S_curr) / 2 * du
    S_prev = S_curr
    h_prev = h_curr
  }
  running_rmst
}

.rmst_dot_psi = function(haz_fn, Q, du) {
  S_mat = .surv_prob_surv_curve(haz_fn, Q, du)
  eval_grid = (0:Q) * du
  cum = matrix(0, nrow(S_mat), Q + 1)
  for (k in 2:(Q + 1))
    cum[, k] = cum[, k - 1] + (S_mat[, k - 1] + S_mat[, k]) / 2 * du
  total = cum[, Q + 1]
  function(u) {
    k = pmax(findInterval(u, eval_grid, rightmost.closed = TRUE), 1)
    -(total - cum[, k])
  }
}

.rmst_psi_discrete = function(h_fn, M, du) {
  n = length(h_fn(1))
  S = rep(1, n)
  rmst = rep(0, n)
  for (k in 1:M) {
    S = S * (1 - h_fn(k))
    rmst = rmst + S * du
  }
  rmst
}

.rmst_dot_psi_discrete = function(h_fn, M, du) {
  n = length(h_fn(1))
  S_vec = matrix(0, n, M)
  h_vals = matrix(0, n, M)
  S = rep(1, n)
  for (k in 1:M) {
    h_vals[, k] = h_fn(k)
    S = S * (1 - h_vals[, k])
    S_vec[, k] = S
  }
  tail = matrix(0, n, M)
  tail[, M] = S_vec[, M] * du
  for (k in (M - 1):1) tail[, k] = tail[, k + 1] + S_vec[, k] * du
  function(k) -tail[, k] / (1 - h_vals[, k])
}

#' RMST TSM: E[int_0^tau S(u | a, X) du].
rmst_tsm = function(arm) {
  list(
    name = paste0("rmst.tsm", arm),
    type = "survival",
    arm = arm,
    contrast = FALSE,

    psi = .rmst_psi,
    psi_discrete = .rmst_psi_discrete,
    dot_psi = .rmst_dot_psi,
    dot_psi_discrete = .rmst_dot_psi_discrete,
    surv_curve = .surv_prob_surv_curve,

    direct_discrete = function(predict_haz, X, M, du) {
      Za = cbind(arm, X)
      h_fn = function(k) predict_haz(k, Za)
      .rmst_psi_discrete(h_fn, M, du)
    },

    dot_psi_Z_discrete = function(predict_haz, X, M, du) {
      Za = cbind(arm, X)
      h_fn = function(k) predict_haz(k, Za)
      .rmst_dot_psi_discrete(h_fn, M, du)
    },

    direct_cts = function(predict_haz, X, horizon, Q) {
      Za = cbind(arm, X)
      du = horizon / Q
      h_fn = function(k) predict_haz(k, Za)
      .rmst_psi(h_fn, Q, du)
    },

    dot_psi_Z_cts = function(predict_haz, X, horizon, Q) {
      Za = cbind(arm, X)
      du = horizon / Q
      h_fn = function(k) predict_haz(k, Za)
      .rmst_dot_psi(h_fn, Q, du)
    },

    grf_target = "RMST"
  )
}

#' RMST ATE: E[int_0^tau S(u|1,X) du] - E[int_0^tau S(u|0,X) du].
rmst_ate = function() {
  list(
    name = "rmst",
    type = "survival",
    contrast = TRUE,

    psi = .rmst_psi,
    psi_discrete = .rmst_psi_discrete,
    dot_psi = .rmst_dot_psi,
    dot_psi_discrete = .rmst_dot_psi_discrete,
    surv_curve = .surv_prob_surv_curve,

    direct_discrete = function(predict_haz, X, M, du) {
      h_fn_1 = function(k) predict_haz(k, cbind(1, X))
      h_fn_0 = function(k) predict_haz(k, cbind(0, X))
      .rmst_psi_discrete(h_fn_1, M, du) - .rmst_psi_discrete(h_fn_0, M, du)
    },

    dot_psi_Z_discrete = function(predict_haz, X, M, du) {
      h_fn_1 = function(k) predict_haz(k, cbind(1, X))
      h_fn_0 = function(k) predict_haz(k, cbind(0, X))
      d1 = .rmst_dot_psi_discrete(h_fn_1, M, du)
      d0 = .rmst_dot_psi_discrete(h_fn_0, M, du)
      function(k) d1(k) - d0(k)
    },

    dot_psi_arm_discrete = function(predict_haz, X, M, du, arm) {
      h_fn = function(k) predict_haz(k, cbind(arm, X))
      .rmst_dot_psi_discrete(h_fn, M, du)
    },

    direct_cts = function(predict_haz, X, horizon, Q) {
      du = horizon / Q
      h_fn_1 = function(k) predict_haz(k, cbind(1, X))
      h_fn_0 = function(k) predict_haz(k, cbind(0, X))
      .rmst_psi(h_fn_1, Q, du) - .rmst_psi(h_fn_0, Q, du)
    },

    dot_psi_Z_cts = function(predict_haz, X, horizon, Q) {
      du = horizon / Q
      h_fn_1 = function(k) predict_haz(k, cbind(1, X))
      h_fn_0 = function(k) predict_haz(k, cbind(0, X))
      d1 = .rmst_dot_psi(h_fn_1, Q, du)
      d0 = .rmst_dot_psi(h_fn_0, Q, du)
      function(u) d1(u) - d0(u)
    },

    grf_target = "RMST"
  )
}

#' Backward-compatible wrapper. Returns rmst_ate().
rmst_estimand = function() rmst_ate()
