## Estimands: the (psi, dot_psi) pair that defines what we estimate.
##
## Single-outcome estimands provide:
##   theta(mu1, mu0)                   — plug-in functional
##   se(mu1, mu0, n)                   — SE of direct term (for Lepski)
##   dr(mu1, mu0, Y, W, gamma1, gamma0) — DR-corrected estimate
##
## Survival estimands provide:
##   psi(lambda_fn, Z, horizon, ...)   — outcome functional (e.g. S(tau|Z))
##   dot_psi(lambda_fn, Z, horizon, ...) — Gateaux derivative w.r.t. lambda
##   plus alpha-based versions for bootstrap
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

#' Survival probability: psi(lambda, Z) = S(tau | Z) = exp(-int_0^tau lambda du).
#' dot_psi(lambda, Z)(u) = -S(tau | Z), constant in u.
#'
#' Time integration (cum_surv, surv_curve) lives here, not on the dispersion.
#' The dispersion says dchis(phi) = exp(phi); the estimand says
#' "integrate exp(-int lambda) to get S(t)."
surv_prob_estimand = function() {
  list(
    name = "survival.probability",
    type = "survival",

    # S(tau | Z) = exp(-int_0^tau lambda du) via trapezoidal rule.
    # haz_fn(k) returns n-vector of hazard rates at eval_grid[k] (1-indexed).
    psi = function(haz_fn, Q, du) {
      running = rep(0, length(haz_fn(1)))
      h_prev = haz_fn(1)
      for (k in 2:(Q + 1)) {
        h_curr = haz_fn(k)
        running = running + (h_prev + h_curr) / 2 * du
        h_prev = h_curr
      }
      exp(-running)
    },

    # dot_psi(u) = -S(tau | Z), constant in u.
    dot_psi = function(haz_fn, Q, du) {
      S = surv_prob_estimand()$psi(haz_fn, Q, du)
      function(u) -S
    },

    # Full survival curve: n x (Q+1) matrix.
    surv_curve = function(haz_fn, Q, du) {
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
    },

    grf_target = "survival.probability"
  )
}

#' RMST: psi(lambda, Z) = int_0^tau S(u | Z) du.
#' dot_psi(lambda, Z)(u) = -int_u^tau S(t | Z) dt (tail RMST), varies with u.
rmst_estimand = function() {
  sp = surv_prob_estimand()
  list(
    name = "rmst",
    type = "survival",

    # RMST = int_0^tau S(u) du via trapezoidal rule on S.
    psi = function(haz_fn, Q, du) {
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
    },

    # dot_psi(u) = -int_u^tau S(t) dt (tail integral).
    dot_psi = function(haz_fn, Q, du) {
      S_mat = sp$surv_curve(haz_fn, Q, du)
      eval_grid = (0:Q) * du
      # Cumulative integral of S
      cum = matrix(0, nrow(S_mat), Q + 1)
      for (k in 2:(Q + 1))
        cum[, k] = cum[, k - 1] + (S_mat[, k - 1] + S_mat[, k]) / 2 * du
      total = cum[, Q + 1]
      function(u) {
        k = pmax(findInterval(u, eval_grid, rightmost.closed = TRUE), 1)
        -(total - cum[, k])
      }
    },

    surv_curve = sp$surv_curve,

    grf_target = "RMST"
  )
}
