#' # Estimand Constructors
#'
#' An estimand is the \eqn{(\\psi, \\dot\\psi_Z)} pair that defines what
#' we estimate and what the balancing weights must balance.
#'
#' **Paper**: the estimand functional \eqn{\\psi(\\mu)} maps a nuisance
#' function \eqn{\\mu} (outcome model or hazard) to a scalar. Its
#' Gateaux derivative \eqn{\\dot\\psi_Z(\\mu)[h]} defines the moment
#' condition:
#' \deqn{\\hat P\\, \\gamma(Z)\\\{Y - \\mu(Z)\\\} \\approx \\dot\\psi_Z(\\mu)[h]
#'   \\quad \\text\{for all \} h \\in B_\\rho}
#'
#' The Riesz representer \eqn{\\gamma_\{\\dot\\psi\}} satisfies
#' \eqn{\\dot\\psi_Z(h) = P\\,\\gamma_\{\\dot\\psi\}(Z)\\,h(Z)} and gives
#' the DR correction:
#' \deqn{\\hat\\psi = \\hat P\\,\\psi_Z(\\hat\\mu) +
#'   \\hat P\\,\\hat\\gamma(Z)\\\{Y - \\hat\\mu(Z)\\\}}
#'
#' ## Interface
#'
#' **Single-outcome** estimands provide:
#'
#' - `theta(mu1, mu0)` — plug-in functional
#' - `se(mu1, mu0, n)` — SE of direct term (for Lepski)
#' - `dr(mu1, mu0, Y, W, gamma1, gamma0)` — DR-corrected estimate + EIF
#'
#' **Survival** estimands provide:
#'
#' - `direct(predict_dLambda, eval_Z, grid)` — plug-in per subject
#' - `dot_psi_Z(predict_dLambda, eval_Z, grid)` — function(k) giving
#'   the Gateaux derivative as an n-vector
#' - `dpsi_grid(predict_dLambda, eval_Z, mesh_u, train_grid)` —
#'   dpsi evaluation grid \eqn{(Z, r)} for gamma embedding
#' - `W_fn(A)` — per-subject sign vector for sign-flip dispersion
#'
#' `predict_dLambda(k, Z)` takes grid index \eqn{k} and covariate
#' matrix \eqn{Z = (W, X)}, returns n-vector of hazard increments
#' \eqn{d\\Lambda_k}. The estimand evaluates at counterfactual arms
#' as needed.
#'
#' Both types carry a `$type` field for dispatch.

# ============================================================
# Single-outcome estimands (spinoff4)
# ============================================================

#' ## Treatment-Specific Mean (TSM)
#'
#' **Paper**: \eqn{\\psi^w(\\mu) = E\\,\\mu(w, X)} with Riesz representer
#' \eqn{\\gamma_\{\\dot\\psi^w\}(W, X) = \\mathbf\{1\}(W = w) / \\pi_w(X)}.
#'
#' **Code**: `theta` = \eqn{\\hat P\\,\\hat\\mu_w(X)}. DR correction adds
#' \eqn{\\hat P\\,\\hat\\gamma_w(Y - \\hat\\mu_w)}. EIF:
#' \eqn{\\hat V_i = \\hat\\mu_w(X_i) - \\hat\\psi + \\hat\\gamma_w(Z_i)
#' \\\{Y_i - \\hat\\mu_w(Z_i)\\\}}.
treatment_specific_mean = function(arm = 1) {
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
  ), class = c("effect_estimand", "estimand"))
}

#' ## Variance of CATE (VCATE)
#'
#' **Paper**: \eqn{\\psi(\\mu) = \\mathrm\{Var\}(\\mu_1(X) - \\mu_0(X))}
#' \eqn{= E[\\tau(X)^2] - (E[\\tau(X)])^2} where \eqn{\\tau(X) = \\mu_1(X) -
#' \\mu_0(X)}. This is a second-order functional (quadratic in
#' \eqn{\\mu}).
#'
#' **Code**: DR correction uses two linear corrections:
#' \eqn{\\hat E[\\tau^2]_\{\\mathrm\{DR\}\} = \\hat P\\,\\hat\\tau^2 +
#' 2\\,\\hat P\\,\\hat\\tau\\,(\\hat\\gamma_1 r_1 - \\hat\\gamma_0 r_0)}
#' where \eqn{r_w = Y - \\hat\\mu_w}. The ATE correction is
#' \eqn{\\hat P\\,\\hat\\tau + \\hat P\\,(\\hat\\gamma_1 r_1 - \\hat\\gamma_0 r_0)}.
#' EIF:
#' \eqn{\\hat V_i = (\\hat\\tau_i - \\hat\\psi_\{\\mathrm\{ATE\}\})^2 +
#' 2(\\hat\\tau_i - \\hat\\psi_\{\\mathrm\{ATE\}\})(\\hat\\gamma_\{1i\} r_\{1i\}
#' - \\hat\\gamma_\{0i\} r_\{0i\}) - \\hat\\psi}.
vcate = function() {
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
  ), class = c("effect_estimand", "estimand"))
}

#' ## Risk Ratio
#'
#' **Paper**: \eqn{\\psi(\\mu) = E[\\mu_1(X)] / E[\\mu_0(X)]}.
#' The EIF uses the delta method on the ratio of two TSMs:
#' \eqn{\\hat V_i^\{\\mathrm\{RR\}\} = \\hat V_i^\{(1)\} / \\hat\\psi_0
#' - (\\hat\\psi_1 / \\hat\\psi_0^2)\\,\\hat V_i^\{(0)\}}
#' where \eqn{\\hat V_i^\{(w)\}} is the TSM EIF for arm \eqn{w}.
#'
#' **Code**: computes both arm-specific TSM DR estimates, takes
#' the ratio, and constructs the delta-method EIF.
risk_ratio = function() {
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
  ), class = c("effect_estimand", "estimand"))
}

# ============================================================
# Survival estimands (spinoff3)
# ============================================================

#' ## Survival Probability TSM
#'
#' **Paper**: \eqn{\\psi^w(\\lambda) = E[S_\{\\bar t\}(w, X)]} where
#'
#' - Discrete: \eqn{S_\{\\bar t\}(w, X) = \\prod_\{u \\in \\mathcal\{U\}\}
#'   \\\{1 - \\lambda_u(w, X)\\\}}
#' - Continuous: \eqn{S_\{\\bar t\}(w, X) = \\exp\\bigl(-\\int_0^\{\\bar t\}
#'   \\lambda_u(w, X)\\,du\\bigr)}
#'
#' The Gateaux derivative evaluates at counterfactual arm \eqn{w}:
#'
#' - Discrete: \eqn{\\dot\\psi^w_Z(\\lambda)[h]_k =
#'   -S_\{\\bar t\}(w, X) / \\\{1 - \\lambda_k(w, X)\\\}}
#' - Continuous: \eqn{\\dot\\psi^w_Z(\\lambda)[h]_k = -S_\{\\bar t\}(w, X)}
#'
#' **Code**: `direct` materializes \eqn{d\\Lambda} at \eqn{(w, X)} and
#' calls `survival_probability`. `dot_psi_Z` calls `survival_probability_dot` which
#' returns a closure over grid index \eqn{k}. `dpsi_grid` builds the
#' counterfactual evaluation grid \eqn{Z_\{\\mathrm\{ctf\}\} = (u_m, w, X_i)}
#' for the gamma embedding vector \eqn{c_\\gamma}.
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
      survival_probability(dL, grid)
    },

    dot_psi_Z = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      Za = cbind(arm, eval_X)
      dL = materialize_dLambda(\(k) predict_dLambda(k, Za), grid$M)
      survival_probability_dot(dL, grid)
    },

    dpsi_grid = function(predict_dLambda, eval_Z, mesh_u, train_grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      Za = cbind(arm, eval_X)
      dL = materialize_dLambda(\(k) predict_dLambda(k, Za), train_grid$M)
      dp = survival_probability_dot(dL, train_grid)
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

#' ## Survival Probability ATE
#'
#' **Paper**: \eqn{\\psi(\\lambda) = E[S_\{\\bar t\}(1, X)] -
#' E[S_\{\\bar t\}(0, X)]} with Gateaux derivative
#' \eqn{\\dot\\psi_Z(\\lambda)[h]_k = \\dot\\psi^1_\{Z,k\} - \\dot\\psi^0_\{Z,k\}}.
#'
#' The Riesz representer (continuous time):
#' \deqn{\\gamma_\{\\dot\\psi\}(u, W, X) = -\\frac\{S_\{\\bar t\}(W, X)
#'   \\cdot W\}\{\\pi(W, X)\\, S_u(W, X)\\, G_u(W, X)\}}
#'
#' **Code**: evaluates both arms' survival and derivatives
#' separately, takes the difference. `dpsi_grid` stacks
#' counterfactual grids for both arms: \eqn{Z_\{\\mathrm\{ctf\}\} =
#' [(u_m, 1, X_i); (u_m, 0, X_i)]} with Riesz targets
#' \eqn{r_1 > 0} and \eqn{r_0 < 0} (negated for ATE = \eqn{\\psi^1 - \\psi^0}).
#' `W_fn(A) = 2A - 1` maps \eqn{\\\{0,1\\\} \\to \\\{-1,+1\\\}} for the
#' sign-flip dispersion.
survival_probability_ate = function() {
  structure(list(
    name = "survival.probability",
    type = "survival",
    W_fn = function(A) 2 * A - 1,

    #' Sigma for the gamma dispersion.  For the ATE Riesz
    #' representer, sign(gamma) = -sign(W): treated weights
    #' are negative, control positive.  Using -W instead of
    #' sign(r_hat) makes the sign constraint robust to
    #' hazard estimates where lambda_hat > 1 could flip
    #' sign(r_hat).
    gamma_disp_sigma = function(A) -(2 * A - 1),

    direct = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), grid$M)
      survival_probability(dL1, grid) - survival_probability(dL0, grid)
    },

    dot_psi_Z = function(predict_dLambda, eval_Z, grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), grid$M)
      d1 = survival_probability_dot(dL1, grid)
      d0 = survival_probability_dot(dL0, grid)
      function(k) d1(k) - d0(k)
    },

    # Per-arm dot_psi (for gamma embedding with ATE sign-flip).
    dot_psi_arm = function(predict_dLambda, eval_Z, grid, arm) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL = materialize_dLambda(\(k) predict_dLambda(k, cbind(arm, eval_X)), grid$M)
      survival_probability_dot(dL, grid)
    },

    dpsi_grid = function(predict_dLambda, eval_Z, mesh_u, train_grid) {
      eval_X = eval_Z[, -1, drop = FALSE]
      dL1 = materialize_dLambda(\(k) predict_dLambda(k, cbind(1, eval_X)), train_grid$M)
      dL0 = materialize_dLambda(\(k) predict_dLambda(k, cbind(0, eval_X)), train_grid$M)
      dp1 = survival_probability_dot(dL1, train_grid)
      dp0 = survival_probability_dot(dL0, train_grid)
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

#' ## RMST TSM
#'
#' **Paper**: \eqn{\\psi^w_\{\\mathrm\{RMST\}\}(\\lambda) =
#' E\\bigl[\\int_0^\{\\bar t\} S_u(w, X)\\,du\\bigr]} (continuous) or
#' \eqn{E\\bigl[\\sum_\{u=1\}^\{\\bar t\} S_u(w, X)\\bigr]} (discrete).
#'
#' The Gateaux derivative (continuous):
#' \eqn{\\dot\\psi^w_\{Z,k\} = -\\int_k^\{\\bar t\} S_u(w, X)\\,du}
#' — the tail integral of the survival curve from grid point \eqn{k}.
#'
#' **Code**: `direct` calls `rmst(dL, grid)` which computes
#' \eqn{\\int S_u\\,du} via the grid's quadrature weights.
#' `dot_psi_Z` calls `rmst_dot` which precomputes the tail
#' integral (discrete: rectangle sum; continuous: trapezoid
#' cumulative) and returns a closure over \eqn{k}.
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

#' ## RMST ATE
#'
#' **Paper**: \eqn{\\psi_\{\\mathrm\{RMST\}\}(\\lambda) =
#' E\\bigl[\\int_0^\{\\bar t\} S_u(1, X)\\,du\\bigr] -
#' E\\bigl[\\int_0^\{\\bar t\} S_u(0, X)\\,du\\bigr]}.
#'
#' Riesz representer (continuous):
#' \deqn{\\gamma^\{\\mathrm\{RMST\}\}(v, W, X) = -\\frac\{W \\cdot
#'   \\int_v^\{\\bar t\} S_u(W, X)\\,du\}\{\\pi(W, X)\\, S_v(W, X)\\,
#'   G_v(W, X)\}}
#'
#' **Code**: same structure as `survival_probability_ate` —
#' evaluates both arms, takes the difference. Uses `rmst_dot`
#' for the tail-integral derivative.
rmst_ate = function() {
  structure(list(
    name = "rmst",
    type = "survival",
    W_fn = function(A) 2 * A - 1,
    gamma_disp_sigma = function(A) -(2 * A - 1),

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
