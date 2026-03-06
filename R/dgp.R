## DGP closures for simulation.
##
## A survival DGP is a list with:
##   name                 — string identifier (used as cache key)
##   hazard(u, Z)         — event hazard rate (cts) or per-step probability (discrete)
##   censor_hazard(u, Z)  — censoring hazard (same convention as event hazard)
##   propensity(X)        — P(A=1|X)
##   horizon              — endpoint: the time at which the estimand is evaluated
##   time                 — "continuous" or "discrete"
##   grid                 — NULL (cts) or c(u_1,...,u_T) (discrete time points)
##   generate(n, p)       — draw a dataset from this DGP
##
## The cumulative hazard Lambda(u) is the bridge between continuous and discrete:
##   Continuous: lambda(u) = Lambda'(u), survival S(u) = exp(-Lambda(u))
##   Discrete:   h_k = 1 - exp(-Delta Lambda_k)
##
## Conversion functions:
##   discretize(dgp, n_steps)  — cts -> discrete (integrate rate over bins)
##   continuize(dgp)           — discrete -> cts (monotone spline of cumhaz)
##
## Cached efficiency bound:
##   efficient_se(dgp, p, estimand)  — semiparametric efficient SD, cached to disk


# ============================================================
# Survival DGP constructor
# ============================================================

#' Build a survival DGP object from its components.
#'
#' @param name String identifier for this DGP (used as cache key).
#' @param hazard function(u, Z) -> vector of length nrow(Z).
#'   Continuous: hazard rate lambda(u, W, X).
#'   Discrete: per-step event probability h_k(W, X).
#' @param censor_hazard Same signature, for censoring.
#' @param propensity function(X) -> P(A=1|X), vector of length nrow(X).
#' @param covariance function(p) -> p x p covariance matrix. Used to generate
#'   X ~ N(0, Sigma). Ignored if draw_X is supplied.
#' @param horizon Numeric, endpoint of follow-up AND estimand target time.
#' @param time "continuous" or "discrete".
#' @param grid For discrete DGPs, the time grid c(u_1, ..., u_T).
#' @param draw_X Optional function(n, p) -> n x p matrix of covariates.
#'   If supplied, overrides covariance-based generation.
#' @return A list of class "survival_dgp".
make_dgp = function(name, hazard, censor_hazard, propensity, covariance = NULL,
                    horizon, time = c("continuous", "discrete"),
                    grid = NULL, draw_X = NULL) {
  time = match.arg(time)
  wt = function(z, a) { z2 = z; z2[, 1] = a; z2 }

  # Draw a dataset of n individuals with p covariates.
  # M_fine: grid resolution for continuous-time event generation (default 1000).
  generate = function(n, p, M_fine = 1000) {
    if (!is.null(draw_X)) {
      X = draw_X(n, p)
    } else {
      Sigma = covariance(p)
      X = matrix(rnorm(n * p), n, p) %*% chol(Sigma)
    }
    pi_x = propensity(X)
    A = rbinom(n, 1, pi_x)
    Z = cbind(A, X)

    if (time == "continuous") {
      du = horizon / M_fine
      fg = (1:M_fine) * du
      T_event = T_censor = rep(Inf, n)
      alive_e = alive_c = rep(TRUE, n)
      for (k in 1:M_fine) {
        he = pmin(hazard(fg[k], Z) * du, 1)
        hc = pmin(censor_hazard(fg[k], Z) * du, 1)
        ev = alive_e & (runif(n) < he); T_event[ev] = fg[k]; alive_e[ev] = FALSE
        ce = alive_c & (runif(n) < hc); T_censor[ce] = fg[k]; alive_c[ce] = FALSE
      }
    } else {
      fg = grid
      M_fine = length(fg)
      T_event = T_censor = rep(Inf, n)
      alive_e = alive_c = rep(TRUE, n)
      for (k in 1:M_fine) {
        he = hazard(fg[k], Z)
        hc = censor_hazard(fg[k], Z)
        ev = alive_e & (runif(n) < he); T_event[ev] = fg[k]; alive_e[ev] = FALSE
        ce = alive_c & (runif(n) < hc); T_censor[ce] = fg[k]; alive_c[ce] = FALSE
      }
    }

    T_obs = pmin(T_event, T_censor)
    D = as.numeric(T_event <= T_censor)

    truth = compute_truth(hazard, Z, horizon, time, if (time == "discrete") grid else NULL,
                          M_fine = if (time == "continuous") M_fine else NULL)

    list(X = X, A = A, Z = Z, T_obs = T_obs, D = D, pi_x = as.vector(pi_x),
         psi1_true = truth$psi1, psi0_true = truth$psi0,
         ate_true = truth$psi1 - truth$psi0,
         rmst1_true = truth$rmst1, rmst0_true = truth$rmst0,
         horizon = horizon, grid = grid, n_steps = if (!is.null(grid)) length(grid) else NULL)
  }

  structure(
    list(name = name, hazard = hazard, censor_hazard = censor_hazard,
         propensity = propensity, covariance = covariance,
         draw_X = draw_X,
         horizon = horizon, time = time, grid = grid,
         generate = generate),
    class = "survival_dgp")
}

#' Compute true E[S_t(w, X)] and E[RMST(w, X)] for w in {0, 1}
#' by numerical integration over a sample Z.
compute_truth = function(hazard, Z, horizon, time, grid = NULL, M_fine = 1000) {
  n = nrow(Z)
  wt = function(z, a) { z2 = z; z2[, 1] = a; z2 }
  Z1 = wt(Z, 1); Z0 = wt(Z, 0)

  if (time == "continuous") {
    du = horizon / M_fine
    fg = (1:M_fine) * du
    ch1 = ch0 = rep(0, n)
    rmst1 = rmst0 = rep(0, n)
    S1_prev = S0_prev = rep(1, n)
    for (k in 1:M_fine) {
      ch1 = ch1 + hazard(fg[k], Z1) * du
      ch0 = ch0 + hazard(fg[k], Z0) * du
      S1_curr = exp(-ch1); S0_curr = exp(-ch0)
      rmst1 = rmst1 + (S1_prev + S1_curr) / 2 * du  # trapezoidal
      rmst0 = rmst0 + (S0_prev + S0_curr) / 2 * du
      S1_prev = S1_curr; S0_prev = S0_curr
    }
    list(psi1 = mean(exp(-ch1)), psi0 = mean(exp(-ch0)),
         rmst1 = mean(rmst1), rmst0 = mean(rmst0))
  } else {
    S1 = S0 = rep(1, n)
    rmst1 = rmst0 = rep(0, n)
    du = if (length(grid) > 1) grid[2] - grid[1] else grid[1]
    for (k in seq_along(grid)) {
      S1 = S1 * (1 - hazard(grid[k], Z1))
      S0 = S0 * (1 - hazard(grid[k], Z0))
      rmst1 = rmst1 + S1 * du
      rmst0 = rmst0 + S0 * du
    }
    list(psi1 = mean(S1), psi0 = mean(S0),
         rmst1 = mean(rmst1), rmst0 = mean(rmst0))
  }
}


# ============================================================
# Continuous -> Discrete
# ============================================================

#' Convert a continuous-time DGP to a discrete-time DGP on an equally-spaced grid.
#'
#' The discrete hazard probability at bin k is exact:
#'   h_k = 1 - exp(-int_{u_{k-1}}^{u_k} lambda(s) ds)
#' where the integral is computed by Simpson's rule with Q_bin quadrature points.
#'
#' For the first bin (lower bound = 0), uses a squared grid t = du * s^2
#' to handle singularities like lambda ~ 1/sqrt(t) (Weibull shape < 1).
discretize = function(dgp, n_steps, Q_bin = 20) {
  stopifnot(dgp$time == "continuous")
  horizon = dgp$horizon
  du = horizon / n_steps
  grid = (1:n_steps) * du

  integrate_bin = function(rate_fn, u, Z) {
    k = which.min(abs(grid - u))
    u_lo = if (k == 1) 0 else grid[k - 1]
    u_hi = grid[k]
    n_z = nrow(Z)

    if (u_lo == 0) {
      # Substitution s in [0, 1], t = u_hi * s^2, dt = 2 * u_hi * s * ds
      # Removes 1/sqrt(t) singularity: 2s * rate(u_hi * s^2) -> finite at s=0
      s_pts = seq(0, 1, length.out = Q_bin + 1)
      ds = s_pts[2] - s_pts[1]
      sw = rep(1, Q_bin + 1); sw[1] = sw[Q_bin + 1] = 1
      for (j in 2:Q_bin) sw[j] = if (j %% 2 == 0) 4 else 2
      integral = rep(0, n_z)
      for (j in seq_along(s_pts)) {
        s = s_pts[j]
        if (s < 1e-15) next
        integral = integral + sw[j] * 2 * s * rate_fn(u_hi * s^2, Z)
      }
      return(integral * u_hi * ds / 3)
    }

    pts = seq(u_lo, u_hi, length.out = Q_bin + 1)
    dp = pts[2] - pts[1]
    sw = rep(1, Q_bin + 1); sw[1] = sw[Q_bin + 1] = 1
    for (j in 2:Q_bin) sw[j] = if (j %% 2 == 0) 4 else 2
    integral = rep(0, n_z)
    for (j in seq_along(pts)) integral = integral + sw[j] * rate_fn(pts[j], Z)
    integral * dp / 3
  }

  disc_hazard = function(u, Z) 1 - exp(-integrate_bin(dgp$hazard, u, Z))
  disc_censor = function(u, Z) 1 - exp(-integrate_bin(dgp$censor_hazard, u, Z))

  disc_dgp = make_dgp(name = paste0(dgp$name, "_disc", n_steps),
                      hazard = disc_hazard, censor_hazard = disc_censor,
                      propensity = dgp$propensity, covariance = dgp$covariance,
                      horizon = horizon, time = "discrete", grid = grid,
                      draw_X = dgp$draw_X)

  # generate_from: take continuous data, bin times to discrete grid.
  # Same subjects, same events, just coarser time resolution.
  disc_dgp$generate_from = function(cts_dat) {
    # Bin continuous T_obs to the discrete grid.
    # Subjects with T_obs > horizon keep T_obs = Inf (admin censored at horizon).
    idx = findInterval(cts_dat$T_obs, c(0, grid), rightmost.closed = TRUE)
    idx = pmin(idx, n_steps)
    T_disc = grid[idx]
    T_disc[cts_dat$T_obs > horizon] = Inf

    # Recompute truth on this sample's Z
    truth = compute_truth(disc_hazard, cts_dat$Z, horizon, "discrete", grid)

    # People who survived past horizon are admin censored
    D_disc = cts_dat$D
    D_disc[cts_dat$T_obs > horizon] = 0

    list(X = cts_dat$X, A = cts_dat$A, Z = cts_dat$Z,
         T_obs = T_disc, D = D_disc, pi_x = cts_dat$pi_x,
         psi1_true = truth$psi1, psi0_true = truth$psi0,
         ate_true = truth$psi1 - truth$psi0,
         rmst1_true = truth$rmst1, rmst0_true = truth$rmst0,
         horizon = horizon, grid = grid, n_steps = n_steps)
  }

  disc_dgp
}


# ============================================================
# Discrete -> Continuous
# ============================================================

#' Convert a discrete-time DGP to a continuous-time DGP via monotone spline
#' interpolation of the cumulative hazard.
#'
#' 1. Compute cumhaz at grid points: Lambda_k = sum_{j<=k} (-log(1 - h_j))
#' 2. Fit monotone cubic Hermite spline (Hyman method) through the cumhaz.
#' 3. Continuous hazard rate is lambda(u) = Lambda'(u).
#'
#' Round-trip exact: int_{k-1}^k Lambda'(s) ds = Lambda_k - Lambda_{k-1} = -log(1-h_k).
continuize = function(dgp) {
  stopifnot(dgp$time == "discrete")
  grid = dgp$grid
  horizon = dgp$horizon
  T_grid = length(grid)
  u_knots = c(0, grid)

  make_spline_fn = function(rate_fn) {
    cache = new.env(parent = emptyenv())
    cache$Z_id = NULL
    cache$splines = NULL

    build = function(Z) {
      n = nrow(Z)
      Lambda = matrix(0, n, T_grid + 1)
      for (k in 1:T_grid) {
        hk = pmin(pmax(rate_fn(grid[k], Z), 1e-15), 1 - 1e-15)
        Lambda[, k + 1] = Lambda[, k] - log(1 - hk)
      }
      lapply(1:n, function(i) splinefun(u_knots, Lambda[i, ], method = "hyman"))
    }

    function(u, Z) {
      zid = attr(Z, "dgp_id")
      if (is.null(zid)) zid = sum(Z[1:min(3, nrow(Z)), ])
      if (!identical(zid, cache$Z_id)) {
        cache$Z_id = zid
        cache$splines = build(Z)
      }
      vapply(seq_along(cache$splines), function(i) {
        max(cache$splines[[i]](u, deriv = 1), 0)
      }, numeric(1))
    }
  }

  make_dgp(name = paste0(dgp$name, "_cts"),
           hazard = make_spline_fn(dgp$hazard),
           censor_hazard = make_spline_fn(dgp$censor_hazard),
           propensity = dgp$propensity, covariance = dgp$covariance,
           horizon = horizon, time = "continuous",
           draw_X = dgp$draw_X)
}


# ============================================================
# Efficient variance (semiparametric efficiency bound)
# ============================================================

#' Compute Var(EIF) at true nuisance parameters via Monte Carlo.
#'
#' The EIF for psi = E[f(1,X) - f(-1,X)] is:
#'   phi_i = [f(1,X_i) - f(-1,X_i) - psi] + int_0^t gamma(u, W_i, X_i) dM_i(u)
#'
#' Discrete time:
#'   gamma(u, W, X) = -W * S_t(W,X) / [(1-h_u(W,X)) * pi(W|X) * G_u(W,X)]
#'
#' Continuous time:
#'   gamma(u, W, X) = -W * S_t(W,X) / [pi(W|X) * S_u(W,X) * G_u^C(W,X)]
#'
#' For RMST: S_t(W,X) replaced by tail sum/integral of S_s from u to t.
#'
#' @return List with eff_var, eff_se, psi, phi, rel_se.
efficient_variance = function(dgp, dat,
                              estimand = c("survival.probability", "rmst"),
                              M_quad = 200) {
  estimand = match.arg(estimand)
  time = dgp$time
  hazard_fn = dgp$hazard
  censor_hazard_fn = dgp$censor_hazard
  n = length(dat$T_obs)
  horizon = dgp$horizon
  Z = dat$Z
  W = 2 * dat$A - 1
  pi_w = ifelse(W == 1, dat$pi_x, 1 - dat$pi_x)

  wt = function(z, a) { z2 = z; z2[, 1] = a; z2 }
  Z1 = wt(Z, 1); Z0 = wt(Z, 0)

  if (time == "discrete") {
    grid = dgp$grid
    M = length(grid)

    S1 = S0 = matrix(0, n, M)
    S1_prev = S0_prev = rep(1, n)
    for (k in 1:M) {
      S1[, k] = S1_prev * (1 - hazard_fn(grid[k], Z1))
      S0[, k] = S0_prev * (1 - hazard_fn(grid[k], Z0))
      S1_prev = S1[, k]; S0_prev = S0[, k]
    }

    h_obs = matrix(0, n, M)
    for (k in 1:M) h_obs[, k] = hazard_fn(grid[k], Z)

    Surv_event = Surv_censor = matrix(1, n, M)
    if (M >= 2) {
      for (k in 2:M) {
        Surv_event[, k] = Surv_event[, k-1] * (1 - h_obs[, k-1])
        Surv_censor[, k] = Surv_censor[, k-1] * (1 - censor_hazard_fn(grid[k-1], Z))
      }
    }
    at_risk = Surv_event * Surv_censor

    du = if (M > 1) grid[2] - grid[1] else grid[1]
    if (estimand == "survival.probability") {
      f1 = S1[, M]; f0 = S0[, M]
    } else {
      f1 = rowSums(S1) * du; f0 = rowSums(S0) * du
    }
    direct = f1 - f0
    psi = mean(direct)

    is1 = W == 1
    S_obs_treat = S1 * is1 + S0 * (!is1)
    if (estimand == "survival.probability") {
      numer = matrix(S_obs_treat[, M], n, M) / (1 - h_obs)
    } else {
      tail_sum = matrix(0, n, M)
      tail_sum[, M] = S_obs_treat[, M]
      for (k in (M-1):1) tail_sum[, k] = tail_sum[, k+1] + S_obs_treat[, k]
      numer = tail_sum / (1 - h_obs)
    }

    gamma_grid = (-W / pi_w) * numer / at_risk

    T_eff = pmin(dat$T_obs, horizon)
    active = outer(T_eff, grid, ">=")
    compensator = rowSums(active * gamma_grid * h_obs)

    dirac = rep(0, n)
    ev_idx = which(dat$D == 1 & dat$T_obs <= horizon)
    if (length(ev_idx) > 0) {
      k_ev = match(dat$T_obs[ev_idx], grid)
      if (any(is.na(k_ev)))
        k_ev[is.na(k_ev)] = findInterval(dat$T_obs[ev_idx[is.na(k_ev)]], grid)
      k_ev = pmax(pmin(k_ev, M), 1)
      dirac[ev_idx] = gamma_grid[cbind(ev_idx, k_ev)]
    }

  } else {
    # --- Continuous time ---
    du = horizon / M_quad
    grid = (1:M_quad) * du
    M = M_quad

    S1 = S0 = matrix(0, n, M)
    ch1 = ch0 = rep(0, n)
    for (k in 1:M) {
      ch1 = ch1 + hazard_fn(grid[k], Z1) * du
      ch0 = ch0 + hazard_fn(grid[k], Z0) * du
      S1[, k] = exp(-ch1); S0[, k] = exp(-ch0)
    }

    lambda_grid = ev_ch = cen_ch = matrix(0, n, M)
    run_ev = run_cen = rep(0, n)
    for (k in 1:M) {
      lk = hazard_fn(grid[k], Z)
      ck = censor_hazard_fn(grid[k], Z)
      run_ev = run_ev + lk * du; run_cen = run_cen + ck * du
      ev_ch[, k] = run_ev; cen_ch[, k] = run_cen
      lambda_grid[, k] = lk
    }

    at_risk = exp(-ev_ch - cen_ch)

    if (estimand == "survival.probability") {
      f1 = S1[, M]; f0 = S0[, M]
    } else {
      f1 = (1 + 2 * rowSums(S1[, -M, drop = FALSE]) + S1[, M]) * du / 2
      f0 = (1 + 2 * rowSums(S0[, -M, drop = FALSE]) + S0[, M]) * du / 2
    }
    direct = f1 - f0
    psi = mean(direct)

    is1 = W == 1
    S_obs = S1 * is1 + S0 * (!is1)
    if (estimand == "survival.probability") {
      numer = matrix(S_obs[, M], n, M)
    } else {
      numer = matrix(0, n, M)
      numer[, M] = S_obs[, M] * du
      for (k in (M-1):1) numer[, k] = numer[, k+1] + S_obs[, k] * du
    }

    gamma_grid = (-W / pi_w) * numer / at_risk

    T_eff = pmin(dat$T_obs, horizon)
    active = outer(T_eff, grid, ">=")
    compensator = rowSums(active * gamma_grid * lambda_grid) * du

    dirac = rep(0, n)
    ev_idx = which(dat$D == 1 & dat$T_obs <= horizon)
    if (length(ev_idx) > 0) {
      k_ev = pmax(pmin(findInterval(dat$T_obs[ev_idx], grid), M), 1)
      dirac[ev_idx] = gamma_grid[cbind(ev_idx, k_ev)]
    }
  }

  phi = direct - psi + dirac - compensator
  eff_var = mean(phi^2)

  # Relative SE of the variance estimate: sqrt(kurt(phi) - 1) / sqrt(n)
  kurtosis = mean(phi^4) / eff_var^2
  rel_se = sqrt(max(kurtosis - 1, 0)) / sqrt(n)

  invisible(list(eff_var = eff_var, eff_se = sqrt(eff_var), psi = psi,
                 phi = phi, mean_phi = mean(phi), rel_se = rel_se,
                 kurtosis = kurtosis, n_mc = n,
                 direct = direct, dirac = dirac, compensator = compensator))
}


# ============================================================
# Cached efficient SE
# ============================================================

#' Compute the semiparametric efficient SD for a DGP, with disk caching
#' and automatic MC sample size selection.
#'
#' @param dgp A DGP object with $name set.
#' @param p Number of covariates.
#' @param estimand "survival.probability" or "rmst".
#' @param target_rel_se Target relative SE of the variance estimate (default 0.01).
#' @param n_pilot Pilot sample size for kurtosis estimation (default 10000).
#' @param max_n_mc Maximum MC sample size (default 500000).
#' @param cache_dir Directory for cached results.
#' @param M_quad Quadrature grid for continuous time (default 200).
#' @param verbose Print progress (default TRUE).
#' @return List with eff_se, eff_var, psi, n_mc, rel_se, kurtosis.
efficient_se = function(dgp, p,
                        estimand = c("survival.probability", "rmst"),
                        target_rel_se = 0.01,
                        n_pilot = 10000,
                        max_n_mc = 500000,
                        cache_dir = NULL,
                        M_quad = 200,
                        verbose = TRUE) {
  estimand = match.arg(estimand)
  if (is.null(dgp$name) || dgp$name == "")
    stop("DGP must have a $name for caching. Use a named DGP constructor.")

  # Default cache dir: next to this file
  if (is.null(cache_dir))
    cache_dir = file.path(dirname(sys.frame(sys.nframe())$ofile %||% "."), "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  cache_key = paste(dgp$name, paste0("p", p), estimand, sep = "_")
  cache_file = file.path(cache_dir, paste0(cache_key, ".rds"))

  # Check cache
  if (file.exists(cache_file)) {
    cached = readRDS(cache_file)
    if (verbose)
      cat(sprintf("[cached] %s: eff_se = %.4f (n_mc=%d, rel_se=%.3f%%)\n",
                  cache_key, cached$eff_se, cached$n_mc, cached$rel_se * 100))
    return(invisible(cached))
  }

  # Pilot: estimate kurtosis to determine MC sample size
  if (verbose) cat(sprintf("[pilot] %s: n=%d ... ", cache_key, n_pilot))
  set.seed(12345)  # reproducible pilot
  dat_pilot = dgp$generate(n_pilot, p, M_fine = 2000)
  res_pilot = efficient_variance(dgp, dat_pilot, estimand = estimand, M_quad = M_quad)
  kurtosis = res_pilot$kurtosis

  # Required n for target relative SE: n = (sqrt(kurtosis - 1) / target)^2
  n_needed = ceiling((sqrt(max(kurtosis - 1, 2)) / target_rel_se)^2)
  n_mc = min(max(n_needed, n_pilot), max_n_mc)

  if (verbose)
    cat(sprintf("kurtosis=%.1f, E[phi]=%.1e, need n=%d\n", kurtosis, res_pilot$mean_phi, n_mc))

  # Full run if pilot wasn't enough
  if (n_mc > n_pilot) {
    if (verbose) cat(sprintf("[full]  %s: n=%d ... ", cache_key, n_mc))
    set.seed(12345)
    dat = dgp$generate(n_mc, p, M_fine = 2000)
    res = efficient_variance(dgp, dat, estimand = estimand, M_quad = M_quad)
    if (verbose) cat(sprintf("done. E[phi]=%.1e\n", res$mean_phi))
  } else {
    res = res_pilot
  }

  out = list(eff_se = res$eff_se, eff_var = res$eff_var, psi = res$psi,
             n_mc = length(res$phi), rel_se = res$rel_se,
             kurtosis = res$kurtosis, mean_phi = res$mean_phi)

  saveRDS(out, cache_file)
  if (verbose)
    cat(sprintf("[saved] %s: eff_se = %.4f (rel_se=%.3f%%)\n",
                cache_key, out$eff_se, out$rel_se * 100))
  invisible(out)
}

# ============================================================
# Named survival DGPs
# ============================================================

# --- Our DGPs (logistic hazard rate) ---

#' Spinoff3 paper, Section 6.2: continuous time on [0, horizon].
#'
#' Event hazard rate: expit(alpha + beta_w * W + beta_x * 1'X + gamma_t * u)
#' Censoring rate:    expit(alpha_c + beta_w_c * W + beta_wx2_c * W * 1'X^2 + gamma_t_c * u)
#' Propensity:        expit(prop_coef * 1'X)
#' Covariates:        X ~ N(0, Sigma), Sigma = 0.8*I + 0.2*11'
paper_cts_dgp = function(horizon = 1, alpha = -0.1, beta_w = 0.4,
                          beta_x = 0.3, gamma_t = 0.5,
                          alpha_c = -3.15, beta_w_c = -0.15,
                          beta_wx2_c = 0.1, gamma_t_c = 0.2,
                          prop_coef = 0.3) {
  make_dgp(
    name = sprintf("paper_cts_h%.2g_p%.2g", horizon, prop_coef),
    hazard = function(u, Z)
      plogis(alpha + beta_w * Z[,1] + beta_x * rowSums(Z[,2:ncol(Z),drop=FALSE]) + gamma_t * u),
    censor_hazard = function(u, Z)
      plogis(alpha_c + beta_w_c * Z[,1] + beta_wx2_c * rowSums(Z[,1] * Z[,2:ncol(Z),drop=FALSE]^2) + gamma_t_c * u),
    propensity = function(X) plogis(X %*% rep(prop_coef, ncol(X))),
    covariance = function(p) 0.8 * diag(p) + 0.2 * matrix(1, p, p),
    horizon = horizon, time = "continuous")
}

#' Pham-style discrete DGP: discrete time on {1, ..., t_max}.
#'
#' Event hazard: expit(alpha + beta_w * A + beta_x * 1'X + gamma_t * u/t_max)
#'   interpreted as a per-step PROBABILITY, not a rate.
pham_discrete_dgp = function(t_max = 8, n_steps = NULL,
                              alpha = -2, beta_w = 0.5,
                              beta_x = 0.2, gamma_t = 0.4,
                              alpha_c = -3, beta_w_c = -0.3,
                              beta_wx2_c = 0.2, gamma_t_c = 0.2) {
  if (is.null(n_steps)) n_steps = t_max
  grid = (1:n_steps) * (t_max / n_steps)
  make_dgp(
    name = sprintf("pham_disc_t%d_s%d", t_max, n_steps),
    hazard = function(u, Z)
      plogis(alpha + beta_w * Z[,1] + beta_x * rowSums(Z[,2:ncol(Z),drop=FALSE]) + gamma_t * u / t_max),
    censor_hazard = function(u, Z)
      plogis(alpha_c + beta_w_c * Z[,1] + beta_wx2_c * rowSums(Z[,1] * Z[,2:ncol(Z),drop=FALSE]^2) + gamma_t_c * u / t_max),
    propensity = function(X) plogis(X %*% rep(0.3, ncol(X))),
    covariance = function(p) 0.8 * diag(p) + 0.2 * matrix(1, p, p),
    horizon = t_max, time = "discrete", grid = grid)
}


# --- GRF DGPs (Cui et al. 2023) ---

#' GRF type2: Cox PH (Weibull shape 1/2) failure, Uniform(0,3) censoring.
#'
#' Failure: T = (-log U / exp(eta))^2, eta = X_1 + (-0.5 + X_2)*W
#'   Hazard: lambda(t) = exp(eta) / (2*sqrt(t))
#' Propensity: e(X) = (1 + dbeta(X_1, 2, 4)) / 4, range ~ [0.25, 0.75].
grf_type2_dgp = function(horizon = 1.2) {
  make_dgp(
    name = sprintf("grf_type2_h%.2g", horizon),
    hazard = function(u, Z) {
      W = Z[, 1]; X1 = Z[, 2]; X2 = Z[, 3]
      eta = X1 + (-0.5 + X2) * W
      exp(eta) / (2 * sqrt(pmax(u, 1e-10)))
    },
    censor_hazard = function(u, Z) {
      rep(1 / (3 - u), nrow(Z))  # Uniform(0,3)
    },
    propensity = function(X) (1 + dbeta(X[, 1], 2, 4)) / 4,
    horizon = horizon, time = "continuous",
    draw_X = function(n, p) matrix(runif(n * p), n, p))
}

#' GRF type1: AFT (log-normal) failure, Weibull (Cox PH) censoring.
#'
#' Failure: T = exp(mu + eps), eps ~ N(0,1)
#'   mu = -1.85 - 0.8*I(X_1<0.5) + 0.7*sqrt(X_2) + 0.2*X_3 + tau*W
#'   tau = 0.7 - 0.4*I(X_1<0.5) - 0.4*sqrt(X_2)
#' Propensity: e(X) = (1 + dbeta(X_1, 2, 4)) / 4.
grf_type1_dgp = function(horizon = 0.8) {
  make_dgp(
    name = sprintf("grf_type1_h%.2g", horizon),
    hazard = function(u, Z) {
      W = Z[, 1]; X1 = Z[, 2]; X2 = Z[, 3]; X3 = Z[, 4]
      I1 = as.numeric(X1 < 0.5)
      tau = 0.7 - 0.4 * I1 - 0.4 * sqrt(X2)
      mu = -1.85 - 0.8 * I1 + 0.7 * sqrt(X2) + 0.2 * X3 + tau * W
      z = log(pmax(u, 1e-10)) - mu
      dnorm(z) / (pmax(u, 1e-10) * pnorm(-z))
    },
    censor_hazard = function(u, Z) {
      W = Z[, 1]; X1 = Z[, 2]; X2 = Z[, 3]; X3 = Z[, 4]
      I1 = as.numeric(X1 < 0.5)
      eta_C = -1.75 - 0.5 * sqrt(X2) + 0.2 * X3 +
              (1.15 + 0.5 * I1 - 0.3 * sqrt(X2)) * W
      exp(eta_C) / (2 * sqrt(pmax(u, 1e-10)))
    },
    propensity = function(X) (1 + dbeta(X[, 1], 2, 4)) / 4,
    horizon = horizon, time = "continuous",
    draw_X = function(n, p) matrix(runif(n * p), n, p))
}


# ============================================================
# Binary-outcome DGPs
# ============================================================

#' Hainmueller design 3 (binary outcome).
#'
#' X ~ MVN(0, R) where R is the Hainmueller correlation matrix.
#' Propensity: probit, pi(X) = Phi(X'beta / sqrt(zeta)).
#' Outcome: mu_1(X) = expit((s^2 - 2)/3), mu_0(X) = expit(s/3),
#'   where s = X_1 + X_2 + X_3.
#'
#' @param n Sample size.
#' @param zeta Propensity scaling (default sqrt(10)).
#' @return List with X, W, Y, pi, mu1, mu0, tau, tsm1, tsm0, ate, rr, vcate.
generate_binary_hain = function(n, zeta = sqrt(10)) {
  Sigma_hain = matrix(c(2, 1, -1, 1, 1, -0.5, -1, -0.5, 1), nrow = 3)
  R_hain = cov2cor(Sigma_hain)

  X = matrix(rnorm(n * 3), n, 3) %*% chol(R_hain)
  linpred = X %*% c(1, 2, -1)
  pi_x = pnorm(as.numeric(linpred) / zeta)
  W = rbinom(n, 1, pi_x)
  s = X[, 1] + X[, 2] + X[, 3]
  mu1 = plogis((s^2 - 2) / 3)
  mu0 = plogis(s / 3)
  Y = rbinom(n, 1, W * mu1 + (1 - W) * mu0)
  list(X = X, W = W, Y = Y, pi = pi_x, mu1 = mu1, mu0 = mu0,
       tau = mu1 - mu0, tsm1 = mean(mu1), tsm0 = mean(mu0),
       ate = mean(mu1 - mu0), rr = mean(mu1) / mean(mu0), vcate = var(mu1 - mu0))
}

#' Hainmueller design 3 (continuous outcome).
#'
#' Same covariates and propensity as binary version.
#' Outcome: mu_1(X) = s^2, mu_0(X) = s, Y = W*mu_1 + (1-W)*mu_0 + N(0,1).
#'
#' @param n Sample size.
#' @param zeta Propensity scaling (default sqrt(10)).
#' @return List with X, W, Y, pi, mu1, mu0, tau, tsm1, ate, vcate.
generate_continuous_hain = function(n, zeta = sqrt(10)) {
  Sigma_hain = matrix(c(2, 1, -1, 1, 1, -0.5, -1, -0.5, 1), nrow = 3)
  R_hain = cov2cor(Sigma_hain)

  X = matrix(rnorm(n * 3), n, 3) %*% chol(R_hain)
  linpred = X %*% c(1, 2, -1)
  pi_x = pnorm(as.numeric(linpred) / zeta)
  W = rbinom(n, 1, pi_x)
  s = X[, 1] + X[, 2] + X[, 3]
  mu1 = s^2
  mu0 = s
  tau = mu1 - mu0
  Y = W * mu1 + (1 - W) * mu0 + rnorm(n)
  list(X = X, W = W, Y = Y, pi = pi_x, mu1 = mu1, mu0 = mu0, tau = tau,
       tsm1 = mean(mu1), ate = mean(tau), vcate = var(tau))
}
