#!/usr/bin/env Rscript
## Report for 50-rep run: sampling distributions, bootstrap comparison, t-stats.
##
## Usage: cd ~/work/balancing && Rscript scratch/report-50rep.R [file1.rds ...]
## If no files given, globs examples/ for coverage-results-paper_cts_2-n200-*-r*.rds

library(ggplot2)
library(svglite)

for (f in list.files("R", full.names = TRUE)) source(f)
source("scratch/plot-traces.R")   # with_annotations, width, whitebg.theme

# ============================================================
# Style
# ============================================================
midnight = '#073b4c'; pink = '#ef476f'; teal = '#118ab2'
green = '#06d6a0'; mustard = '#ffd166'
method_colors = c(
  "1S"            = midnight,
  "Pham_quad"     = pink,
  "Pham_normed"   = teal,
  "Pham_logistic" = green,
  "KBW_disc"      = mustard,
  "PI"            = "purple"
)
method_order = c("1S", "Pham_quad", "Pham_normed", "Pham_logistic", "KBW_disc", "PI")

# ============================================================
# Load data
# ============================================================
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  files = list.files("examples",
    pattern = "coverage-results-paper_cts_2-n200-.*-r[0-9]+-[0-9]+\\.rds$",
    full.names = TRUE)
  if (length(files) == 0) stop("No matching files found in examples/")
} else {
  files = args
}
cat(sprintf("Loading %d file(s)...\n", length(files)))

rlist  = lapply(files, readRDS)
r      = do.call(rbind, lapply(rlist, function(x) x$results))
params = rlist[[1]]$params
truths = params$truths

cat(sprintf("Loaded %d rows, %d unique reps\n", nrow(r), length(unique(r$rep))))
cat(sprintf("Truths: surv_prob=%.5f  rmst=%.5f  tsm1=%.5f\n",
            truths$surv_prob, truths$rmst, truths$tsm1))

# Fixed tuning only
rf = r[r$tuning %in% c("fixed", "none"), ]
rf$bias = rf$est - sapply(rf$estimand, function(e) truths[[e]])

# Summary table
for (en in c("surv_prob", "rmst", "tsm1")) {
  sub = rf[rf$estimand == en, ]
  if (nrow(sub) == 0) next
  cat(sprintf("\n=== %s (truth=%.5f) ===\n", en, truths[[en]]))
  for (m in method_order) {
    ms = sub[sub$method == m & is.finite(sub$est), ]
    if (nrow(ms) == 0) next
    cat(sprintf("  %-15s n=%d  bias=%+.4f  sd=%.4f  mean_se_eif=%.4f  cov_eif=%.0f%%  cov_boot=%.0f%%\n",
        m, nrow(ms), mean(ms$bias), sd(ms$est),
        mean(ms$se_eif, na.rm=TRUE),
        100 * mean(ms$covers_eif, na.rm=TRUE),
        100 * mean(ms$covers_boot_t, na.rm=TRUE)))
  }
}

# ============================================================
# Plot 1: Sampling distribution per estimand
# Faceted histogram, one panel per method, truth vline.
# ============================================================
cat("\n--- Sampling distribution plots ---\n")

for (en in c("surv_prob", "rmst", "tsm1")) {
  sub = rf[rf$estimand == en & is.finite(rf$est), ]
  if (nrow(sub) < 4) next
  tr  = truths[[en]]
  sub$method = factor(sub$method, levels = intersect(method_order, unique(sub$method)))

  # Efficient reference: N(truth, se_eff^2), se_eff = mean EIF SE across all methods
  se_eff = mean(sub$se_eif, na.rm = TRUE)
  xr  = range(sub$est)
  pad = diff(xr) * 0.3
  xs  = seq(xr[1] - pad, xr[2] + pad, length.out = 300)
  goal_df = data.frame(x = xs, y = dnorm(xs, tr, se_eff))

  p = ggplot(sub, aes(x = est, fill = method, color = method)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.2, bins = 20,
                   position = "identity") +
    geom_line(data = goal_df, aes(x = x, y = y), inherit.aes = FALSE,
              color = 'green', linewidth = 0.6, alpha = 0.7) +
    geom_vline(xintercept = tr, color = 'green', linewidth = 1.5, alpha = 0.8) +
    facet_wrap(~ method, ncol = 3) +
    scale_fill_manual(values = method_colors, guide = "none") +
    scale_color_manual(values = method_colors, guide = "none") +
    whitebg.theme(labs = TRUE)

  fname = sprintf("scratch/rep50-dist-%s.svg", en)
  svglite(fname, width = 9, height = 5)
  print(p)
  dev.off()
  cat(sprintf("  %s -> %s\n", en, fname))
}

# ============================================================
# Plot 2: Bootstrap comparison — per-rep EIF CI bars (purple, low)
# vs per-rep boot CI bars (orange, high), on top of histogram.
# ============================================================
cat("\n--- Bootstrap comparison plots ---\n")

boot_methods = intersect(method_order, c("Pham_quad", "Pham_normed", "KBW_disc"))

for (m in boot_methods) {
  sub = rf[rf$estimand == "surv_prob" & rf$method == m & is.finite(rf$est), ]
  if (nrow(sub) < 4) next
  tr = truths$surv_prob

  w_eif  = mean(sub$ci_width_eif,  na.rm = TRUE)
  w_boot = mean(sub$ci_width_boot_t, na.rm = TRUE)

  estimates = sub$est
  n   = length(estimates)
  ht  = default_height(estimates)
  delta = ht * 0.003

  two_arm_fn = local({
    w_e = w_eif; w_b = w_boot; d = delta
    function(dot.data) {
      eif_d  = transform(dot.data, y = y - d)
      boot_d = transform(dot.data, y = y + d)
      lw = 0.6
      list(
        geom_segment(aes(x=x-w_e/2, xend=x+w_e/2, y=y, yend=y, alpha=I(alpha)),
                     color='purple', linewidth=lw, data=eif_d),
        geom_segment(aes(x=x-w_b/2, xend=x+w_b/2, y=y, yend=y, alpha=I(alpha)),
                     color='darkorange', linewidth=lw, data=boot_d),
        geom_point(aes(y=y, x=x, alpha=I(alpha)), color='black',
                   size=2, data=dot.data)
      )
    }
  })

  p = ggplot() +
    geom_histogram(aes(x = estimates, y = after_stat(density)),
                   fill = method_colors[[m]], alpha = 0.15, bins = 20) +
    with_annotations(estimates, estimand = tr, dots = n,
                     arm.fn = two_arm_fn, highlight.dots = 0,
                     unhighlighted.dot.alpha = 0.5,
                     height = ht * 0.7) +
    whitebg.theme(labs = TRUE)

  fname = sprintf("scratch/rep50-boot-surv_prob-%s.svg", m)
  svglite(fname, width = 7, height = 4)
  print(p)
  dev.off()
  cat(sprintf("  %s -> %s\n", m, fname))
}

# ============================================================
# Plot 3: T-stat distribution per estimand
# Faceted histogram vs N(0,1). Data-driven xlim.
# ============================================================
cat("\n--- T-stat plots ---\n")

for (en in c("surv_prob", "rmst", "tsm1")) {
  sub = rf[rf$estimand == en & is.finite(rf$est), ]
  if (nrow(sub) < 4) next
  tr = truths[[en]]
  sub$method = factor(sub$method, levels = intersect(method_order, unique(sub$method)))

  # EIF t-stats
  sub_t = sub[is.finite(sub$se_eif) & sub$se_eif > 0, ]
  sub_t$tstat = (sub_t$est - tr) / sub_t$se_eif

  if (nrow(sub_t) >= 4) {
    xrng = quantile(sub_t$tstat, c(0.01, 0.99), na.rm = TRUE)
    xpad = diff(xrng) * 0.2
    xlims = c(min(-3.5, xrng[1] - xpad), max(3.5, xrng[2] + xpad))
    xs = seq(xlims[1], xlims[2], length.out = 300)
    norm_df = data.frame(x = xs, y = dnorm(xs))

    p = ggplot(sub_t, aes(x = tstat, fill = method, color = method)) +
      geom_histogram(aes(y = after_stat(density)), alpha = 0.2, bins = 20,
                     position = "identity") +
      geom_line(data = norm_df, aes(x = x, y = y), inherit.aes = FALSE,
                color = 'green', linewidth = 1.5, alpha = 0.8) +
      facet_wrap(~ method, ncol = 3) +
      coord_cartesian(xlim = xlims) +
      scale_fill_manual(values = method_colors, guide = "none") +
      scale_color_manual(values = method_colors, guide = "none") +
      whitebg.theme(labs = TRUE)

    fname = sprintf("scratch/rep50-tstat-%s.svg", en)
    svglite(fname, width = 9, height = 5)
    print(p)
    dev.off()
    cat(sprintf("  %s (EIF) -> %s\n", en, fname))
  }

  # Boot t-stats
  sub_b = sub[is.finite(sub$se_boot) & sub$se_boot > 0, ]
  sub_b$tstat = (sub_b$est - tr) / sub_b$se_boot

  if (nrow(sub_b) >= 4) {
    xrng = quantile(sub_b$tstat, c(0.01, 0.99), na.rm = TRUE)
    xpad = diff(xrng) * 0.2
    xlims = c(min(-3.5, xrng[1] - xpad), max(3.5, xrng[2] + xpad))
    xs = seq(xlims[1], xlims[2], length.out = 300)
    norm_df = data.frame(x = xs, y = dnorm(xs))

    p2 = ggplot(sub_b, aes(x = tstat, fill = method, color = method)) +
      geom_histogram(aes(y = after_stat(density)), alpha = 0.2, bins = 20,
                     position = "identity") +
      geom_line(data = norm_df, aes(x = x, y = y), inherit.aes = FALSE,
                color = 'green', linewidth = 1.5, alpha = 0.8) +
      facet_wrap(~ method, ncol = 3) +
      coord_cartesian(xlim = xlims) +
      scale_fill_manual(values = method_colors, guide = "none") +
      scale_color_manual(values = method_colors, guide = "none") +
      whitebg.theme(labs = TRUE)

    fname2 = sprintf("scratch/rep50-tstat-boot-%s.svg", en)
    svglite(fname2, width = 9, height = 5)
    print(p2)
    dev.off()
    cat(sprintf("  %s (boot) -> %s\n", en, fname2))
  }
}

cat("\nDone.\n")
