#!/usr/bin/env Rscript
## Plot optimizer traces from a saved log_data object.
##
## Usage:
##   x = readRDS("examples/simple-n200-rep001.rds")
##   source("scratch/plot-traces.R")
##   plot_delta(x$results$surv_prob$log_data, "gamma")
##   plot_delta(x$results$surv_prob$log_data, "lambda")

library(ggplot2)
library(svglite)

# From qtm285-1/shared-code.qmd
default.seed = 0
default.dot.size = 2
default.highlight.colors = c('#073b4c', '#ef476f', '#118ab2', '#06d6a0')
default_height = function(samples) { dnorm(0, mean=0, sd=sd(samples)) }

with_seed = function(thunk, seed) {
  set.seed(seed)
  thunk()
}

annotation_dot_heights = function(dots, highlight.dots, height, seed=default.seed) {
  (\() {
    highlight.dot.heights = sample(1:highlight.dots) / (highlight.dots + 1)
    unhighlighted.dot.heights = setdiff(sample(1:dots), 1:highlight.dots) / (dots + 1)
    height * c(highlight.dot.heights[1:min(dots,highlight.dots)],
               unhighlighted.dot.heights)
  }) |> with_seed(seed)
}

with_annotations = function(samples, estimand = NA,
                             shade.center = mean(samples),
                             shade.width = NA, shade.color = 'green', shade.alpha = .15,
                             dots = 100, interval.width = NA,
                             arm.fn = NULL,
                             dot.size = default.dot.size,
                             highlight.dots = 3,
                             highlighted.dot.alpha = 1,
                             unhighlighted.dot.alpha = .15,
                             height = default_height(samples), vline.alpha = .8,
                             seed = default.seed) {
  highlight.colors = if (highlight.dots) default.highlight.colors[1:min(dots,highlight.dots)] else c()
  summary.data = data.frame(mus = mean(samples), mu = shade.center,
                            w66 = width(samples, alpha=1/3),
                            w95 = width(samples, alpha=.05),
                            shade.width = shade.width, estimand = estimand)
  estimand.geom = if (is.na(summary.data$estimand)) geom_blank() else
    geom_vline(aes(xintercept=estimand), color='green', linewidth=1.5, alpha=.8, data=summary.data)
  shade.geom = if (is.na(summary.data$shade.width)) geom_blank() else
    geom_rect(aes(xmin=mu-shade.width/2, xmax=mu+shade.width/2, ymin=0, ymax=Inf),
              fill=shade.color, alpha=shade.alpha, data=summary.data)
  dot.geom = if (dots == 0) geom_blank() else {
    heights = annotation_dot_heights(dots, highlight.dots, height, seed)
    dot.data = data.frame(
      x = samples[1:dots], y = heights, w = interval.width,
      color = c(highlight.colors, rep('purple', dots - length(highlight.colors))),
      alpha = c(rep(highlighted.dot.alpha,   length(highlight.colors)),
                rep(unhighlighted.dot.alpha, dots - length(highlight.colors))))
    if (!is.null(arm.fn))
      arm.fn(dot.data)
    else if (any(is.na(interval.width)))
      geom_point(aes(y=y, x=x, color=I(color), alpha=I(alpha)), size=dot.size, data=dot.data)
    else
      geom_pointrange(aes(y=y, x=x, xmin=x-w/2, xmax=x+w/2, color=I(color), alpha=I(alpha)),
                      size=dot.size/4, linewidth=dot.size/4, data=dot.data)
  }
  vlines = if (is.na(vline.alpha)) geom_blank() else {
    list(geom_vline(aes(xintercept=mus), color="blue", data=summary.data, alpha=vline.alpha),
         geom_vline(aes(xintercept=mus-w66/2), color="blue", linetype='dashed',  data=summary.data, alpha=vline.alpha),
         geom_vline(aes(xintercept=mus+w66/2), color="blue", linetype='dashed',  data=summary.data, alpha=vline.alpha),
         geom_vline(aes(xintercept=mus-w95/2), color="blue", linetype='dotted',  data=summary.data, alpha=vline.alpha),
         geom_vline(aes(xintercept=mus+w95/2), color="blue", linetype='dotted',  data=summary.data, alpha=vline.alpha))
  }
  list(estimand.geom, shade.geom, dot.geom, vlines)
}

whitebg.theme = function(lightgray="#cccccc", gridcolor=rgb(0,0,0,.07,maxColorValue=1), labs=FALSE) {
  theme(
    plot.background        = element_rect(fill="transparent", colour=NA),
    panel.background       = element_rect(fill="transparent", colour=NA),
    legend.background      = element_rect(fill="transparent", colour=NA),
    legend.box.background  = element_rect(fill="transparent", colour=NA),
    legend.key             = element_rect(fill="transparent", colour=NA),
    axis.ticks.x           = element_blank(),
    axis.ticks.y           = element_blank(),
    axis.text.x            = element_text(colour=lightgray),
    axis.text.y            = element_text(colour=lightgray),
    axis.title.x           = if (labs) element_text(colour=lightgray) else element_blank(),
    axis.title.y           = if (labs) element_text(colour=lightgray, angle=90) else element_blank(),
    panel.grid.major       = element_line(colour=gridcolor),
    panel.grid.minor       = element_line(colour=gridcolor, linewidth=0.25))
}

width = function(samples, center=mean(samples), alpha=.05) {
  limits = c(0, 2*max(abs(samples - center)))
  uniroot(function(w) mean(center - w/2 <= samples & samples <= center + w/2) - (1-alpha),
          interval = limits)$root
}

# Extract all ::trace entries from a saved log$all() list.
# Returns a tidy data frame: key, solver, context, iter, gnorm, delta, accepted.
extract_traces = function(log_data) {
  trace_keys = grep("::trace$", names(log_data), value = TRUE)
  if (length(trace_keys) == 0) { message("No traces found."); return(NULL) }

  rows = lapply(trace_keys, function(k) {
    tr = log_data[[k]]
    if (is.null(tr) || nrow(tr) == 0) return(NULL)
    ctx     = sub("::trace$", "", k)
    solver  = if (grepl("::gamma::", k)) "gamma" else "lambda"
    boot_rep = {
      m = regmatches(ctx, regexpr("boot \\d+/\\d+", ctx))
      if (length(m)) m else "init"
    }
    fold = {
      m = regmatches(ctx, regexpr("fold \\d+/\\d+", ctx))
      if (length(m)) m else NA_character_
    }
    cbind(data.frame(key = k, solver = solver, boot_rep = boot_rep,
                     fold = fold, stringsAsFactors = FALSE),
          tr)
  })
  do.call(rbind, rows[!sapply(rows, is.null)])
}

# gnorm trajectories: main fit + bootstrap on same axes, faceted by solver.
# rep = "main" for the initial surv_dr solve, "boot k" for bootstrap reps.
plot_gnorm = function(log_data, out = NULL) {
  df = extract_traces(log_data)
  if (is.null(df)) return(invisible(NULL))
  df$rep = ifelse(df$boot_rep == "init", "main", df$boot_rep)

  df$alpha = ifelse(df$rep == "main", 1.0, 0.25)

  p = ggplot(df, aes(x = iter, y = gnorm, group = key,
                     color = rep, linetype = fold, alpha = I(alpha))) +
    geom_line() +
    geom_point(aes(shape = accepted), size = 0.8) +
    scale_y_log10() +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4), guide = "none") +
    facet_wrap(~ solver, scales = "free_y") +
    labs(x = "Newton iteration", y = "gradient norm",
         color = NULL, linetype = NULL) +
    theme_minimal(base_size = 11)

  if (!is.null(out)) { svglite(out, width = 8, height = 4); print(p); dev.off() }
  p
}

# Trust-region delta vs gnorm for one solver.
plot_delta = function(log_data, solver = "gamma", out = NULL) {
  df = extract_traces(log_data)
  if (is.null(df)) return(invisible(NULL))
  df = df[df$solver == solver, ]
  if (nrow(df) == 0) { message("No traces for solver: ", solver); return(NULL) }
  df$rep = ifelse(df$boot_rep == "init", "main", df$boot_rep)

  base = df[, c("key","rep","fold","iter","accepted")]

  rows = list(
    cbind(base, metric = "gnorm",   value = df$gnorm),
    cbind(base, metric = "delta",   value = df$delta))
  if (!all(is.na(df$bregman)))
    rows = c(rows, list(cbind(base, metric = "bregman", value = abs(df$bregman))))
  if (!all(is.na(df$wdelta)))
    rows = c(rows, list(cbind(base, metric = "wdelta",  value = df$wdelta)))
  long_diag = do.call(rbind, rows)
  long_diag$alpha = ifelse(long_diag$rep == "main", 1.0, 0.25)

  p_diag = ggplot(long_diag,
    aes(x = iter, y = value, group = interaction(key, metric),
        color = rep, linetype = metric, alpha = I(alpha))) +
    geom_line() +
    scale_y_log10() +
    scale_linetype_manual(values = c(gnorm = "solid", delta = "dashed",
                                     bregman = "dotted", wdelta = "longdash")) +
    labs(x = "iteration", y = "value", title = paste(solver, "— gnorm / delta"),
         color = NULL, linetype = NULL) +
    theme_minimal(base_size = 11)

  # val shifted so minimum = 1 (log-plottable), per key
  val_shifted = ave(df$val, df$key, FUN = function(v) v - min(v) + 1e-10)
  val_df = cbind(base, value = val_shifted)
  val_df$alpha = ifelse(val_df$rep == "main", 1.0, 0.25)

  p_val = ggplot(val_df,
    aes(x = iter, y = value, group = key, color = rep, alpha = I(alpha))) +
    geom_line() +
    scale_y_log10() +
    labs(x = "iteration", y = "val - min(val)", title = paste(solver, "— objective"),
         color = NULL) +
    theme_minimal(base_size = 11)

  if (!is.null(out)) { svglite(out, width = 10, height = 4); gridExtra::grid.arrange(p_diag, p_val, ncol = 2); dev.off() }
  invisible(list(diag = p_diag, val = p_val))
}

# Bootstrap sampling distribution.
# boot_ates: vector of bootstrap estimates.
# est:       point estimate (blue vline).
# truth:     true estimand (green vline + green normal overlay).
# eff_se:    semiparametric efficient SD (per-obs). Divide by sqrt(n) for estimator SD.
# n:         sample size.
plot_boot = function(boot_ates, est, truth, eff_se, n, out = NULL) {
  estimator_sd = eff_se / sqrt(n)
  df = data.frame(x = boot_ates)
  h = default_height(boot_ates)

  p = ggplot(df, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40) +
    stat_function(fun = dnorm, args = list(mean = truth, sd = estimator_sd),
                  color = "green", linewidth = 0.5) +
    stat_function(fun = dnorm, args = list(mean = est, sd = estimator_sd),
                  color = "#073b4c", linewidth = 0.5) +
    with_annotations(boot_ates, estimand = truth,
                     dots = length(boot_ates), highlight.dots = 0) +
    geom_vline(xintercept = est, color = "#073b4c", linewidth = 1, alpha = 0.8) +
    scale_x_continuous(limits = c(min(truth, est) - 3 * estimator_sd,
                                  max(truth, est) + 3 * estimator_sd)) +
    whitebg.theme(labs = TRUE)

  if (!is.null(out)) { svglite(out, width = 5, height = 4); print(p); dev.off() }
  p
}

# Bootstrap t-statistic sampling distribution.
# t_b = (boot_ates[b] - est) / boot_ses[b]. Should be N(0,1) if pivotal.
plot_tstat = function(boot_ates, boot_ses, est, out = NULL) {
  good = !is.na(boot_ates) & !is.na(boot_ses) & boot_ses > 0
  tstats = (boot_ates[good] - est) / boot_ses[good]
  df = data.frame(x = tstats)
  h = default_height(tstats)

  p = ggplot(df, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  color = "green", linewidth = 0.5) +
    with_annotations(tstats, estimand = 0, dots = length(tstats), highlight.dots = 0) +
    scale_x_continuous(limits = c(-4, 4)) +
    whitebg.theme(labs = TRUE)

  if (!is.null(out)) { svglite(out, width = 5, height = 4); print(p); dev.off() }
  p
}
