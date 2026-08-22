#' DP Kolmogorov-Smirnov test
#'
#' Releases the ECDFs of two bounded samples on a common grid under DP
#' (ECDF sensitivity \eqn{1/n} per sample; budget split evenly between the two
#' samples) and computes the KS statistic from the privatized ECDFs.
#'
#' @inheritParams dp_t_test
#' @param n_grid Number of grid points for ECDF evaluation (default 100).
#' @return An object of class \code{c("dp_htest", "htest")}.
#' @examples
#' set.seed(8)
#' dp_ks_test(rnorm(150), rnorm(150, 1), 1.0, bounds = c(-6, 6))
#' @export
dp_ks_test <- function(x, y, epsilon,
                       bounds = c(min(c(x, y), na.rm = TRUE),
                                  max(c(x, y), na.rm = TRUE)),
                       n_grid = 100) {
  validate_bounds(bounds)
  x <- clip_numeric(stats::na.omit(x), bounds)
  y <- clip_numeric(stats::na.omit(y), bounds)
  nx <- length(x); ny <- length(y)
  if (nx < 1 || ny < 1) {
    stop("Both samples must contain observations.", call. = FALSE)
  }
  grid <- seq(bounds[1], bounds[2], length.out = max(n_grid, 10L))
  fx_true <- vapply(grid, function(t) mean(x <= t), numeric(1))
  fy_true <- vapply(grid, function(t) mean(y <= t), numeric(1))
  fx_priv <- pmax(pmin(laplace_mechanism(identity, fx_true, epsilon / 2,
                                         1 / nx), 1), 0)
  fy_priv <- pmax(pmin(laplace_mechanism(identity, fy_true, epsilon / 2,
                                         1 / ny), 1), 0)
  d_stat <- max(abs(fx_priv - fy_priv))
  ne <- nx * ny / (nx + ny)
  p_value <- asymptotic_ks_pvalue(d_stat * sqrt(ne))
  structure(
    list(
      statistic = c(D = d_stat),
      p.value = p_value,
      method = "Differentially Private Two-Sample Kolmogorov-Smirnov Test",
      epsilon = epsilon,
      delta = 0
    ),
    class = c("dp_htest", "htest")
  )
}

#' DP one-way ANOVA F-test
#'
#' Privatizes group counts, means and variances (equal budget splits) and
#' computes the one-way ANOVA F statistic from privatized summaries only.
#'
#' @param formula Formula of the form \code{response ~ group}.
#' @param data Data frame.
#' @inheritParams dp_t_test
#' @return An object of class \code{c("dp_htest", "htest")}.
#' @examples
#' dat <- data.frame(y = c(rnorm(60), rnorm(60, 1)),
#'                   g = factor(rep(c("A", "B"), each = 60)))
#' dp_anova(y ~ g, dat, 2.0, bounds = c(-5, 5))
#' @export
dp_anova <- function(formula, data, epsilon, bounds, delta = NULL) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula: response ~ group.", call. = FALSE)
  }
  mf <- stats::model.frame(formula, data)
  ok <- stats::complete.cases(mf)
  validate_bounds(bounds)
  y_all <- as.numeric(mf[[1]][ok])
  grp_all <- droplevels(factor(mf[[2]][ok]))
  groups <- levels(grp_all)
  k <- length(groups)
  if (k < 2) stop("Need at least 2 groups.", call. = FALSE)

  eps_each <- epsilon / (3 * k)
  mech <- if (!is.null(delta)) "gaussian" else "laplace"
  n_priv <- numeric(k); m_priv <- numeric(k); v_priv <- numeric(k)
  for (i in seq_len(k)) {
    xi <- y_all[grp_all == groups[i]]
    n_priv[i] <- max(laplace_mechanism(length, xi, eps_each, 1), 2)
    m_priv[i] <- dp_mean(xi, eps_each, bounds, mech, delta)$estimate
    v_priv[i] <- dp_variance(xi, eps_each, bounds, mech, delta)$estimate
  }
  grand_m <- sum(n_priv * m_priv) / sum(n_priv)
  ss_between <- sum(n_priv * (m_priv - grand_m)^2)
  ss_within <- sum((n_priv - 1) * v_priv)
  df_b <- k - 1
  df_w <- sum(n_priv) - k
  f_stat <- (ss_between / df_b) /
    max(ss_within / df_w, .Machine$double.eps)
  p_value <- stats::pf(f_stat, df_b, df_w, lower.tail = FALSE)
  structure(
    list(
      statistic = c(F = f_stat),
      parameter = c(df1 = df_b, df2 = df_w),
      p.value = p_value,
      method = "Differentially Private One-Way ANOVA",
      group_means_private = stats::setNames(m_priv, groups),
      epsilon = epsilon,
      delta = if (!is.null(delta)) delta else 0
    ),
    class = c("dp_htest", "htest")
  )
}