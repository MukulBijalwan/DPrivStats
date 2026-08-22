#' DP two-sample t-test
#'
#' Differentially private two-sample t-test. The privacy budget is split
#' evenly across the DP group counts, means and variances; the test statistic
#' is computed entirely by post-processing of privatized quantities.
#'
#' @param x,y Numeric vectors (NAs ignored).
#' @param epsilon Total privacy parameter epsilon (> 0).
#' @param bounds Common bounds c(L, U) assumed known for both samples.
#' @param alternative Character, \code{"two.sided"}, \code{"less"} or
#'   \code{"greater"}.
#' @param delta Optional Gaussian-mechanism delta; if NULL Laplace is used.
#' @return An object of class \code{c("dp_htest", "htest")} with elements
#'   \code{statistic}, \code{estimate}, \code{std.error}, \code{p.value},
#'   \code{alternative}, and \code{epsilon}.
#' @examples
#' set.seed(7)
#' dp_t_test(rnorm(200), rnorm(200, 0.5), 1.0, bounds = c(-5, 5))
#' @export
dp_t_test <- function(x, y, epsilon,
                      bounds = c(min(c(x, y), na.rm = TRUE),
                                 max(c(x, y), na.rm = TRUE)),
                      alternative = c("two.sided", "less", "greater"),
                      delta = NULL) {
  validate_bounds(bounds)
  alternative <- match.arg(alternative)
  eps_each <- epsilon / 3
  x <- as.numeric(stats::na.omit(x))
  y <- as.numeric(stats::na.omit(y))
  if (length(x) < 2 || length(y) < 2) {
    stop("Both samples must contain at least 2 observations.", call. = FALSE)
  }
  use_gauss <- !is.null(delta)

  n_x_priv <- max(laplace_mechanism(length, x, eps_each, 1), 2)
  n_y_priv <- max(laplace_mechanism(length, y, eps_each, 1), 2)
  # post-processing: counts floored at 2 to avoid division blow-ups

  mech <- if (use_gauss) "gaussian" else "laplace"
  mean_x_priv <- dp_mean(x, eps_each, bounds, mech, delta)$estimate
  mean_y_priv <- dp_mean(y, eps_each, bounds, mech, delta)$estimate
  var_x_priv <- dp_variance(x, eps_each, bounds, mech, delta)$estimate
  var_y_priv <- dp_variance(y, eps_each, bounds, mech, delta)$estimate

  diff_priv <- mean_x_priv - mean_y_priv
  se_priv <- sqrt(var_x_priv / n_x_priv + var_y_priv / n_y_priv)
  t_stat <- if (se_priv > 0) diff_priv / se_priv else 0

  df <- n_x_priv + n_y_priv - 2
  p_value <- switch(
    alternative,
    two.sided = 2 * stats::pt(-abs(t_stat), df = df),
    less = stats::pt(t_stat, df = df),
    greater = stats::pt(-t_stat, df = df, lower.tail = FALSE)
  )

  structure(
    list(
      statistic = c(t = t_stat),
      estimate = c(mean_difference = diff_priv,
                   mean_of_x = mean_x_priv,
                   mean_of_y = mean_y_priv),
      std.error = se_priv,
      p.value = p_value,
      alternative = alternative,
      parameter = c(df = df),
      method = "Differentially Private Two Sample t-test",
      epsilon = epsilon,
      delta = if (use_gauss) delta else 0
    ),
    class = c("dp_htest", "htest")
  )
}

#' DP chi-square test of independence
#'
#' Privatizes each cell of a contingency table with Laplace noise
#' (sensitivity 1 per cell), truncates negative counts (post-processing),
#' then computes the standard chi-square test on the privatized table.
#'
#' @param table A numeric matrix or table of non-negative counts.
#' @param epsilon Privacy parameter epsilon (> 0).
#' @return An object of class \code{c("dp_htest", "htest")} wrapping a
#'   standard chi-square test on the privatized table.
#' @examples
#' tab <- matrix(c(30, 20, 10, 40), nrow = 2)
#' dp_chisq_test(tab, 1.0)
#' @export
dp_chisq_test <- function(table, epsilon) {
  table <- as.matrix(as.data.frame(table))
  storage.mode(table) <- "double"
  if (!is.numeric(table) || any(table < 0)) {
    stop("`table` must be a numeric matrix/table of non-negative counts.",
         call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  priv_counts <- laplace_mechanism(function(t) as.vector(t), table, epsilon, 1)
  priv_counts <- pmax(priv_counts, 0) # post-processing: non-negative counts
  priv_table <- matrix(priv_counts, nrow = nrow(table),
                       dimnames = dimnames(table))
  res <- stats::chisq.test(round(priv_table))
  structure(
    list(
      statistic = res$statistic,
      parameter = res$parameter,
      p.value = res$p.value,
      method =
        "Chi-squared test on differentially private contingency table",
      data.name = deparse(substitute(table)),
      observed_private = priv_table,
      epsilon = epsilon,
      delta = 0
    ),
    class = c("dp_htest", "htest")
  )
}