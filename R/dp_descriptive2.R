#' DP variance
#'
#' Differentially private variance of a bounded variable, released as the DP
#' second moment minus the square of the DP mean under a split budget
#' (\eqn{\epsilon/2} each), with post-processing to ensure non-negativity.
#'
#' @inheritParams dp_mean
#' @param delta Privacy parameter for Gaussian mechanisms (if chosen).
#' @return An object of class \code{dp_estimate}.
#' @examples
#' set.seed(2)
#' dp_variance(runif(300), 1.0, bounds = c(0, 1))$estimate
#' @export
dp_variance <- function(x, epsilon, bounds = range(x, na.rm = TRUE),
                        mechanism = c("laplace", "gaussian"), delta = NULL) {
  validate_bounds(bounds)
  mechanism <- match.arg(mechanism)
  eps_half <- epsilon / 2
  xc <- clip_numeric(stats::na.omit(x), bounds)
  n <- length(xc)
  if (n == 0) stop("`x` contains no usable observations.", call. = FALSE)
  mean_obj <- dp_mean(xc, eps_half, bounds, mechanism, delta)
  sens_second <- diff(bounds)^2 / n # L1 sensitivity of mean(x^2)
  second_moment <- switch(
    mechanism,
    laplace = laplace_mechanism(function(z) z, mean(xc^2), eps_half,
                                sens_second),
    gaussian = {
      check_delta(delta)
      gaussian_mechanism(function(z) z, mean(xc^2), eps_half, delta,
                         sens_second)
    }
  )
  var_est <- max(second_moment - mean_obj$estimate^2, 0) # post-processing
  structure(
    list(
      estimate = var_est,
      private = TRUE,
      epsilon = epsilon,
      delta = if (mechanism == "laplace") 0 else delta,
      mechanism = mechanism,
      statistic = "variance",
      n = n
    ),
    class = "dp_estimate"
  )
}

#' DP quantile via the exponential mechanism
#'
#' Releases arbitrary quantiles \code{probs} of a bounded variable using the
#' exponential mechanism with score based on distance from the target order
#' statistic. Each quantile receives budget \eqn{\epsilon / |\text{probs}|}.
#'
#' @inheritParams dp_median
#' @param probs Numeric vector of probabilities in (0, 1).
#' @return An object of class \code{dp_estimate}; \code{estimate} is a named
#'   vector with one entry per probability in \code{probs}.
#' @examples
#' set.seed(3)
#' dp_quantile(rnorm(200), 1.0, bounds = c(-5, 5), probs = c(.25, .5, .75))
#' @export
dp_quantile <- function(x, epsilon, bounds, probs = 0.5, n_bins = 1000) {
  validate_bounds(bounds)
  if (!is.numeric(probs) || any(probs <= 0) || any(probs >= 1)) {
    stop("`probs` must be numeric values in (0, 1).", call. = FALSE)
  }
  x <- as.numeric(stats::na.omit(x))
  n <- length(x)
  if (n == 0) stop("`x` contains no usable observations.", call. = FALSE)
  t_grid <- seq(bounds[1], bounds[2], length.out = max(n_bins, 2L))
  ranks <- vapply(t_grid, function(t) sum(x <= t), numeric(1))
  eps_each <- epsilon / length(probs)
  estimates <- vapply(probs, function(p) {
    scores <- -(ranks - p * n)^2
    exponential_mechanism(t_grid, scores, eps_each, score_sensitivity = n)
  }, numeric(1))
  names(estimates) <- paste0(probs * 100, "%")
  structure(
    list(
      estimate = estimates,
      private = TRUE,
      epsilon = epsilon,
      delta = 0,
      mechanism = "exponential",
      statistic = "quantile",
      n = n,
      bounds = bounds
    ),
    class = "dp_estimate"
  )
}

#' DP histogram
#'
#' Releases bin counts of a bounded variable using Laplace noise with
#' sensitivity 1 per bin, followed by non-negativity post-processing and
#' optional normalization to probabilities.
#'
#' @param x Numeric vector (NAs ignored).
#' @param epsilon Privacy parameter epsilon (> 0).
#' @param breaks Numeric vector of cut points (monotone). Defaults to deciles
#'   of the observed data.
#' @param normalize Logical; return relative frequencies instead of counts.
#' @return An object of class \code{dp_estimate} with element \code{estimate}
#'   (named numeric vector of privatized bin counts/proportions).
#' @examples
#' set.seed(4)
#' fit <- dp_histogram(rnorm(500), 1.0, breaks = seq(-3, 3, by = 1))
#' fit$estimate
#' @export
dp_histogram <- function(x, epsilon, breaks = NULL, normalize = FALSE) {
  x <- as.numeric(stats::na.omit(x))
  if (is.null(breaks)) breaks <- seq(min(x), max(x), length.out = 11)
  if (!is.numeric(breaks) || length(breaks) < 2 || any(diff(breaks) <= 0)) {
    stop("`breaks` must be an increasing numeric vector of >= 2 values.",
         call. = FALSE)
  }
  bins <- cut(x, breaks = breaks, include.lowest = TRUE)
  counts <- as.numeric(table(bins))
  noisy <- laplace_mechanism(identity, counts, epsilon, l1_sensitivity = 1)
  noisy <- pmax(noisy, 0) # post-processing: non-negative counts
  if (normalize) noisy <- noisy / sum(noisy)
  structure(
    list(
      estimate = stats::setNames(noisy, levels(bins)),
      private = TRUE,
      epsilon = epsilon,
      delta = 0,
      mechanism = "laplace",
      statistic = "histogram",
      n = length(x),
      breaks = breaks
    ),
    class = "dp_estimate"
  )
}