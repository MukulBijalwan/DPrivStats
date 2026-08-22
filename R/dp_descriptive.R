#' DP mean
#'
#' Differentially private sample mean of a bounded variable. Data are clipped
#' to \code{bounds} before the statistic is computed.
#'
#' @param x Numeric vector (NAs ignored).
#' @param epsilon Privacy parameter epsilon (> 0).
#' @param bounds Bounds c(L, U) assumed known; data are clipped to them.
#'   Defaults to the observed range (use a priori bounds in practice).
#' @param mechanism One of \code{"laplace"}, \code{"gaussian"},
#'   \code{"analytic_gaussian"}.
#' @param delta Required (in (0,1)) when a Gaussian mechanism is used.
#' @return An object of class \code{dp_estimate}.
#' @examples
#' set.seed(1)
#' dp_mean(rnorm(500), 1.0, bounds = c(-5, 5))$estimate
#' @export
dp_mean <- function(x, epsilon, bounds = range(x, na.rm = TRUE),
                    mechanism = c("laplace", "gaussian",
                                  "analytic_gaussian"),
                    delta = NULL) {
  validate_bounds(bounds)
  x <- as.numeric(stats::na.omit(x))
  n <- length(x)
  if (n == 0) stop("`x` contains no usable observations.", call. = FALSE)
  l1_sens <- diff(bounds) / n
  true_mean <- mean(clip_numeric(x, bounds))
  mech <- match.arg(mechanism)
  est <- switch(
    mech,
    laplace = laplace_mechanism(function(z) z, true_mean, epsilon, l1_sens),
    gaussian = {
      check_delta(delta)
      gaussian_mechanism(function(z) z, true_mean, epsilon, delta, l1_sens)
    },
    analytic_gaussian = {
      check_delta(delta)
      analytic_gaussian_mechanism(function(z) z, true_mean, epsilon, delta,
                                  l1_sens)
    }
  )
  sigma2 <- if (mech == "laplace") (l1_sens / epsilon)^2 else {
    check_delta(delta)
    (analytic_gaussian_sigma(epsilon, delta, l1_sens))^2
  }
  structure(
    list(
      estimate = unname(est),
      private = TRUE,
      epsilon = epsilon,
      delta = if (mech == "laplace") 0 else delta,
      mechanism = mech,
      statistic = "mean",
      n = n,
      sensitivity = list(l1 = unname(l1_sens), l2 = unname(l1_sens)),
      noise_variance = sigma2
    ),
    class = "dp_estimate"
  )
}

check_delta <- function(delta) {
  if (!is.numeric(delta) || length(delta) != 1L || delta <= 0 || delta >= 1) {
    stop("`delta` must be a numeric scalar in (0, 1) for Gaussian mechanisms.",
         call. = FALSE)
  }
  invisible(TRUE)
}