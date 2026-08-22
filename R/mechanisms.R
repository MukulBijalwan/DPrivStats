# ---------------------------------------------------------------- Laplace ----

#' Sample Laplace noise
#'
#' Draws \code{n} independent Laplace(0, scale) deviates using inverse-CDF
#' sampling with base R's uniform generator.
#'
#' @param n Number of draws.
#' @param location Location parameter (default 0).
#' @param scale Scale parameter; must be positive.
#' @return Numeric vector of length \code{n}.
#' @examples
#' rlaplace(5, 0, 1)
#' @export
rlaplace <- function(n, location = 0, scale = 1) {
  if (!is.numeric(scale) || length(scale) != 1L || scale <= 0) {
    stop("`scale` must be a positive scalar.", call. = FALSE)
  }
  u <- stats::runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

#' Laplace mechanism
#'
#' Releases \code{f(x)} plus Laplace noise calibrated to the L1 sensitivity,
#' providing pure epsilon-DP.
#'
#' @param f Function computing the target statistic on data \code{x}.
#' @param x The private dataset.
#' @param epsilon Privacy parameter epsilon (> 0).
#' @param l1_sensitivity L1 sensitivity of \code{f}.
#' @return The privatized value of \code{f(x)} (scalar or vector).
#' @examples
#' laplace_mechanism(mean, c(1:100), 1.0, l1_sensitivity = 100 / 99)
#' @export
laplace_mechanism <- function(f, x, epsilon, l1_sensitivity) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(l1_sensitivity) || length(l1_sensitivity) != 1L ||
      l1_sensitivity < 0) {
    stop("`l1_sensitivity` must be a non-negative scalar.", call. = FALSE)
  }
  true_value <- f(x)
  noise <- rlaplace(length(true_value), 0, l1_sensitivity / epsilon)
  true_value + noise
}

# --------------------------------------------------------------- Gaussian ----

#' Gaussian mechanism
#'
#' Releases \code{f(x)} plus Gaussian noise using the classic
#' \eqn{\sigma = \Delta_2 \sqrt{2 \log(1.25/\delta)} / \epsilon}
#' calibration for (epsilon, delta)-DP.
#'
#' @param f Function computing the target statistic on data \code{x}.
#' @param x The private dataset.
#' @param epsilon Privacy parameter epsilon (> 0).
#' @param delta Privacy parameter delta in (0, 1).
#' @param l2_sensitivity L2 sensitivity of \code{f}.
#' @return The privatized value of \code{f(x)} (scalar or vector).
#' @examples
#' gaussian_mechanism(sum, c(1:100), 1.0, 1e-6, l2_sensitivity = 100)
#' @export
gaussian_mechanism <- function(f, x, epsilon, delta, l2_sensitivity) {
  sigma <- gaussian_sigma(epsilon, delta, l2_sensitivity)
  true_value <- f(x)
  true_value + stats::rnorm(length(true_value), 0, sigma)
}

#' Classic Gaussian noise scale
#'
#' Computes the standard deviation used by the classic Gaussian mechanism.
#'
#' @inheritParams gaussian_mechanism
#' @return Scalar noise standard deviation.
#' @export
gaussian_sigma <- function(epsilon, delta, l2_sensitivity) {
  check_epsilon_delta(epsilon, delta)
  if (!is.numeric(l2_sensitivity) || length(l2_sensitivity) != 1L ||
      l2_sensitivity < 0) {
    stop("`l2_sensitivity` must be a non-negative scalar.", call. = FALSE)
  }
  l2_sensitivity * sqrt(2 * log(1.25 / delta)) / epsilon
}