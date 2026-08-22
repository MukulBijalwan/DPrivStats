#' Analytic Gaussian calibration (Balle & Wang, 2018)
#'
#' Solves for the smallest sigma satisfying the analytic Gaussian mechanism
#' privacy guarantee:
#' \deqn{\Phi(\Delta/(2\sigma) - \epsilon\sigma/\Delta)
#'   - e^{\epsilon}\Phi(-\Delta/(2\sigma) - \epsilon\sigma/\Delta) \le \delta}
#' via bisection. Always at least as tight as the classic calibration.
#'
#' @inheritParams gaussian_mechanism
#' @return Scalar noise standard deviation.
#' @references
#' Balle, B., & Wang, Y. (2018). Improving the Gaussian Mechanism for
#' Differential Privacy via Analytic Concentration of the Gaussian
#' Distribution. \emph{NeurIPS 2018}.
#' @examples
#' analytic_gaussian_sigma(1.0, 1e-6, 1.0)
#' gaussian_sigma(1.0, 1e-6, 1.0) # looser
#' @export
analytic_gaussian_sigma <- function(epsilon, delta, l2_sensitivity) {
  check_epsilon_delta(epsilon, delta)
  if (!is.numeric(l2_sensitivity) || length(l2_sensitivity) != 1L ||
      l2_sensitivity < 0) {
    stop("`l2_sensitivity` must be a non-negative scalar.", call. = FALSE)
  }
  d <- l2_sensitivity
  if (d == 0) {
    return(0)
  }
  target <- function(sigma) {
    a <- d / (2 * sigma)
    b <- epsilon * sigma / d
    stats::pnorm(a - b) - exp(epsilon) * stats::pnorm(-a - b)
  }
  # Bracket: the classic sigma is always sufficient (Balle & Wang, 2018).
  lo <- 1e-12
  hi <- max(gaussian_sigma(epsilon, delta, d), .Machine$double.xmin) * 2
  while (target(hi) > delta && hi < 1e300) hi <- hi * 2
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (target(mid) <= delta) hi <- mid else lo <- mid
  }
  (lo + hi) / 2
}

#' Analytic Gaussian mechanism
#'
#' Releases \code{f(x)} plus Gaussian noise calibrated with the analytic
#' Gaussian mechanism (Balle & Wang, 2018).
#'
#' @inheritParams gaussian_mechanism
#' @return The privatized value of \code{f(x)} (scalar or vector).
#' @export
analytic_gaussian_mechanism <- function(f, x, epsilon, delta, l2_sensitivity) {
  sigma <- analytic_gaussian_sigma(epsilon, delta, l2_sensitivity)
  true_value <- f(x)
  true_value + stats::rnorm(length(true_value), 0, sigma)
}

# ---------------------------------------------------- Exponential mechanism --

#' Exponential mechanism sampler
#'
#' Samples a candidate from \code{candidates} with probability proportional to
#' \eqn{\exp(\epsilon \, s(c) / (2 \Delta s))}, where \eqn{\Delta s} is the
#' sensitivity of the score function.
#'
#' @param candidates Numeric vector of candidate outputs.
#' @param scores Numeric vector of scores (higher is better).
#' @param epsilon Privacy parameter epsilon (> 0).
#' @param score_sensitivity L1 sensitivity of the score function.
#' @return One sampled candidate.
#' @examples
#' exponential_mechanism(1:10, -(1:10 - 5)^2, 1.0, 1)
#' @export
exponential_mechanism <- function(candidates, scores, epsilon,
                                  score_sensitivity = 1) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  if (length(candidates) != length(scores)) {
    stop("`candidates` and `scores` must have equal length.", call. = FALSE)
  }
  log_w <- epsilon * scores / (2 * score_sensitivity)
  log_w <- log_w - max(log_w)
  w <- exp(log_w)
  sample(candidates, size = 1, prob = w)
}

# ------------------------------------------------------------ Validation -----

check_epsilon_delta <- function(epsilon, delta) {
  if (!is.numeric(delta) || length(delta) != 1L ||
      delta <= 0 || delta >= 1) {
    stop("`delta` must be a numeric scalar in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_bounds <- function(bounds, n = 2L) {
  if (!is.numeric(bounds) || length(bounds) != n || anyNA(bounds) ||
      bounds[1] >= bounds[n]) {
    stop("`bounds` must be a numeric vector of length ", n,
         " with lower < upper.", call. = FALSE)
  }
  invisible(TRUE)
}

clip_numeric <- function(x, bounds) {
  pmax(bounds[1], pmin(as.numeric(x), bounds[2]))
}