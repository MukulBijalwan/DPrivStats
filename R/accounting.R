#' Privacy loss random variables for the Laplace mechanism
#'
#' The privacy loss of one Laplace release with scale \eqn{s} is
#' \eqn{L = |d|/s} with \eqn{d \sim \mathrm{Laplace}(0, s)}, i.e.
#' \eqn{L \sim \mathrm{Exp}(\epsilon)} when \eqn{s = \Delta_1/\epsilon}.
#' These helpers expose its tail probabilities and quantiles for PLRV-based
#' accounting and simulation studies.
#'
#' @param t Numeric vector of thresholds on the privacy loss.
#' @param epsilon Privacy parameter (> 0).
#' @return Tail probabilities \eqn{P(L > t)} (\code{laplace_plr_tail}) or
#'   quantiles at probabilities \code{p} (\code{laplace_plr_quantile}).
#' @examples
#' laplace_plr_tail(3, 1.0) # exp(-3)
#' @export
laplace_plr_tail <- function(t, epsilon) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  exp(-epsilon * pmax(t, 0))
}

#' @rdname laplace_plr_tail
#' @param p Numeric vector of probabilities in [0, 1].
#' @export
laplace_plr_quantile <- function(p, epsilon) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  -log(1 - p) / epsilon
}

#' Simulate privacy losses of Gaussian releases
#'
#' Draws the realized privacy loss \eqn{L = z^2/(2\sigma^2) + ... } approximated
#' by \eqn{L = |z|/\sigma + 1/(2\sigma^2)} for two neighbouring datasets under
#' the Gaussian mechanism, useful for empirical delta estimation.
#'
#' @param n Number of draws.
#' @param sigma Noise standard deviation.
#' @return Numeric vector of simulated losses.
#' @examples
#' l <- simulate_gaussian_losses(1000, 5); mean(l)
#' @export
simulate_gaussian_losses <- function(n, sigma) {
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("`sigma` must be a positive scalar.", call. = FALSE)
  }
  z <- stats::rnorm(n, 0, sigma)
  abs(z) / sigma + 1 / (2 * sigma^2)
}

#' Empirical delta from simulated privacy losses
#'
#' Estimates the delta needed so that \eqn{P(L > \epsilon) \le \delta} given a
#' sample of simulated losses.
#'
#' @param losses Numeric vector of simulated privacy losses.
#' @param epsilon Privacy parameter threshold.
#' @return Empirical tail probability.
#' @examples
#' empirical_delta_from_losses(simulate_gaussian_losses(1000, 5), 2.0)
#' @export
empirical_delta_from_losses <- function(losses, epsilon) {
  mean(losses > epsilon)
}
