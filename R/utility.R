#' Utility diagnostics for DP estimators
#'
#' Computes bias, mean squared error, RMSE, and (optionally) CI coverage of a
#' DP estimator relative to its non-private counterpart or a known truth,
#' across Monte Carlo replicates.
#'
#' @param estimates Numeric vector of estimates from repeated runs.
#' @param truth Scalar true value (or the non-private estimate).
#' @param ci_matrix Optional B x 2 matrix of interval endpoints
#'   (\code{lower}, \code{upper}) aligned with \code{estimates}, for coverage.
#' @return List with \code{bias}, \code{mse}, \code{rmse} and optional
#'   \code{coverage}.
#' @examples
#' ests <- replicate(100, dp_mean(rnorm(200), 1.0, c(-5, 5))$estimate)
#' dp_utility_diagnostics(ests, 0)
#' @export
dp_utility_diagnostics <- function(estimates, truth,
                                   ci_matrix = NULL) {
  if (!is.numeric(estimates) || length(estimates) == 0) {
    stop("`estimates` must be a non-empty numeric vector.", call. = FALSE)
  }
  bias <- mean(estimates - truth)
  mse <- mean((estimates - truth)^2)
  out <- list(bias = bias, mse = mse, rmse = sqrt(mse), n_reps =
                length(estimates))
  if (!is.null(ci_matrix)) {
    if (!is.matrix(ci_matrix) || ncol(ci_matrix) != 2 ||
        nrow(ci_matrix) != length(estimates)) {
      stop("`ci_matrix` must be a length(estimates) x 2 matrix.",
           call. = FALSE)
    }
    cov <- mean(ci_matrix[, 1] <= truth & truth <= ci_matrix[, 2])
    width <- mean(ci_matrix[, 2] - ci_matrix[, 1])
    out$coverage <- cov
    out$mean_ci_width <- width
  }
  out
}

#' Compare utility of non-private and private fits
#'
#' Fits \code{\link{dp_lm}} at each epsilon on an epsilon grid and reports MSE
#' of coefficients relative to the ordinary least squares fit.
#'
#' @param formula Model formula.
#' @param data Data frame.
#' @param epsilon_grid Numeric vector of epsilon values to compare.
#' @param delta Delta for Gaussian mechanism.
#' @param bounds Bounds list passed to \code{\link{dp_lm}}.
#' @param n_reps Number of repetitions per epsilon (default 10).
#' @return Data frame with columns \code{epsilon}, \code{rep},
#'   \code{mse_coef}, \code{max_abs_error}.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x = rnorm(200)); d$y <- 1 + d$x + rnorm(200)
#' compare_utility(y ~ x, d, c(0.5, 2), 1e-6, list(y = c(-15, 15)),
#'                 n_reps = 3)
#' }
#' @export
compare_utility <- function(formula, data, epsilon_grid, delta = 1e-6,
                            bounds = NULL, n_reps = 10) {
  non_dp <- stats::lm(formula, data)
  beta0 <- coef(non_dp)
  rows <- lapply(epsilon_grid, function(eps) {
    do.call(rbind, lapply(seq_len(n_reps), function(r) {
      fit <- dp_lm(formula, data, epsilon = eps, delta = delta,
                   bounds = bounds)
      err <- fit$coefficients - beta0[names(fit$coefficients)]
      data.frame(
        epsilon = eps,
        rep = r,
        mse_coef = mean(err^2),
        max_abs_error = max(abs(err))
      )
    }))
  })
  do.call(rbind, rows)
}