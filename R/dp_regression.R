#' Differentially private linear regression
#'
#' Fits OLS on clipped data and releases the coefficient vector under the
#' Gaussian mechanism, calibrated to the L2 sensitivity of
#' \eqn{\hat\beta = (X^\top X)^{-1} X^\top y} given response bounds
#' \code{bounds$y} (and optionally a bound \eqn{C} on row norms of X; see
#' \code{\link{dp_lm_sensitivity}}).
#'
#' @param formula Model formula.
#' @param data Data frame.
#' @param epsilon Privacy parameter epsilon (> 0).
#' @param delta Privacy parameter delta in (0, 1); required.
#' @param bounds Named list: \code{bounds$y} = c(L_y, U_y) for the response
#'   (required); optional \code{bounds$x_norm} giving a priori bound C on row
#'   norms of the design matrix when the design itself is private.
#' @param assume_public_design Logical; see \code{\link{dp_lm_sensitivity}}.
#' @return An object of class \code{c("dp_lm", "lm")} with privatized
#'   coefficients, residuals and fitted values computed from the private fit,
#'   plus elements \code{epsilon}, \code{delta}, \code{noise_variance}
#'   (privacy noise variance per coefficient), and \code{sampling_vcov}
#'   (the non-private sampling covariance estimate used by
#'   \code{\link{dp_confint}}).
#' @examples
#' set.seed(11)
#' d <- data.frame(x = rnorm(300), z = rnorm(300))
#' d$y <- 1 + 2 * d$x - d$z + rnorm(300)
#' fit <- dp_lm(y ~ x + z, d, epsilon = 2.0, delta = 1e-6,
#'              bounds = list(y = c(-30, 30)))
#' coef(fit)
#' @export
dp_lm <- function(formula, data, epsilon, delta,
                  bounds = NULL, assume_public_design = TRUE) {
  if (is.null(bounds) || is.null(bounds$y)) {
    stop("`bounds` must be a list containing `y` = c(L_y, U_y).",
         call. = FALSE)
  }
  check_delta(delta)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  if (!is.numeric(y)) stop("Response must be numeric.", call. = FALSE)
  X <- stats::model.matrix(formula, mf)

  y_clipped <- clip_numeric(y, bounds$y)
  # Clip numeric predictors column-wise to their observed range so that no
  # single record can move the released fit arbitrarily.
  X_clipped <- X
  num_cols <- which(vapply(seq_len(ncol(X)), function(j)
    is.numeric(X[, j]) && !all(X[, j] %in% c(0, 1)), logical(1)))
  for (j in num_cols) {
    b <- range(X[, j])
    if (diff(b) > 0) X_clipped[, j] <- clip_numeric(X[, j], b)
  }

  XtX <- crossprod(X_clipped)
  XtX_inv <- MASS::ginv(XtX)
  beta_hat <- as.vector(XtX_inv %*% crossprod(X_clipped, y_clipped))

  l2_sens <- dp_lm_sensitivity(
    X_clipped, bounds$y,
    x_norm_bound = bounds$x_norm,
    assume_public_design = assume_public_design
  )
  sigma <- analytic_gaussian_sigma(epsilon, delta, l2_sens)
  noise <- stats::rnorm(length(beta_hat), 0, sigma)
  beta_priv <- beta_hat + noise

  names(beta_priv) <- colnames(X)
  fitted <- as.vector(X_clipped %*% beta_priv)
  resid <- y_clipped - fitted

  n <- nrow(X_clipped); p <- ncol(X_clipped)
  sigma2_hat <- sum(resid^2) / max(n - p, 1)
  sampling_vcov <- sigma2_hat * XtX_inv

  structure(
    list(
      coefficients = beta_priv,
      residuals = resid,
      fitted.values = fitted,
      df.residual = max(n - p, 1),
      rank = p,
      call = match.call(),
      model = mf,
      epsilon = epsilon,
      delta = delta,
      l2_sensitivity = l2_sens,
      noise_variance = rep(sigma^2, p),
      sampling_vcov = sampling_vcov,
      bounds = bounds
    ),
    class = c("dp_lm", "lm")
  )
}

#' Print method for dp_lm objects
#'
#' @param x A \code{dp_lm} object.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x}, invisibly.
#' @export
print.dp_lm <- function(x, ...) {
  cat("\nDifferentially Private Linear Regression\n")
  cat(sprintf("epsilon = %.3f, delta = %.1e\n\n", x$epsilon, x$delta))
  print(round(coef(x), 4))
  invisible(x)
}