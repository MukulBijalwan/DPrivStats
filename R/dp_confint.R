#' Privacy-aware confidence intervals
#'
#' Confidence intervals for \code{\link{dp_lm}} coefficients that account for
#' both sampling variability and privacy noise.
#'
#' Three methods:
#' \describe{
#'   \item{\code{"analytical"}}{Normal intervals using
#'     \eqn{\sqrt{\text{diag}(\Sigma_{\text{sampling}}) +
#'     v_{\text{privacy}}}}, where privacy noise variance is known exactly
#'     from the mechanism calibration.}
#'   \item{\code{"parametric_bootstrap"}}{Simulates
#'     \eqn{\beta^* \sim N(\hat\beta_{\text{priv}}, \Sigma_{\text{sampling}})}
#'     and adds fresh privacy noise; returns percentile intervals.}
#'   \item{\code{"privacy_aware_bootstrap"}}{Resamples the data and refits the
#'     DP estimator on each resample with per-resample budget
#'     \eqn{\epsilon / \sqrt{B}} (advanced composition heuristic); returns
#'     percentile intervals. Most conservative but computationally heavy.}
#' }
#'
#' @param object A \code{dp_lm} object.
#' @param parm Coefficient indices or names (default all).
#' @param level Confidence level (default 0.95).
#' @param method One of \code{"analytical"}, \code{"parametric_bootstrap"},
#'   \code{"privacy_aware_bootstrap"}.
#' @param B Number of bootstrap replicates for bootstrap methods.
#' @param ... Additional arguments (ignored).
#' @return Numeric matrix of class \code{c("dp_confint", "matrix")} with
#'   columns \code{estimate}, \code{lower}, \code{upper}.
#' @examples
#' set.seed(11)
#' d <- data.frame(x = rnorm(300))
#' d$y <- 1 + 2 * d$x + rnorm(300)
#' fit <- dp_lm(y ~ x, d, epsilon = 2.0, delta = 1e-6,
#'              bounds = list(y = c(-20, 20)))
#' dp_confint(fit)
#' @export
dp_confint <- function(object, parm, level = 0.95,
                       method = c("analytical", "parametric_bootstrap",
                                  "privacy_aware_bootstrap"),
                       B = 200, ...) {
  if (!inherits(object, "dp_lm")) {
    stop("`object` must be a dp_lm object.", call. = FALSE)
  }
  method <- match.arg(method)
  beta <- object$coefficients
  if (missing(parm)) parm <- seq_along(beta)
  idx <- if (is.character(parm)) match(parm, names(beta)) else parm
  if (anyNA(idx)) stop("Unknown coefficient(s) in `parm`.", call. = FALSE)

  alpha <- 1 - level
  est <- beta[idx]
  lo <- hi <- numeric(length(idx))

  if (method == "analytical") {
    se_sampling <- sqrt(diag(object$sampling_vcov))
    se_privacy <- sqrt(object$noise_variance)
    total_se <- sqrt(se_sampling[idx]^2 + se_privacy[idx]^2)
    z <- stats::qnorm(1 - alpha / 2)
    lo <- est - z * total_se
    hi <- est + z * total_se
  } else if (method == "parametric_bootstrap") {
    draws <- MASS::mvrnorm(B, beta, object$sampling_vcov) +
      matrix(stats::rnorm(B * length(beta), 0,
                          rep(sqrt(object$noise_variance),
                              each = B)), B, length(beta), byrow = TRUE)
    lo <- apply(draws[, idx, drop = FALSE], 2, stats::quantile,
                probs = alpha / 2)
    hi <- apply(draws[, idx, drop = FALSE], 2, stats::quantile,
                probs = 1 - alpha / 2)
  } else {
    call_obj <- object$call
    data_env <- eval(call_obj$data, environment(formula(eval(
      call_obj$formula, parent.frame()))))
    eps_b <- object$epsilon / sqrt(B)
    boot <- matrix(NA_real_, B, length(beta))
    for (b in seq_len(B)) {
      rows <- sample(nrow(data_env), replace = TRUE)
      refit <- tryCatch(
        dp_lm(call_obj$formula, data_env[rows, , drop = FALSE],
              epsilon = eps_b, delta = object$delta,
              bounds = object$bounds, assume_public_design = TRUE),
        error = function(e) NULL
      )
      if (!is.null(refit)) boot[b, ] <- refit$coefficients
    }
    ok <- stats::complete.cases(boot)
    if (sum(ok) < max(B %/% 10, 5)) {
      warning("Too many failed bootstrap refits; falling back to analytical.",
              call. = FALSE)
      return(dp_confint(object, parm, level, "analytical"))
    }
    lo <- apply(boot[ok, idx, drop = FALSE], 2, stats::quantile,
                probs = alpha / 2)
    hi <- apply(boot[ok, idx, drop = FALSE], 2, stats::quantile,
                probs = 1 - alpha / 2)
  }

  out <- cbind(estimate = est, lower = lo, upper = hi)
  rownames(out) <- names(beta)[idx]
  class(out) <- c("dp_confint", "matrix")
  out
}

#' Confint interface for dp_lm objects
#'
#' Convenience wrapper calling \code{\link{dp_confint}} with default settings.
#'
#' @inheritParams dp_confint
#' @return Numeric matrix with columns \code{estimate}, \code{lower},
#'   \code{upper}.
#' @export
confint.dp_lm <- function(object, parm, level = 0.95, method = "analytical",
                          B = 200, ...) {
  dp_confint(object, parm = parm, level = level, method = method, B = B)
}