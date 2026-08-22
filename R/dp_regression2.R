#' DP-SGD for generalized linear models
#'
#' Fits a GLM with differentially private stochastic gradient descent:
#' per-sample gradients are clipped to \code{max_grad_norm}, averaged, and
#' Gaussian noise is added at each iteration. Privacy is accounted with the
#' advanced-composition heuristic over \code{n_iter} iterations.
#'
#' @inheritParams dp_lm
#' @param family A GLM family (e.g. \code{\link[stats]{binomial}()}).
#' @param max_grad_norm Per-sample gradient clipping norm C (> 0).
#' @param n_iter Number of SGD iterations.
#' @param lr Learning rate.
#' @param batch_size Mini-batch size (default full batch).
#' @return A list of class \code{"dp_glm"} with \code{coefficients},
#'   \code{epsilon}, \code{delta}, and convergence details.
#' @examples
#' set.seed(12)
#' d <- data.frame(x = rnorm(400))
#' d$y <- rbinom(400, 1, plogis(0.5 * d$x))
#' fit <- dp_glm(y ~ x, d, binomial(), epsilon = 2.0, delta = 1e-6,
#'               bounds = list(y = c(0, 1)))
#' fit$coefficients
#' @export
dp_glm <- function(formula, data, family = stats::gaussian(),
                   epsilon, delta, max_grad_norm = 1,
                   n_iter = 1000, lr = 0.01, batch_size = NULL,
                   bounds = NULL) {
  check_delta(delta)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, mf)
  n <- nrow(X); p <- ncol(X)
  if (is.null(batch_size)) batch_size <- n

  if (!is.null(bounds$y)) y <- clip_numeric(y, bounds$y)
  y <- as.numeric(y)

  eps_iter <- epsilon / sqrt(n_iter * log(1 / delta))

  beta <- rep(0, p)
  trace <- numeric(n_iter)
  for (t in seq_len(n_iter)) {
    idx <- if (batch_size >= n) seq_len(n) else sample.int(n, batch_size)
    eta <- as.vector(X[idx, , drop = FALSE] %*% beta)
    mu <- switch(
      family$family,
      binomial = stats::plogis(eta),
      poisson = exp(pmin(eta, 700)),
      gaussian = eta
    )
    w <- switch(
      family$family,
      binomial = mu * (1 - mu),
      poisson = pmax(mu, .Machine$double.eps),
      gaussian = rep(1, length(mu))
    )
    r <- y[idx] - mu
    # Negative log-likelihood gradient per sample: -x_i * w_i * r_i
    # (raw per-sample gradients are clipped, then averaged)
    grads <- -sweep(X[idx, , drop = FALSE], 1, w, `*`) * r
    cl <- clip_gradients(grads, max_grad_norm)
    avg_grad <- colMeans(cl$clipped)
    sigma <- max_grad_norm * sqrt(2 * log(1.25 / delta)) /
      (batch_size * eps_iter)
    noisy_grad <- avg_grad + stats::rnorm(p, 0, sigma)
    beta <- beta - lr * noisy_grad
    trace[t] <- mean(r^2)
  }

  structure(
    list(
      coefficients = stats::setNames(as.vector(beta), colnames(X)),
      epsilon = epsilon,
      delta = delta,
      max_grad_norm = max_grad_norm,
      n_iter = n_iter,
      lr = lr,
      family = family$family,
      loss_trace = trace
    ),
    class = "dp_glm"
  )
}

#' Print method for dp_glm objects
#'
#' @param x A \code{dp_glm} object.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x}, invisibly.
#' @export
print.dp_glm <- function(x, ...) {
  cat("\nDifferentially Private GLM (DP-SGD)\n")
  cat(sprintf("family = %s, iterations = %d\n", x$family, x$n_iter))
  cat(sprintf("epsilon = %.3f, delta = %.1e\n\n", x$epsilon, x$delta))
  print(round(x$coefficients, 4))
  invisible(x)
}