# ------------------------------------------------------- Sensitivity rules ---

#' Sensitivity of common bounded-data statistics
#'
#' Computes the L1/L2 sensitivities under bounded-DP (neighbouring datasets
#' differ in one record, same size n) for a record \eqn{x_i \in [L, U]}.
#'
#' @param bounds Numeric vector c(L, U) giving the data bounds.
#' @param n Sample size (positive integer).
#' @return A list with elements \code{l1} and \code{l2}, each a named numeric
#'   vector of sensitivities for \code{"mean"}, \code{"sum"}, \code{"count"},
#'   \code{"variance"} and \code{"ecdf"}.
#' @examples
#' sens <- stat_sensitivities(c(0, 500000), 1000)
#' sens$l1["mean"]
#' @export
stat_sensitivities <- function(bounds, n) {
  validate_bounds(bounds)
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != floor(n)) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  range <- bounds[2] - bounds[1]
  list(
    l1 = c(
      mean     = range / n,
      sum      = range,
      count    = 1,
      variance = range^2 / n,
      ecdf     = 1 / n
    ),
    l2 = c(
      mean     = range / n,
      sum      = range,
      count    = 1,
      variance = range^2 / n,
      ecdf     = 1 / n
    )
  )
}

#' Global L2 sensitivity of OLS coefficients
#'
#' Under the assumption that the design matrix \eqn{X} is public/fixed and the
#' response is bounded (\eqn{y_i \in [L_y, U_y]}), neighbouring datasets differ
#' in one response value only. With \eqn{C = \max_i \|x_i\|_2} and
#' \eqn{D = U_y - L_y},
#' \deqn{\Delta_2(\hat\beta) = C \cdot D / \lambda_{\min}(X^\top X).}
#'
#' If the design matrix itself is bounded per-row by \code{x_norm_bound} and
#' treated as private, a conservative multiplicative factor of 3 is applied to
#' account for simultaneous changes in \eqn{X'X} and \eqn{X'y} (documented
#' heuristic; see vignette "regression-guide" for discussion).
#'
#' @param X Numeric model matrix (n x p).
#' @param y_bounds Bounds c(L_y, U_y) on the response.
#' @param x_norm_bound Optional bound \eqn{C = \max_i \|x_i\|_2}. If NULL it is
#'   computed from X (use this when X is public). When the design is private,
#'   supply an a priori bound instead.
#' @param assume_public_design Logical; if TRUE (default) the sensitivity
#'   accounts only for changes in y.
#' @return Scalar L2 sensitivity of the OLS coefficient vector.
#' @examples
#' X <- cbind(1, 1:10)
#' dp_lm_sensitivity(X, y_bounds = c(0, 10))
#' @export
dp_lm_sensitivity <- function(X, y_bounds, x_norm_bound = NULL,
                              assume_public_design = TRUE) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }
  validate_bounds(y_bounds)
  XtX <- crossprod(X)
  ev <- eigen(XtX, symmetric = TRUE, only.values = TRUE)$values
  ev_pos <- ev[ev >= 0]
  tol <- sqrt(.Machine$double.eps) * max(abs(ev))
  if (length(ev_pos) < ncol(X) || min(ev_pos) <= tol) {
    stop("Design matrix is rank-deficient; cannot compute sensitivity.",
         call. = FALSE)
  }
  lambda_min <- min(ev_pos)
  C <- if (!is.null(x_norm_bound)) x_norm_bound else max(sqrt(rowSums(X^2)))
  D <- diff(y_bounds)
  delta <- C * D / lambda_min
  if (!assume_public_design) delta <- 3 * delta
  as.numeric(delta)
}

#' Per-sample gradient clipping
#'
#' Clips each row of a gradient matrix to a maximum L2 norm:
#' \eqn{\tilde g_i = g_i \min(1, C/\|g_i\|_2)}. Used by DP-SGD.
#'
#' @param grads Numeric matrix (n x p) of per-sample gradients.
#' @param max_grad_norm Positive clipping constant C.
#' @return List with \code{clipped} (the clipped gradients) and
#'   \code{norms} (original norms).
#' @examples
#' g <- matrix(rnorm(30), 10, 3)
#' clip_gradients(g, 1)$norms
#' @export
clip_gradients <- function(grads, max_grad_norm) {
  if (!is.matrix(grads) || !is.numeric(grads)) {
    stop("`grads` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(max_grad_norm) || length(max_grad_norm) != 1L ||
      max_grad_norm <= 0) {
    stop("`max_grad_norm` must be a positive scalar.", call. = FALSE)
  }
  norms <- sqrt(rowSums(grads^2))
  factors <- pmin(1, max_grad_norm / pmax(norms, .Machine$double.eps))
  list(clipped = grads * factors, norms = norms)
}

#' DP median via the exponential mechanism
#'
#' Releases a differentially private estimate of the median using the
#' exponential mechanism over a grid of candidate values, with score
#' \eqn{s(t; x) = -|\#\{i : x_i \le t\}| - n/2|} which has L1
#' sensitivity 1.
#'
#' @param x Numeric vector (may contain NAs, ignored).
#' @param epsilon Privacy parameter epsilon (> 0).
#' @param bounds Bounds c(L, U) assumed known for the data.
#' @param n_bins Number of candidate grid points (default 1000).
#' @return A list of class \code{dp_estimate} with element \code{estimate}.
#' @examples
#' set.seed(42)
#' fit <- dp_median(rnorm(200), 1.0, bounds = c(-5, 5))
#' fit$estimate
#' @export
dp_median <- function(x, epsilon, bounds, n_bins = 1000) {
  validate_bounds(bounds)
  x <- as.numeric(stats::na.omit(x))
  n <- length(x)
  if (n == 0) stop("`x` contains no usable observations.", call. = FALSE)
  t_grid <- seq(bounds[1], bounds[2], length.out = max(n_bins, 2L))
  ranks <- vapply(t_grid, function(t) sum(x <= t), numeric(1))
  scores <- -abs(ranks - n / 2)
  t_hat <- exponential_mechanism(t_grid, scores, epsilon, score_sensitivity = 1)
  structure(
    list(
      estimate = unname(t_hat),
      private = TRUE,
      epsilon = epsilon,
      delta = 0,
      mechanism = "exponential",
      statistic = "median",
      n = n,
      bounds = bounds
    ),
    class = "dp_estimate"
  )
}