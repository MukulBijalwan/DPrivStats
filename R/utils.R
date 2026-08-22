#' Fast Laplace sampling (C++)
#'
#' Vectorized inverse-CDF Laplace sampler implemented in C++ via Rcpp.
#' Equivalent to \code{\link{rlaplace}} but faster for large draws.
#'
#' @param n Number of draws.
#' @param scale Scale parameter (positive).
#' @return Numeric vector of length \code{n}.
#' @examples
#' cpp_rlaplace(5, 1.0)
#' @useDynLib DPrivStats, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
cpp_rlaplace <- function(n, scale = 1) {
  .Call(`_DPrivStats_rlaplace_cpp`, as.integer(n), as.numeric(scale))
}

#' Fast per-sample gradient clipping (C++)
#'
#' Clips rows of a gradient matrix to a maximum L2 norm; equivalent to
#' \code{\link{clip_gradients}} but implemented in C++.
#'
#' @param grads Numeric matrix of per-sample gradients.
#' @param max_grad_norm Positive clipping constant.
#' @return Numeric matrix of clipped gradients.
#' @examples
#' g <- matrix(rnorm(30), 10, 3)
#' cpp_clip_gradients(g, 1)
#' @export
cpp_clip_gradients <- function(grads, max_grad_norm) {
  if (!is.matrix(grads)) grads <- matrix(grads, nrow = 1)
  .Call(`_DPrivStats_clip_gradients_cpp`,
        matrix(as.numeric(grads), nrow = nrow(grads)),
        as.numeric(max_grad_norm))
}
