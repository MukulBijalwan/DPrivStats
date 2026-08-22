#' Asymptotic Kolmogorov distribution survival function
#'
#' Computes \eqn{P(D_n \sqrt{n} > t)} under the Kolmogorov-Smirnov null via
#' the series \eqn{2 \sum_{k\ge1} (-1)^{k-1} e^{-2 k^2 t^2}}.
#'
#' @param q Numeric vector of scaled KS statistics (>= 0).
#' @return Tail probabilities in [0, 1].
#' @keywords internal
asymptotic_ks_pvalue <- function(q) {
  q <- pmax(q, 0)
  vapply(q, function(tq) {
    k <- seq_len(max(50L, ceiling(tq)))
    s <- sum((-1)^{seq_along(k) - 1} * exp(-2 * k^2 * tq^2))
    max(pmin(2 * s, 1), 0)
  }, numeric(1))
}