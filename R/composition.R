#' Check whether a release fits in the remaining budget
#'
#' @param budget A \code{"privacy_budget"} object.
#' @param epsilon Proposed epsilon spend.
#' @param delta Proposed delta spend (default 0).
#' @return Logical; TRUE if the proposed release fits.
#' @export
can_spend <- function(budget, epsilon, delta = 0) {
  if (!inherits(budget, "privacy_budget")) {
    stop("`budget` must be a privacy_budget object.", call. = FALSE)
  }
  probe <- suppressWarnings(spend(new_privacy_budget(budget$epsilon_total,
                                                     budget$delta_total,
                                                     budget$composition),
                                  epsilon, delta))
  probe$epsilon_spent <= budget$epsilon_total &&
    probe$delta_spent <= budget$delta_total
}

#' Print method for privacy_budget objects
#'
#' @param x A \code{"privacy_budget"} object.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x}, invisibly.
#' @export
print.privacy_budget <- function(x, ...) {
  cat("\nPrivacy Budget (composition:", x$composition, ")\n")
  cat(sprintf("  total epsilon:  %.3f\n", x$epsilon_total))
  cat(sprintf("  spent epsilon:  %.3f\n", x$epsilon_spent))
  if (x$composition == "rdp") {
    cat(sprintf("  accumulated rho: %.4f\n", x$rho_spent))
  }
  cat(sprintf("  remaining:      %.3f\n",
              x$epsilon_total - x$epsilon_spent))
  invisible(x)
}

# ------------------------------------------------- Composition rules -------

#' Basic (sequential) composition
#'
#' @param epsilons Numeric vector of per-release epsilons.
#' @param deltas Numeric vector of per-release deltas.
#' @return List with total \code{epsilon} and \code{delta}.
#' @export
basic_composition <- function(epsilons, deltas = rep(0, length(epsilons))) {
  list(epsilon = sum(epsilons), delta = sum(deltas))
}

#' Advanced composition
#'
#' Computes the advanced-composition bound for \eqn{k} homogeneous releases:
#' \deqn{\varepsilon' = \sqrt{2 k \log(1/\delta')} \varepsilon
#'   + k \varepsilon (e^{\varepsilon} - 1).}
#'
#' @param prior_epsilons Numeric vector of epsilons already spent.
#' @param new_epsilon Epsilon of the new release.
#' @param delta_prime Target overall failure parameter \eqn{\delta'} in (0,1).
#' @return Total epsilon under advanced composition.
#' @export
advanced_composition_epsilon <- function(prior_epsilons, new_epsilon,
                                         delta_prime = 1e-6) {
  eps_all <- c(prior_epsilons, new_epsilon)
  k <- length(eps_all)
  if (k == 0) return(0)
  eps <- max(eps_all) # per-release epsilon (assumed homogeneous)
  sqrt(2 * k * log(1 / delta_prime)) * eps +
    k * eps * (exp(eps) - 1)
}

#' Convert epsilon-DP to zCDP (RDP) parameter
#'
#' An (epsilon, 0)-DP mechanism satisfies (epsilon^2 / 2)-zCDP
#' (Bun & Steinke, 2016).
#'
#' @param epsilon Pure-DP epsilon.
#' @return zCDP rho parameter.
#' @export
epsilon_to_rdp <- function(epsilon) {
  epsilon^2 / 2
}

#' Convert zCDP (RDP) parameter to (epsilon, delta)-DP
#'
#' Uses \eqn{\varepsilon = \rho + \sqrt{2 \rho \log(1/\delta)}}.
#'
#' @param rho zCDP parameter (>= 0).
#' @param delta Target delta in (0, 1).
#' @return Equivalent epsilon.
#' @export
rdp_to_epsilon <- function(rho, delta) {
  if (rho <= 0) return(0)
  rho + sqrt(2 * rho * log(1 / delta))
}

#' RDP composition
#'
#' Composes a sequence of pure-DP releases via the zCDP/RDP accumulator and
#' converts the result back to (epsilon, delta)-DP. Typically much tighter
#' than basic composition for many small releases.
#'
#' @param epsilons Numeric vector of per-release epsilons.
#' @param delta Target delta in (0, 1).
#' @return List with elements \code{rho} and \code{epsilon}.
#' @examples
#' rdp_composition(rep(0.1, 20), 1e-6)$epsilon
#' basic_composition(rep(0.1, 20))$epsilon
#' @export
rdp_composition <- function(epsilons, delta = 1e-6) {
  rho <- sum(vapply(epsilons, epsilon_to_rdp, numeric(1)))
  list(rho = rho, epsilon = rdp_to_epsilon(rho, delta))
}