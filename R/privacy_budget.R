# ------------------------------------------------------- Budget class -------

#' Create a privacy budget object
#'
#' An S3 object tracking cumulative privacy loss across multiple DP releases
#' under basic, advanced, or Renyi DP (RDP) composition.
#'
#' @param epsilon Total epsilon available (> 0).
#' @param delta Total delta available (default 1e-6; must be in (0,1)).
#' @param composition One of \code{"basic"}, \code{"advanced"}, \code{"rdp"}.
#' @return An object of class \code{"privacy_budget"}.
#' @examples
#' b <- new_privacy_budget(3.0, 1e-6)
#' b <- spend(b, 1.0, description = "DP mean")
#' b$epsilon_spent
#' @export
new_privacy_budget <- function(epsilon, delta = 1e-6,
                               composition = c("basic", "advanced", "rdp")) {
  composition <- match.arg(composition)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("`epsilon` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || delta <= 0 ||
      delta >= 1) {
    stop("`delta` must be a numeric scalar in (0, 1).", call. = FALSE)
  }
  structure(
    list(
      epsilon_total = epsilon,
      delta_total = delta,
      epsilon_spent = 0,
      delta_spent = 0,
      rho_spent = 0, # zCDP/RDP accumulator for "rdp"
      composition = composition,
      history = data.frame(
        epsilon = numeric(0), delta = numeric(0),
        description = character(0), time = as.POSIXct(character(0))
      )
    ),
    class = "privacy_budget"
  )
}

#' Spend privacy budget
#'
#' Records one DP release against the budget and updates the cumulative
#' privacy loss under the configured composition rule.
#'
#' @param budget A \code{"privacy_budget"} object.
#' @param epsilon Epsilon spent by this release (> 0).
#' @param delta Delta spent by this release (default 0).
#' @param description Optional text label.
#' @return The updated \code{"privacy_budget"} object (invisibly also sets
#'   \code{remaining_epsilon}).
#' @seealso \code{\link{new_privacy_budget}}
#' @export
spend <- function(budget, epsilon, delta = 0, description = "") {
  if (!inherits(budget, "privacy_budget")) {
    stop("`budget` must be a privacy_budget object.", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon < 0) {
    stop("`epsilon` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || delta < 0) {
    stop("`delta` must be a non-negative scalar.", call. = FALSE)
  }
  k <- nrow(budget$history) + 1

  if (budget$composition == "basic") {
    budget$epsilon_spent <- budget$epsilon_spent + epsilon
    budget$delta_spent <- budget$delta_spent + delta
  } else if (budget$composition == "advanced") {
    budget$epsilon_spent <- advanced_composition_epsilon(
      budget$history$epsilon, epsilon, budget$delta_total)
    budget$delta_spent <- budget$delta_spent + delta
  } else { # rdp / zCDP
    budget$rho_spent <- budget$rho_spent + epsilon^2 / 2
    budget$epsilon_spent <- rdp_to_epsilon(budget$rho_spent,
                                           budget$delta_total)
    budget$delta_spent <- budget$delta_spent + delta
  }

  budget$history <- rbind(
    budget$history,
    data.frame(epsilon = epsilon, delta = delta,
               description = as.character(description),
               time = Sys.time())
  )
  budget$remaining_epsilon <- budget$epsilon_total - budget$epsilon_spent
  if (budget$remaining_epsilon < 0) {
    warning("Privacy budget exceeded!", call. = FALSE)
  }
  invisible(budget)
}