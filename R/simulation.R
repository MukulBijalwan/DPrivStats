#' Simulate synthetic microdata
#'
#' Generates a census-like dataset with correlated predictors used in examples
#' and vignettes: income as a function of education, age, and hours worked.
#'
#' @param n Sample size (default 1000).
#' @return A data frame with columns \code{education} (years, 0-20),
#'   \code{age} (18-80), \code{hours} (weekly hours, 0-80),
#'   \code{income} (annual income, >= 0), and \code{region} (factor).
#' @examples
#' head(simulate_data(5))
#' @export
simulate_data <- function(n = 1000) {
  if (!is.numeric(n) || length(n) != 1L || n < 4) {
    stop("`n` must be a positive integer >= 4.", call. = FALSE)
  }
  education <- pmin(pmax(round(stats::rnorm(n, 13, 3)), 0), 20)
  age <- pmin(pmax(round(stats::rnorm(n, 42, 12)), 18), 80)
  region <- factor(sample(c("North", "South", "East", "West"), n,
                          replace = TRUE))
  hours <- pmin(pmax(stats::rnorm(n, 38, 8), 0), 80)
  income <- pmax(
    -30000 + 2500 * education + 350 * age + 400 * hours +
      5000 * (region == "West") +
      stats::rnorm(n, 0, 12000),
    0
  )
  data.frame(education = education, age = age, hours = hours,
             income = income, region = region)
}

#' Validate confidence-interval coverage by Monte Carlo
#'
#' Repeatedly simulates data, fits \code{\link{dp_lm}}, builds privacy-aware
#' CIs with \code{\link{dp_confint}}, and reports empirical coverage and mean
#' interval width against the true generating coefficients.
#'
#' @param formula Model formula matching \code{true_betas} names.
#' @param true_betas Named numeric vector of true coefficients.
#' @param sigma Error standard deviation of the simulated model.
#' @param data_gen Function returning a data frame with columns needed by
#'   \code{formula}; defaults to a Gaussian design built from \code{n}.
#' @param n Sample size per replicate (used by the default generator).
#' @param n_sims Number of Monte Carlo replicates.
#' @param epsilon Privacy parameter epsilon.
#' @param delta Privacy parameter delta.
#' @param y_bounds Response bounds passed to \code{\link{dp_lm}}.
#' @param confint_method Method passed to \code{\link{dp_confint}}.
#' @param level Confidence level.
#' @return List with \code{coverage_rate}, \code{mean_ci_width},
#'   \code{target_coverage}, \code{n_sims}.
#' @examples
#' set.seed(1)
#' validate_coverage(y ~ x, c(`(Intercept)` = 1, x = 2), sigma = 1,
#'                   n = 200, n_sims = 30, epsilon = 3, delta = 1e-6,
#'                   y_bounds = c(-15, 15))
#' @export
validate_coverage <- function(formula, true_betas, sigma = 1,
                              data_gen = NULL, n = 500, n_sims = 200,
                              epsilon = 1.0, delta = 1e-6,
                              y_bounds = c(-50, 50),
                              confint_method = "analytical",
                              level = 0.95) {
  if (is.null(data_gen)) {
    vars <- setdiff(all.vars(formula[[3]]), character(0))
    data_gen <- function(m) {
      base <- as.data.frame(matrix(stats::rnorm(m * length(vars), 2, 1),
                                   ncol = length(vars)))
      names(base) <- vars
      base
    }
  }
  resp <- all.vars(formula[[2]])[1]
  covered <- logical(n_sims)
  widths <- numeric(n_sims)
  for (i in seq_len(n_sims)) {
    d <- data_gen(n)
    X <- stats::model.matrix(formula[-2], data = d)
    betas <- as.numeric(true_betas[colnames(X)])
    if (length(betas) != ncol(X) || anyNA(betas)) {
      stop("`true_betas` must be named to match the model matrix columns.",
           call. = FALSE)
    }
    d[[resp]] <- as.vector(X %*% betas) + stats::rnorm(nrow(d), 0, sigma)
    fit <- tryCatch(
      dp_lm(formula, d, epsilon = epsilon, delta = delta,
            bounds = list(y = y_bounds)),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      covered[i] <- NA; widths[i] <- NA; next
    }
    ci <- dp_confint(fit, level = level, method = confint_method)
    b <- ci[, "estimate"] # names align with true_betas
    nm <- rownames(ci)
    covered[i] <- all(ci[, "lower"] <= true_betas[nm] &
                        true_betas[nm] <= ci[, "upper"])
    widths[i] <- mean(ci[, "upper"] - ci[, "lower"])
  }
  ok <- !is.na(covered)
  list(
    coverage_rate = mean(covered[ok]),
    mean_ci_width = mean(widths[ok]),
    target_coverage = level,
    n_sims = n_sims
  )
}

#' Compare composition rules over a workflow
#'
#' Given a sequence of per-release epsilons, returns total privacy cost under
#' basic, advanced, and RDP composition — useful for the research question of
#' which rule dominates in practice.
#'
#' @param epsilons Numeric vector of per-release epsilons.
#' @param delta Target delta in (0, 1).
#' @return Data frame with one row per composition rule and total epsilon.
#' @examples
#' compare_composition(c(1, 0.5, 0.5, 0.25), 1e-6)
#' @export
compare_composition <- function(epsilons, delta = 1e-6) {
  basic <- basic_composition(epsilons)$epsilon
  k <- length(epsilons)
  adv <- vapply(seq_len(k), function(j)
    advanced_composition_epsilon(epsilons[seq_len(j - 1)], epsilons[j],
                                 delta), numeric(1))[k]
  rdpc <- rdp_composition(epsilons, delta)$epsilon
  data.frame(
    composition = c("basic", "advanced", "rdp"),
    epsilon = c(basic, adv, rdpc)
  )
}