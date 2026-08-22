#' Print method for dp_estimate objects
#'
#' @param x A \code{dp_estimate} object.
#' @param digits Number of digits to print.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x}, invisibly.
#' @export
print.dp_estimate <- function(x, digits = 4, ...) {
  cat(sprintf("\nDP %s (%s mechanism)\n", x$statistic, x$mechanism))
  cat(sprintf("epsilon = %.3f, delta = %s\n",
              x$epsilon, format(x$delta, digits = 2)))
  cat("estimate:\n")
  print(round(x$estimate, digits))
  invisible(x)
}

#' Print method for dp_htest objects
#'
#' @param x A \code{dp_htest} object.
#' @param digits Number of digits to print.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x}, invisibly.
#' @export
print.dp_htest <- function(x, digits = 4, ...) {
  cat(sprintf("\n\t%s\n\n", x$method))
  if (!is.null(x$data.name)) cat("data:", x$data.name, "\n\n")
  stat <- format(round(as.numeric(x$statistic), digits))
  nm <- paste(names(x$statistic), collapse = ", ")
  pv <- if (!is.null(x$p.value)) format.pval(x$p.value, digits = digits) else NA
  cat(sprintf("%s = %s, p-value %s\n", nm, stat, pv))
  if (!is.null(x$parameter)) {
    par_txt <- paste(sprintf("%s = %.1f", names(x$parameter),
                             as.numeric(x$parameter)), collapse = ", ")
    cat(par_txt, "\n")
  }
  if (!is.null(x$alternative)) {
    cat("alternative hypothesis:", x$alternative, "\n")
  }
  if (!is.null(x$std.error)) {
    cat(sprintf("effect size = %.4f (SE = %.4f)\n",
                as.numeric(x$estimate)[1], x$std.error))
  }
  cat(sprintf("\nprivacy cost: epsilon = %.3f, delta = %s\n",
              x$epsilon, format(x$delta, digits = 2)))
  invisible(x)
}

#' Print method for dp_confint objects
#'
#' @param x A \code{dp_confint} matrix.
#' @param digits Number of digits to print.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x}, invisibly.
#' @export
print.dp_confint <- function(x, digits = 4, ...) {
  print(round(unclass(x), digits))
  invisible(x)
}

#' Summarize a dp_lm fit
#'
#' Prints privatized coefficients with analytical privacy-aware standard
#' errors and confidence intervals.
#'
#' @param object A \code{dp_lm} object.
#' @param level Confidence level (default 0.95).
#' @param ... Additional arguments (ignored).
#' @return Invisibly, the list with coefficients table.
#' @method summary dp_lm
#' @export
summary.dp_lm <- function(object, level = 0.95, ...) {
  ci <- dp_confint(object, level = level)
  se_total <- sqrt(diag(object$sampling_vcov) + object$noise_variance)
  tbl <- cbind(object$coefficients, se_total, ci[, "lower"], ci[, "upper"])
  colnames(tbl) <- c("Estimate", "Std. Error", "CI Lower", "CI Upper")
  out <- list(coefficients = tbl, call = object$call,
              epsilon = object$epsilon, delta = object$delta)
  class(out) <- "summary.dp_lm"
  print.summary.dp_lm(out)
  invisible(out)
}

#' Print method for summary.dp_lm objects
#'
#' @param x A \code{summary.dp_lm} object.
#' @param digits Number of digits to print.
#' @param ... Additional arguments (ignored).
#' @return The object \code{x}, invisibly.
#' @export
print.summary.dp_lm <- function(x, digits = 4, ...) {
  cat("\nDifferentially Private Linear Regression Summary\n")
  cat(sprintf("epsilon = %.3f, delta = %.1e\n\n", x$epsilon, x$delta))
  print(round(x$coefficients, digits))
  invisible(x)
}