#' Census-like example microdata
#'
#' A synthetic census-style microdata set of 2000 individuals generated for
#' use in examples and vignettes. Income is a linear function of education,
#' age, hours worked, region, plus noise; it is not real personal data.
#'
#' @format A data frame with 2000 rows and 5 variables:
#' \describe{
#'   \item{education}{Years of education (integer, 0-20).}
#'   \item{age}{Age in years (integer, 18-80).}
#'   \item{hours}{Weekly working hours (numeric, 0-80).}
#'   \item{region}{Factor with levels North, South, East, West.}
#'   \item{income}{Annual income (numeric, non-negative).}
#' }
#' @source Simulated; see \code{data-raw/generate_data.R}.
"example_microdata"