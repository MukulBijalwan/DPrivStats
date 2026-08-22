## data-raw/generate_data.R
##
## Regenerates data/example_microdata.rda.
## Run from the package root with devtools installed:
##   Rscript data-raw/generate_data.R

set.seed(20260822)
n <- 2000

education <- pmin(pmax(round(rnorm(n, 13, 3)), 0), 20)
age <- pmin(pmax(round(rnorm(n, 42, 12)), 18), 80)
hours <- pmin(pmax(rnorm(n, 38, 8), 0), 80)
region <- factor(sample(c("North", "South", "East", "West"), n,
                        replace = TRUE))
income <- pmax(
  -30000 + 2500 * education + 350 * age + 400 * hours +
    5000 * (region == "West") + rnorm(n, 0, 12000),
  0
)

example_microdata <- data.frame(
  education = education,
  age = age,
  hours = hours,
  region = region,
  income = round(income)
)

save(example_microdata, file = "data/example_microdata.rda",
     compress = "bzip2", version = 3)

