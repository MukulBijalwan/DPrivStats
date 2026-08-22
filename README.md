# DPrivStats — Differentially Private Classical Inference

[![R-CMD-check](https://github.com/MukulBijalwan/DPrivStats/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MukulBijalwan/DPrivStats/actions)

An R package implementing **differentially private (DP)** versions of common
classical statistical procedures, with privacy budget accounting and
privacy-aware confidence intervals.

## Features

- **Mechanisms** — Laplace (pure ε-DP), Gaussian (classic calibration), analytic
  Gaussian calibration (Balle & Wang 2018), exponential mechanism.
- **Descriptive statistics** — DP mean, variance, quantiles, median, histogram,
  all with post-processing to keep outputs valid.
- **Hypothesis tests** — DP two-sample t-test, chi-square test of independence,
  Kolmogorov–Smirnov test, one-way ANOVA.
- **Regression** — closed-form DP linear regression (`dp_lm`) calibrated to the
  L2 sensitivity `Δ₂ = C·D / λmin(XᵀX)`, and DP-SGD for GLMs (`dp_glm`).
- **Privacy-aware confidence intervals** — analytical intervals combining
  sampling and privacy noise variance, parametric bootstrap, and a
  privacy-aware bootstrap that re-privatizes each resample.
- **Privacy budget management** — an S3 budget tracker supporting basic,
  advanced, and Rényi (zCDP) composition rules.
- **Diagnostics** — Monte Carlo coverage validation, utility comparison across
  ε-grids, composition comparison studies.

## Installation

```r
# install.packages("remotes")
remotes::install_github("MukulBijalwan/DPrivStats")
```

Or from a local clone:

```r
install.packages("path/to/DPrivStats", repos = NULL, type = "source")
```

## Quick start

```r
library(DPrivStats)

data(example_microdata)

# DP mean income under a Laplace mechanism
m <- dp_mean(example_microdata$income, epsilon = 0.5,
             bounds = c(0, 500000))
print(m)

# DP linear regression with privacy-aware confidence intervals
fit <- dp_lm(income ~ education + age + hours, example_microdata,
             epsilon = 2.0, delta = 1e-6,
             bounds = list(y = c(0, 500000)))
dp_confint(fit)
```

See the vignettes (`browseVignettes("DPrivStats")`) for full workflows,
including official-statistics tabulation under a shared budget.

## Author

Mukul Bijalwan — <mukulbijalwan555@gmail.com>

## License

MIT © Mukul Bijalwan
