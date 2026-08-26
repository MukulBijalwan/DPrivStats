# DPrivStats: An R Package for Differentially Private Classical Statistical Inference

**Mukul Bijalwan**
*mukulbijalwan555@gmail.com*

**Package version:** 0.1.0 · **License:** MIT · **URL:** https://github.com/MukulBijalwan/DPrivStats

---

## Abstract

We present **DPrivStats**, an R package implementing differentially private (DP) versions of common classical statistical procedures with first-class support for privacy budget accounting and privacy-aware uncertainty quantification. The package provides the Laplace mechanism (pure ε-DP), the classic and analytically calibrated Gaussian mechanisms, and the exponential mechanism; DP descriptive statistics (mean, variance, quantiles, median, histogram); DP hypothesis tests (two-sample t-test, chi-square test of independence, Kolmogorov–Smirnov test, one-way ANOVA); and regression via closed-form DP linear regression (`dp_lm`) and DP-SGD for generalized linear models (`dp_glm`). A distinguishing feature is its suite of confidence-interval procedures that jointly account for sampling variability and privacy noise—analytical intervals, parametric bootstrap, and a privacy-aware bootstrap that re-privatizes each resample—as well as an S3 budget tracker supporting basic, advanced, and Rényi/zCDP composition rules. Diagnostics for Monte Carlo coverage validation and utility–privacy trade-off studies are included. We describe the design, mathematical foundations, usage, and validation workflows of the package.

**Keywords:** differential privacy, statistical inference, R, privacy budgets, confidence intervals, official statistics.

---

## 1. Introduction

Statistical agencies and data analysts increasingly must publish summaries of sensitive microdata while providing formal guarantees that individual records are protected. Differential privacy (DP) (Dwork et al., 2006) has become the de facto standard: a randomized mechanism M satisfies (ε, δ)-DP if for all neighbouring datasets x ≅ x′ differing in one record,

  P[M(x) ∈ S] ≤ e^ε · P[M(x′) ∈ S] + δ,  for all S ⊆ Range(M).

While the theory of DP statistics is mature, practical software for *classical* inferential tasks—means, variances, quantiles, hypothesis tests, linear models—with coherent privacy accounting remains scarce in R. DPrivStats fills this gap by offering:

1. **Mechanisms** with rigorous sensitivity-based calibration;
2. **DP estimators and tests** whose outputs remain valid through post-processing;
3. **Privacy-aware confidence intervals** combining sampling and privacy noise;
4. **Budget management** under basic, advanced, and Rényi (zCDP) composition;
5. **Diagnostics** for empirical coverage and utility across ε-grids.

Section 2 reviews the mathematical machinery. Section 3 describes package architecture. Sections 4–7 cover mechanisms, descriptive statistics, inference (tests and regression), and confidence intervals. Section 8 presents budget accounting, Section 9 diagnostics, and Section 10 a worked workflow for official-statistics tabulation under a shared budget.

---

## 2. Mathematical Foundations

### 2.1 Mechanisms

**Laplace mechanism.** For a statistic f with L1 sensitivity Δ₁ = max_{x≅x′} ‖f(x) − f(x′)‖₁, releasing f(x) + Lap(0, Δ₁/ε) yields pure ε-DP. Noise is drawn by inverse-CDF:

  Z = − sgn(u) · s · log(1 − 2|u|), u ~ Unif(−½, ½), s = Δ₁/ε.

**Gaussian mechanism (classic calibration).** With L2 sensitivity Δ₂, adding N(0, σ²) noise where

  σ = Δ₂ √(2 log(1.25/δ)) / ε

gives (ε, δ)-DP (Dwork & Roth, 2014).

**Analytic Gaussian calibration.** DPrivStats also implements the tight analytic Gaussian mechanism of Balle & Wang (2018), which solves for the minimal σ satisfying the DP constraint directly rather than relying on the sufficient classic bound, typically reducing noise for a given (ε, δ).

**Exponential mechanism.** Used for quantiles/median selection with a utility score u; output probability ∝ exp(ε u(x, r)/2Δu).

### 2.2 Sensitivity of bounded statistics

All descriptive estimators clip inputs to known bounds [L, U]:

- mean: Δ₁ = (U − L)/n
- variance: bounded by (U − L)²/4 per pair-difference argument
- counts (histogram cells): Δ₁ = 1
- quantiles (exponential mechanism): score sensitivity Δu/n

For linear regression, the closed-form estimator β̂ = (XᵀX)⁻¹Xᵀy on clipped data has L2 sensitivity

  Δ₂ = C·D / λ_min(XᵀX),

where C bounds the row norms ‖(x_i, y_i)‖ (or only ‖y_i‖ ≤ U_y when the design is assumed public), D is the scale of the release, and λ_min(XᵀX) is the smallest eigenvalue of the Gram matrix (`dp_lm_sensitivity`). Predictors are clipped column-wise so no single record can move the released fit arbitrarily.

### 2.3 Post-processing invariance

Because post-processing cannot degrade DP, all downstream quantities—t-statistics, p-values, chi-square statistics, fitted values—are computed from privatized primitives only, inheriting the privacy guarantee without additional budget.

---

## 3. Package Architecture

```
DPrivStats/
├── R/
│   ├── mechanisms.R, mechanisms2.R   # Laplace, Gaussian, analytic Gaussian,
│   │                                 #   exponential mechanism
│   ├── dp_descriptive.R/.R2          # dp_mean, dp_variance, dp_quantile,
│   │                                 #   dp_median, dp_histogram
│   ├── dp_tests.R, dp_tests2.R       # dp_t_test, dp_chisq_test,
│   │                                 #   dp_ks_test, dp_anova
│   ├── dp_regression.R/.R2           # dp_lm, dp_glm (DP-SGD)
│   ├── dp_confint.R                  # privacy-aware confidence intervals
│   ├── accounting.R, composition.R   # epsilon<->RDP conversion, composition rules
│   ├── privacy_budget.R              # S3 "privacy_budget" tracker (spend/can_spend)
│   ├── sensitivity.R                 # stat_sensitivities, dp_lm_sensitivity
│   ├── simulation.R, utility.R       # validate_coverage, compare_utility,
│   │                                 #   compare_composition
│   └── print_methods.R, utils.R      # S3 print methods, clipping helpers
├── src/mechanisms.cpp                # Rcpp-accelerated Laplace draws & gradient clipping
├── data/example_microdata.rda        # synthetic public-use example dataset
├── tests/testthat/                   # unit + coverage/regression tests
└── vignettes/                        # getting-started, regression-guide,
                                      #   official-statistics
```

Key S3 classes:

| Class | Producer | Purpose |
|---|---|---|
| `dp_estimate` | `dp_mean`, `dp_variance`, … | estimate + metadata (ε, δ, mechanism, sensitivity, noise_variance) |
| `dp_htest` / `htest` | `dp_t_test`, `dp_chisq_test`, … | test results mirroring base-`htest` structure |
| `dp_lm` / `lm` | `dp_lm` | privatized coefficients plus `noise_variance`, `sampling_vcov` |
| `dp_confint` | `dp_confint` | interval matrix with estimate/lower/upper |
| `privacy_budget` | `new_privacy_budget` | cumulative-loss tracker with audit history |

Every result object carries its privacy parameters and noise variance, enabling downstream inference without re-specifying the analysis.

---

## 4. Mechanisms in Practice

```r
library(DPrivStats)

# Pure epsilon-DP via Laplace
laplace_mechanism(mean, 1:100, epsilon = 1.0, l1_sensitivity = 100/99)

# Classic Gaussian calibration
gaussian_sigma(epsilon = 1.0, delta = 1e-6, l2_sensitivity = 10)
gaussian_mechanism(sum, 1:100, 1.0, 1e-6, l2_sensitivity = 100)

# Analytic Gaussian (tighter sigma than classic)
analytic_gaussian_sigma(epsilon = 1.0, delta = 1e-6, l2_sensitivity = 10)
```

The analytic Gaussian implementation follows Balle & Wang (2018): it solves the calibration equation for σ given (ε, δ, Δ₂), guaranteeing (ε, δ)-DP with strictly less noise than the classic bound when δ < e^(−ε²/8)/2.

---

## 5. DP Descriptive Statistics

```r
data(example_microdata)
m <- dp_mean(example_microdata$income, epsilon = 0.5, bounds = c(0, 500000))
m$estimate; m$noise_variance
```

Design points:

- **Clipping** to `bounds` bounds sensitivity; defaults fall back to the observed range for exploration, but *a priori* bounds should be used in production.
- **Mechanism choice**: `"laplace"` (pure DP), `"gaussian"`, or `"analytic_gaussian"` (requires `delta`).
- **Validity by construction**: histograms truncate negative noisy counts; quantiles use the exponential mechanism so outputs always lie in the supported range.
- Each returned `dp_estimate` records `noise_variance` exactly from the mechanism calibration, which powers the analytical intervals of Section 7.

---

## 6. Hypothesis Testing and Regression

### 6.1 DP two-sample t-test

```r
set.seed(7)
tt <- dp_t_test(rnorm(200), rnorm(200, 0.5), epsilon = 1.0,
                bounds = c(-5, 5))
print(tt)
```

The total ε is split evenly (ε/3 each) across group counts, means, and variances; private counts are floored at 2 (post-processing) to prevent division blow-ups, and the t-statistic, degrees of freedom, and p-value are computed entirely from privatized quantities.

### 6.2 Chi-square test of independence

Each cell receives Laplace(1/ε) noise; negative counts are truncated and rounded before calling `stats::chisq.test`. The privatized table is retained as `observed_private`.

### 6.3 Kolmogorov–Smirnov and ANOVA

`dp_ks_test` privatizes ECDF evaluations and computes p-values against the asymptotic KS distribution (`asymptotic_ks_pvalue`); `dp_anova` privatizes group means/variances/counts and reconstructs the F statistic by post-processing.

### 6.4 Closed-form DP linear regression (`dp_lm`)

```r
fit <- dp_lm(income ~ education + age + hours, example_microdata,
             epsilon = 2.0, delta = 1e-6,
             bounds = list(y = c(0, 500000)))
summary(fit); coef(fit)
```

OLS is fit on clipped data via a pseudo-inverse of XᵀX (`MASS::ginv` for numerical stability), and coefficients are released under the analytic Gaussian mechanism calibrated to Δ₂ = C·D/λ_min(XᵀX). The fit stores both the exact privacy noise variance per coefficient and the classical sampling covariance estimate σ̂²(XᵀX)⁻¹, consumed by `dp_confint`.

### 6.5 DP-SGD for GLMs (`dp_glm`)

For logistic/Poisson-type responses, `dp_glm` performs differentially private stochastic gradient descent with per-example gradient clipping (`clip_gradients`, with an Rcpp fast path `cpp_clip_gradients`) and Gaussian noise on gradients; privacy loss accumulates over iterations via the package's RDP accountant.

---

## 7. Privacy-Aware Confidence Intervals

Standard CIs ignore privacy noise; DPrivStats offers three methods in `dp_confint(fit, method = ...)`:

1. **`"analytical"`** — Normal intervals with standard error
   √( diag(Σ_sampling) + v_privacy ), where v_privacy is the exact mechanism noise variance:
   ```r
   dp_confint(fit)                          # default: analytical, 95%
   confint.dp_lm(fit, level = 0.90)
   ```
2. **`"parametric_bootstrap"`** — draws β* ~ N(β̂_priv, Σ_sampling), adds fresh privacy noise, returns percentile intervals.
3. **`"privacy_aware_bootstrap"`** — resamples rows and refits the full DP estimator per replicate with per-resample budget ε/√B (advanced-composition heuristic). Most conservative and computationally heaviest; falls back to analytical if too many refits fail.

Empirical coverage of these intervals can be validated with `validate_coverage()` (Section 9).

---

## 8. Privacy Budget Accounting

An S3 tracker records every release with a description and timestamp:

```r
b <- new_privacy_budget(epsilon = 3.0, delta = 1e-6, composition = "rdp")
b <- spend(b, 0.5, description = "DP mean income")
b <- spend(b, 1.0, description = "DP lm")
can_spend(b, 2.0)          # FALSE: proposed spend does not fit
print(b)                    # totals, rho accumulator, remaining budget
```

Composition rules implemented:

| Rule | Formula | Function |
|---|---|---|
| Basic (sequential) | Σεᵢ, Σδᵢ | `basic_composition()` |
| Advanced (homogeneous k releases) | √(2k log(1/δ′))·ε + kε(e^ε−1) | `advanced_composition_epsilon()` |
| Rényi / zCDP | ρ = Σεᵢ²/2; ε = ρ + √(2ρ log(1/δ)) | `epsilon_to_rdp()`, `rdp_to_epsilon()`, `rdp_composition()` |

RDP/zCDP composition is typically much tighter for many small releases—for example, twenty ε = 0.1 releases cost ε_total = 2.0 under basic composition but far less under `rdp_composition(rep(0.1, 20), 1e-6)`.

---

## 9. Diagnostics

- **Coverage validation** — `validate_coverage()` runs Monte Carlo replications of a DP procedure and reports empirical coverage of the analytical/bootstrap intervals, exposing under-coverage from privacy noise.
- **Utility comparison** — `compare_utility()` sweeps an ε-grid and summarizes estimation error (e.g., MSE/bias) of competing mechanisms or estimators.
- **Composition comparison** — `compare_composition()` contrasts basic vs advanced vs RDP total-loss curves for planned release schedules, aiding budget planning.

These functions double as regression tests in `tests/testthat/`, guarding statistical behaviour (sensitivity formulas, monotonicity of loss, CI containment).

---

## 10. Worked Workflow: Official Statistics Under a Shared Budget

```r
library(DPrivStats)
data(example_microdata)

# 1. Plan a total privacy budget for the release cycle
budget <- new_privacy_budget(epsilon = 3.0, delta = 1e-6, composition = "rdp")

# 2. Allocate pieces to each published statistic
stopifnot(can_spend(budget, 0.5)); budget <- spend(budget, 0.5, "mean income")
stopifnot(can_spend(budget, 0.5)); budget <- spend(budget, 0.5, "income histogram")

h <- dp_histogram(example_microdata$income, epsilon = 0.5,
                  breaks = seq(0, 500000, by = 50000))

# 3. Regression with privacy-aware intervals
stopifnot(can_spend(budget, 2.0)); budget <- spend(budget, 2.0, "wage regression")
fit <- dp_lm(income ~ education + age + hours, example_microdata,
             epsilon = 2.0, delta = 1e-6, bounds = list(y = c(0, 500000)))
dp_confint(fit)

print(budget)  # auditable ledger of the entire release
```

Because every call consumes declared budget up front and every output object carries its own (ε, δ) provenance, the workflow yields a reproducible, auditable privacy ledger suitable for disclosure-control review.

---

## 11. Installation and Reproducibility

```r
# install.packages("remotes")
remotes::install_github("MukulBijalwan/DPrivStats")
```

Requirements: R ≥ 4.0.0; imports stats, MASS, Rcpp (LinkingTo Rcpp); suggests testthat (≥ 3.0.0), knitr, rmarkdown, ggplot2. Continuous integration runs `R CMD check` via GitHub Actions (`.github/workflows/R-CMD-check.yaml`). Vignettes: `getting-started`, `regression-guide`, `official-statistics` (`browseVignettes("DPrivStats")`). The bundled `example_microdata` is synthetically generated by `data-raw/generate_data.R` for reproducibility.

## 12. Limitations

- Bounds defaulting to the observed range leak range information; analysts should supply *a priori* bounds in practice.
- Advanced composition assumes homogeneous per-release ε.
- The privacy-aware bootstrap's ε/√B allocation is a heuristic; formal adaptive-composition accounting is left to future work.
- `dp_lm` sensitivity assumes well-conditioned designs; near-singular Gram matrices inflate Δ₂ through small λ_min.

## References

1. Dwork, C., McSherry, F., Nissim, K., & Smith, A. (2006). Calibrating noise to sensitivity in private data analysis. *TCC 2006*.
2. Dwork, C., & Roth, A. (2014). The algorithmic foundations of differential privacy. *Foundations and Trends in Theoretical Computer Science*, 9(3–4).
3. Balle, B., & Wang, Y.-X. (2018). Improving the Gaussian mechanism for differential privacy: Analytical calibration and optimal denoising. *ICML 2018*.
4. Bun, M., & Steinke, T. (2016). Concentrated differential privacy: Simplifications, extensions, and lower bounds. *TCC 2016-B*.
5. Abadi, M. et al. (2016). Deep learning with differential privacy. *CCS 2016*. (DP-SGD)

---

*© Mukul Bijalwan, MIT License.*
