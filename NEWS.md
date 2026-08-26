# DPrivStats 0.1.0

* Initial CRAN submission.
* Privacy mechanisms: Laplace (pure ε-DP), Gaussian (classic calibration),
  analytic Gaussian calibration (Balle & Wang, 2018), exponential mechanism.
* DP descriptive statistics: mean, variance, quantiles, median, histogram.
* DP hypothesis tests: two-sample t-test, chi-square test of independence,
  Kolmogorov–Smirnov test, one-way ANOVA.
* Regression: closed-form DP linear regression (`dp_lm`) and DP-SGD for
  GLMs (`dp_glm`).
* Privacy-aware confidence intervals: analytical, parametric bootstrap,
  and privacy-aware bootstrap.
* Privacy budget tracker with basic, advanced, and Rényi/zCDP composition.
* Diagnostics: coverage validation, utility comparison, composition comparison.
