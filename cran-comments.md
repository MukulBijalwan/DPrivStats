# cran-comments.md

## Submission

Initial submission of **DPrivStats 0.1.0**.

## R CMD check results

- 0 errors | 0 warnings | 0 notes (to be confirmed on the CRAN maintainers'
  check; local checks were run on Windows, R 4.x with testthat edition 3)

## Reviewer notes

- This package implements differentially private statistical procedures.
  All randomness is seeded by the user; examples use `set.seed()` for
  reproducibility and run quickly (< 5 s each).
- `src/` contains a small Rcpp layer (`mechanisms.cpp`) for fast Laplace
  sampling and gradient clipping; compiled with C++11 defaults via Rcpp.
- The bundled dataset `example_microdata` is fully synthetic (generated in
  `data-raw/generate_data.R`); it contains no real personal data.
- `paper/` contains documentation artifacts only and is excluded from the
  build via `.Rbuildignore`.
- Vignettes use `rmarkdown::html_vignette` and are precompiled-friendly;
  no internet access is required by any example, test, or vignette.

## Test environments

- GitHub Actions CI: windows-latest, macOS-latest, ubuntu-latest (release
  and devel), see `.github/workflows/R-CMD-check.yaml`.
