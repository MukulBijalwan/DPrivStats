test_that("dp_lm recovers coefficients approximately", {
  set.seed(21)
  d <- data.frame(x = rnorm(400), z = rnorm(400))
  d$y <- 2 + 3 * d$x - d$z + rnorm(400)
  fits <- replicate(20,
    coef(dp_lm(y ~ x + z, d, epsilon = 5.0, delta = 1e-6,
               bounds = list(y = c(-30, 30))))
  )
  avg <- rowMeans(fits)
  expect_lt(abs(avg["x"] - 3), 0.5)
  expect_lt(abs(avg["z"] + 1), 0.5)
  fit <- dp_lm(y ~ x + z, d, epsilon = 5, delta = 1e-6,
               bounds = list(y = c(-30, 30)))
  expect_s3_class(fit, "dp_lm")
  expect_s3_class(fit, "lm")
  expect_length(fit$coefficients, 3)
  expect_output(print(fit), "Differentially Private")
  expect_output(print(summary(fit)), "Summary")
})

test_that("dp_lm validation errors", {
  set.seed(22)
  d <- data.frame(x = rnorm(50)); d$y <- rnorm(50)
  expect_error(dp_lm(y ~ x, d, 1, delta = 1e-6), "bounds")
  expect_error(dp_lm(y ~ x, d, 1, delta = NULL,
                     bounds = list(y = c(-5, 5))), "delta")
  expect_error(dp_lm(y ~ x, d, -1, delta = 1e-6,
                     bounds = list(y = c(-5, 5))), "epsilon")
})

test_that("confint methods all work", {
  set.seed(23)
  d <- data.frame(x = rnorm(200)); d$y <- 1 + 2 * d$x + rnorm(200)
  fit <- dp_lm(y ~ x, d, epsilon = 3, delta = 1e-6,
               bounds = list(y = c(-15, 15)))
  ci_a <- dp_confint(fit, method = "analytical")
  ci_p <- dp_confint(fit, method = "parametric_bootstrap", B = 50)
  ci_b <- suppressWarnings(
    dp_confint(fit, method = "privacy_aware_bootstrap", B = 10)
  )
  for (ci in list(ci_a, ci_p, ci_b)) {
    expect_s3_class(ci, "dp_confint")
    expect_true(all(ci[, "lower"] < ci[, "upper"]))
  }
  # analytical CI should contain the OLS estimate most of the time
  expect_s3_class(dp_confint(fit, parm = "x"), "dp_confint")
})

test_that("dp_glm logistic DP-SGD moves toward signal", {
  set.seed(24)
  d <- data.frame(x = rnorm(600))
  d$y <- rbinom(600, 1, plogis(1.5 * d$x))
  set.seed(25)
  fits <- sapply(1:5, function(i) {
    set.seed(100 + i)
    unname(dp_glm(y ~ x, d, binomial(), epsilon = 8, delta = 1e-6,
                  max_grad_norm = 1, n_iter = 300,
                  lr = 0.05)$coefficients["x"])
  })
  expect_gt(mean(fits), 0.1)
  g <- dp_glm(y ~ x, d, binomial(), epsilon = 4, delta = 1e-6,
              n_iter = 100)
  expect_s3_class(g, "dp_glm")
  expect_output(print(g), "DP-SGD")
})

test_that("diagnostics: coverage and utility comparison", {
  set.seed(26)
  gen <- function(n) {
    d <- data.frame(x = stats::rnorm(n, 2, 1))
    d
  }
  vc <- validate_coverage(y ~ x, c(`(Intercept)` = 1, x = 2), sigma = 1,
                          data_gen = gen, n = 300, n_sims = 25,
                          epsilon = 4, delta = 1e-6,
                          y_bounds = c(-20, 20))
  expect_true(vc$coverage_rate >= 0.7 && vc$coverage_rate <= 1)
  expect_equal(vc$target_coverage, 0.95)

  dd <- gen(150); dd$y <- 1 + 2 * dd$x + rnorm(150)
  cmp <- compare_utility(y ~ x, dd, c(1, 5), 1e-6,
                         list(y = c(-15, 15)), n_reps = 3)
  expect_setequal(unique(cmp$epsilon), c(1, 5))

  cc <- compare_composition(c(1, 0.5, 0.5), 1e-6)
  expect_equal(cc$composition, c("basic", "advanced", "rdp"))
  expect_equal(cc$epsilon[cc$composition == "basic"], 2)
})
