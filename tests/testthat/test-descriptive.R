test_that("dp_mean / dp_variance / dp_histogram are unbiased-ish and valid", {
  set.seed(31)
  x <- rnorm(2000)
  ests <- replicate(200, dp_mean(x, 1.0, bounds = c(-5, 5),
                                 mechanism = "laplace")$estimate)
  expect_lt(abs(mean(ests) - mean(x)), 0.05)
  m <- dp_mean(x, 1.0, c(-5, 5))
  expect_s3_class(m, "dp_estimate")
  expect_equal(m$sensitivity$l1, 10 / 2000)
  expect_output(print(m), "mean")

  v <- dp_variance(runif(500), 1.0, c(0, 1), "laplace")
  expect_true(v$estimate >= 0)
  expect_error(dp_variance(x, 1.0, c(0, 1), "gaussian", delta = NULL),
               "delta")

  h <- dp_histogram(rnorm(300), 2.0, breaks = seq(-3, 3, by = 1))
  expect_true(all(h$estimate >= 0))
})

test_that("analytic gaussian dp_mean works", {
  set.seed(32)
  est <- dp_mean(rnorm(500), 1.0, c(-5, 5), "analytic_gaussian",
                 delta = 1e-6)
  expect_s3_class(est, "dp_estimate")
  expect_equal(est$delta, 1e-6)
})

test_that("dp_median and dp_quantile land near truth", {
  set.seed(33)
  x <- rcauchy(1000) # heavy tail but bounded grid used
  med <- dp_median(x, 3.0, bounds = c(-10, 10))
  expect_lt(abs(med$estimate - stats::median(x)), 1.5)
  q <- dp_quantile(x, 4.0, c(-10, 10), probs = c(.25, .75))
  expect_length(q$estimate, 2)
  expect_true(all(names(q$estimate) == c("25%", "75%")))
})
