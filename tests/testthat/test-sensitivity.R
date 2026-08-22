test_that("stat_sensitivities match closed forms", {
  sens <- stat_sensitivities(c(0, 10), 100)
  expect_equal(unname(sens$l1["mean"]), 0.1)
  expect_equal(unname(sens$l1["sum"]), 10)
  expect_equal(unname(sens$l1["count"]), 1)
  expect_equal(unname(sens$l2["variance"]), 1)
  expect_error(stat_sensitivities(c(5, 5), 10), "bounds")
})

test_that("dp_lm_sensitivity matches C*D/lambda_min", {
  X <- cbind(1, 1:10)
  XtX <- crossprod(X)
  lm_ev <- min(eigen(XtX)$values)
  expected <- max(sqrt(rowSums(X^2))) * diff(c(0, 10)) / lm_ev
  expect_equal(dp_lm_sensitivity(X, c(0, 10)), expected)
  expect_error(dp_lm_sensitivity(cbind(1, c(1, 1, 1)), c(0, 1)),
               "rank-deficient")
})

test_that("clip_gradients clips to max norm", {
  g <- matrix(rnorm(60, sd = 10), 20, 3)
  cl <- clip_gradients(g, 1)
  norms <- sqrt(rowSums(cl$clipped^2))
  expect_true(all(norms <= 1 + 1e-12))
  big <- g > 1
  # rows originally under the norm are unchanged
  orig <- sqrt(rowSums(g^2))
  same <- orig <= 1
  expect_equal(cl$clipped[same, , drop = FALSE], g[same, , drop = FALSE])
})

test_that("cpp_clip_gradients agrees with R implementation", {
  skip_if_not_installed("Rcpp")
  g <- matrix(rnorm(30), 10, 3)
  expect_equal(cpp_clip_gradients(g, 1),
               clip_gradients(g, 1)$clipped,
               tolerance = 1e-12)
})
