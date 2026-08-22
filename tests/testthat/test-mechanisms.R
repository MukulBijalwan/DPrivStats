test_that("rlaplace has correct moments and support", {
  set.seed(1)
  z <- rlaplace(200000, 0, 2)
  expect_length(z, 200000)
  expect_lt(abs(mean(z)), 0.1)          # mean ~ 0
  expect_lt(abs(stats::var(z) - 2 * 4), 0.3) # var = 2*scale^2 = 8
})

test_that("mechanism input validation", {
  expect_error(laplace_mechanism(mean, 1:10, -1, 1), "epsilon")
  expect_error(laplace_mechanism(mean, 1:10, 1, -1), "sensitivity")
  expect_error(gaussian_mechanism(mean, 1:10, 1, 0, 1), "delta")
  expect_error(gaussian_mechanism(mean, 1:10, 1, 1.5, 1), "delta")
})

test_that("analytic sigma is tighter than classic sigma", {
  expect_true(analytic_gaussian_sigma(1, 1e-6, 1) <=
                gaussian_sigma(1, 1e-6, 1))
  expect_equal(analytic_gaussian_sigma(1, 1e-6, 0), 0)
})

test_that("exponential mechanism concentrates on best score", {
  set.seed(2)
  candidates <- 1:100
  scores <- -(candidates - 42)^2
  draws <- replicate(500, exponential_mechanism(candidates, scores, 2, 1))
  # mode should be near 42
  expect_lt(mean(abs(draws - 42)), 15)
})
