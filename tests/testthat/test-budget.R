test_that("budget spend and composition rules", {
  b <- new_privacy_budget(3, 1e-6, "basic")
  b <- spend(b, 1, description = "mean")
  b <- spend(b, 1, description = "regression")
  expect_equal(b$epsilon_spent, 2)
  expect_equal(b$remaining_epsilon, 1)
  expect_s3_class(b, "privacy_budget")

  expect_warning(spend(b, 5), "exceeded")

  # RDP is tighter than basic for many small releases
  eps_seq <- rep(0.1, 20)
  expect_true(rdp_composition(eps_seq)$epsilon <
                basic_composition(eps_seq)$epsilon)
})

test_that("advanced composition formula", {
  k <- 10; eps <- 0.1; dp <- 1e-6
  expected <- sqrt(2 * k * log(1 / dp)) * eps + k * eps * (exp(eps) - 1)
  expect_equal(advanced_composition_epsilon(rep(eps, k - 1), eps, dp),
               expected)
})

test_that("can_spend probes correctly", {
  b <- new_privacy_budget(2, 1e-6, "basic")
  expect_true(can_spend(b, 1))
  expect_false(can_spend(b, 3))
})

test_that("rdp conversions are self-consistent", {
  expect_equal(epsilon_to_rdp(1), 0.5)
  expect_equal(rdp_to_epsilon(0, 1e-6), 0)
  expect_true(rdp_to_epsilon(0.5, 1e-6) > 0.5)
})

test_that("PLRV helpers", {
  expect_equal(laplace_plr_tail(3, 1), exp(-3))
  expect_equal(laplace_plr_quantile(c(0, 0.5, 1), 2),
               c(0, -log(0.5) / 2, Inf))
  set.seed(3)
  l <- simulate_gaussian_losses(1000, 5)
  expect_length(l, 1000)
  expect_true(all(l >= 0))
})
