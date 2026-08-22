test_that("dp_t_test detects shift with enough epsilon", {
  set.seed(41)
  pvals <- replicate(50, {
    tt <- dp_t_test(rnorm(150), rnorm(150, 0.6), 2.0, bounds = c(-5, 5))
    tt$p.value
  })
  expect_lt(mean(pvals < 0.05), 0.9)
  expect_gt(mean(pvals < 0.05), 0.1)
  tt <- dp_t_test(rnorm(50), rnorm(50), 2.0, c(-4, 4),
                  alternative = "greater")
  expect_s3_class(tt, "dp_htest")
  expect_output(print(tt), "t-test")
})

test_that("dp_chisq_test runs on a table", {
  set.seed(42)
  tab <- matrix(c(80, 20, 30, 70), nrow = 2)
  res <- dp_chisq_test(tab, 3.0)
  expect_s3_class(res, "dp_htest")
  expect_true(all(as.vector(res$observed_private) >= 0))
  expect_error(dp_chisq_test(matrix(c(-1, 2)), 1), "non-negative")
})

test_that("dp_ks_test and dp_anova return htest objects", {
  set.seed(43)
  ks <- dp_ks_test(rnorm(200), rnorm(200, 1.5), 2.0, c(-6, 6))
  expect_s3_class(ks, c("dp_htest", "htest"))
  expect_true(ks$statistic >= 0)

  dat <- data.frame(y = c(rnorm(90), rnorm(90, 2)),
                    g = factor(rep(c("A", "B"), each = 90)))
  an <- dp_anova(y ~ g, dat, 3.0, c(-8, 8))
  expect_s3_class(an, c("dp_htest", "htest"))
  expect_named(an$group_means_private, c("A", "B"))
  expect_error(dp_anova(y ~ g, dat[dat$g == "A", ], 1, c(-5, 5)),
               "at least 2 groups")
})
