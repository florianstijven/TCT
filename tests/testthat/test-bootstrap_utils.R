test_that("bootstrap intervals are correct in a simple setting", {
  # Some arbitary data vector with no NAs or NaNs
  bootstrap_replicates = rep(c(1.2, 2.3, 1.4, 1.5, 1.6, 1.8, 3.9, 2.0), 20)
  # Compute the bootstrap CIs
  BC_percentile_limits = BC_percentile_CI(boot_replicates = bootstrap_replicates, estimate = 2)
  percentile_limits = percentile_CI(bootstrap_replicates)
  # Check that the limits are correct
  expect_equal(BC_percentile_limits, c(1.4, 3.9), ignore_attr = TRUE)
  expect_equal(percentile_limits, c(1.2, 3.9), ignore_attr = TRUE)
  })
