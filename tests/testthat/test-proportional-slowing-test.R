test_that("Test for proportional slowing works for test based on GLS criterion", {
  # Test for proportional slowing using measurement occasion 1 to 4.
  test1_4 = proportional_slowing_gls_test(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    interpolation = "spline",
    vcov = vcov_mmrm,
    j = 1:4,
    start_gamma = 1
  )
  # Test for proportional slowing using measurement occasion 2 to 4.
  test2_4 = proportional_slowing_gls_test(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    interpolation = "spline",
    vcov = vcov_mmrm,
    j = 2:4,
    start_gamma = 1
  )
  # Test for proportional slowing using measurement occasion 3 to 4.
  test3_4 = proportional_slowing_gls_test(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    interpolation = "spline",
    vcov = vcov_mmrm,
    j = 3:4,
    start_gamma = 1
  )
  expect_equal(
    test1_4,
    c(chi_squared = 3.412201487001, df = 3, p_value = 0.332329086671)
  )
  expect_equal(
    test2_4,
    c(chi_squared = 3.408077341064, df = 2, p_value = 0.181947213346)
  )
  expect_equal(
    test3_4,
    c(chi_squared = 0.750142563188, df = 1, p_value = 0.386431098178)
  )
})
