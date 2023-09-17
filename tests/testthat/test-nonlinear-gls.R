test_that("nonlinear_gls_estimator() works for estimating the common acceleration factor.", {
  # Run the nonlinear GLS estimator.
  nl_gls_spline = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline"
  )
  nl_gls_linear = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "linear"
  )
  check_object_spline = list(
    estimates = c(22.5404426, 22.9390485, 25.6086292, 28.3469626, 30.2499492, 0.8550425),
    criterion = 3.42888792
  )
  check_object_linear = list(
    estimates = c(22.5414633, 22.8546460, 25.5366020, 28.4245472, 30.2297674, 0.8596865),
    criterion = 1.8174775
  )
  expect_equal(nl_gls_spline, check_object_spline, ignore_attr = "names")
  expect_equal(nl_gls_linear, check_object_linear, ignore_attr = "names")
})


test_that("nonlinear_gls_conf_int_common() works.", {
  # Run the nonlinear GLS estimator.
  conf_int_spline = nonlinear_gls_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline"
  )
  conf_int_linear = nonlinear_gls_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "linear"
  )
  check_objects_spline = c(0.70102076, 1.02806351)
  check_objects_linear = c(0.73293483, 0.98855244)
  expect_equal(conf_int_spline, check_objects_spline, ignore_attr = "names")
  expect_equal(conf_int_linear, check_objects_linear, ignore_attr = "names")
})
