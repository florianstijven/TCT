test_that("nonlinear_gls_estimator() works for estimating the common acceleration factor.", {
  # Run the nonlinear GLS estimator.
  nl_gls_spline = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline"
  )
  nl_gls_linear = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "linear"
  )
  check_object_spline = list(
    estimates = c(22.54044773, 22.93901666, 25.60866976, 28.34691764, 30.24994326, 0.85504336),
    criterion = 3.4122015
  )
  check_object_linear = list(
    estimates = c(22.54145662, 22.85460212, 25.53659717, 28.42443216, 30.22969846, 0.85970034),
    criterion = 1.8085288
  )
  expect_equal(nl_gls_spline, check_object_spline, ignore_attr = "names")
  expect_equal(nl_gls_linear, check_object_linear, ignore_attr = "names")
})

test_that("nonlinear_gls_conf_int_common() works.", {
  # Run the nonlinear GLS estimator.
  conf_int_spline = nonlinear_gls_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline"
  )
  conf_int_linear = nonlinear_gls_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "linear"
  )
  check_objects_spline = c(0.70072621, 1.02848644)
  check_objects_linear = c(0.73257342, 0.98887667)
  expect_equal(conf_int_spline, check_objects_spline, ignore_attr = "names")
  expect_equal(conf_int_linear, check_objects_linear, ignore_attr = "names")
})

test_that("nonlinear_gls_estimator_se() works for estimating the SE of the common acceleration factor.", {
  # Run the nonlinear GLS estimator.
  nl_gls_spline = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline"
  )
  nl_gls_linear = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "linear"
  )
  se_spline = nonlinear_gls_estimator_se(
    time_points = 0:4,
    interpolation = "spline",
    vcov = vcov_mmrm,
    j = 1:4,
    gamma_est = nl_gls_spline$estimates[6],
    alpha_est = nl_gls_spline$estimates[1:5]
  )
  se_linear = nonlinear_gls_estimator_se(
    time_points = 0:4,
    interpolation = "linear",
    vcov = vcov_mmrm,
    j = 1:4,
    gamma_est = nl_gls_linear$estimates[6],
    alpha_est = nl_gls_linear$estimates[1:5]
  )
  print(se_spline, digits = 8)
  print(se_linear, digits = 8)
  expect_equal(se_spline, 0.080210856)
  expect_equal(se_linear, 0.064338448)
})


test_that(
  "nonlinear_gls_estimator_vcov() works for estimating the variance-covariance matrix.",
  {
    # Run the nonlinear GLS estimator.
    nl_gls_spline = nonlinear_gls_estimator(
      time_points = 0:4,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      vcov = vcov_mmrm,
      interpolation = "spline"
    )
    nl_gls_linear = nonlinear_gls_estimator(
      time_points = 0:4,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      vcov = vcov_mmrm,
      interpolation = "linear"
    )
    vcov_spline = nonlinear_gls_estimator_vcov(
      time_points = 0:4,
      interpolation = "spline",
      vcov = vcov_mmrm,
      j = 1:4,
      gamma_est = nl_gls_spline$estimates[6],
      alpha_est = nl_gls_spline$estimates[1:5]
    )
    vcov_linear = nonlinear_gls_estimator_vcov(
      time_points = 0:4,
      interpolation = "linear",
      vcov = vcov_mmrm,
      j = 1:4,
      gamma_est = nl_gls_linear$estimates[6],
      alpha_est = nl_gls_linear$estimates[1:5]
    )
    expect_equal(
      vcov_spline[1,],
      c(
        1.4675637e-01,
        1.3563979e-01,
        1.4569601e-01,
        1.4919334e-01,
        1.5865166e-01,
        9.9603518e-05
      ),
      tolerance = 1e-6
    )
    expect_equal(
      vcov_linear[2,],
      c(
        0.13485206382,
        0.19213263889,
        0.19468419248,
        0.22334610970,
        0.24495854755,
        -0.00092349316
      ),
      tolerance = 1e-6
    )
  }
)


test_that("nonlinear_gls_test() works.", {
  # Run the nonlinear GLS estimator.
  nl_gls_spline = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline"
  )
  nl_gls_linear = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "linear"
  )
  # Compute test statistics for gamma_0 = 1
  test_spline = nonlinear_gls_test(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    j = 1:4,
    gamma_0 = 1
  )
  test_linear = nonlinear_gls_test(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "linear",
    j = 1:4,
    gamma_0 = 1
  )
  expect_equal(test_spline,
               c("chi-squared" = 2.826856557, "p-value" = 0.092699643),
               ignore_attr = "names")
  expect_equal(test_linear,
               c("chi-squared" = 4.43052924, "p-value" = 0.03530156),
               ignore_attr = "names")
  # Compute test statistics for gamma_0 = 0.7
  test_spline = nonlinear_gls_test(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    j = 1:4,
    gamma_0 = 0.7
  )
  test_linear = nonlinear_gls_test(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "linear",
    j = 1:4,
    gamma_0 = 0.7
  )
  expect_equal(test_spline,
               c("chi-squared" = 3.879983935, "p-value" = 0.048865025),
               ignore_attr = "names")
  expect_equal(test_linear,
               c("chi-squared" = 5.526636925, "p-value" = 0.018729071),
               ignore_attr = "names")
})
