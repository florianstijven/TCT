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
    estimates = c(22.54145917, 22.85460550, 25.53660298, 28.42443752, 30.22970660, 0.85969996),
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
  expect_equal(se_spline, 0.080210856)
  expect_equal(se_linear, 0.064338500)
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


test_that("nonlinear_gls_conf_int_common() works correctly in previously failing setting", {
  # We start by setting up the setting where the function failed. Specifically,
  # the CI did not contain the estimate.
  time_points = c(0, 6, 12, 18, 24, 36)
  ctrl_estimates = c(19.32327, 20.03718, 20.30573, 22.01748, 23.35910, 26.72632)
  exp_estimates = c(20.23691, 20.71274, 20.78565, 20.98616, 22.53080)
  vcov = matrix(
    data = c(0.04664180,  0.04117288, 0.04649988, 0.05743188, 0.05469797, 0.06234498,  0.04117288, 0.04649988, 0.05743188, 0.05469797,  0.06234498,
             0.04117288,  0.08381315, 0.07175819, 0.08765594, 0.08199624, 0.09586312,  0.03634521, 0.04104760, 0.05069778, 0.04828444,  0.05503481,
             0.04649988,  0.07175819, 0.10232254, 0.10905066, 0.10505128, 0.12495016,  0.04104760, 0.04635839, 0.05725712, 0.05453153,  0.06215528,
             0.05743188,  0.08765594, 0.10905066, 0.15706612, 0.13673738, 0.17093999,  0.05069778, 0.05725712, 0.07071812, 0.06735175,  0.07676782,
             0.05469797,  0.08199624, 0.10505128, 0.13673738, 0.16097673, 0.18180898,  0.04828444, 0.05453153, 0.06735175, 0.06414564,  0.07311347,
             0.06234498,  0.09586312, 0.12495016, 0.17093999, 0.18180898, 0.29281802,  0.05503481, 0.06215528, 0.07676782, 0.07311347,  0.08333505,
             0.04117288,  0.03634521, 0.04104760, 0.05069778, 0.04828444, 0.05503481,  0.08381315, 0.07175819, 0.08765594, 0.08199624,  0.09586312,
             0.04649988,  0.04104760, 0.04635839, 0.05725712, 0.05453153, 0.06215528,  0.07175819, 0.10232254, 0.10905066, 0.10505128,  0.12495016,
             0.05743188,  0.05069778, 0.05725712, 0.07071812, 0.06735175, 0.07676782,  0.08765594, 0.10905066, 0.15706612, 0.13673738,  0.17093999,
             0.05469797,  0.04828444, 0.05453153, 0.06735175, 0.06414564, 0.07311347,  0.08199624, 0.10505128, 0.13673738, 0.16097673,  0.18180898,
             0.06234498,  0.05503481, 0.06215528, 0.07676782, 0.07311347, 0.08333505,  0.09586312, 0.12495016, 0.17093999, 0.18180898,  0.29281802),
          nrow = 11, ncol = 11
  )
  TCT_fit = TCT_meta(
    time_points = time_points,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    interpolation = "spline",
    vcov = vcov
  )
  TCT_common = TCT_meta_common(TCT_fit, inference = "least-squares", start_gamma = 0.5)
  sum_TCT = summary(TCT_common)
  conf_int_estimated = sum_TCT$gamma_common_ci
  conf_int_expected = matrix(c(0.42446364, 0.57704044), nrow = 1)
  expect_equal(conf_int_estimated, conf_int_expected)
})
