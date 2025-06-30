test_that("score TCT test is equivalent to comparison of means for gamma = 1", {
  # z-value for a direct comparison of means
  z_mean = as.numeric(mmrm::df_1d(mmrm_fit, c(-1, 0, 0, 0, 1, 0, 0, 0, 0))$t)
  ref_fun = ref_fun_constructor(0:4,
                                coef(mmrm_fit)[c(9, 1:4)],
                                "spline")
  # z-value for TCT score test
  z_score = score_test(time_points = 0:4,
             ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
             exp_estimates = coef(mmrm_fit)[5:8],
             vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
             interpolation = "spline",
             ref_fun = ref_fun,
             j = 1,
             gamma_0 = 1)[1]
  expect_equal(z_score, c("z" = z_mean))
})

test_that("score TCT CI is correct", {
  ref_fun = ref_fun_constructor(0:4,
                                ctrl_estimates,
                                "spline")
  # z-value for TCT score test
  conf_int = score_conf_int(time_points = 0:4,
                       ctrl_estimates = ctrl_estimates,
                       exp_estimates = exp_estimates,
                       vcov = vcov_mmrm,
                       interpolation = "spline",
                       ref_fun = ref_fun,
                       j = 1)
  expect_equal(conf_int, c(-0.28764977, 1.35489195))
})

test_that("score TCT CI does not fail for very wide first CI", {
  ref_fun = ref_fun_constructor(0:4, ctrl_estimates, "spline")
  # Make the SE for the second time point very wide.
  vcov_mmrm_modified = vcov_mmrm
  vcov_mmrm_modified[6, 6] = vcov_mmrm[6, 6] * 50000
  # z-value for TCT score test
  suppressWarnings(
    expect_warning(
      conf_int <- score_conf_int(
        time_points = 0:4,
        ctrl_estimates = ctrl_estimates,
        exp_estimates = exp_estimates,
        vcov = vcov_mmrm_modified,
        interpolation = "spline",
        ref_fun = ref_fun,
        j = 1
      ))
  )

  expect_equal(conf_int, c(-Inf, +Inf))
})

test_that("omnibus score TCT test for common treatment effect is correct", {
  ref_fun = ref_fun_constructor(0:4,
                                ctrl_estimates,
                                "spline")
  # t-value for TCT score test
  t_sq = score_test_common(time_points = 0:4,
                       ctrl_estimates = ctrl_estimates,
                       exp_estimates = exp_estimates,
                       vcov = vcov_mmrm,
                       interpolation = "spline",
                       ref_fun = ref_fun,
                       gamma_0 = 1)[1]
  expect_equal(as.numeric(t_sq), 6.2390581)
})

test_that("omnibus score TCT test for common treatment effect is equivalent to to linear hypothesis test for gamma = 1", {
  ref_fun = ref_fun_constructor(0:4,
                                coef(mmrm_fit)[c(9, 1:4)],
                                "spline")
  contrast_matrix = matrix(c(-1, 0, 0, 0, 1, 0, 0, 0, 0,
                             0, -1, 0, 0, 0, 1, 0, 0, 0,
                             0, 0, -1, 0, 0, 0, 1, 0, 0,
                             0, 0, 0, -1, 0, 0, 0, 1, 0), nrow = 4, byrow = TRUE)
  # Linear hypothesis test, chi-squared test statistic.
  chi_linear_test = as.numeric(mmrm::df_md(mmrm_fit, contrast_matrix)$f_stat) * 4
  # t-value for TCT score test
  t_sq = score_test_common(time_points = 0:4,
                           ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
                           exp_estimates = coef(mmrm_fit)[5:8],
                           vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
                           interpolation = "spline",
                           ref_fun = ref_fun,
                           gamma_0 = 1)[1]
  expect_equal(as.numeric(t_sq), chi_linear_test)
})

test_that("all type of score TCT test for common treatment effect are correct", {
  ref_fun = ref_fun_constructor(0:4,
                                ctrl_estimates,
                                "spline")

  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0
  )
  # z-value for TCT score test
  z_omnibus = score_test_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "omnibus",
    j = 1:4
  )[1]
  z_directional = score_test_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "directional",
    j = 1:4
  )[1]
  z_inv_var = score_test_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "inverse variance",
    j = 1:4
  )[1]
  z_custom = score_test_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "custom",
    j = 1:4,
    weights = c(0.2, 1, 1, 2)
  )[1]
  output_vector = c(z_omnibus, z_directional, z_inv_var, z_custom)
  expect_vector = c(6.2390580536, 0.0023884661, -1.4091152443, -1.7817867581)
  expect_equal(output_vector, expect_vector, ignore_attr = "names")
})

test_that("one-dimensional score TCT tests for common treatment effect are equivalent", {
  ref_fun = ref_fun_constructor(0:4,
                                ctrl_estimates,
                                "spline")

  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0
  )
  # z-value for TCT score test
  z_omnibus = score_test_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "omnibus",
    j = 4
  )[1]
  z_directional = score_test_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "directional",
    j = 4
  )[1]
  z_inv_var = score_test_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "inverse variance",
    j = 4
  )[1]
  z_custom = score_test_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "custom",
    j = 4,
    weights = c(10)
  )[1]
  output_vector = c(z_omnibus, z_directional, z_inv_var ** 2, z_custom ** 2)
  expect_equal(output_vector, rep(z_omnibus, 4), ignore_attr = "names")
})

test_that("all type of mutlivariate score TCT confidence intervals for common treatment effect are correct", {
  ref_fun = ref_fun_constructor(0:4,
                                ctrl_estimates,
                                "spline")

  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0
  )


  # z-value for TCT score test
  conf_int_omnibus = score_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "omnibus",
    j = 1:4,
    gamma_est = 0.8,
  )
  conf_int_directional = score_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "directional",
    j = 1:4,
    gamma_est = 0.8,
  )
  conf_int_inv_var = score_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "inverse variance",
    j = 1:4,
    gamma_est = 0.8,
  )
  conf_int_custom = score_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "custom",
    j = 1:4,
    gamma_est = 0.8,
    weights = c(0, 1, 1, 2)
  )
  output_vector = c(conf_int_omnibus, conf_int_directional, conf_int_inv_var, conf_int_custom)
  expect_vector = c(0.66415999, 1.08546660, 0.36736952, 1.51726499, 0.57410875, 1.08242731, 0.62474983, 1.02097324)
  expect_equal(output_vector, expect_vector)
})

test_that("all type of multivariate score TCT confidence intervals for common treatment effect are equivalent to their univariate counterparts", {
  ref_fun = ref_fun_constructor(0:4,
                                coef(mmrm_fit)[c(9, 1:4)],
                                "spline")

  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0
  )

  # z-value for TCT score test
  conf_int_omnibus = score_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "omnibus",
    j = 2,
    gamma_est = 0.8,
  )
  conf_int_directional = score_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "directional",
    j = 2,
    gamma_est = 0.8,
  )
  conf_int_inv_var = score_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "inverse variance",
    j = 2,
    gamma_est = 0.8,
  )
  conf_int_custom = score_conf_int_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "custom",
    j = 2,
    gamma_est = 0.8,
    weights = 10
  )
  conf_int_univariate = score_conf_int(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    j = 2
  )
  output_vector = c(conf_int_omnibus, conf_int_directional, conf_int_inv_var, conf_int_custom)
  expect_vector = rep(conf_int_univariate, 4)
  expect_equal(output_vector, expect_vector)
})

# Score based estimators
test_that("all type of multivariate score TCT estimators are correct", {
  ref_fun = ref_fun_constructor(0:4,
                                ctrl_estimates,
                                "spline")

  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0
  )

  # z-value for TCT score test
  set.seed(1)
  gamma_omnibus = score_estimate_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "omnibus",
    j = 1:4
  )
  gamma_directional = score_estimate_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "directional",
    j = 1:4
  )
  gamma_inv_var = score_estimate_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "inverse variance",
    j = 1:4
  )
  gamma_custom = score_estimate_common(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "custom",
    j = 1:4,
    weights = c(0, 1, 2, 2)
  )
  output_vector = c(gamma_omnibus, gamma_directional, gamma_inv_var, gamma_custom)
  expect_vector = c(0.85504339, 0.98628002, 0.80792870, 0.78630815)
  expect_equal(output_vector, expect_vector, tolerance = 1e-5)
})

