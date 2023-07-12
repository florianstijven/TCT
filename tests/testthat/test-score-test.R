test_that("score TCT test is equivalent to comparison of means for gamma = 1", {
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
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
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
  # z-value for a direct comparison of means
  z_mean = as.numeric(mmrm::df_1d(mmrm_fit, c(-1, 0, 0, 0, 1, 0, 0, 0, 0))$t)
  ref_fun = ref_fun_constructor(0:4,
                                coef(mmrm_fit)[c(9, 1:4)],
                                "spline")
  # z-value for TCT score test
  conf_int = score_conf_int(time_points = 0:4,
                       ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
                       exp_estimates = coef(mmrm_fit)[5:8],
                       vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
                       interpolation = "spline",
                       ref_fun = ref_fun,
                       j = 1)
  expect_equal(conf_int, c(-0.2869702204, 1.3540556235))
})

test_that("omnibus score TCT test for common treatment effect is correct", {
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
  ref_fun = ref_fun_constructor(0:4,
                                coef(mmrm_fit)[c(9, 1:4)],
                                "spline")
  # t-value for TCT score test
  t_sq = score_test_common(time_points = 0:4,
                       ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
                       exp_estimates = coef(mmrm_fit)[5:8],
                       vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
                       interpolation = "spline",
                       ref_fun = ref_fun,
                       gamma_0 = 1)[1]
  expect_equal(as.numeric(t_sq), 6.2673327)
})

test_that("omnibus score TCT test for common treatment effect is equivalent to to linear hypothesis test for gamma = 1", {
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
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
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
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
  z_omnibus = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "omnibus",
    j = 1:4
  )[1]
  z_directional = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "directional",
    j = 1:4
  )[1]
  z_inv_var = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "inverse variance",
    j = 1:4
  )[1]
  z_custom = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "custom",
    j = 1:4,
    weights = c(0.2, 1, 1, 2)
  )[1]
  output_vector = c(z_omnibus, z_directional, z_inv_var, z_custom)
  expect_vector = c(6.26733272, 0.00236691, -1.41203381, -1.78553411)
  expect_equal(output_vector, expect_vector, ignore_attr = "names")
})

test_that("one-dimensional score TCT tests for common treatment effect are equivalent", {
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
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
  z_omnibus = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "omnibus",
    j = 4
  )[1]
  z_directional = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "directional",
    j = 4
  )[1]
  z_inv_var = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "inverse variance",
    j = 4
  )[1]
  z_custom = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
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
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
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
    j = 1:4,
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
    j = 1:4,
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
    j = 1:4,
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
    j = 1:4,
    gamma_est = 0.8,
    weights = c(0, 1, 1, 2)
  )
  output_vector = c(conf_int_omnibus, conf_int_directional, conf_int_inv_var, conf_int_custom)
  expect_vector = c(0.66475028, 1.08446326, 0.36763236, 1.5162254, 0.57461092, 1.08179608, 0.62504568, 1.02042558)
  expect_equal(output_vector, expect_vector)
})

test_that("all type of multivariate score TCT confidence intervals for common treatment effect are equivalent to their univariate counterparts", {
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
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
  library(dplyr)
  # Example data set transformed to format required by TCT_meta()
  data_test = simulated_test_trial %>%
    mutate(
      time_int = (Week %/% 25) + 1,
      arm_time = ifelse(time_int == 1L,
                        "baseline",
                        paste0(arm, ":", time_int))
    )
  # Fit a MMRM model to the data set. The parameter estimates of this model form
  # the basis to perform the time component tests.
  mmrm_fit = analyze_mmrm(data_test)
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
  set.seed(1)
  gamma_omnibus = score_estimate_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "omnibus",
    j = 1:4
  )
  gamma_directional = score_estimate_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "directional",
    j = 1:4
  )
  gamma_inv_var = score_estimate_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "inverse variance",
    j = 1:4
  )
  gamma_custom = score_estimate_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    type = "custom",
    j = 1:4,
    weights = c(0, 1, 2, 2)
  )
  output_vector = c(gamma_omnibus, gamma_directional, gamma_inv_var, gamma_custom)
  expect_vector = c(0.85504244, 0.98637208, 0.80793279, 0.78630034)
  expect_equal(output_vector, expect_vector, tolerance = 1e-5)
})


# v_function(
#   time_points = 0:4,
#   ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#   exp_estimates = coef(mmrm_fit)[5:8],
#   vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#   interpolation = "spline",
#   ref_fun = ref_fun,
#   gamma_0 = 0.78630034,
#   j = 1:4,
#   weights = c(0, 1, 2, 2)
# )
#
# gr_gamma_w = gradient_gamma_w(
#   time_points = 0:4,
#   ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#   exp_estimates = coef(mmrm_fit)[5:8],
#   vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#   interpolation = "spline",
#   gamma_0 = 0.8743958,
#   j = 1:4,
#   weights = c(w_opt, 1 - sum(w_opt))
#
#
# gr_gamma_w = gradient_gamma_w(
#   time_points = 0:4,
#   ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#   exp_estimates = coef(mmrm_fit)[5:8],
#   vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#   interpolation = "spline",
#   gamma_0 = 0.8743958,
#   j = 1:4,
#   weights = c(0, 0, 1, 0))
#
#
#
# matrix(gr_gamma_w, nrow = 1) %*%
#   vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)] %*%
#   matrix(gr_gamma_w, ncol = 1)
#
#
# score_estimate_common(
#   time_points = 0:4,
#   ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#   exp_estimates = coef(mmrm_fit)[5:8],
#   vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#   interpolation = "spline",
#   ref_fun = ref_fun,
#   type = "custom",
#   j = 1:4,
#   weights = c(w_opt, 1 - sum(w_opt))
# )
#
# w_opt = optimize_weights(
#   time_points = 0:4,
#   ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#   exp_estimates = coef(mmrm_fit)[5:8],
#   vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#   interpolation = "spline",
#   ref_fun = ref_fun,
#   weights = c(0.2, 1, 3, 1),
#   j = 1:4
# )
