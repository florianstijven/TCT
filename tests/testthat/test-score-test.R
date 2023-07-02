test_that("score TCT test is equivalent to comparison of means for gamma = 1", {
  library(dplyr)
  # Example data set transformed to format required by TCT()
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
             gamma_0 = 1)
  expect_equal(z_score, z_mean)
})

test_that("score TCT CI is correct", {
  library(dplyr)
  # Example data set transformed to format required by TCT()
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
  # Example data set transformed to format required by TCT()
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
                       gamma_0 = 1)
  expect_equal(as.numeric(t_sq), 6.2673327)
})

test_that("omnibus score TCT test for common treatment effect is equivalent to to linear hypothesis test for gamma = 1", {
  library(dplyr)
  # Example data set transformed to format required by TCT()
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
                           gamma_0 = 1)
  expect_equal(as.numeric(t_sq), chi_linear_test)
})

test_that("all type of score TCT test for common treatment effect are correct", {
  library(dplyr)
  # Example data set transformed to format required by TCT()
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

  TCT_Fit = TCT(
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
  )
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
  )
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
  )
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
  )
  output_vector = c(z_omnibus, z_directional, z_inv_var, z_custom)
  expect_vector = c(6.26733272, -0.04865092, -1.41203381, -1.78553411)
  expect_equal(output_vector, expect_vector)
})

test_that("one-dimensional score TCT tests for common treatment effect are equivalent", {
  library(dplyr)
  # Example data set transformed to format required by TCT()
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

  TCT_Fit = TCT(
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
  )
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
  )
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
  )
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
  )
  output_vector = c(z_omnibus, z_directional ** 2, z_inv_var ** 2, z_custom ** 2)
  expect_equal(output_vector, rep(z_omnibus, 4))
})

# Score based estimates

# Score based CI for common gamma
