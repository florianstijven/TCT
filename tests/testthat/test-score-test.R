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

# WIP
test_that("directional score TCT test for common treatment effect is correct", {
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
  summary.TCT(TCT_Fit)
  set.seed(1)
  TCT_common_fit = TCT_common(TCT_Fit = TCT_Fit, B = 1e3, bs_fix_vcov = FALSE)
  summary(TCT_common_fit)
  # z-value for TCT score test
  z = score_test_common(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun,
    gamma_0 = 1,
    type = "custom",
    j = 1:4,
    weights = c(0, 0.5, 1, 1)
  )
  expect_equal(z, 6.2673327)
})
