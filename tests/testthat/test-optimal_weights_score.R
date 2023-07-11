test_that("Optimisation for weights vector works", {
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
  weights_opt1 = optimize_weights(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    ref_fun = ref_fun,
    interpolation = "spline",
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    j = 1:4
  )
  weights_opt2 = optimize_weights(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    ref_fun = ref_fun,
    interpolation = "spline",
    vcov = diag(1, nrow = 9),
    j = 1:4
  )
  output_vector = c(weights_opt1, weights_opt2)
  check_vector = c(
    0,
    0,
    1,
    0,
    0.0240284838615804,
    0.238064736941874,
    0.385870163560325,
    0.35203661563622
  )
  expect_equal(output_vector,
               check_vector)
})
