test_that("nonlinear_gls_estimator() works for estimating the common acceleration factor.", {
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
  # Run the nonlinear GLS estimator.
  nl_gls = nonlinear_gls_estimator(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    ref_fun = ref_fun
  )
  check_object = list(
    estimates = c(22.2663151, 22.7008012, 25.3371983, 28.4786173, 30.3118762, 0.7805827),
    criterion = 47.131846
  )
  expect_equal(nl_gls, check_object, ignore_attr = "names")
})
