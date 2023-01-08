test_that("TCT() function works", {
  data = simulated_test_trial %>%
    mutate(time_int = (Week %/% 25)) %>%
    arrange(trial_number, SubjId, time_int) %>%
    mutate(time_int = as.integer(time_int) + 1L) %>%
    mutate(arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)))
    mmrm_fit = analyze_mmrm(data)
    set.seed(1)
    TCT_Fit = TCT(
      time_points = 0:4,
      ctrl_estimates = x$coefficients[c(9, 1:4)],
      exp_estimates = x$coefficients[5:8],
      vcov = x$varBeta[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
      interpolation = "spline",
      B = 1e3
    )
    TCT_output_vctr = c(TCT_Fit$coefficients,
                        TCT_Fit$bootstrap_estimates[1:2],
                        TCT_Fit$vcov[1, 2])
    check_vctr = c(0.7327657260, 0.8734756221, 0.7936971718, 0.7510185731,
                   0.4256797100, 1.0165886470, 0.0411074803)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_common() function works", {
  data = simulated_test_trial %>%
    mutate(time_int = (Week %/% 25)) %>%
    arrange(trial_number, SubjId, time_int) %>%
    mutate(time_int = as.integer(time_int) + 1L) %>%
    mutate(arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = x$coefficients[c(9, 1:4)],
    exp_estimates = x$coefficients[5:8],
    vcov = x$varBeta[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = TRUE)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.8027742481, 0.8488277428, 0.8906516826, 0.0046274862)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})


