test_that("TCT() function works with cubic spline interpolation", {
  library(mmrm)
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)))
    mmrm_fit = analyze_mmrm(data)
    set.seed(1)
    TCT_Fit = TCT(
      time_points = 0:4,
      ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
      exp_estimates = coef(mmrm_fit)[5:8],
      vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
      inference = "wald",
      interpolation = "spline",
      B = 1e3
    )
    TCT_output_vctr = c(TCT_Fit$coefficients,
                        TCT_Fit$bootstrap_estimates[1:2],
                        TCT_Fit$vcov[1, 2])
    check_vctr = c(0.7327657260, 0.8734756221, 0.7936971718, 0.7510185731,
                   0.4256797100, 1.0165886470, 0.0411074803)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names", tolerance = 1e-3)
})

test_that("TCT_common() function works with cubic spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = TRUE)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.8027742481, 0.8488277428, 0.8906516826, 0.0046274862)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names", tolerance = 1e-3)
})

test_that("TCT() function works with monoH.FC spline interpolation", {
  library(mmrm)
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "monoH.FC",
    B = 1e3
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(0.68776822, 0.86289314, 0.79971832, 0.75102503,
                 0.34939350, 1.01529410, 0.06732386)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names", tolerance = 1e-3)
})

test_that("TCT_common() function works with monoH.FC spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "monoH.FC",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_common(
    TCT_Fit = TCT_Fit,
    B = 1e1,
    bs_fix_vcov = FALSE,
    null_bs = TRUE
  )
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov,
                      TCT_common_fit$bootstrap_estimates_null$estimates_bootstrap[1:2])
  check_vctr = c(0.812935193, 0.824029215, 0.840300615, 0.003978323, 0.775201628, 1.213386432)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names", tolerance = 1e-5)
})

test_that("TCT() function works with linear spline interpolation", {
  library(mmrm)
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "linear",
    B = 1e3
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(0.2746882, 0.8443806, 0.8024704, 0.7513858,
                 0.2452254, 1.0094905, 0.3594113)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names", tolerance = 1e-3)
})

test_that("TCT_common() function works with linear spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "linear",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = FALSE)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.8144168, 0.8091335, 0.8858162, 0.2172113)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names", tolerance = 1e-5)
})


test_that("TCT() function works with fourPL interpolation", {
  library(mmrm)
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "fourPL",
    B = 1e3
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(0.4888728, 0.8896603, 0.8004555, 0.7428344,
                 0.4196664, 0.9506719, 0.1495659)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names", tolerance = 1e-3)
})

test_that("TCT_common() function works with fourPL interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "fourPL",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = FALSE)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.77765914, 0.78482953, 1.01569215, 0.08541886)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names", tolerance = 1e-5)
})


