test_that("TCT_meta() function works with cubic spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)))
    mmrm_fit = analyze_mmrm(data)
    set.seed(1)
    TCT_Fit = TCT_meta(
      time_points = 0:4,
      ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
      exp_estimates = coef(mmrm_fit)[5:8],
      vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
      inference = "wald",
      interpolation = "spline",
      B = 1e1
    )
    TCT_output_vctr = c(coef(TCT_Fit),
                        TCT_Fit$bootstrap_estimates[1:2],
                        TCT_Fit$vcov[1, 2])
    check_vctr = c(0.73268086,
                   0.87350717,
                   0.79369280,
                   0.75100870,
                   0.42588694,
                   1.01608204,
                   0.04095813)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() function works with cubic spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_meta_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = TRUE)
  summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.802753988, 0.848726277, 0.890420727, 0.004609946)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta() function works with monoH.FC spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "monoH.FC",
    B = 1e1
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(0.68776822,
                 0.86289314,
                 0.79971832,
                 0.75102503,
                 0.34939350,
                 1.01529410,
                 0.06732386)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() function works with monoH.FC spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "monoH.FC",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    B = 1e1,
    bs_fix_vcov = FALSE
  )
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.812935193, 0.824029215, 0.840300615, 0.003978323)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta() function works with linear spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "linear",
    B = 1e1
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(
    0.274688245,
    0.844380586,
    0.802470387,
    0.751385771,
    0.245225389,
    1.009490511,
    0.064547695
  )
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() function works with linear spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "linear",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_meta_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = FALSE)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.847130645, 0.753917984, 0.840351151, 0.003785658)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

# test_that("TCT_meta() function works with fourPL interpolation", {
#   library(mmrm)
#   data = simulated_test_trial %>%
#     dplyr::mutate(time_int = (Week %/% 25)) %>%
#     dplyr::arrange(trial_number, SubjId, time_int) %>%
#     dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
#     dplyr::mutate(arm_time = ifelse(time_int == 1L,
#                                     "baseline",
#                                     paste0(arm, ":", time_int)))
#   mmrm_fit = analyze_mmrm(data)
#   set.seed(1)
#   TCT_Fit = TCT_meta(
#     time_points = 0:4,
#     ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#     exp_estimates = coef(mmrm_fit)[5:8],
#     vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#     interpolation = "fourPL",
#     B = 1e3
#   )
#   TCT_output_vctr = c(TCT_Fit$coefficients,
#                       TCT_Fit$bootstrap_estimates[1:2],
#                       TCT_Fit$vcov[1, 2])
#   check_vctr = c(0.4888728, 0.8896603, 0.8004555, 0.7428344,
#                  0.4196664, 0.9506719, 0.1495659)
#   expect_equal(TCT_output_vctr, check_vctr,
#                ignore_attr = "names", tolerance = 1e-3)
# })

# test_that("TCT_meta_common() function works with fourPL interpolation", {
#   data = simulated_test_trial %>%
#     dplyr::mutate(time_int = (Week %/% 25)) %>%
#     dplyr::arrange(trial_number, SubjId, time_int) %>%
#     dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
#     dplyr::mutate(arm_time = ifelse(time_int == 1L,
#                                     "baseline",
#                                     paste0(arm, ":", time_int)))
#   mmrm_fit = analyze_mmrm(data)
#   set.seed(1)
#   TCT_Fit = TCT_meta(
#     time_points = 0:4,
#     ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#     exp_estimates = coef(mmrm_fit)[5:8],
#     vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#     interpolation = "fourPL",
#     B = 0
#   )
#   set.seed(1)
#   TCT_common_fit = TCT_meta_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = FALSE)
#   TCT_output_vctr = c(TCT_common_fit$coefficients,
#                       TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
#                       TCT_common_fit$vcov)
#   check_vctr = c(0.77765914, 0.78482953, 1.01569215, 0.08541886)
#   expect_equal(TCT_output_vctr, check_vctr,
#                ignore_attr = "names", tolerance = 1e-5)
# })

test_that("TCT_meta() function works with cubic spline interpolation", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    inference = "wald",
    interpolation = "spline",
    B = 1e1
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(0.73268086,
                 0.87350717,
                 0.79369280,
                 0.75100870,
                 0.42588694,
                 1.01608204,
                 0.04095813)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta() and its summary work with cubic spline interpolation and score-based inference", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0,
    inference = "score"
  )
  TCT_Fit_summary = summary(TCT_Fit)
  TCT_output_vctr = c(TCT_Fit$coefficients[1:2],
                      TCT_Fit_summary$ci_matrix[1, 1],
                      TCT_Fit_summary$ci_matrix[2, 2])
  check_vctr = c(
    0.732680858024503,
    0.873507167779775,
    -0.286970220428998,
    1.07477707220273
  )
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() and its summary work with cubic spline interpolation and score-based inference", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0,
    inference = "score"
  )
  TCT_Fit_summary = summary(TCT_Fit)
  set.seed(1)
  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    B = 10,
    bs_fix_vcov = TRUE,
    inference = "score",
    type = "custom"
  )
  TCT_common_summary = summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_summary$gamma_common_ci[1, 1:2])
  check_vctr = c(0.793693254456409,
                 0.643696602719346,
                 0.978881018916929)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common(inference = score) can be combine with TCT_meta(inference = wald)", {
  data = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25)) %>%
    dplyr::arrange(trial_number, SubjId, time_int) %>%
    dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
    dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data)
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0,
    inference = "wald"
  )
  TCT_Fit_summary = summary(TCT_Fit)
  set.seed(1)
  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    B = 10,
    bs_fix_vcov = TRUE,
    inference = "score",
    type = "custom"
  )
  TCT_common_summary = summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_summary$gamma_common_ci[1, 1:2])
  check_vctr = c(0.793693254456409,
                 0.643696602719346,
                 0.978881018916929)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})


