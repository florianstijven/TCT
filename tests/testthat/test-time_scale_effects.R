test_that("TCT_meta() function works with cubic spline interpolation", {
    set.seed(1)
    TCT_Fit = TCT_meta(
      time_points = 0:4,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      vcov = vcov_mmrm,
      inference = "wald",
      interpolation = "spline",
      B = 1e1
    )
    TCT_output_vctr = c(coef(TCT_Fit),
                        TCT_Fit$bootstrap_estimates[1:2],
                        TCT_Fit$vcov[1, 2])
    check_vctr = c(0.73272539,
                   0.87347422,
                   0.79370231,
                   0.75102487,
                   0.42556428,
                   1.01657371,
                   0.04112038)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() function works with cubic spline interpolation", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_meta_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = TRUE)
  summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.802783786, 0.848841763, 0.890680251, 0.004629814)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta() function works with monoH.FC spline interpolation", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "monoH.FC",
    B = 1e1
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(0.68781306,
                 0.86286059,
                 0.79972841,
                 0.75104144,
                 0.34923210,
                 1.01576094,
                 0.06758570)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() function works with monoH.FC spline interpolation", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
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
  check_vctr = c(0.8129680985, 0.8240954613, 0.8403840360, 0.0039954344)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta() function works with linear spline interpolation", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "linear",
    B = 1e1
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(
    0.274709954,
    0.844347957,
    0.802482431,
    0.751407810,
    0.245249935,
    1.009786220,
    0.064813433
  )
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() function works with linear spline interpolation", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "linear",
    B = 0
  )
  set.seed(1)
  TCT_common_fit = TCT_meta_common(TCT_Fit = TCT_Fit, B = 1e1, bs_fix_vcov = FALSE)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.8471565808, 0.7539046524, 0.8404511275, 0.0038022224)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta() function works with cubic spline interpolation", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    inference = "wald",
    interpolation = "spline",
    B = 1e1
  )
  TCT_output_vctr = c(TCT_Fit$coefficients,
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  check_vctr = c(0.732725393, 0.873474215, 0.793702312, 0.751024869, 0.425564278, 1.016573713, 0.041120375)
  expect_equal(TCT_output_vctr, check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta() and its summary work with cubic spline interpolation and score-based inference", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0,
    inference = "score"
  )
  TCT_Fit_summary = summary(TCT_Fit)
  TCT_output_vctr = c(TCT_Fit$coefficients[1:2],
                      TCT_Fit_summary$ci_matrix[1, 1],
                      TCT_Fit_summary$ci_matrix[2, 2])
  check_vctr = c(0.73272539, 0.87347422,-0.28764977, 1.07517890)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() and its summary work with cubic spline interpolation and score-based inference", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
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
  check_vctr = c(0.79370276, 0.64340024, 0.97939229)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common(inference = score) can be combine with TCT_meta(inference = wald)", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
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
  check_vctr = c(0.79370276, 0.64340024, 0.97939229)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() and its summary work with the nonlinear GLS estimator", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
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
    inference = "least-squares"
  )
  TCT_common_summary = summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_summary$gamma_common_ci[1, 1:2])
  check_vctr = c(0.85504336, 0.70072621, 1.02848644)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common() works with subset of estimates in exp_estimates", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0,
    inference = "wald"
  )
  TCT_Fit_summary = summary(TCT_Fit)
  set.seed(1)
  # non-linear least squares method
  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    inference = "least-squares",
    select_coef = 3:4,
    B = 30
  )
  TCT_common_summary = summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_summary$gamma_common_ci[1, 1:2],
                      TCT_common_summary$se_bootstrap)
  check_vctr = c(0.83071747, 0.64949799, 0.99620621, 0.088907886)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})
