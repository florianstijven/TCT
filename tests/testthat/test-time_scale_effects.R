test_that("TCT_meta() function works with cubic spline interpolation", {
    set.seed(1)
    TCT_Fit = TCT_meta(
      time_points = 0:4,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      vcov = vcov_mmrm,
      inference = "delta-method",
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
  TCT_common_fit = TCT_meta_common(TCT_Fit = TCT_Fit, B = 1e1, inference = "delta-method")
  summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_fit$bootstrap_estimates$estimates_bootstrap[1:2],
                      TCT_common_fit$vcov)
  check_vctr = c(0.802783786, 0.80783893, 0.85520306, 0.004629814)
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
    inference = "delta-method"
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
  TCT_common_fit = TCT_meta_common(TCT_Fit = TCT_Fit, B = 1e1, inference = "delta-method")
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
    inference = "delta-method",
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

test_that("TCT_meta() and TCT_meta_common() function works with cubic spline interpolation and known origin mean", {
  set.seed(1)
  # Set first row and column of vcov_mmrm to zeroes, reflecting the fact that
  # the mean at the origin (i.e., time zero) is known. This reflects, for
  # instance, the setting where the endpoint is change from baseline.
  vcov_mmrm_new = vcov_mmrm
  vcov_mmrm_new[1, ] = 0
  vcov_mmrm_new[, 1] = 0

  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm_new,
    inference = "contrast",
    interpolation = "spline",
    B = 1e1
  )
  TCT_Fit_summary = summary(TCT_Fit)
  # Expect warning if the bounds for the search interval for contrast-based CIs are bad.
  suppressWarnings(
    expect_warning(
      summary(TCT_Fit, bounds = c(5, 6))
    )
  )

  # Expect no warning if the bounds are good.
  expect_no_warning(
    TCT_Fit_summary <- summary(TCT_Fit, bounds = c(-1, 2))
  )

  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    B = 10,
    inference = "least-squares"
  )
  TCT_common_summary = summary(TCT_common_fit)

  # Extract the computed values.
  TCT_output_vctr = c(coef(TCT_Fit),
                      TCT_Fit$bootstrap_estimates[1:2],
                      TCT_Fit$vcov[1, 2])
  TCT_common_output_vctr = c(TCT_common_fit$coefficients,
                             TCT_common_summary$gamma_common_ci[1, 1:2])
  # Set expected values.
  check_vctr_TCT_meta = c(0.73272539, 0.87347422, 0.79370231, 0.75102487, 0.43709139, 1.01568520, 0.04100701)
  check_vctr_TCT_common = c(0.85419054, 0.69587063, 1.02827121)
  # Check that the computed values match the expected values.
  expect_equal(TCT_output_vctr, check_vctr_TCT_meta,
               ignore_attr = "names")
  expect_equal(TCT_common_output_vctr, check_vctr_TCT_common,
               ignore_attr = "names")
})


test_that("TCT_meta() and TCT_meta_common() yield similar results with known and quasi-known parameter", {
  set.seed(1)
  # Set first row and column of vcov_mmrm to zeroes, reflecting the fact that
  # the mean at the origin (i.e., time zero) is known. This reflects, for
  # instance, the setting where the endpoint is change from baseline.
  vcov_mmrm_known = vcov_mmrm
  vcov_mmrm_known[1, ] = 0
  vcov_mmrm_known[, 1] = 0

  # Set first row and column of `vcov_mmrm` to very small number. This ensures
  # that the standard algorithm is run, but with a quasi-known parameter (i.e.,
  # the first time point). The result should be very similar when the
  # corresponding variance is set to zero.
  vcov_mmrm_quasi_known = vcov_mmrm
  vcov_mmrm_quasi_known[1, ] = 1e-8 # Set the first row to a small value.
  vcov_mmrm_quasi_known[, 1] = 1e-8 # Set the first column to a small value.

  TCT_Fit_known = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm_known,
    inference = "contrast",
    interpolation = "spline",
    B = 1e1
  )
  TCT_Fit_known_summary = summary(TCT_Fit_known)
  TCT_common_fit_known = TCT_meta_common(
    TCT_Fit = TCT_Fit_known,
    B = 10,
    inference = "least-squares"
  )
  TCT_common_summary_known = summary(TCT_common_fit_known)

  TCT_Fit_quasi_known = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm_quasi_known,
    inference = "contrast",
    interpolation = "spline",
    B = 1e1
  )
  TCT_Fit_quasi_known_summary = summary(TCT_Fit_quasi_known)
  TCT_common_fit_quasi_known = TCT_meta_common(
    TCT_Fit = TCT_Fit_quasi_known,
    B = 10,
    inference = "least-squares"
  )
  TCT_common_summary_quasi_known = summary(TCT_common_fit_quasi_known)


  # Check whether the results from TCT_meta() and TCT_meta_common() are similar
  # for both known and quasi-known parameter settings.
  expect_equal(
    c(TCT_Fit_known$coefficients,
      TCT_Fit_known_summary$ci_matrix[1, 1],
      TCT_Fit_known_summary$ci_matrix[2, 2]),
    c(TCT_Fit_quasi_known$coefficients,
      TCT_Fit_quasi_known_summary$ci_matrix[1, 1],
      TCT_Fit_quasi_known_summary$ci_matrix[2, 2])
  )
  expect_equal(
    c(TCT_common_fit_known$coefficients,
      TCT_common_summary_known$gamma_common_ci[1, 1:2]),
    c(TCT_common_fit_quasi_known$coefficients,
      TCT_common_summary_quasi_known$gamma_common_ci[1, 1:2])
  )
})

test_that("TCT_meta() and its summary work with cubic spline interpolation and contrast-based inference", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0,
    inference = "contrast"
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

test_that("TCT_meta_common() and its summary work with cubic spline interpolation and contrast-based inference", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0,
    inference = "contrast"
  )
  TCT_Fit_summary = summary(TCT_Fit)
  set.seed(1)
  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    B = 10,
    inference = "contrast"
  )
  TCT_common_summary = summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_summary$gamma_common_ci[1, 1:2])
  check_vctr = c(0.79370279, 0.64340024, 0.97939229)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta_common(inference = contrast) can be combine with TCT_meta(inference = delta-method)", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 0,
    inference = "delta-method"
  )
  TCT_Fit_summary = summary(TCT_Fit)
  set.seed(1)
  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    B = 10,
    inference = "contrast"
  )
  TCT_common_summary = summary(TCT_common_fit)
  TCT_output_vctr = c(TCT_common_fit$coefficients,
                      TCT_common_summary$gamma_common_ci[1, 1:2])
  check_vctr = c(0.79370279, 0.64340024, 0.97939229)
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
    inference = "delta-method"
  )
  TCT_Fit_summary = summary(TCT_Fit)
  set.seed(1)
  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    B = 10,
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
    inference = "delta-method"
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
  check_vctr = c(0.830717202, 0.649497989, 0.996206209, 0.088907909)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})

test_that("TCT_meta() and TCT_meta_common() work with BC percentile interval", {
  set.seed(1)
  TCT_Fit = TCT_meta(
    time_points = 0:4,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov_mmrm,
    interpolation = "spline",
    B = 50,
    inference = "delta-method"
  )
  TCT_Fit_summary = summary(TCT_Fit, bootstrap_type = "BC percentile")
  TCT_common_fit = TCT_meta_common(
    TCT_Fit = TCT_Fit,
    inference = "least-squares",
    B = 50
  )
  TCT_common_summary = summary(TCT_common_fit, bootstrap_type = "BC percentile")
  TCT_output_vctr = c(TCT_Fit_summary$ci_bootstrap[1, 1:2],
                      TCT_common_summary$ci_bootstrap)
  check_vctr = c( -0.118015572370, 1.377234650547, 0.712604553192, 1.275732088628)
  expect_equal(TCT_output_vctr,
               check_vctr,
               ignore_attr = "names")
})
