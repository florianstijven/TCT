test_that("Optimisation for weights vector works", {
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
