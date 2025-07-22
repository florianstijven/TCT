test_that("contrast_vector_jacobian() works as expected for gamma_0 = 1", {
  J =  contrast_vector_jacobian(
    time_points = 0:4,
    ctrl_estimates = c(22, 22.9, 25.6, 28.7, 30.1),
    interpolation = "spline",
    gamma_0 = 1,
    j = 1:4
  )
  expect_equal(
    J[1:4, 2:5],
    diag(-1, nrow = 4)
  )
})

test_that("contrast_vector_jacobian() works as expected for gamma_0 != 1", {
  J = contrast_vector_jacobian(
    time_points = 0:4,
    ctrl_estimates = c(22, 22.9, 25.6, 28.7, 30.1),
    interpolation = "spline",
    gamma_0 = .55,
    j = 1:4
  )
  expect_equal(
    J[1, 3],
    0.164410714060068
  )
})

test_that("v_deriv_alpha_beta() works as expected for gamma_0 != 1", {
  v_deriv = v_deriv_alpha_beta(
    time_points = 0:4,
    ctrl_estimates = c(22, 22.9, 25.6, 28.7, 30.1),
    interpolation = "spline",
    gamma_0 = .55,
    j = 1:4,
    weights = 1:4
  )
  expect_equal(
    v_deriv[c(2, 9)],
    c(-3.42105356366353, 4)
  )
})

test_that("f0_gradient_t_jacobian() works as expected", {
  J = f0_gradient_t_jacobian(
    t_vec = -1:6,
    ctrl_estimates = c(22, 22.9, 25.6, 28.7, 30.1),
    exp_estimates = 1:4,
    interpolation = "spline",
    time_points = 0:4
  )
  expect_equal(
    J[5, 3],
    -0.857142854947597
  )
})

test_that("f0_gradient_t_jacobian_numerical() and f0_gradient_t_jacobian_analytical() are equivalent", {
  grad_analytical = gradient_gamma_w_analytical(
    time_points = 0:4,
    ctrl_estimates = c(22, 22.9, 25.6, 28.7, 30.1),
    exp_estimates = c(22.9, 25.6, 28.7, 30.1) - 0.5,
    interpolation = "spline",
    gamma_0 = 0.55,
    j = 1:4,
    weights = 1:4
  )

  grad_numerical = gradient_gamma_w_numerical(
    time_points = 0:4,
    ctrl_estimates = c(22, 22.9, 25.6, 28.7, 30.1),
    exp_estimates = c(22.9, 25.6, 28.7, 30.1) - 0.5,
    interpolation = "spline",
    gamma_0 = 0.55,
    j = 1:4,
    weights = 1:4
  )

  expect_equal(
    grad_analytical,
    grad_numerical,
    tolerance = 1e-5
  )
})


