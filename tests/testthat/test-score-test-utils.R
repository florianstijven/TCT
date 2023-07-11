test_that("score_vector_jacobian() works as expected for gamma_0 = 1", {
  J =  score_vector_jacobian(
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
