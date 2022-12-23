test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

testthat::test_that("time mapping with cubic splines works",
                    {
                      x_ref = 0:5
                      y_ref = x_ref**3
                      y_obs = c(y_ref,
                                stats::spline(x_ref,
                                              y_ref,
                                              method = "natural",
                                              xout = c(-50, 23.5, 50))$y
                                )
                      testthat::expect_equal(
                        get_new_time(y_ref, x_ref, y_obs, method = "spline"),
                        c(x_ref, c(-50, 23.5, 50))
                      )
                    })

testthat::test_that("time mapping with linear interpolation works",
                    {
                      x_ref = 1:5
                      y_ref = x_ref**3
                      y_obs = c(y_ref, 186, -6)
                      testthat::expect_equal(
                        get_new_time(y_ref, x_ref, y_obs, method = "linear"),
                        c(x_ref, c(6, 0))
                      )
                    })
