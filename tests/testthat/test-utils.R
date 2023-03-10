testthat::test_that("time mapping with cubic splines works",
                    {
                      x_ref = 0:5
                      y_ref = x_ref**3
                      y_obs = c(y_ref, 250, -125)
                      testthat::expect_equal(
                        get_new_time(y_ref, x_ref, y_obs, method = "spline"),
                        c(x_ref, c(10, -5))
                      )
                    })

testthat::test_that("time mapping with linear interpolation works",
                    {
                      x_ref = 0:5
                      y_ref = x_ref**3
                      y_obs = c(y_ref, 250, -125)
                      testthat::expect_equal(
                        get_new_time(y_ref, x_ref, y_obs, method = "linear"),
                        c(x_ref, c(10, -5))
                      )
                    })

testthat::test_that("time mapping with monoH.FC interpolation works",
                    {
                      x_ref = 0:5
                      y_ref = x_ref**3
                      y_obs = c(y_ref, 250, -125)
                      testthat::expect_equal(
                        get_new_time(y_ref, x_ref, y_obs, method = "monoH.FC"),
                        c(x_ref, c(10, -5))
                      )
                    })
