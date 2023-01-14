test_that("mmrm is correctly fitted for simulated_test_trial with ML", {
  data_trial = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25) + 1,
                  arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data_trial)
  expect_equal(c(coef(mmrm_fit)[1],
                 mmrm_fit$neg_log_lik),
               c(22.94564411,
                 7955.777072), ignore_attr = "names")
})

test_that("mmrm is correctly fitted for simulated_test_trial with REML", {
  data_trial = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25) + 1,
                  arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data_trial, method = "REML")
  expect_equal(c(coef(mmrm_fit)[1],
                 mmrm_fit$neg_log_lik),
               c(22.94563423 ,
                 7957.205066), ignore_attr = "names")
})

test_that("mmrm (null model) is correctly fitted for simulated_test_trial with ML", {
  data_trial = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25) + 1,
                  arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data_trial, type = "null")
  expect_equal(c(coef(mmrm_fit)[1],
                 mmrm_fit$neg_log_lik),
               c(22.544,
                 7958.8597402), ignore_attr = "names")
})

test_that("mmrm (null model) is correctly fitted for simulated_test_trial with REML", {
  data_trial = simulated_test_trial %>%
    dplyr::mutate(time_int = (Week %/% 25) + 1,
                  arm_time = ifelse(time_int == 1L,
                                    "baseline",
                                    paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm(data_trial, method = "REML", type = "null")
  expect_equal(c(coef(mmrm_fit)[1],
                 mmrm_fit$neg_log_lik),
               c(22.544 ,
                 7960.975908), ignore_attr = "names")
})
