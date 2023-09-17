library(dplyr)
# Example data set transformed to format required by TCT_meta()
data_test = simulated_test_trial %>%
  dplyr::mutate(
    time_int = (Week %/% 25) + 1,
    arm_time = ifelse(time_int == 1L,
                      "baseline",
                      paste0(arm, ":", time_int))
  )
# Fit a MMRM model to the data set. The parameter estimates of this model form
# the basis to perform the time component tests.
mmrm_fit = analyze_mmrm(data_test)
# The mmrm_fit object is used throughout the tests.
