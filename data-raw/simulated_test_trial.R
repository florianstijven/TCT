## code to prepare `simulated_test_trial` dataset goes here
set.seed(1)
simulated_test_trial = simulate_single_trial(
  control = "placebo",
  delta = 1.3,
  total_time = 100,
  n_measurements = 5,
  sample_size = 500,
  trial_number = 1,
  details = FALSE
)

usethis::use_data(simulated_test_trial, overwrite = TRUE)
