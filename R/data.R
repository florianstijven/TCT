#' Simulated clinical trial data
#'
#' Data that was simulated with the [simulate_single_trial()] function with the
#' following argument:
#' * `control = "placebo"`
#' * `delta = 1.3`
#' * `total_time = 100`
#' * `n_measurements = 5`
#' * `sample_size = 1000`
#' * `trial_number = 1`
#' * `details = FALSE`
#' @details
#' This is corresponds to a trial with 5 equally spaced measurements per patient
#' over a period of 100 weeks. The total sample size is 1000, and there are
#' patient groups: (i) placebo and (ii) disease-modifying treatment. The
#' treatment effect corresponds to a deceleration factor of 1.3.
#'
#' @format ## `simulated_test_trial` A data frame with 5,000 rows and 9 columns:
#' \describe{
#'   \item{Week}{The number of weeks since randomization the measurement is taken.}
#'   \item{SubjId}{Integer that uniquely identifies patients}
#'   \item{Apo}{Number of ApoE4 alleles}
#'   \item{Gender}{Gender of the patient}
#'   \item{ADAScog_integer}{Measured score on the ADAScog scale}
#'   \item{Age_integer}{Age at time of randomization}
#'   \item{BMMSE}{Measured Mini Mental State Examination score at time of randomization}
#'   \item{arm}{Integer that indicates treatment arm: 0 for control, 1 for active treatment}
#'   \item{trial_number}{Integer that uniquely identifies trials. Not relevant since there is data for only a single trial.}
#' }
"simulated_test_trial"
