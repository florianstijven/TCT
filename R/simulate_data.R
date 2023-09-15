#' Simulates N clinical trials under the supplied model parameters.
#'
#' @param control Type of control group, can be:
#' * `"placebo"`: Control group receives placebo
#' * `"symptomatic"`: control group receives a symptomatic treatment
#' @param delta deceleration factor. This parameterizes the disease-modifying
#'   treatment effect.
#' @param total_time Total trial duration in weeks.
#' @param n_measurements Number of equally spaced measurement occasions for each
#'   patient.
#' @param sample_size Total sample size in each trial.
#' @param N_trials Total number of trial to simulate independently.
#' @param details Return detailed data frame:
#' * `TRUE`: the returned data frame contains more variables such as
#' the underlying random effect values and numeric values for sampled covariates
#' and outcome values.
#' *  `FALSE` (default): the returned data frame only contains variables that
#' are available in real clinical trials.
#'
#' @return Data frame with the same structure as [simulated_test_trial].
#' @export
#'
#' @examples
#' # simulate a single trial with default arguments
#' simulate_N_trials()
simulate_N_trials = function(control = "placebo",
                             delta = 1.3,
                             total_time = 50,
                             n_measurements = 5,
                             sample_size = 1000,
                             N_trials = 1,
                             details = FALSE) {
  list_of_trials = lapply(
    X = 1:N_trials,
    FUN = simulate_single_trial,
    control = control,
    delta = delta,
    total_time = total_time,
    n_measurements = n_measurements,
    sample_size = sample_size,
    details = FALSE
  )
  dplyr::bind_rows(list_of_trials)
}


#' Simulates one clinical trial under the supplied model parameters.
#'
#' @param control Type of control group, can be:
#' * `"placebo"`: Control group receives placebo
#' * `"symptomatic"`: control group receives a symptomatic treatment
#' @param delta deceleration factor. This parameterizes the disease-modifying
#'   treatment effect.
#' @param total_time Total trial duration in weeks.
#' @param n_measurements Number of equally spaced measurement occasions for each
#'   patient.
#' @param sample_size Total sample size in each trial.
#' @param trial_number Defaults to 1. This is used by [simulate_N_trial()] to
#' distinguish different trials, but is not relevant when
#' [simulate_single_trial()] is used by the user directly to simulate a single
#' trial.
#' @param details Return detailed data frame:
#' * `TRUE`: the returned data frame contains more variables such as
#' the underlying random effect values and numeric values for sampled covariates
#' and outcome values.
#' *  `FALSE` (default): the returned data frame only contains variables that
#' are available in real clinical trials.
#'
#' @return Data frame with the same structure as [simulated_test_trial].
#' @export
#'
#' @examples
#' # simulate a single trial with default parameters
#' simulate_single_trial()
#'
simulate_single_trial = function(control = "placebo",
                                 delta = 1.3,
                                 total_time = 50,
                                 n_measurements = 5,
                                 sample_size = 1000,
                                 trial_number = 1,
                                 details = FALSE) {
  # Simulate data for both treatment arms.
  data_control = simulate_single_arm(
    control = control,
    arm_sample_size = ceiling(sample_size / 2),
    delta = 1,
    total_time = total_time,
    n_measurements = n_measurements,
    details = details
  )
  data_control = dplyr::mutate(data_control,
                               arm = 0L)

  data_active = simulate_single_arm(
    control = control,
    arm_sample_size = ceiling(sample_size / 2),
    delta = delta,
    total_time = total_time,
    n_measurements = n_measurements,
    details = details
  )
  data_active = dplyr::mutate(data_active,
                              arm = 1L,
                              SubjId = SubjId + as.integer(ceiling(sample_size / 2)))

  # Returned combined data set with added variable to indicate trial number.
  dplyr::mutate(dplyr::bind_rows(data_active, data_control),
                trial_number = trial_number)
}

# The following three functions are reproduced from the adsim package to
# sample from the covariate distribution.
simApoE = function (n, p, ...)
{
  pApo <- c(p["pApo[1]"], p["pApo[2]"], p["pApo[3]"])
  simdat <- sample(0:2, n, replace = TRUE, prob = pApo)
  simdat
}

simAge = function (n,
                   p,
                   ApoE,
                   lbAge = -Inf,
                   ubAge = Inf,
                   ...)
{
  Eage <- c(p["Eage[1]"], p["Eage[2]"], p["Eage[3]"])
  Vage <- p["Vage"]
  age <- rnorm.trunc(n, Eage[ApoE + 1], sqrt(Vage), lower = lbAge,
                     upper = ubAge)
  age
}

simGen = function (n, p, ...)
{
  stats::rbinom(n, 1, p["pGen"])
}

rnorm.trunc = function (n,
                        mean = 0,
                        sd = 1,
                        lower = -Inf,
                        upper = Inf)
{
  qnorm.trunc(stats::runif(n), mean, sd, lower, upper)
}

qnorm.trunc = function (p,
                        mean = 0,
                        sd = 1,
                        lower = -Inf,
                        upper = Inf)
{
  stats::qnorm(p * stats::pnorm(upper, mean, sd) + (1 - p) * stats::pnorm(lower,
                                                     mean, sd),
        mean,
        sd)
}

# Function to generate data from a single arm in a single clinical trial.
simulate_single_arm = function(control,
                               arm_sample_size,
                               delta,
                               total_time,
                               n_measurements,
                               details = FALSE) {
  # Posterior medians of estimated model in Rogers et al. (2012).
  popPars = popPars_function()

  # Time point at which measurements are taken
  times <- seq(0, total_time, length.out = n_measurements)
  # Sample covariate value for each patient.
  bmmse <- stats::runif(arm_sample_size, 16, 26)
  apo <- simApoE(arm_sample_size, popPars)
  age <- simAge(arm_sample_size, popPars, apo)
  gender <- simGen(arm_sample_size, popPars)

  # Extract the relevant parameter estimates. For tauEta, tauAlpha, and tauResid
  # the mean of the corresponding inter-study random effects distribution is used.
  # We are thus considering the "average" study. In principle, we could also
  # "sample a trial" from the estimated trial-level random effects distribution.
  muEta <- popPars['nuEta']
  muAlpha <- popPars['nuAlpha']
  tauEta <- 1 / popPars['phiEta']^2
  sigmaEta <- 1 / sqrt(tauEta)
  tauAlpha <- 1 / popPars['phiAlpha']^2
  sigmaAlpha <- 1 / sqrt(tauAlpha)
  tauResid <- 1 / popPars['phiEpsilon']^2
  muEtaAdj <- muEta + popPars['lambdaEtaBMMSE'] * (bmmse - 21)
  muAlphaAdj <- muAlpha + popPars['lambdaAlphaBMMSE'] * (bmmse - 21) +
    popPars['lambdaAlphaApo1'] * (apo==1) + popPars['lambdaAlphaApo2'] * (apo==2) +
    popPars['lambdaAlphaAge'] * (age -75) + popPars['lambdaAlphaGen'] * gender
  # Sample patient-level random effects.
  eta <- stats::rnorm(arm_sample_size, muEtaAdj, sigmaEta)
  alpha <- stats::rnorm(arm_sample_size, muAlphaAdj, sigmaAlpha)

  # Vector of subject id's.
  subjId <- as.integer(1:arm_sample_size)
  # For each subject, all time points in the times vector are considered. This
  # gives us a dataframe in long format without covariates.
  dat <- expand.grid(times, subjId)
  names(dat) <- c("Week", "SubjId")
  # The covariates that were sampled before, are added to the data frame.
  dat$BMMSE <- bmmse[dat$SubjId]
  dat$Apo <- apo[dat$SubjId]
  dat$Age <- age[dat$SubjId]
  dat$Gender <- gender[dat$SubjId]
  # Eta_i is added.
  dat$Eta <- eta[dat$SubjId] # (not available in a real data set)
  # alpha_i is added.
  dat$Alpha <- alpha[dat$SubjId] # (not available in a real data set )

  # The placebo effect is added. The placebo effect is common for all patient
  # since no random effects are involved here.
  kel <- popPars['kel']
  keq <- kel + popPars['keqMinusKel']
  beta <- - popPars['aucPlacebo'] / (1 / kel - 1 / keq)
  # This is the placebo effect on the logit scale for each time point in the times
  # vector.
  ePlacebo <- beta * ( exp( -kel * times) - exp( -keq * times ) )

  if (control == "symptomatic"){
    # The symptomatic treatment effect is added. The symptomatic treatment effect is
    # common for all patient since no random effects are involved here. First, the
    # posterior medians are transformed to the actual model parameters.
    eStar <- popPars['eStar[1]']
    et50 <- popPars['et50[1]']
    gamma <- popPars['gamma[1]']
    b <- 12 / et50
    eDelta <- - (1 + b) * eStar / b
    # This is the symptomatic drug effect on the logit scale for each time point in
    # the times vector.
    eDon10 <- eDelta * times / (et50 + times)
    refDose <- 10
    eDon10 <- eDon10 * (10/refDose)^gamma
  }
  else {
    # no symptomatic effect
    eDon10 = 0 * times
  }

  # Add all effect to the logit scale/systematic component.
  dat = dplyr::mutate(dat,
                      LogitCondExp =
                        Eta + (1 / delta) * Alpha * Week + rep(ePlacebo + eDon10, arm_sample_size))
  # Transform the patient-level conditional expectations to the original ADASCog
  # scale.
  dat = dplyr::mutate(dat,
                      theta = exp(LogitCondExp) / (1 + exp(LogitCondExp)))
  # Add residual variability. We have two versions of the final outcome: (1) the
  # exactly simulated value and (2) the rounded value. In practice, ADAScog
  # scores are always integers, so it is advised to use the rounded values for
  # further analyses.
  dat = dplyr::mutate(dat,
                      ADAScog_numeric = 70 * stats::rbeta(
                        n = nrow(dat),
                        shape1 = theta * tauResid,
                        shape2 = (1 - theta) * tauResid
                      ))
  dat = dplyr::mutate(dat,
                      ADAScog_integer = as.integer(round(ADAScog_numeric, 0)))

  # Round of values for the covariates. In any realist setting, we will not have
  # baseline age more precise than years. The MMSE score can only take integer
  # values, so non-integer MMSE scores cannot occur in practice.
  dat = dplyr::mutate(dat,
                      Age_integer = as.integer(round(Age, digits = 0)),
                      BMMSE_integer = as.integer(round(Age, digits = 0)))


  if (!details) {
    dat$theta = NULL
    dat$LogitCondExp = NULL
    dat$Alpha = NULL
    dat$Eta = NULL
    dat$ADAScog_numeric = NULL
    dat$Age = NULL
    dat$BMMSE = NULL
  }
  return(dat)
}

# The only use of the following function is to store  and return the posterior
# medians of the finally estimated model in Rogers et al. (2012).
popPars_function = function() {
  popPars = c(
    0.005074,
    -0.7952,
    -0.1251,
    -0.0006761,
    -0.0001143,
    -2.4915e-05,
    0.0001254,
    -0.00065305,
    5.2945,
    0.027475,
    0.05362,
    0.1052,
    75.55,
    0.1684,
    0.3243,
    0.4138,
    48.505,
    0.0077795,
    2.7125,
    0.79345,
    1.0115,
    0.008905,
    0.4561,
    0.09536,
    0.001637,
    0.4219,
    0.435,
    0.1434,
    75.27,
    75.64,
    66.48,
    86.02,
    0.5454,
    0.24145,
    0.20285,
    0.11415,
    0.1382,
    0.1553,
    0.12615,
    1.6185,
    8.978,
    7.9065,
    -0.2159,
    0.080835,
    -23515,
    -4.883,
    0.016255,
    -0.07367,
    0.93725,
    0.3736,
    0.5199,
    0.18615,
    -0.36675,
    0.1921,
    -0.3476,
    -0.00226,
    0.05049,
    0.17905,
    -0.431,
    6699
  )
  names(popPars) = c(
    'nuAlpha',
    'nuEta',
    'lambdaEtaBMMSE',
    'lambdaAlphaBMMSE',
    'lambdaAlphaAge',
    'lambdaAlphaApo1',
    'lambdaAlphaApo2',
    'lambdaAlphaGen',
    'aucPlacebo',
    'kel',
    'keqMinusKel',
    'phiEpsilon',
    'kappaEpsilon',
    'phiEpsilonM',
    'kappaEpsilonM',
    'phiEta',
    'kappaEta',
    'phiAlpha',
    'kappaAlpha',
    'phiEtaM',
    'kappaEtaM',
    'phiAlphaM',
    'kappaAlphaM',
    'psiEta',
    'psiAlpha',
    'pApo[1]',
    'pApo[2]',
    'pApo[3]',
    'Eage[1]',
    'Eage[2]',
    'Eage[3]',
    'Vage',
    'pGen',
    'gamma[1]',
    'gamma[2]',
    'gamma[3]',
    'eStar[1]',
    'eStar[2]',
    'eStar[3]',
    'et50[1]',
    'et50[2]',
    'et50[3]',
    'beta',
    'keq',
    'deviance',
    'beta0',
    'betaAge',
    'betaBmmse',
    'v',
    'sdStudy',
    'omegaStudy[1]',
    'omegaStudy[2]',
    'omegaStudy[3]',
    'omegaStudy[4]',
    'omegaStudy[5]',
    'omegaStudy[6]',
    'omegaStudy[7]',
    'omegaStudy[8]',
    'omegaStudy[9]',
    'deviance'
  )
  return(popPars)
}
