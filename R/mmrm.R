#' Analyze trial data with MMRM
#'
#' The [analyze_mmrm()] function analyzes data from a single trial with the
#' mixed model for repeated measures (See Details).
#'
#' @param data_trial A data frame from a single clinical trial. The following
#'   variables should be present:
#'   * `ADAScog_integer`: integer-valued ADAScog score
#'   * `arm_time`: factor that characterizes the treatment arm-time combination.
#'   * `time_int`: integer that indicates the measurement occasion. Should start
#'   with `1` for the first occasion, values for subsequent occasions should be
#'   consecutive integers.
#' @param method Which method should be used for fitting the MMRM? Should be one
#'   of the following strings:
#'   * `"ML"` (default): maximum likelihood
#'   * `"REML"`: restricted maximum likelihood
#'
#' @return fitted-model object from [nlme::gls()]
#' @export
#'
#' @examples
#'
analyze_mmrm = function(data_trial, method = "ML") {
  # The following options could be elaborated on in the future.
  baseline = FALSE
  covariates = FALSE
  if (baseline) {
    # Transform data to "baseline format".
    data_trial = baseline_data_format(data_trial)
    formula = formula(ADAScog_integer~arm_time + b_ADAScog_integer + 0)
  }
  else {
    formula = formula(ADAScog_integer~arm_time + 0)
  }

  if (covariates) {
    formula = update.formula(old = formula, .~. +
                               BMMSE_integer + BMMSE_integer:as.factor(time_int)
    )
  }

  nlme::gls(
    formula,
    data = data_trial,
    correlation = nlme::corSymm(form = ~ time_int | SubjId),
    weights = nlme::varIdent(form = ~ 1 | time_int),
    method = method
  )
}
