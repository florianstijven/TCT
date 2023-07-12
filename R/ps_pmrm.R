
#' Analyze trial data with the proportion slowing PMRM
#'
#'
#'
#' @param data_trial A data frame from a single clinical trial. The following
#'   variables should be present:
#'   * `ADAScog_integer`: integer-valued ADAScog score
#'   * `arm_time`: factor that characterizes the treatment arm-time combination.
#'   * `time_int`: integer that indicates the measurement occasion. Should start
#'   with `1` for the first occasion, values for subsequent occasions should be
#'   consecutive integers.
#' @param b_start Starting value for the proportional slowing parameter.
#'   Defaults to 0.
#'
#' @return
analyze_ps_pmrm = function(data_trial,
                           b_start = 0.20) {

  baseline = FALSE
  time_points <<- unique(data_trial$Week)
  TPMRM_generalized_helper <- function(t, v, b) {
    TPMRM_generalized(t, v, b, time_points)
  }
  if (baseline) {
    time_points = time_points[-1]
    data_trial = baseline_data_format(data_trial)
    formula = formula(ADAScog_integer ~
                        TPMRM_generalized_helper(
                          Week,
                          v = cbind(v1, v2, v3, v4),
                          b = cbind(b, b, b, b)
                        ))
    params = list(v1 + v2 + v3 + v4 ~ 1 + b_ADAScog_integer,
                  b ~ arm + 0)
  }
  else {
    formula = ADAScog_integer ~
      TPMRM_generalized_helper(Week,
                               v = cbind(v0, v1, v2, v3, v4),
                               b = cbind(0, b, b, b, b))
    params = list(v0 + v1 + v2 + v3 + v4 ~ 1,
                  b ~ arm - 1)
  }
  #compute starting values
  start_vec = c(aggregate(data_trial$ADAScog_integer ~ data_trial$Week,
                          FUN = mean)[, 2],
                b_start)
  nlme::gnls(model = formula,
             data = data_trial,
             params = params,
             correlation = nlme::corSymm(form = ~ time_int | SubjId),
             weights = nlme::varIdent(form = ~ 1 | time_int),
             start = start_vec,
             control = nlme::gnlsControl(nlsTol = 500))
}
