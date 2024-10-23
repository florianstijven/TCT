#' Test for proportional slowing based on GLS criterion
#'
#' [proportional_slowing_gls_test()] performs a test for proportional slowing
#' based on the minimized GLS criterion.
#'
#' @inheritParams nonlinear_gls_estimator
#'
#'
#' @return A named vector with the following elements:
#' * `chi_squared`: The chi-squared test statistic that asymptotically follows
#' a chi-squared distribution with `df` degrees of freedom under proportional
#' slowing.
#' * `df`: Degrees of freedom of the chi-squared distribution under proportional
#' slowing.
#' * `p_value`: The p-value based on the asymptotic chi-squared distribution
#' with `df` degrees of freedom.
proportional_slowing_gls_test = function(time_points,
                                         ctrl_estimates,
                                         exp_estimates,
                                         interpolation,
                                         vcov,
                                         j,
                                         start_gamma) {
  minimized_criterion_restricted = nonlinear_gls_estimator(
    time_points = time_points,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    interpolation = interpolation,
    vcov = vcov,
    j = j,
    start_gamma = start_gamma
  )$criterion

  return(c(
    chi_squared = minimized_criterion_restricted,
    df = length(j) - 1,
    p_value = stats::pchisq(
      q = minimized_criterion_restricted,
      df = length(j) - 1,
      lower.tail = FALSE
    )
  ))
}
