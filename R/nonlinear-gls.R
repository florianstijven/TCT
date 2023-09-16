#' TCT Estimation by Nonlinear GLS
#'
#' @inheritParams score_estimate_common
#'
#' @return A list with the following elements:
#' * `estimates`: Estimated parameters by minimizing the generalized least
#'  squares criterion. The first `length(ctrl_estimates)` values are the estimated
#'  means in the control group. The last element is the estimated common acceleration
#'  factor.
#'  * `criterion`: Value of the criterion function in the parameter estimates.

nonlinear_gls_estimator = function(time_points,
                                   ctrl_estimates,
                                   exp_estimates,
                                   ref_fun,
                                   interpolation,
                                   vcov,
                                   j = 1:length(exp_estimates),
                                   ...) {
  # The inverted variance-covariance matrix.
  Sigma_inv = solve(vcov)
  # Define a function that returns that generalized least squares criterion.
  objective_function = function(param) {
    # Vector of mean parameters in the control group.
    alpha_vec = param[1:length(ctrl_estimates)]
    # Acceleration factor parameter.
    gamma = param[length(ctrl_estimates) + 1]
    # Predicted means. This is the alpha vector for the means in the control
    # group, and the trajectory function with accelerated time for the
    # experimental group.
    predicted_means = matrix(data = c(alpha_vec, ref_fun(time_points[j] * gamma)),
                             ncol = 1)
    # The "observed data" vector. In this line of thought, the estimated means
    # represent the observed data.
    data_vec = matrix(data = c(ctrl_estimates, exp_estimates), ncol = 1)
    # Compute and return the GLS criterion function.
    criterion = t(data_vec - predicted_means) %*% Sigma_inv %*% (data_vec - predicted_means)
  }
  # Minimize the GLS criterion in the parameters.
  optim_object = stats::optim(
    par = c(ctrl_estimates, 1),
    fn = objective_function
  )
  # Return the estimated parameters and the value of the minimized object function.
  return(
    list(
      estimates = optim_object$par,
      criterion = optim_object$value
    )
  )
}
