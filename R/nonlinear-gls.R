#' TCT Estimation by Nonlinear GLS
#'
#' @param gamma_0 Defaults to `NULL`. If a non-null value is given to this argument,
#' the common acceleration factor is held constant at that value.
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
                                   interpolation,
                                   vcov,
                                   j = 1:length(exp_estimates),
                                   gamma_0 = NULL,
                                   ...) {
  objective_function = nonlinear_gls_criterion_constructor(
    time_points,
    ctrl_estimates,
    exp_estimates,
    interpolation,
    vcov,
    j = 1:length(exp_estimates),
    gamma_0 = NULL
  )
  # Minimize the GLS criterion in the parameters.
  if (is.null(gamma_0)) {
    optim_object = stats::optim(
      par = c(ctrl_estimates, 1),
      fn = objective_function,
      method = "BFGS",
      control = list(abstol = 1e-7)
    )
  }
  else {
    optim_object = stats::optim(
      par = c(ctrl_estimates),
      fn = objective_function,
      method = "BFGS",
      control = list(abstol = 1e-7)
    )
  }
  # Return the estimated parameters and the value of the minimized object function.
  return(
    list(
      estimates = optim_object$par,
      criterion = optim_object$value
    )
  )
}


nonlinear_gls_conf_int_common = function(time_points,
                                 ctrl_estimates,
                                 exp_estimates,
                                 interpolation,
                                 vcov,
                                 j = 1:length(exp_estimates),
                                 alpha = 0.05){
  # Force argument values. This is required because we're using these arguments
  # in a function factory.
  force(time_points); force(ctrl_estimates); force(exp_estimates)
  force(interpolation); force(vcov)
  force(j); force(alpha)

  # Extract the criterion value in the full model.
  nl_gls_object = nonlinear_gls_estimator(
    time_points = time_points,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov,
    interpolation = interpolation,
    j = 1:length(exp_estimates)
  )
  criterion_full = nl_gls_object$criterion
  gamma_est = nl_gls_object$estimates[length(nl_gls_object$estimates)]
  force(criterion_full)
  # Define function that computes the difference in criterion values of the full
  # model and the constrained model. In the latter model, the value for gamma is
  # held fixed.
    t_sq_value = function(gamma) {
      criterion_reduced = nonlinear_gls_estimator(
        time_points = time_points,
        ctrl_estimates = ctrl_estimates,
        exp_estimates = exp_estimates,
        vcov = vcov,
        interpolation = interpolation,
        j = 1:length(exp_estimates),
        gamma_0 = gamma
      )$criterion
      return(criterion_reduced - criterion_full)
    }
  # Compute critical value for the test-statistic.
  t_sq_critical = stats::pchisq(alpha, df = 1)
  # If the right limit of the test statistic, as a function of gamma, does not
  # cross the critical value, the upper confidence limit is infinity. The same
  # principle applies to the lower limit.
  if (t_sq_value(10) < t_sq_critical) {
    upper_limit = +Inf
  }
  else {
    upper_limit = stats::uniroot(
      f = function(gamma)
        sqrt(t_sq_value(gamma)) - sqrt(t_sq_critical),
      interval = c(gamma_est,
                   10),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }

  # Find lower limit
  if (t_sq_value(-10) < t_sq_critical) {
    lower_limit = -Inf
  }
  else {
    lower_limit = stats::uniroot(
      f = function(gamma)
        sqrt(t_sq_value(gamma)) - sqrt(t_sq_critical),
      interval = c(-10,
                   gamma_est),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }

  # Return estimated confidence interval.
  return(
    c(lower_limit, upper_limit)
  )
}


