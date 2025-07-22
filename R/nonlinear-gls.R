#' TCT Estimation by Nonlinear GLS
#'
#' @param gamma_0 Defaults to `NULL`. If a non-null value is given to this argument,
#' the common acceleration factor is held constant at that value.
#' @param start_gamma Starting value for the acceleration factor used for for the
#' generalized least squares estimator.
#' @inheritParams contrast_estimate_common
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
                                   start_gamma = 0.75,
                                   ...) {
  objective_function = nonlinear_gls_criterion_constructor(
    time_points,
    ctrl_estimates,
    exp_estimates,
    interpolation,
    vcov,
    j = j,
    gamma_0 = gamma_0
  )
  gradient_function = gradient_gls_criterion_constructor(
    time_points,
    ctrl_estimates,
    exp_estimates,
    interpolation,
    vcov,
    j = j,
    gamma_0 = gamma_0
  )
  # If there is a row/column of zeroes in vcov, this signifies a known mean at a
  # certain time point. We need to remove these rows/columns from the vcov
  # matrix and further treat the correspond mean "estimate" as known.
  if (any(rowSums(vcov) == 0)) {
    # Compute the indices of the known means.
    known_means_index = which(rowSums(vcov) == 0)
    # Compute the indices of the known means in the ctrl_estimates.
    known_ctrl_means_index = known_means_index[known_means_index <= length(ctrl_estimates)]
    # Compute the indices of the estimated means in the ctrl_estimates.
    est_ctrl_means_index = (1:length(ctrl_estimates))[-known_ctrl_means_index]
  } else {
    known_means_index = c()
    known_ctrl_means_index = c()
    est_ctrl_means_index = 1:length(ctrl_estimates)
  }
  # Minimize the GLS criterion in the parameters.
  if (is.null(gamma_0)) {
    optim_object = stats::optim(
      par = c(ctrl_estimates[est_ctrl_means_index], start_gamma),
      fn = objective_function,
      gr = gradient_function,
      method = "BFGS",
      control = list(abstol = 1e-7, reltol = 1e-8)
    )
  }
  else {
    optim_object = stats::optim(
      par = c(ctrl_estimates[est_ctrl_means_index]),
      fn = objective_function,
      gr = gradient_function,
      method = "BFGS",
      control = list(abstol = 1e-7, reltol = 1e-8)
    )
  }
  # If there are known means, we need to add them to the estimated parameters.
  # The `estimates` vector contains both estimated and known means.
  if (length(known_means_index) > 0) {
    estimates = ctrl_estimates
    if (is.null(gamma_0)) {
      estimates[c(est_ctrl_means_index, length(ctrl_estimates) + 1)] = optim_object$par
    } else {
      estimates[est_ctrl_means_index] = optim_object$par
    }
  } else {
    estimates = optim_object$par
  }

  # Return the estimated parameters and the value of the minimized object function.
  return(
    list(
      estimates = estimates,
      criterion = optim_object$value
    )
  )
}

#' Compute test statistic and p-value for the nonlinear GLS estimator
#'
#' @inheritParams contrast_test
#' @inheritParams nonlinear_gls_estimator
#'
#'
#' @inherit contrast_test_common return
nonlinear_gls_test = function(time_points,
                              ctrl_estimates,
                              exp_estimates,
                              interpolation,
                              vcov,
                              j = 1:length(exp_estimates),
                              gamma_0 = 1,
                              start_gamma = 0.75,
                              ...) {
  # Compute criterion function in optimum.
  criterion_full = nonlinear_gls_estimator(time_points,
                                           ctrl_estimates,
                                           exp_estimates,
                                           interpolation,
                                           vcov,
                                           j = j,
                                           gamma_0 = NULL,
                                           start_gamma = start_gamma,
                                           ...)$criterion
  # Compute criterion function under the null.
  criterion_reduced = nonlinear_gls_estimator(time_points,
                                           ctrl_estimates,
                                           exp_estimates,
                                           interpolation,
                                           vcov,
                                           j = j,
                                           gamma_0 = gamma_0,
                                           start_gamma = start_gamma,
                                           ...)$criterion
  test_statistic = criterion_reduced - criterion_full
  return(c(
    "chi-squared" = test_statistic,
    "p-value" = stats::pchisq(test_statistic, df = 1, lower.tail = FALSE)
  ))

}

#' Standard Error of nonlinear GLS estimator of the common acceleration factor
#'
#' @inheritParams nonlinear_gls_estimator
#' @param gamma_est Estimated value for the common acceleration factor.
#' @param alpha_est Estimated value of the mean vector in the control group.
#'   Note that this should be the estimate as returned by
#'   [nonlinear_gls_estimator()], and not the original estimates.
#'
#' @return (numeric) Estimated SE of the estimator.
nonlinear_gls_estimator_se = function(time_points,
                                      interpolation,
                                      vcov,
                                      j,
                                      gamma_est,
                                      alpha_est) {
  vcov_matrix = nonlinear_gls_estimator_vcov(time_points,
                                             interpolation,
                                             vcov,
                                             j,
                                             gamma_est,
                                             alpha_est)
  se = sqrt(vcov_matrix[nrow(vcov_matrix), nrow(vcov_matrix)])
  return(se)
}

#' Variance-Covariance matrix of nonlinear GLS estimator
#'
#' @inheritParams nonlinear_gls_estimator_se
#'
#' @return (numeric) variance-covariance matrix where the rows and columns
#' correspond to
#' \eqn{(\hat{\alpha}_0, \hat{\alpha}_1, ..., \hat{\alpha}_K, \hat{\gamma})'}.
nonlinear_gls_estimator_vcov = function(time_points,
                                        interpolation,
                                        vcov,
                                        j,
                                        gamma_est,
                                        alpha_est){
  # If there is a row/column of zeroes in vcov, this signifies a known mean at a
  # certain time point. We need to remove these rows/columns from the vcov
  # matrix and further treat the correspond mean "estimate" as known.
  if (any(rowSums(vcov) == 0)) {
    # Compute the indices of the known means.
    known_means_index = which(rowSums(vcov) == 0)
    # Compute the indices of the known means in the alpha_est.
    known_ctrl_means_index = known_means_index[known_means_index <= length(alpha_est)]
    # Compute the indices of the estimated means in the ctrl_estimates.
    est_ctrl_means_index = (1:length(alpha_est))[-known_ctrl_means_index]
  } else {
    known_means_index = c()
    known_ctrl_means_index = c()
    est_ctrl_means_index = 1:length(alpha_est)
  }

  # The truly estimated mean parameters in the control group.
  alpha_vec_est = alpha_est[est_ctrl_means_index]

  # The mean parameters in the experimental group corresponding to `j` may not
  # be known. If this is the case, an error is raised.
  if (!is.null(known_means_index)) {
    if (known_means_index %in% (j + length(alpha_est))) {
      stop("The mean parameters in the experimental group corresponding to `j` must have a non-zero variance.")
    }
  }

  # Compute subsetting vector.
  subset_vec = c(1:length(alpha_est), j + length(alpha_est))
  # Remove mean "estimates" which are known. If all means are estimated, no
  # elements are removed.
  if (length(known_means_index) > 0) {
    subset_vec = subset_vec[!(subset_vec %in% known_means_index)]
  }
  vcov = vcov[subset_vec, subset_vec]

  # Compute the reference trajectory.
  ref_fun = ref_fun_constructor(x_ref = time_points,
                                y_ref = alpha_est,
                                method = interpolation)
  # Compute the Jacobian for the trajectory function with respect to the
  # alpha-parameters.
  A = deriv_f0_alpha(
    t_m = gamma_est * time_points[j + 1],
    x_ref = time_points,
    y_ref = alpha_est,
    method = interpolation
  )
  # The "known" alpha parameters have to be excluded.
  A = A[, est_ctrl_means_index]

  D_t = diag(time_points[j + 1])
  B = cbind(diag(1, nrow = length(alpha_vec_est)), matrix(0, nrow = length(alpha_vec_est), ncol = 1))
  C = cbind(A, D_t %*% f0_gradient_t(ref_fun, time_points[j + 1] * gamma_est))
  J = rbind(B, C)
  # Compute and return the variance-covariance matrix.
  # Ensure that vcov is a symmetric matrix
  vcov[lower.tri(vcov)] = t(vcov)[lower.tri(vcov)]
  # Ensure that the following matrix is a symmetric matrix.
  Z = t(J) %*% mnormt::pd.solve(vcov) %*% J
  Z[lower.tri(Z)] = t(Z)[lower.tri(Z)]
  vcov_matrix = mnormt::pd.solve(Z)
  return(vcov_matrix)
}


nonlinear_gls_conf_int_common = function(time_points,
                                 ctrl_estimates,
                                 exp_estimates,
                                 interpolation,
                                 vcov,
                                 j = 1:length(exp_estimates),
                                 alpha = 0.05,
                                 start_gamma = 0.75){
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
    j = j,
    start_gamma = start_gamma
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
      j = j,
      gamma_0 = gamma,
      start_gamma = start_gamma
    )$criterion
    # The absolute value of the difference is returned because in a epsilon
    # neighborhood around gamma, the difference could become negative.
    return(abs(criterion_reduced - criterion_full))
  }
  # Compute critical value for the test-statistic.
  t_sq_critical = stats::qchisq(1 - alpha, df = 1)

  # Compute limits of search interval for computing the confidence intervals.
  # This is based on the wald confidence interval, but adding additional
  # tolerance.
  gamma_se = nonlinear_gls_estimator_se(
    time_points = time_points,
    interpolation = interpolation,
    vcov = vcov,
    j = j,
    gamma_est = gamma_est,
    alpha_est = nl_gls_object$estimates[-1 * length(nl_gls_object$estimates)]
  )
  start_upper = gamma_est + stats::qnorm(1 - alpha / 10) * gamma_se
  start_lower = gamma_est - stats::qnorm(1 - alpha / 10) * gamma_se

  # If the right limit of the test statistic, as a function of gamma, does not
  # cross the critical value, the upper confidence limit is infinity. The same
  # principle applies to the lower limit.
  if (t_sq_value(start_upper) < t_sq_critical) {
    upper_limit = +Inf
  }
  else {
    upper_limit = stats::uniroot(
      f = function(gamma)
        sqrt(t_sq_value(gamma)) - sqrt(t_sq_critical),
      interval = c(gamma_est,
                   start_upper),
      tol = .Machine$double.eps ^ 0.5,
      extendInt = "upX",
      f.lower = - sqrt(t_sq_critical),
      maxiter = 1e3
    )$root
  }

  # Find lower limit
  if (t_sq_value(start_lower) < t_sq_critical) {
    lower_limit = -Inf
  }
  else {
    lower_limit = stats::uniroot(
      f = function(gamma)
        sqrt(t_sq_value(gamma)) - sqrt(t_sq_critical),
      interval = c(start_lower,
                   gamma_est),
      tol = .Machine$double.eps ^ 0.5,
      extendInt = "downX",
      f.upper = - sqrt(t_sq_critical),
      maxiter = 1e3
    )$root
  }

  # Return estimated confidence interval.
  return(
    c(lower_limit, upper_limit)
  )
}


