nonlinear_gls_criterion_constructor = function(time_points,
                                               ctrl_estimates,
                                               exp_estimates,
                                               interpolation,
                                               vcov,
                                               j = 1:length(exp_estimates),
                                               gamma_0 = NULL) {
  # Compute subsetting vector.
  subset_vec = c(1:length(ctrl_estimates), j + length(ctrl_estimates))
  # The inverted variance-covariance matrix.
  Sigma_inv = solve(vcov[subset_vec, subset_vec])
  if (is.null(gamma_0)) {
    # Define a function that returns that generalized least squares criterion for
    # the full model. That is the model where gamma is estimated.
    objective_function = function(param) {
      # Vector of mean parameters in the control group.
      alpha_vec = param[1:length(ctrl_estimates)]
      # Compute the reference trajectory.
      ref_fun = ref_fun_constructor(x_ref = time_points,
                                    y_ref = alpha_vec,
                                    method = interpolation)
      # Acceleration factor parameter.
      gamma = param[length(ctrl_estimates) + 1]
      # Predicted means. This is the alpha vector for the means in the control
      # group, and the trajectory function with accelerated time for the
      # experimental group.
      predicted_means = matrix(data = c(alpha_vec, ref_fun(time_points[j + 1] * gamma)),
                               ncol = 1)
      # The "observed data" vector. In this line of thought, the estimated means
      # represent the observed data.
      data_vec = matrix(data = c(ctrl_estimates, exp_estimates[j]), ncol = 1)
      # Compute and return the GLS criterion function.
      criterion = t(data_vec - predicted_means) %*% Sigma_inv %*% (data_vec - predicted_means)
      return(criterion)
    }
  }
  else {
    # Define a function that returns that generalized least squares criterion for
    # the reduced model. That is the model where gamma is held fixed at a certain
    # value.
    objective_function= function(param) {
      # Vector of mean parameters in the control group.
      alpha_vec = param[1:length(ctrl_estimates)]
      # Compute the reference trajectory.
      ref_fun = ref_fun_constructor(x_ref = time_points,
                                    y_ref = alpha_vec,
                                    method = interpolation)
      # Acceleration factor parameter is held fixed.
      gamma = gamma_0
      # Predicted means. This is the alpha vector for the means in the control
      # group, and the trajectory function with accelerated time for the
      # experimental group.
      predicted_means = matrix(data = c(alpha_vec, ref_fun(time_points[j + 1] * gamma)),
                               ncol = 1)
      # The "observed data" vector. In this line of thought, the estimated means
      # represent the observed data.
      data_vec = matrix(data = c(ctrl_estimates, exp_estimates[j]), ncol = 1)
      # Compute and return the GLS criterion function.
      criterion = t(data_vec - predicted_means) %*% Sigma_inv %*% (data_vec - predicted_means)
      return(criterion)
    }
  }
  return(objective_function)
}

gradient_gls_criterion_constructor = function(time_points,
                                              ctrl_estimates,
                                              exp_estimates,
                                              interpolation,
                                              vcov,
                                              j = 1:length(exp_estimates),
                                              gamma_0 = NULL) {
  # Compute subsetting vector.
  subset_vec = c(1:length(ctrl_estimates), j + length(ctrl_estimates))
  # The inverted variance-covariance matrix.
  Sigma_inv = solve(vcov[subset_vec, subset_vec])
  # Compute values that will be re-used. This improves efficiency.
  k = length(time_points)
  X = splines::ns(
    x = time_points,
    knots = time_points[2:(k - 1)],
    Boundary.knots = time_points[c(1, k)],
    intercept = TRUE
  )
  A_helper = solve(t(X) %*% X) %*% t(X)
  if (is.null(gamma_0)) {
    # Define a function that returns that gradient of the generalized least
    # squares criterion for the full model. That is the model where gamma is
    # estimated.
    gradient_function = function(param) {
      # Vector of mean parameters in the control group.
      alpha_vec = param[1:length(ctrl_estimates)]
      # Compute the reference trajectory.
      ref_fun = ref_fun_constructor(x_ref = time_points,
                                    y_ref = alpha_vec,
                                    method = interpolation)
      # Acceleration factor parameter.
      gamma = param[length(ctrl_estimates) + 1]
      # Predicted means. This is the alpha vector for the means in the control
      # group, and the trajectory function with accelerated time for the
      # experimental group.
      predicted_means = matrix(data = c(alpha_vec, ref_fun(time_points[j + 1] * gamma)),
                               ncol = 1)
      # The "observed data" vector. In this line of thought, the estimated means
      # represent the observed data.
      data_vec = matrix(data = c(ctrl_estimates, exp_estimates[j]), ncol = 1)
      # Compute the Jacobian for the trajectory function with respect to the
      # alpha-parameters.
      A = deriv_f0_alpha(
        t_m = gamma * time_points[j + 1],
        x_ref = time_points,
        y_ref = alpha_vec,
        method = interpolation,
        A = A_helper
      )
      D_t = diag(time_points[j + 1])
      B = cbind(diag(1, nrow = length(alpha_vec)), matrix(0, nrow = length(alpha_vec), ncol = 1))
      C = cbind(A, D_t %*% f0_gradient_t(ref_fun, time_points[j + 1] * gamma))
      # Compute and return the GLS criterion function.
      gradient = -2 * t(data_vec - predicted_means) %*%
        Sigma_inv %*%
        rbind(B, C)
      return(gradient)
    }
  }
  else {
    # Define a function that returns that gradient of the generalized least
    # squares criterion for the reduced model. That is the model where gamma is
    # fixed at some value.
    gradient_function = function(param) {
      # Vector of mean parameters in the control group.
      alpha_vec = param[1:length(ctrl_estimates)]
      # Compute the reference trajectory.
      ref_fun = ref_fun_constructor(x_ref = time_points,
                                    y_ref = alpha_vec,
                                    method = interpolation)
      # Acceleration factor parameter.
      gamma = gamma_0
      # Predicted means. This is the alpha vector for the means in the control
      # group, and the trajectory function with accelerated time for the
      # experimental group.
      predicted_means = matrix(data = c(alpha_vec, ref_fun(time_points[j + 1] * gamma)),
                               ncol = 1)
      # The "observed data" vector. In this line of thought, the estimated means
      # represent the observed data.
      data_vec = matrix(data = c(ctrl_estimates, exp_estimates[j]), ncol = 1)
      # Compute the Jacobian for the trajectory function with respect to the
      # alpha-parameters.
      A = deriv_f0_alpha(
        t_m = gamma * time_points[j + 1],
        x_ref = time_points,
        y_ref = alpha_vec,
        method = interpolation,
        A = A_helper
      )
      B = diag(1, nrow = length(alpha_vec))
      C = A
      # Compute and return the GLS criterion function.
      gradient = -2 * t(data_vec - predicted_means) %*%
        Sigma_inv %*%
        rbind(B, C)
      return(gradient)
    }
  }
  return(gradient_function)
}
