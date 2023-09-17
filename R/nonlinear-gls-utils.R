nonlinear_gls_criterion_constructor = function(time_points,
                                               ctrl_estimates,
                                               exp_estimates,
                                               interpolation,
                                               vcov,
                                               j = 1:length(exp_estimates),
                                               gamma_0 = NULL) {
  # The inverted variance-covariance matrix.
  Sigma_inv = solve(vcov)
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
      data_vec = matrix(data = c(ctrl_estimates, exp_estimates), ncol = 1)
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
      data_vec = matrix(data = c(ctrl_estimates, exp_estimates), ncol = 1)
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
  # The inverted variance-covariance matrix.
  Sigma_inv = solve(vcov)
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
      data_vec = matrix(data = c(ctrl_estimates, exp_estimates), ncol = 1)
      # Compute the Jacobian for the trajectory function with respect to the
      # alpha-parameters.
      A = attr(
        deriv_f0_alpha(
          t_m = gamma * time_points[j + 1],
          x_ref = time_points,
          y_ref = alpha_vec,
          method = interpolation
        ),
        "gradient"
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
      data_vec = matrix(data = c(ctrl_estimates, exp_estimates), ncol = 1)
      # Compute the Jacobian for the trajectory function with respect to the
      # alpha-parameters.
      A = attr(
        deriv_f0_alpha(
          t_m = gamma * time_points[j + 1],
          x_ref = time_points,
          y_ref = alpha_vec,
          method = interpolation
        ),
        "gradient"
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
}
