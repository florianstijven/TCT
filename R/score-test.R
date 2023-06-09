#' Compute TCT score test z-value
#'
#' The [score_test()] function computes the z-value for the score test at
#' measurement `j`. This corresponds to `exp_estimates[j]`.
#'
#' @param j Measurement occasion to test acceleration factor for.
#' @param gamma_0 Value under the null hypothesis for the acceleration factor.
#' @inheritParams DeltaMethod
#'
#'
#' @return (numeric) z-value. This test-statistic follows a standard normal
#'   distribution under the null hypothesis.
score_test = function(time_points,
                      ctrl_estimates,
                      exp_estimates,
                      ref_fun,
                      interpolation,
                      vcov,
                      j,
                      gamma_0 = 1) {
  # Number of measurement occasions in the experimental group.
  K = length(exp_estimates)
  # Compute gradient of g_{\gamma}. The first length(time_points) elements of
  # the gradient are computed.
  deriv_f0_alpha_t_j  = -1 * attr(
    deriv_f0_alpha(
      t_m = gamma_0 * time_points[j + 1],
      x_ref = time_points,
      y_ref = ctrl_estimates,
      method = interpolation
    ),
    "gradient"
  )
  # The remaining elements are all zeros and 1 one.
  grad_g = c(deriv_f0_alpha_t_j,
             rep(0, j - 1),
             1,
             rep(0, K - j))
  grad_g = matrix(grad_g, ncol = 1)
  # Compute the variance of the test statistic under the null.
  sigma_sq = t(grad_g) %*% vcov %*% grad_g
  # Compute z-value
  z = (exp_estimates[j] - ref_fun(gamma_0 * time_points[j + 1])) / sqrt(sigma_sq)
  # Return z-value
  return(as.numeric(z))
}

#' Compute Confidence Interval Based on Score Test
#'
#' The [score_conf_int()] function computes the confidence interval for the
#' score test at measurement `j`.
#'
#' @param alpha `1 - alpha` represent the two-sided confidence level. Defaults
#'   to `0.05`.
#' @inheritParams score_test
#'
#' @return (numeric) vector with two elements. The first element is the lower
#'   confidence limit, the second element is the upper confidence limit.
score_conf_int = function(time_points,
                          ctrl_estimates,
                          exp_estimates,
                          ref_fun,
                          interpolation,
                          vcov,
                          j,
                          alpha = 0.05){
  # Construct function of gamma that return the z-value.
  z_value = function(gamma) {
    return(score_test(time_points,
                      ctrl_estimates,
                      exp_estimates,
                      ref_fun,
                      interpolation,
                      vcov,
                      j,
                      gamma_0 = gamma))
  }
  # Find upper limit
  z_critical = qnorm(p = 1 - alpha / 2)
  upper_limit = stats::uniroot(
    f = function(gamma)
      z_value(gamma) + z_critical,
    interval = c(-5,
                 5),
    tol = .Machine$double.eps ^ 0.5,
    maxiter = 1e3
  )$root
  # Find lower limit
  lower_limit = stats::uniroot(
    f = function(gamma)
      z_value(gamma) - z_critical,
    interval = c(-5,
                 5),
    tol = .Machine$double.eps ^ 0.5,
    maxiter = 1e3
  )$root
  # Return estimated confidence interval.
  return(
    c(lower_limit, upper_limit)
  )
}

score_test_common = function(time_points,
                             ctrl_estimates,
                             exp_estimates,
                             ref_fun,
                             interpolation,
                             vcov,
                             gamma_0 = 1,
                             type = "omnibus"){
  # Number of measurement occasions in the experimental group.
  K = length(exp_estimates)
  # The Jacobian matrix is computed in two parts. First, the (K x K + 1) part is
  # computed.
  J_f0_alpha_t  = -1 * attr(
    deriv_f0_alpha(
      t_m = gamma_0 * time_points[-1],
      x_ref = time_points,
      y_ref = ctrl_estimates,
      method = interpolation
    ),
    "gradient"
  )
  # Next, the second part of the Jacobian. This corresponds to the identify
  # matrix. We join both parts two get the Jacobian.
  J = t(cbind(J_f0_alpha_t, diag(x = 1, nrow = K)))
  # Compute the variance of the test statistic under the null.
  Sigma_g = t(J) %*% vcov %*% J
  g = (exp_estimates - ref_fun(gamma_0 * time_points[-1]))
  if (type == "omnibus") {
    # Compute test-statistic
    t_sq =  g %*% solve(Sigma_g) %*% g
    return(t_sq)
  }
  else if (type == "directional") {
    # SHOULD BE WEIGHTED IN SOME OTHER WAY!!!!!!!!!!
    ones = matrix(1, nrow = K)
    sigma_sum_g = as.numeric(t(ones) %*% Sigma_g %*% ones)
    z_value = sum(g) / sqrt(sigma_sum_g)
    return(z_value)
  }
}

score_conf_int_common = function(time_points,
                                 ctrl_estimates,
                                 exp_estimates,
                                 ref_fun,
                                 interpolation,
                                 vcov,
                                 gamma_est, type = "omnibus"){
  # Number of measurements after treatment.
  K = length(exp_estimates)
  if (type == "omnibus") {
    # Construct function of gamma that return the z-value.
    t_sq_value = function(gamma) {
      return(score_test_common(time_points,
                               ctrl_estimates,
                               exp_estimates,
                               ref_fun,
                               interpolation,
                               vcov,
                               gamma_0 = gamma))
    }
    # Find upper limit
    t_sq_critical = qchisq(p = 1 - alpha / 2, df = K)
    upper_limit = stats::uniroot(
      f = function(gamma)
        t_sq_value(gamma) - t_sq_critical,
      interval = c(gamma_est,
                   5),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
    # Find lower limit
    lower_limit = stats::uniroot(
      f = function(gamma)
        t_sq_value(gamma) - t_sq_critical,
      interval = c(-5,
                   gamma_est),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }
  else if (type == "directional") {
    # Construct function of gamma that return the z-value.
    z_value = function(gamma) {
      return(score_test_common(time_points,
                        ctrl_estimates,
                        exp_estimates,
                        ref_fun,
                        interpolation,
                        vcov,
                        gamma_0 = gamma,
                        type = "directional"))
    }
    # Find upper limit
    z_critical = qnorm(p = 1 - alpha / 2)
    upper_limit = stats::uniroot(
      f = function(gamma)
        z_value(gamma) + z_critical,
      interval = c(-5,
                   5),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
    # Find lower limit
    lower_limit = stats::uniroot(
      f = function(gamma)
        z_value(gamma) - z_critical,
      interval = c(-5,
                   5),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }

  # Return estimated confidence interval.
  return(
    c(lower_limit, upper_limit)
  )
}

