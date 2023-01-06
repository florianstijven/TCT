#' Delta method with numerical derivatives
#'
#'
#'
#' @param par vector of parameter estimates
#' @param fct function that transforms the parameter values to the appropriate
#'   function
#' @param vcov The variance-covariance matrix of the parameter estimates in the
#'   `par` vector argument
#' @param ...
#'
#' @return A list with three element:
#'  * `estimate`: the transformed parameter estimates
#'  * `variance`: the variance-covariance matrix for the transformed parameter
#'   estimates
#'  * `partial`: the Jacobian matrix
#'
#' @examples
DeltaMethod = function (par, fct, vcov,  ...) {
  theta <- function(par)
    fct(par, ...)
  partial <- t(numDeriv::jacobian(theta, x = par))
  variance <- t(partial) %*% vcov %*% partial
  return(list(estimate = theta(par), variance = variance, partial = partial))
}

#' Transform vertical treatment effect estimates to time scale
#'
#' This function transform the vertical parameter estimates to parameter
#' estimates on the time scale.
#'
#' @details
#'
#' @param time_points Ordered vector that contains the times corresponding to
#'   the estimated means in the `ctrl_means` vector.
#' @param ctrl_means Estimated mean outcome in the control group at fixed
#'   occasions.
#' @param exp_means Estimated mean outcomes in the experimental group at fixed
#'   occasions.
#' @param vcov The variance-covariance matrix for the means. In order to map to
#'   the correct estimates, this matrix should be the variance-covariance matrix
#'   of `c(ctrl_means, exp_means)`.
#' @param interpolation Which interpolation method to use?
#'  * `"linear"`: linear interpolation
#'  * `"spline"`: natural cubic spline interpolation
#'  * `"monoH.FC`: monotone Hermite spline according to the method of Fritsch
#'  and Carlson
#'
#' @return A list with three element:
#'  * `estimate`: the transformed parameter estimates
#'  * `variance`: the variance-covariance matrix for the transformed parameter
#'   estimates
#'  * `partial`: the Jacobian matrix
#' @export
#'
#' @examples
delta_vertical_to_horizontal = function(time_points,
                                        ctrl_means,
                                        exp_means,
                                        vcov,
                                        interpolation = "spline"
                                        ) {
  # number of time points
  n_points = length(time_points)
  # rewrite get_new_time() for the DeltaMethod() function
  g_Delta = function(par, x_ref, method = interpolation) {
    y_ref = par[1:n_points]
    y_obs = par[(n_points + 1):length(par)]
    t_mapped = sapply(
      X = y_obs,
      FUN = get_new_time,
      y_ref = y_ref,
      x_ref = x_ref,
      method = method
    )
    return(t_mapped/time_points[-1])
  }
  se_delta = DeltaMethod(
    par = c(ctrl_means, exp_means),
    fct = g_Delta,
    vcov = vcov,
    x_ref = time_points,
    method = interpolation
  )
  return(se_delta)
}

pm_bootstrap_vertical_to_horizontal = function(time_points,
                                               ctrl_means,
                                               exp_means,
                                               vcov,
                                               interpolation = "spline",
                                               B = 100) {
  # number of time points
  n_points = length(time_points)
  # rewrite get_new_time() for the DeltaMethod() function
  g_Delta = function(par, x_ref, method = interpolation) {
    y_ref = par[1:n_points]
    y_obs = par[(n_points + 1):length(par)]
    t_mapped = sapply(
      X = y_obs,
      FUN = get_new_time,
      y_ref = y_ref,
      x_ref = x_ref,
      method = method
    )
    return(t_mapped/time_points[-1])
  }

  par_sampled = mvtnorm::rmvnorm(n = B,
                                 mean = c(ctrl_means, exp_means),
                                 sigma = vcov
                                 )
  estimates = matrix(0, nrow = B, ncol = 4)
  for (i in 1:B) {
    estimates[i, ] = g_Delta(par_sampled[i, ], time_points, interpolation)
  }

  return(estimates)
}

new_time_scale_estimates = function(coefficients,
                                    vcov) {

  structure(list(coefficients = coefficients,
                 vcov = vcov),
            class = "time_scale_estimates")
}



common_deceleration = function(estimates, vcov) {
  p = length(estimates)
  vec_1 = matrix(1, nrow = p, ncol = 1)
  common_estimate = (t(vec_1) %*% solve(vcov) %*% matrix(estimates, ncol = 1) ) /
    (t(vec_1) %*% solve(vcov) %*% vec_1)
  common_variance = (t(vec_1) %*% solve(vcov) %*% vec_1)**(-1)
  return(list(common_estimate = common_estimate,
              common_variance = common_variance))
}

bootstrap_c_decel = function(time_points,
                             ctrl_means,
                             exp_means,
                             vcov,
                             interpolation = "spline",
                             B = 100) {
  # number of time points
  n_points = length(time_points)
  # rewrite get_new_time() for the DeltaMethod() function
  g_Delta = function(par, x_ref, method = interpolation) {
    y_ref = par[1:n_points]
    y_obs = par[(n_points + 1):length(par)]
    t_mapped = sapply(
      X = y_obs,
      FUN = get_new_time,
      y_ref = y_ref,
      x_ref = x_ref,
      method = method
    )
    return(t_mapped/time_points[-1])
  }


  estimates = 1:B
  par_sampled = mvtnorm::rmvnorm(n = B,
                                 mean = c(ctrl_means, exp_means),
                                 sigma = vcov)
  for (i in 1:B) {

    tct_results = delta_vertical_to_horizontal(time_points,
                                 par_sampled[i, 1:n_points],
                                 par_sampled[i, (n_points + 1):length(par_sampled[1,])],
                                 vcov,
                                 interpolation = "spline")
    estimates[i] = common_deceleration(estimates = tct_results$estimate,
                                       vcov = tct_results$variance)$common_estimate
  }

  return(estimates)
}


full_tct_analysis = function(time_points,
                             ctrl_means,
                             exp_means,
                             vcov,
                             interpolation = "spline") {
  tct_results = delta_vertical_to_horizontal(time_points,
                                             ctrl_means,
                                             exp_means,
                                             vcov,
                                             interpolation = "spline")
  print(car::linearHypothesis(model = t,
                        vcov. = tct_results$variance,
                        coef. = tct_results$estimate,
                        rhs = c(1, 1, 1, 1),
                        hypothesis.matrix = matrix(c(1, 0, 0, 0,
                                                     0, 1, 0, 0,
                                                     0, 0, 1, 0,
                                                     0, 0, 0, 1), nrow = 4, byrow = TRUE)))
  print(car::linearHypothesis(model = t,
                              vcov. = tct_results$variance,
                              coef. = tct_results$estimate,
                              rhs = c(0, 0, 0),
                              hypothesis.matrix = matrix(c(1, -1, 0, 0,
                                                           1, 0, -1, 0,
                                                           1, 0, 0, -1), nrow = 3, byrow = TRUE)))
  tct_common = common_deceleration(tct_results$estimate, tct_results$variance)
  print(
    paste0("gamma common: ", tct_common$common_estimate, "; se = ", sqrt(tct_common$common_variance))
  )
  print(
    paste0("[", -1.96 * sqrt(tct_common$common_variance) + tct_common$common_estimate,
           ", ", 1.96 * sqrt(tct_common$common_variance) + tct_common$common_estimate)
  )
  print(2 * (1 - pnorm(abs(tct_common$common_estimate) / sqrt(tct_common$common_variance))))

}
