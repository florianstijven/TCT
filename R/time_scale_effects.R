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
