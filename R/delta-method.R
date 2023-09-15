#' Delta method
#'
#' The [DeltaMethod()] function applies the delta-method in the context of the
#' time-component tests. This function is based on a combination of analytical
#' (where possible) and numerical derivatives.
#'
#' @param ref_fun Intra- and extrapolation function that is returned by
#'   `ref_fun_constructor()`.
#' @inheritParams TCT_meta
#'
#' @details
#'
#' We can write the slowing factor at \eqn{t_j} as a function of the the mean
#' parameters,
#' \deqn{
#' \boldsymbol{\gamma} = \boldsymbol{\gamma}(\boldsymbol{t}; \boldsymbol{\alpha}, \boldsymbol{\beta}).
#' }
#' Further assume that \eqn{\hat{\boldsymbol{\alpha}}} and \eqn{\hat{\boldsymbol{\beta}}}
#'are normally distributed,
#'\deqn{
#' (\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})' \sim N\left( (\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})', D \right)
#' }
#' where \eqn{(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})'} is the true
#' parameter vector. We can apply the delta method to obtain that \eqn{\hat{\boldsymbol{\gamma}}}
#' is approximately distributed as
#' \deqn{
#' \hat{\boldsymbol{\gamma}} \sim N\left( \boldsymbol{\gamma_0},  J(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0}) \cdot D \cdot J(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})^t \right)
#' }
#' where \eqn{J(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})} is the Jacobian
#' matrix of \eqn{\boldsymbol{\gamma}( \boldsymbol{t}; \boldsymbol{\alpha}, \boldsymbol{\beta})}
#' evaluated in \eqn{(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})'}.
#'
#' To obtain the approximate sampling distribution of \eqn{\hat{\boldsymbol{\gamma}}}
#' in practice, we replace the true parameters by the estimated counterparts,
#' \deqn{(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})' \text{ and } \hat{D}.}
#' We thus obtain the following quantities,
#' \deqn{
#' \hat{\boldsymbol{\gamma}} = \boldsymbol{\gamma}(\boldsymbol{t}; \hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})
#' }
#' and
#' \deqn{
#' \Sigma_{\boldsymbol{\gamma}} = J(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}}) \cdot \hat{D} \cdot J(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})^t.
#' }
#'
#' @return A list with three element:
#'  * `estimate`: the transformed parameter estimates, \eqn{\hat{\boldsymbol{\gamma}}}.
#'  * `variance`: the (estimated) variance-covariance matrix for the transformed parameter
#'   estimates, \eqn{\Sigma_{\boldsymbol{\gamma}}}.
#'  * `partial`: the Jacobian matrix, \eqn{J(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})}.
#'
DeltaMethod = function (time_points,
                        ctrl_estimates,
                        exp_estimates,
                        ref_fun,
                        interpolation,
                        vcov)
{
  # Vector of estimated acceleration factors
  gamma_est = g_Delta_bis(par = c(ctrl_estimates, exp_estimates),
                          interpolation = interpolation,
                          time_points = time_points)
  # Jacobian matrix
  # Condition handler is used here since derivatives may not exist in special
  # situations such as when no unique time mapping is available.
  partial <- variance <- NA
  tryCatch(expr = {
    partial <- t(
      jacobian_tct(
        t_m = gamma_est * time_points[-1],
        t_j = time_points[-1],
        x_ref = time_points,
        y_ref = ctrl_estimates,
        ref_fun = ref_fun,
        method = interpolation
      )
    )
    # Apply delta method to obtain the variance-covariance matrix of the estimated
    # acceleration factors.
    variance <- t(partial) %*% vcov %*% partial
  },
  error = function(e) {
    warning("Jacobian matrix could not be computed.")
    return(NA)
  })
  return(list(
    estimate = gamma_est,
    variance = variance,
    partial = partial
  ))
}


#' Compute acceleration factor
#'
#' The [g_Delta_bis()] function computes the time-specific acceleration factor
#' from the estimated trajectories in both treatment groups.
#'
#' @param par Vector with estimated points of the trajectories at `time_points`.
#'   This vector is a combinations of the estimates in both treatment groups.
#'   The `length(time_points)` first elements correspond to the reference group.
#'   The `length(time_points) - 1` remaining elements correspond to the active
#'   treatment group.
#' @param time_points Time points to which the first `length(time_points)`
#'   elements in `par` correspond.
#' @inheritParams TCT_meta
#'
#' @details
#' We can write the slowing factor at \eqn{t_j} as a function of the the mean
#' parameters,
#' \deqn{
#' \boldsymbol{\gamma} = \boldsymbol{\gamma}(\boldsymbol{t}; \boldsymbol{\alpha}, \boldsymbol{\beta}).
#' }
#' The [g_Delta_bis()] function computes \eqn{\boldsymbol{\gamma}} where `par`
#' corresponds to \eqn{(\boldsymbol{\alpha}, \boldsymbol{\beta})'} and
#' `time_points` corresponds to \eqn{\boldsymbol{t}}.
#'
#'
#'
#' @return A (numeric) vector of length `length(time_points) - 1` that
#' corresponds to \eqn{\boldsymbol{\gamma}}.
g_Delta_bis = function(par,
                       interpolation,
                       time_points) {
  n_points = length(time_points)
  ctrl_estimates = par[1:n_points]
  exp_estimates = par[(n_points + 1):length(par)]

  t_mapped = sapply(
    X = exp_estimates,
    FUN = get_new_time,
    y_ref = ctrl_estimates,
    x_ref = time_points,
    method = interpolation
  )

  return(t_mapped/time_points[(n_points - length(exp_estimates) + 1):length(time_points)])
}
