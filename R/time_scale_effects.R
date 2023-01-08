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
# delta_vertical_to_horizontal = function(time_points,
#                                         ctrl_means,
#                                         exp_means,
#                                         vcov,
#                                         interpolation = "spline"
#                                         ) {
#   # number of time points
#   n_points = length(time_points)
#   # rewrite get_new_time() for the DeltaMethod() function
#   g_Delta = function(par, x_ref, method = interpolation) {
#     y_ref = par[1:n_points]
#     y_obs = par[(n_points + 1):length(par)]
#     t_mapped = sapply(
#       X = y_obs,
#       FUN = get_new_time,
#       y_ref = y_ref,
#       x_ref = x_ref,
#       method = method
#     )
#     return(t_mapped/time_points[-1])
#   }
#   se_delta = DeltaMethod(
#     par = c(ctrl_means, exp_means),
#     fct = g_Delta,
#     vcov = vcov,
#     x_ref = time_points,
#     method = interpolation
#   )
#   return(se_delta)
# }

g_Delta_bis = function(par,
                       method,
                       time_points) {
  n_points = length(time_points)
  y_ref = par[1:n_points]
  y_obs = par[(n_points + 1):length(par)]
  t_mapped = sapply(
    X = y_obs,
    FUN = get_new_time,
    y_ref = y_ref,
    x_ref = time_points,
    method = method
  )
  return(t_mapped/time_points[-1])
}

#' Title
#'
#' @param time_points
#' @param ctrl_means
#' @param exp_means
#' @param vcov
#' @param interpolation
#' @param B
#'
#' @return
#'
#' @examples
pm_bootstrap_vertical_to_horizontal = function(time_points,
                                               ctrl_estimates,
                                               exp_estimates,
                                               vcov,
                                               interpolation = "spline",
                                               B = 100) {
  if (B == 0)
    return(NULL)

  par_sampled = mvtnorm::rmvnorm(
    n = B,
    mean = c(ctrl_estimates, exp_estimates),
    sigma = vcov
  )
  estimates = matrix(0, nrow = B, ncol = 4)
  for (i in 1:B) {
    estimates[i, ] = g_Delta_bis(par = par_sampled[i, ],
                                 method = interpolation,
                                 time_points = time_points)
  }
  return(estimates)
}

new_TCT = function(coefficients,
                   vcov,
                   bootstrap_estimates,
                   interpolation,
                   type,
                   vertical_model) {

  structure(list(coefficients = coefficients,
                 vcov = vcov,
                 bootstrap_estimates = bootstrap_estimates,
                 interpolation = interpolation,
                 type = type,
                 vertical_model = vertical_model),
            class = "TCT")
}

TCT = function(time_points,
               ctrl_estimates,
               exp_estimates,
               vcov,
               interpolation = "spline",
               B = 0) {
  se_delta = DeltaMethod(
    par = c(ctrl_estimates, exp_estimates),
    fct = g_Delta_bis,
    vcov = vcov,
    time_points = time_points,
    method = interpolation
  )

  bootstrap_estimates = pm_bootstrap_vertical_to_horizontal(
    time_points = time_points,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov,
    interpolation = interpolation,
    B = B
  )

  return(new_TCT(coefficients = se_delta$estimate,
                 vcov = se_delta$variance,
                 bootstrap_estimates = bootstrap_estimates,
                 interpolation = interpolation,
                 type = "time-based treatment effects",
                 vertical_model = list(time_points = time_points,
                                       ctrl_estimates = ctrl_estimates,
                                       exp_estimates = exp_estimates,
                                       vcov = vcov)
                 ))
}

print.TCT = function(x) {
  cat(
    paste0(
      "Time Component Test: ",
      x$type,
      "\n\n"
    )
  )
  cat("Coefficients: \n")
  print(coef(x))
}

summary.TCT = function(x,
                       alpha = 0.05) {
  # inference based on delta method
  se_delta = sqrt(diag(x$vcov))
  z_delta = (1 - coef(x)) / se_delta
  ci_delta_lower = coef(x) - qnorm(1 - (alpha / 2)) * se_delta
  ci_delta_upper = coef(x) + qnorm(1 - (alpha / 2)) * se_delta
  ci_delta = matrix(
    data = c(ci_delta_lower, ci_delta_upper),
    ncol = 2,
    byrow = FALSE
  )
  lht_delta = car::linearHypothesis(
    model = x,
    vcov. = x$vcov,
    coef. = coef(x),
    rhs = rep(1, length(coef(x))),
    hypothesis.matrix = diag(1, nrow = length(coef(x)), ncol = length(coef(x)))
  )
  p_delta =  (1 - pnorm(abs(z_delta))) * 2

  # inference based on parametric bootstrap
  vcov_bootstrap = var(x$bootstrap_estimates)
  se_bootstrap = sqrt(diag(vcov_bootstrap))
  ci_bootstrap = t(apply(
    X = x$bootstrap_estimates,
    MARGIN = 2,
    FUN = quantile,
    probs = c(alpha / 2, 1 - alpha / 2)
  ))
  p_bootstrap = apply(
    X = x$bootstrap_estimates,
    MARGIN = 2,
    FUN = function(x) {
      prop = mean(x > 1)
      return(min(prop * 2, (1 - prop) * 2))
    }
  )

  new_summary.TCT(
    x = x,
    se_delta = se_delta,
    z_delta = z_delta,
    ci_delta = ci_delta,
    p_delta = p_delta,
    lht_delta = lht_delta,
    vcov_bootstrap = vcov_bootstrap,
    se_bootstrap = se_bootstrap,
    ci_bootstrap = ci_bootstrap,
    p_bootstrap = p_bootstrap,
    alpha = alpha
  )
}

new_summary.TCT = function(
    x,
    se_delta,
    z_delta,
    ci_delta,
    p_delta,
    lht_delta,
    vcov_bootstrap,
    se_bootstrap,
    ci_bootstrap,
    p_bootstrap,
    alpha
) {
  structure(append(
    x,
    list(
      se_delta = se_delta,
      z_delta = z_delta,
      ci_delta = ci_delta,
      p_delta = p_delta,
      lht_delta = lht_delta,
      vcov_bootstrap = vcov_bootstrap,
      se_bootstrap = se_bootstrap,
      ci_bootstrap = ci_bootstrap,
      p_bootstrap = p_bootstrap,
      alpha = alpha
    )
  ),
  class = "summary.TCT")
}

print.summary.TCT = function(x) {
  cat(
    paste0(
      "Time Component Test: ",
      x$type,
      "\n\n"
    )
  )
  cat("Coefficients: \n")
  coefficients_df = data.frame(
    Value = coef(x),
    `Std. Error (delta)` = x$se_delta,
    `z-value (delta)` = x$z_delta,
    `p-value (delta)` = x$p_delta,
    `CI lower (delta)` = x$ci_delta[, 1],
    `CI upper (delta)` = x$ci_delta[, 2],
    `CI lower (bootstrap)` = x$ci_bootstrap[, 1],
    `CI upper (bootstrap)` = x$ci_bootstrap[, 2]
  )
  print(coefficients_df)
  cat(paste0("alpha = ", x$alpha))
}

pm_bootstrap_vertical_to_common = function(time_points,
                                           ctrl_estimates,
                                           exp_estimates,
                                           vcov,
                                           TCT_vcov,
                                           interpolation = "spline",
                                           B = 100,
                                           bs_fix_vcov = TRUE) {
  if (B == 0)
    return(NULL)

  p = length(exp_estimates)
  vec_1 = matrix(1, nrow = p, ncol = 1)
  estimates_bootstrap = 1:B
  ctrl_estimates = ctrl_estimates
  exp_estimates = exp_estimates
  time_points = time_points
  par_sampled = mvtnorm::rmvnorm(n = B,
                                 mean = c(ctrl_estimates, exp_estimates),
                                 sigma = vcov)
  for (i in 1:B) {
    if (bs_fix_vcov) {
      vcov_gls = TCT_vcov
      coef_gls = g_Delta_bis(par = par_sampled[i, ],
                             time_points = time_points,
                             method = interpolation)
    }
    else {
      tct_results = TCT(
        time_points = time_points,
        ctrl_estimates = par_sampled[i, 1:n_points],
        exp_estimates = par_sampled[i, (n_points + 1):length(par_sampled[1, ])],
        vcov = TCT_Fit$vertical_model$vcov,
        interpolation = interpolation,
        B = 0
      )
      vcov_gls = tct_results$vcov
      coef_gls = coef(tct_results)
    }

    est_bs = (t(vec_1) %*% solve(vcov_gls) %*% matrix(coef_gls, ncol = 1) ) /
      (t(vec_1) %*% solve(vcov_gls) %*% vec_1)
    estimates_bootstrap[i] = est_bs
  }
  return(estimates_bootstrap)
}

TCT_common = function(TCT_Fit,
                      B = 0,
                      bs_fix_vcov = FALSE) {
  estimates = coef(TCT_Fit)
  vcov = TCT_Fit$vcov
  n_points = length(TCT_Fit$vertical_model$time_points)

  # delta method
  p = length(estimates)
  vec_1 = matrix(1, nrow = p, ncol = 1)
  est_delta = (t(vec_1) %*% solve(vcov) %*% matrix(estimates, ncol = 1) ) /
    (t(vec_1) %*% solve(vcov) %*% vec_1)
  vcov_delta = (t(vec_1) %*% solve(vcov) %*% vec_1)**(-1)

  # bootstrap
  # if (B > 0) {
  #   estimates_bootstrap = 1:B
  #   ctrl_estimates = TCT_Fit$vertical_model$ctrl_estimates
  #   exp_estimates = TCT_Fit$vertical_model$exp_estimates
  #   time_points = TCT_Fit$vertical_model$time_points
  #   par_sampled = mvtnorm::rmvnorm(n = B,
  #                                  mean = c(ctrl_estimates, exp_estimates),
  #                                  sigma = TCT_Fit$vertical_model$vcov)
  #   for (i in 1:B) {
  #     if (bs_fix_vcov) {
  #       vcov_gls = TCT_Fit$vcov
  #       coef_gls = g_Delta_bis(par = par_sampled[i, ],
  #                              time_points = time_points,
  #                              method = "spline")
  #     }
  #     else {
  #       tct_results = TCT(
  #         time_points = time_points,
  #         ctrl_estimates = par_sampled[i, 1:n_points],
  #         exp_estimates = par_sampled[i, (n_points + 1):length(par_sampled[1, ])],
  #         vcov = TCT_Fit$vertical_model$vcov,
  #         interpolation = "spline",
  #         B = 0
  #       )
  #       vcov_gls = tct_results$vcov
  #       coef_gls = coef(tct_results)
  #     }
  #
  #     est_bs = (t(vec_1) %*% solve(vcov_gls) %*% matrix(coef_gls, ncol = 1) ) /
  #       (t(vec_1) %*% solve(vcov_gls) %*% vec_1)
  #     estimates[i] = est_bs
  #   }
  # }
  # else {
  #   estimates = NULL
  # }
  estimates = pm_bootstrap_vertical_to_common(time_points = TCT_Fit$vertical_model$time_points,
                                              ctrl_estimates = TCT_Fit$vertical_model$ctrl_estimates,
                                              exp_estimates = TCT_Fit$vertical_model$exp_estimates,
                                              vcov = TCT_Fit$vertical_model$vcov,
                                              TCT_vcov = TCT_Fit$vcov,
                                              interpolation = TCT_Fit$interpolation,
                                              B = B
                                              )


  new_TCT_common(
    coefficients = est_delta,
    vcov = vcov_delta,
    bootstrap_estimates = estimates
  )
}

new_TCT_common = function(coefficients,
                          vcov,
                          bootstrap_estimates
                          ) {
  structure(
    list(
      coefficients = coefficients,
      vcov = vcov,
      bootstrap_estimates = bootstrap_estimates
    ),
    class = "TCT_common"
  )
}

print.TCT_common = function(x) {
  cat(
    paste0(
      "Time Component Test: ",
      "proportional slowing",
      "\n\n"
    )
  )
  cat("Coefficients: \n")
  print(x$coefficients)
}

summary.TCT_common = function(x,
                              alpha = 0.05) {
  # inference based on delta method
  se_delta = sqrt(diag(x$vcov))
  z_delta = (1 - coef(x)) / se_delta
  ci_delta_lower = coef(x) - qnorm(1 - (alpha / 2)) * se_delta
  ci_delta_upper = coef(x) + qnorm(1 - (alpha / 2)) * se_delta
  ci_delta = matrix(
    data = c(ci_delta_lower, ci_delta_upper),
    ncol = 2,
    byrow = FALSE
  )
  p_delta =  (1 - pnorm(abs(z_delta))) * 2

  # inference based on parametric bootstrap
  if (!(is.null(x$bootstrap_estimates))) {
    vcov_bootstrap = var(x$bootstrap_estimates)
    se_bootstrap = sqrt(vcov_bootstrap)
    ci_bootstrap = quantile(
      x = x$bootstrap_estimates,
      probs = c(alpha / 2, 1 - alpha / 2)
    )
    p_bootstrap = min(mean(x$bootstrap_estimates > 1) * 2,
                      (1 - mean(x$bootstrap_estimates > 1)) * 2)
  }
  else {
    vcov_bootstrap = NULL
    se_bootstrap = NULL
    ci_bootstrap = NULL
    p_bootstrap = NULL
  }


  new_summary.TCT_common(
    x = x,
    se_delta = se_delta,
    z_delta = z_delta,
    ci_delta = ci_delta,
    p_delta = p_delta,
    vcov_bootstrap = vcov_bootstrap,
    se_bootstrap = se_bootstrap,
    ci_bootstrap = ci_bootstrap,
    p_bootstrap = p_bootstrap,
    alpha = alpha
  )
}

new_summary.TCT_common = function(
    x,
    se_delta,
    z_delta,
    ci_delta,
    p_delta,
    vcov_bootstrap,
    se_bootstrap,
    ci_bootstrap,
    p_bootstrap,
    alpha
) {
  structure(append(
    x,
    list(
      se_delta = se_delta,
      z_delta = z_delta,
      ci_delta = ci_delta,
      p_delta = p_delta,
      vcov_bootstrap = vcov_bootstrap,
      se_bootstrap = se_bootstrap,
      ci_bootstrap = ci_bootstrap,
      p_bootstrap = p_bootstrap,
      alpha = alpha
    )
  ),
  class = "summary.TCT_common")
}

print.summary.TCT_common = function(x) {
  cat(
    paste0(
      "Time Component Test: ",
      "proportional slowing",
      "\n\n"
    )
  )
  cat("Coefficients: \n")
  coefficients_df = data.frame(
    Value = coef(x),
    `Std. Error (delta)` = x$se_delta,
    `z-value (delta)` = x$z_delta,
    `p-value (delta)` = x$p_delta,
    `CI lower (delta)` = x$ci_delta[, 1],
    `CI upper (delta)` = x$ci_delta[, 2],
    `CI lower (bootstrap)` = x$ci_bootstrap[1],
    `CI upper (bootstrap)` = x$ci_bootstrap[2],
    `p-value (bootstrap)` = x$p_bootstrap[1]
  )
  print(coefficients_df)
  cat(paste0("alpha = ", x$alpha))
}




