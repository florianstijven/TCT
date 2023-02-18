#' Delta method
#'
#' The DeltaMethod() function applies the delta-method in the context of the
#' time-component tests. This function is based on a combination of analytical
#' (where possible) and numerical derivatives.
#'
#' @param time_points Time points to which the estimates in `y_ref` correspond.
#' @param y_ref Vector with estimated values for reference trajectory at
#'   `time_points`. This vector should thus have the same length as
#'   `time_points`.
#' @param y_obs Vector with estimated values for the active treatment group at
#'   `time_points`. This vector should thus have length `length(time_points) -
#'   1` since the baseline values are assumed to be equal in both treatment
#'   groups.
#' @param ref_fun Intra- and extrapolation function that is returned by
#'   `ref_fun_constructor()`.
#' @param method Interpolation method; see [TCT()].
#' @param vcov The (estimated) variance-covariance matrix of the parameter
#'   estimates in the `c(y_ref, y_obs)` vector.
#'
#' @return A list with three element:
#'  * `estimate`: the transformed parameter estimates
#'  * `variance`: the (estimated) variance-covariance matrix for the transformed parameter
#'   estimates
#'  * `partial`: the Jacobian matrix
#'
DeltaMethod = function (time_points,
                        y_ref,
                        y_obs,
                        ref_fun,
                        method,
                        vcov)
{
  # Vector of estimated acceleration factors
  gamma_est = g_Delta_bis(par = c(y_ref, y_obs),
                          method = method,
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
        y_ref = y_ref,
        ref_fun = ref_fun,
        method = method
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
#' @param method Interpolation method; see [TCT()].
#' @param time_points Time points to which the first `length(time_points)`
#'   elements in `par` correspond.
#'
#' @return A vector of length `length(time_points) - 1` containing the estimated
#'   acceleration factors for the estimates in the active treatment group.
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
  return(t_mapped/time_points[(n_points - length(y_obs) + 1):length(time_points)])
}

#' Parametric Bootstrap for Time-Specific Acceleration Factors
#'
#'
#' @param time_points Time points to which the estimates in `ctrl_estimates`
#'   correspond.
#' @param ctrl_estimates Vector with estimated values for reference trajectory
#'   at `time_points`. This vector should thus have the same length as
#'   `time_points`.
#' @param exp_estimates Vector with estimated values for the active treatment
#'   group at `time_points`. This vector should thus have length
#'   `length(time_points) - 1` since the baseline values are assumed to be equal
#'   in both treatment groups.
#' @param vcov The (estimated) variance-covariance matrix of the parameter
#'   estimates in the `c(ctrl_estimates, exp_estimates)` vector.
#' @param interpolation Interpolation method; see [TCT()].
#' @param B number of bootstrap replications.
#' @param null (boolean): conduct the bootstrap under the null hypothesis of no
#'   treatment effect?
#'
#' @return Matrix where each row corresponds to the bootstrap replicates of the
#'   acceleration factors.
pm_bootstrap_vertical_to_horizontal = function(time_points,
                                               ctrl_estimates,
                                               exp_estimates,
                                               vcov,
                                               interpolation = "spline",
                                               B = 100,
                                               null = FALSE) {
  if (B == 0)
    return(NULL)

  if (null) {
    # estimated covariance matrix of arm-specific vertical estimates
    vcov_gls = vcov[-1, -1]
    # "covariate" matrix for gls that maps the corresponding vertical estimates
    # in both treatment arms to the same value.
    X_gls = rbind(matrix(),
                  diag(x = 1, nrow = length(exp_estimates)),
                  diag(x = 1, nrow = length(exp_estimates)))
    # Generalized Least Squares estimate of common vertical parameters. First,
    # do a precomputation to prevent doing computations twice.
    vcov_gls[lower.tri(vcov_gls)] = t(vcov_gls)[lower.tri(vcov_gls)]
    A = t(X_gls) %*% mnormt::pd.solve(vcov_gls)
    alpha_gls = solve( A %*% X_gls ) %*% A %*% matrix(c(ctrl_estimates[-1], exp_estimates), ncol = 1)
    # The first element in ctrl_estimates is always the estimated common vertical
    # parameter at time of randomization.
    alpha_gls = c(
      ctrl_estimates[1],
      alpha_gls,
      alpha_gls
    )
    par_sampled = mvtnorm::rmvnorm(
      n = B,
      mean = alpha_gls,
      sigma = vcov
    )
  }
  else {
    par_sampled = mvtnorm::rmvnorm(
      n = B,
      mean = c(ctrl_estimates, exp_estimates),
      sigma = vcov
    )
  }

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


#' Transform vertical Treatment Effect Estimates to the Time Scale
#'
#' The [TCT()] function transforms the vertical parameter estimates to parameter
#' estimates on the time scale. These treatment effect estimates are
#' acceleration factors. This is an application of the Time-Component Test (TCT)
#' methodology.
#'
#' @param time_points Ordered vector that contains the times corresponding to
#'   the estimated means in the `ctrl_estimates` vector.
#' @param ctrl_estimates Estimated mean outcome in the control group at fixed
#'   occasions.
#' @param exp_estimates Estimated mean outcomes in the experimental group at
#'   fixed occasions.
#' @param vcov The variance-covariance matrix for the means. In order to map to
#'   the correct estimates, this matrix should be the variance-covariance matrix
#'   of `c(ctrl_means, exp_means)`.
#' @param interpolation Which interpolation method to use?
#'  * `"linear"`: linear interpolation
#'  * `"spline"`: natural cubic spline interpolation
#'  * `"monoH.FC`: monotone Hermite spline according to the method of Fritsch
#'   and Carlson
#' @param B Number of bootstrap replications. If `B = 0`, no bootstrap is
#'   performed (default).
#'
#' @return An object from the TCT-class
#' @export
#'
#' @details
#' # Time-Component Test
#' Explanation of TCT
#'
#' # Inference
#'
#' @examples
#' # transform example data set to desired format
#' library(dplyr)
#' data = simulated_test_trial %>%
#' mutate(time_int = (Week %/% 25)) %>%
#'   arrange(trial_number, SubjId, time_int) %>%
#'   mutate(time_int = as.integer(time_int) + 1L) %>%
#'   mutate(arm_time = ifelse(time_int == 1L,
#'                            "baseline",
#'                            paste0(arm, ":", time_int)))
#' # fit e.g. MMRM model to obtain estimates of profiles
#' mmrm_fit = analyze_mmrm(data)
#' set.seed(1)
#' TCT_Fit = TCT(
#'   time_points = 0:4,
#'   ctrl_estimates = mmrm_fit$coefficients[c(9, 1:4)],
#'   exp_estimates = mmrm_fit$coefficients[5:8],
#'   vcov = mmrm_fit$varBeta[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#'   interpolation = "spline",
#'   B = 1e3
#' )
TCT = function(time_points,
               ctrl_estimates,
               exp_estimates,
               vcov,
               interpolation = "spline",
               B = 0) {
  ref_fun = ref_fun_constructor(time_points,
                                ctrl_estimates,
                                interpolation)
  se_delta = DeltaMethod(
    time_points = time_points,
    y_ref = ctrl_estimates,
    y_obs = exp_estimates,
    ref_fun = ref_fun,
    method = interpolation,
    vcov = vcov
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


print.TCT = function(x, ...) {
  cat(
    paste0(
      "Time Component Test: ",
      x$type,
      "\n\n"
    )
  )
  cat("Coefficients: \n")
  print(coef(x))
  cat("\n Interpolation Method: ")
  cat(x$interpolation)
}


summary.TCT = function(x,
                       alpha = 0.05,
                       delta_transformation = "identity") {
  # inference based on delta method
  if (delta_transformation == "identity") {
    se_delta = sqrt(diag(x$vcov))
    z_delta = (1 - coef(x)) / se_delta
    ci_delta_lower = coef(x) - qnorm(1 - (alpha / 2)) * se_delta
    ci_delta_upper = coef(x) + qnorm(1 - (alpha / 2)) * se_delta
  }
  else if (delta_transformation == "log") {
    se_delta = (1 / coef(x)) * sqrt(diag(x$vcov))
    z_delta = (log(coef(x))) / se_delta
    ci_delta_lower = exp(log(coef(x)) - qnorm(1 - (alpha / 2)) * se_delta)
    ci_delta_upper = exp(log(coef(x)) + qnorm(1 - (alpha / 2)) * se_delta)
  }
  else if (delta_transformation == "log10") {
    se_delta = log10(exp(1)) * (1 / coef(x)) * sqrt(diag(x$vcov))
    z_delta = (log10(coef(x))) / se_delta
    ci_delta_lower = 10**(log10(coef(x)) - qnorm(1 - (alpha / 2)) * se_delta)
    ci_delta_upper = 10**(log10(coef(x)) + qnorm(1 - (alpha / 2)) * se_delta)
  }
  ci_delta = matrix(
    data = c(ci_delta_lower, ci_delta_upper),
    ncol = 2,
    byrow = FALSE
  )

  lht_delta = linearHypothesis.default(
    model = x,
    vcov. = x$vcov,
    coef. = coef(x),
    rhs = rep(1, length(coef(x))),
    hypothesis.matrix = diag(1, nrow = length(coef(x)), ncol = length(coef(x)))
  )
  p_delta =  (1 - pnorm(abs(z_delta))) * 2

  # inference based on parametric bootstrap
  if (!(is.null(x$bootstrap_estimates))) {
    vcov_bootstrap = var(x$bootstrap_estimates, na.rm = TRUE)
    se_bootstrap = sqrt(diag(vcov_bootstrap))
    ci_bootstrap = t(apply(
      X = x$bootstrap_estimates,
      MARGIN = 2,
      FUN = quantile,
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE
    ))
    p_bootstrap = apply(
      X = x$bootstrap_estimates,
      MARGIN = 2,
      FUN = function(x) {
        prop = mean(x > 1, na.rm = TRUE)
        return(min(prop * 2, (1 - prop) * 2))
      }
    )
  }
  else {
    vcov_bootstrap = NULL
    se_bootstrap = NULL
    ci_bootstrap = NULL
    p_bootstrap = NULL
  }

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
  if (is.null(x$ci_bootstrap)) {
    coefficients_df = data.frame(
      Value = coef(x),
      "Std. Error (delta)" = x$se_delta,
      "z-value (delta)" = x$z_delta,
      "p-value (delta)" = x$p_delta,
      "CI (delta)" = paste0("(",
                            format(x$ci_delta[, 1], digits = 5),
                            ", ",
                            format(x$ci_delta[, 2], digits = 5),
                            ")"),
      check.names = FALSE
    )
  }
  else {
    coefficients_df = data.frame(
      Value = coef(x),
      `Std. Error (delta)` = x$se_delta,
      `z-value (delta)` = x$z_delta,
      `p-value (delta)` = x$p_delta,
      `CI (delta)` = paste0("(",
                            format(x$ci_delta[, 1], digits = 5),
                            ", ",
                            format(x$ci_delta[, 2], digits = 5),
                            ")"),
      `CI (bootstrap)` = paste0("(",
                                format(x$ci_bootstrap[, 1], digits = 5),
                                ", ",
                                format(x$ci_bootstrap[, 2], digits = 5),
                                ")"),
      check.names = FALSE
    )
  }

  print(coefficients_df, digits = 5)
  cat(paste0("alpha = ", x$alpha))
  cat("\n Interpolation Method: ")
  cat(x$interpolation)
  cat("\n")
}

#' Parametric Bootstrap for a Common Acceleration Factor
#'
#' @param time_points
#' @param ctrl_estimates
#' @param exp_estimates
#' @param vcov
#' @param TCT_vcov
#' @param interpolation
#' @param B
#' @param bs_fix_vcov
#' @param return_se
#'
#' @return
pm_bootstrap_vertical_to_common = function(time_points,
                                           ctrl_estimates,
                                           exp_estimates,
                                           vcov,
                                           TCT_vcov,
                                           interpolation = "spline",
                                           B = 100,
                                           bs_fix_vcov = TRUE,
                                           return_se = TRUE,
                                           null = FALSE) {
  if (B == 0)
    return(NULL)

  n_points = length(time_points)
  p = length(exp_estimates)
  vec_1 = matrix(1, nrow = p, ncol = 1)
  estimates_bootstrap = 1:B
  se_bootstrap = rep(NA, B)
  if (null) {
    # estimated covariance matrix of arm-specific vertical estimates
    vcov_gls = vcov[-1, -1]
    # "covariate" matrix for gls that maps the corresponding vertical estimates
    # in both treatment arms to the same value.
    X_gls = rbind(diag(x = 1, nrow = length(exp_estimates)),
                  diag(x = 1, nrow = length(exp_estimates)))
    # Generalized Least Squares estimate of common vertical parameters. First,
    # do a precomputation to prevent doing computations twice.
    vcov_gls[lower.tri(vcov_gls)] = t(vcov_gls)[lower.tri(vcov_gls)]
    A = t(X_gls) %*% mnormt::pd.solve(vcov_gls)
    alpha_gls = solve( A %*% X_gls ) %*% A %*% matrix(c(ctrl_estimates[-1], exp_estimates), ncol = 1)
    # The first element in ctrl_estimates is always the estimated common vertical
    # parameter at time of randomization.
    alpha_gls = c(
      ctrl_estimates[1],
      alpha_gls,
      alpha_gls
    )
    par_sampled = mvtnorm::rmvnorm(
      n = B,
      mean = alpha_gls,
      sigma = vcov
    )
  }
  else {
    par_sampled = mvtnorm::rmvnorm(
      n = B,
      mean = c(ctrl_estimates, exp_estimates),
      sigma = vcov
    )
  }
  for (i in seq_along(estimates_bootstrap)) {
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
        vcov = vcov,
        interpolation = interpolation,
        B = 0
      )
      vcov_gls = tct_results$vcov
      coef_gls = stats::coef(tct_results)
    }
    # If vcov cannot be computed (eg, derivative is infinite) NA is returned.
    if (is.na(vcov_gls)[[1]]) {
      estimates_bootstrap[i] = NA
      if (return_se) {
        se_bootstrap[i] = sqrt((t(vec_1) %*% inv_vcov_gls %*% vec_1)**(-1))
      }
    }
    else {
      vcov_gls[lower.tri(vcov_gls)] = t(vcov_gls)[lower.tri(vcov_gls)]
      inv_vcov_gls = mnormt::pd.solve(vcov_gls)
      est_bs = (t(vec_1) %*% inv_vcov_gls %*% matrix(coef_gls, ncol = 1) ) /
        (t(vec_1) %*% inv_vcov_gls %*% vec_1)
      estimates_bootstrap[i] = est_bs
      if (return_se) {
        se_bootstrap[i] = sqrt((t(vec_1) %*% inv_vcov_gls %*% vec_1)**(-1))
      }
    }
  }
  return(list(estimates_bootstrap  = estimates_bootstrap,
              se_bootstrap = se_bootstrap))
}

#' Title
#'
#' @param TCT_Fit
#' @param B
#' @param bs_fix_vcov
#'
#' @return
#' @export
#'
#' @examples
TCT_common = function(TCT_Fit,
                      B = 0,
                      bs_fix_vcov = FALSE,
                      select_coef = 1:length(coef(TCT_Fit)),
                      null_bs = FALSE) {

  estimates = coef(TCT_Fit)[select_coef]
  vcov = TCT_Fit$vcov[select_coef, select_coef]
  n_points = length(TCT_Fit$vertical_model$time_points)

  # delta method
  p = length(estimates)
  vec_1 = matrix(1, nrow = p, ncol = 1)
  est_delta = (t(vec_1) %*% solve(vcov) %*% matrix(estimates, ncol = 1) ) /
    (t(vec_1) %*% solve(vcov) %*% vec_1)
  vcov_delta = (t(vec_1) %*% solve(vcov) %*% vec_1)**(-1)


  bs_estimates = pm_bootstrap_vertical_to_common(
    time_points = TCT_Fit$vertical_model$time_points,
    ctrl_estimates = TCT_Fit$vertical_model$ctrl_estimates,
    exp_estimates = TCT_Fit$vertical_model$exp_estimates[select_coef],
    vcov = TCT_Fit$vertical_model$vcov[c(1:n_points, n_points + select_coef), c(1:n_points, n_points + select_coef)],
    TCT_vcov = vcov,
    interpolation = TCT_Fit$interpolation,
    B = B,
    bs_fix_vcov = bs_fix_vcov,
    return_se = TRUE
  )

  bs_estimates_null = NULL
  if (null_bs) {
    bs_estimates_null = pm_bootstrap_vertical_to_common(
      time_points = TCT_Fit$vertical_model$time_points,
      ctrl_estimates = TCT_Fit$vertical_model$ctrl_estimates,
      exp_estimates = TCT_Fit$vertical_model$exp_estimates[select_coef],
      vcov = TCT_Fit$vertical_model$vcov[c(1:n_points, n_points + select_coef), c(1:n_points, n_points + select_coef)],
      TCT_vcov = vcov,
      interpolation = TCT_Fit$interpolation,
      B = B,
      bs_fix_vcov = bs_fix_vcov,
      return_se = TRUE,
      null = TRUE
    )
  }

  # Test for common slowing parameter
  lht_matrix = matrix(0, nrow = length(estimates) - 1, ncol = length(estimates) - 1)
  diag(lht_matrix) = -1
  lht_matrix = cbind(1, lht_matrix)

  lht_common = linearHypothesis.default(
    model = TCT_Fit,
    vcov. = vcov,
    coef. = estimates,
    rhs = rep(0, length(estimates) - 1),
    hypothesis.matrix = lht_matrix
  )


  new_TCT_common(
    coefficients = est_delta,
    vcov = vcov_delta,
    bootstrap_estimates = bs_estimates,
    bootstrap_estimates_null = bs_estimates_null,
    interpolation = TCT_Fit$interpolation,
    lht_common = lht_common
  )
}

new_TCT_common = function(coefficients,
                          vcov,
                          bootstrap_estimates,
                          bootstrap_estimates_null,
                          interpolation,
                          lht_common
                          ) {
  structure(
    list(
      coefficients = coefficients,
      vcov = vcov,
      bootstrap_estimates = bootstrap_estimates,
      bootstrap_estimates_null = bootstrap_estimates_null,
      interpolation = interpolation,
      lht_common = lht_common
    ),
    class = "TCT_common"
  )
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
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
  cat("\n Interpolation Method: ")
  cat(x$interpolation)
  cat("\n")

}

#' Title
#'
#' @param x
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
summary.TCT_common = function(x,
                              alpha = 0.05,
                              delta_transformation = "identity") {
  # inference based on delta method

  if (delta_transformation == "identity") {
    se_delta = sqrt(diag(x$vcov))
    z_delta = (1 - coef(x)) / se_delta
    ci_delta_lower = coef(x) - qnorm(1 - (alpha / 2)) * se_delta
    ci_delta_upper = coef(x) + qnorm(1 - (alpha / 2)) * se_delta
  }
  else if (delta_transformation == "log") {
    se_delta = (1 / coef(x)) * sqrt(diag(x$vcov))
    z_delta = (log(coef(x))) / se_delta
    ci_delta_lower = exp(log(coef(x)) - qnorm(1 - (alpha / 2)) * se_delta)
    ci_delta_upper = exp(log(coef(x)) + qnorm(1 - (alpha / 2)) * se_delta)
  }
  else if (delta_transformation == "log10") {
    se_delta = log10(exp(1)) * (1 / coef(x)) * sqrt(diag(x$vcov))
    z_delta = log10(coef(x)) / se_delta
    ci_delta_lower = 10**(log10(coef(x)) - qnorm(1 - (alpha / 2)) * se_delta)
    ci_delta_upper = 10**(log10(coef(x)) + qnorm(1 - (alpha / 2)) * se_delta)
  }
  ci_delta = matrix(
    data = c(ci_delta_lower, ci_delta_upper),
    ncol = 2,
    byrow = FALSE
  )
  p_delta =  (1 - pnorm(abs(z_delta))) * 2

  # inference based on parametric bootstrap
  if (!(is.null(x$bootstrap_estimates))) {
    vcov_bootstrap = var(x$bootstrap_estimates[[1]], na.rm = TRUE)
    se_bootstrap = sqrt(vcov_bootstrap)
    ci_bootstrap = quantile(
      x = x$bootstrap_estimates[[1]],
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE
    )
    p_bootstrap = min(mean(x$bootstrap_estimates[[1]] > 1, na.rm = TRUE) * 2,
                      (1 - mean(x$bootstrap_estimates[[1]] > 1, na.rm = TRUE)) * 2)
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

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
print.summary.TCT_common = function(x) {
  cat(
    paste0(
      "Time Component Test: ",
      "proportional slowing",
      "\n\n"
    )
  )
  cat("Coefficients: \n")
  if (is.null(x$vcov_bootstrap)) {
    coefficients_df = data.frame(
      Value = coef(x),
      `Std. Error (delta)` = x$se_delta,
      `z-value (delta)` = x$z_delta,
      `p-value (delta)` = x$p_delta,
      `CI (delta)` = paste0("(",
                            format(x$ci_delta[1], digits = 5),
                            ", ",
                            format(x$ci_delta[2], digits = 5),
                            ")"),
      check.names = FALSE
    )
  }
  else {
    coefficients_df = data.frame(
      Value = coef(x),
      `Std. Error (delta)` = x$se_delta,
      `z-value (delta)` = x$z_delta,
      `p-value (delta)` = x$p_delta,
      `CI (delta)` = paste0("(",
                            format(x$ci_delta[1], digits = 5),
                            ", ",
                            format(x$ci_delta[2], digits = 5),
                            ")"),
      `CI (bootstrap)` = paste0("(",
                                format(x$ci_bootstrap[1], digits = 5),
                                ", ",
                                format(x$ci_bootstrap[2], digits = 5),
                                ")"),
      `p-value (bootstrap)` = x$p_bootstrap[1],
      check.names = FALSE
    )
  }
  print(coefficients_df, digits = 5)
  cat(paste0("alpha = ", x$alpha))
  cat("\n\n")
  cat("Test for proportional slowing factor:\n")
  print(
    data.frame("Df" = x$lht_common$Df[2],
               "Chisq" = x$lht_common$Chisq[2],
               "p-value" = x$lht_common$`Pr(>Chisq)`[2]),
    digits = 5
  )
}





