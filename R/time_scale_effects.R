new_TCT = function(coefficients,
                   vcov,
                   inference,
                   bootstrap_estimates,
                   interpolation,
                   type,
                   vertical_model) {

  structure(list(coefficients = coefficients,
                 vcov = vcov,
                 bootstrap_estimates = bootstrap_estimates,
                 interpolation = interpolation,
                 type = type,
                 vertical_model = vertical_model,
                 inference = inference),
            class = "TCT")
}


#' Transform vertical Treatment Effect Estimates to the Time Scale
#'
#' The [TCT()] function transforms so-called vertical parameter estimates to
#' parameter estimates on the time scale. These treatment effect estimates on
#' the time scale are so-called acceleration factors. This is an application of
#' the Time-Component Test (TCT) methodology.
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
#' @param inference Which approach is used for estimation and inference? Should
#'   be `"wald"` or `"score"`.
#' @param interpolation Which interpolation method to use?
#'  * `"linear"`: linear interpolation
#'  * `"spline"`: natural cubic spline interpolation
#'  * `"monoH.FC`: monotone Hermite spline according to the method of Fritsch
#'   and Carlson
#' @param B Number of bootstrap replications. If `B = 0`, no bootstrap is
#'   performed (default).
#' @param constraints Use the constrained generalized least squares estimator
#'   for the vertical treatment effects.
#' @return An object from the TCT-class
#' @export
#'
#' @details
#'
#' # Time-Component Tests
#'
#' Time-component tests (TCT) constitutes a general methodology to evaluating
#' treatment effects with longitudinal data on the time scale. Conventional
#' treatment effects with longitudinal data are so-called vertical treatment
#' effects; these are comparisons of group means at fixed measurement occasions.
#'
#' Let \eqn{\boldsymbol{t} = (t_0 = 0, t_1, ..., t_K)'} be the fixed measurement
#' occasions (`timepoints` in this function). Let \eqn{\boldsymbol{\alpha} =
#' (\alpha_0, \alpha_1, ..., \alpha_K)'} be the corresponding means in the
#' control group. Let \eqn{\boldsymbol{\beta} = (\beta_1, ..., \beta_{K})} be
#' the corresponding means in the experimental group. Note that the index starts
#' here at 1, i.e., the first measurement *after* start of the treatment. Let
#' the mean trajectory in the control and experimental group be, respectively,
#' \eqn{f_0(t; \boldsymbol{\alpha})} and \eqn{f_1(t; \boldsymbol{\beta})}
#'
#' The treatment effects on the time scale are acceleration factors, analogous
#' to accelerated failure time models. These are defined as follows at
#' \eqn{t_j}, \deqn{f_1(t; \boldsymbol{\beta}) = f_0(\gamma_j \cdot t;
#' \boldsymbol{\alpha})}
#' where \eqn{\gamma_j} is the so-called acceleration factor at \eqn{t_j}, i.e.,
#' treatment causes a acceleration of \eqn{\gamma_j}. For example, if
#' \eqn{\gamma_j = 0.5}, patients in the active treatment group progress half as
#' slow as patients in the control group.
#'
#' # Inference
#'
#' Following options for inference are available:
#' * Delta method
#' * Parametric bootstrap
#' * Score test: More information in [score_test()] and [score_test_common()]
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
               inference = "wald",
               B = 0,
               constraints = FALSE) {
  # Apply constraint GLS estimator if asked.
  if (constraints) {
    contrained_estimates = constrained_vertical_estimator(ctrl_estimates,
                                                          exp_estimates,
                                                          vcov)
    ctrl_estimates = contrained_estimates[1:length(ctrl_estimates)]
    exp_estimates = contrained_estimates[(length(ctrl_estimates) + 1):(length(ctrl_estimates) + length(exp_estimates))]
  }

  ref_fun = ref_fun_constructor(time_points,
                                ctrl_estimates,
                                interpolation)

  se_delta = DeltaMethod(
    time_points = time_points,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    ref_fun = ref_fun,
    interpolation = interpolation,
    vcov = vcov
  )
  coefficients = se_delta$estimate
  vcov_matrix = se_delta$variance



  bootstrap_estimates = pm_bootstrap_vertical_to_horizontal(
    time_points = time_points,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    vcov = vcov,
    interpolation = interpolation,
    B = B,
    constraints = constraints
  )

  return(new_TCT(coefficients = coefficients,
                 vcov = vcov_matrix,
                 inference = inference,
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
      "Time Component Test ",
      "(", x$inference, ")",
      x$type,
      "\n\n"
    )
  )
  cat("Coefficients: \n")
  print(coef(x))
  cat("\n Interpolation Method: ")
  cat(x$interpolation)
}


#' Summarize fitted Time Component Test model
#'
#' @param x Object returned by [TCT()].
#' @param alpha Two-sided confidence level for confidence intervals.
#' @param delta_transformation Transformation when applying the delta-method to
#'   obtain confidence intervals.
#'
#' @return
#' @export
#' @import stats
summary.TCT = function(x,
                       alpha = 0.05,
                       delta_transformation = "identity") {
  # Wald-based inference if the inference option is equal to "wald".
  if (x$inference == "wald") {
    if (delta_transformation == "identity") {
      se_delta = sqrt(diag(x$vcov))
      z_values = (1 - coef(x)) / se_delta
      ci_delta_lower = coef(x) - stats::qnorm(1 - (alpha / 2)) * se_delta
      ci_delta_upper = coef(x) + stats::qnorm(1 - (alpha / 2)) * se_delta
    }
    else if (delta_transformation == "log") {
      se_delta = (1 / coef(x)) * sqrt(diag(x$vcov))
      z_values = (log(coef(x))) / se_delta
      ci_delta_lower = exp(log(coef(x)) - stats::qnorm(1 - (alpha / 2)) * se_delta)
      ci_delta_upper = exp(log(coef(x)) + stats::qnorm(1 - (alpha / 2)) * se_delta)
    }
    else if (delta_transformation == "log10") {
      se_delta = log10(exp(1)) * (1 / coef(x)) * sqrt(diag(x$vcov))
      z_values = (log10(coef(x))) / se_delta
      ci_delta_lower = 10**(log10(coef(x)) - stats::qnorm(1 - (alpha / 2)) * se_delta)
      ci_delta_upper = 10**(log10(coef(x)) + stats::qnorm(1 - (alpha / 2)) * se_delta)
    }
    ci_matrix = matrix(
      data = c(ci_delta_lower, ci_delta_upper),
      ncol = 2,
      byrow = FALSE
    )
  }
  # Score-based inference if the inference option is equal to "score".
  if (x$inference == "score") {
    se_delta = sqrt(diag(x$vcov))
    # (Re)construct reference trajectory.
    ref_fun = ref_fun_constructor(
      x$vertical_model$time_points,
      x$vertical_model$ctrl_estimates,
      x$interpolation
    )
    # Compute confidence intervals for the measurement occasion-specific
    # acceleration factors based on the score test.
    ci_list = lapply(
      X = 1:length(x$vertical_model$exp_estimates),
      FUN = function(j) {
        score_conf_int(
          time_points = x$vertical_model$time_points,
          ctrl_estimates = x$vertical_model$ctrl_estimates,
          exp_estimates = x$vertical_model$exp_estimates,
          ref_fun = ref_fun,
          interpolation = x$interpolation,
          vcov = x$vertical_model$vcov,
          j = j,
          alpha = alpha
        )
      }
    )
    # Convert list of confidence intervals to a matrix. The first column
    # contains the lower limits, the second column contains the upper limits.
    ci_matrix = matrix(unlist(ci_list), ncol = 2, byrow = TRUE)
    # Compute the z-statistics for the measurement occasion-specific
    # acceleration factors based on the score test.
    z_values = vapply(
      X = 1:length(x$vertical_model$exp_estimates),
      FUN = function(j) {
        score_test(
          time_points = x$vertical_model$time_points,
          ctrl_estimates = x$vertical_model$ctrl_estimates,
          exp_estimates = x$vertical_model$exp_estimates,
          ref_fun = ref_fun,
          interpolation = x$interpolation,
          vcov = x$vertical_model$vcov,
          j = j,
          gamma_0 = 1
        )[1]
      },
      FUN.VALUE = 1.1
    )
  }

  lht_delta = linearHypothesis.default(
    model = x,
    vcov. = x$vcov,
    coef. = coef(x),
    rhs = rep(1, length(coef(x))),
    hypothesis.matrix = diag(1, nrow = length(coef(x)), ncol = length(coef(x)))
  )
  p_values =  (1 - pnorm(abs(z_values))) * 2

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


  return(
    new_summary.TCT(
      x = x,
      se_delta = se_delta,
      z_values = z_values,
      ci_matrix = ci_matrix,
      p_values = p_values,
      lht_delta = lht_delta,
      vcov_bootstrap = vcov_bootstrap,
      se_bootstrap = se_bootstrap,
      ci_bootstrap = ci_bootstrap,
      p_bootstrap = p_bootstrap,
      alpha = alpha
    )
  )
}

new_summary.TCT = function(
    x,
    se_delta,
    z_values,
    ci_matrix,
    p_values,
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
      z_values = z_values,
      ci_matrix = ci_matrix,
      p_values = p_values,
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

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
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
      "z-value (delta)" = x$z_values,
      "p-value (delta)" = x$p_values,
      "CI (delta)" = paste0("(",
                            format(x$ci_matrix[, 1], digits = 5),
                            ", ",
                            format(x$ci_matrix[, 2], digits = 5),
                            ")"),
      check.names = FALSE
    )
  }
  else {
    coefficients_df = data.frame(
      Value = coef(x),
      `Std. Error (delta)` = x$se_delta,
      `z-value (delta)` = x$z_values,
      `p-value (delta)` = x$p_values,
      `CI (delta)` = paste0("(",
                            format(x$ci_matrix[, 1], digits = 5),
                            ", ",
                            format(x$ci_matrix[, 2], digits = 5),
                            ")"),
      `CI (bootstrap)` = paste0("(",
                                format(x$ci_bootstrap[, 1], digits = 5),
                                ", ",
                                format(x$ci_bootstrap[, 2], digits = 5),
                                ")"),
      check.names = FALSE
    )

  }
  # Change column names if inference is based on the score test.
  if (x$inference == "score") {
    colnames(coefficients_df)[2:5] = c("Std. Error (delta)",
                                  "z-value (score)",
                                  "p-value (score)",
                                  "CI (score)")
  }

  print(coefficients_df, digits = 5)
  cat(paste0("alpha = ", x$alpha))
  cat("\n Interpolation Method: ")
  cat(x$interpolation)
  cat("\n")
}

#' Estimate Common Acceleration Factor
#'
#' The [TCT_common()] function estimates the common acceleration factor.
#'
#' @param TCT_Fit
#' @param bs_fix_vcov
#'
#' @return
#' @export
#'
#' @examples
TCT_common = function(TCT_Fit,
                      inference = "wald",
                      B = 0,
                      bs_fix_vcov = FALSE,
                      select_coef = 1:length(coef(TCT_Fit)),
                      null_bs = FALSE,
                      gls_est = TRUE,
                      constraints = FALSE,
                      type = NULL,
                      weights = NULL) {
  # (Re)construct reference trajectory.
  ref_fun = ref_fun_constructor(
    TCT_Fit$vertical_model$time_points,
    TCT_Fit$vertical_model$ctrl_estimates,
    TCT_Fit$interpolation
  )
  # Use wald-based inference if this is asked by the user.
  estimates = coef(TCT_Fit)[select_coef]
  if (gls_est) {
    vcov = TCT_Fit$vcov[select_coef, select_coef]
  }
  else {
    vcov = diag(diag(TCT_Fit$vcov[select_coef, select_coef]))
  }
  if (inference == "wald") {
    n_points = length(TCT_Fit$vertical_model$time_points)

    # delta method
    p = length(estimates)
    vec_1 = matrix(1, nrow = p, ncol = 1)
    gamma_common_estimate = (t(vec_1) %*% solve(vcov) %*% matrix(estimates, ncol = 1) ) /
      (t(vec_1) %*% solve(vcov) %*% vec_1)
    gamma_common_vcov = (t(vec_1) %*% solve(vcov) %*% vec_1)**(-1)
  }
  else if (inference == "score") {
    gamma_common_estimate = score_estimate_common(
      time_points = TCT_Fit$vertical_model$time_points,
      ctrl_estimates = TCT_Fit$vertical_model$ctrl_estimates,
      exp_estimates = TCT_Fit$vertical_model$exp_estimates,
      ref_fun = ref_fun,
      interpolation = TCT_Fit$interpolation,
      vcov = TCT_Fit$vertical_model$vcov,
      type = type,
      j = select_coef,
      weights = weights
    )
    gamma_common_vcov = NA
  }




  bs_estimates = pm_bootstrap_vertical_to_common(
    time_points = TCT_Fit$vertical_model$time_points,
    ctrl_estimates = TCT_Fit$vertical_model$ctrl_estimates,
    exp_estimates = TCT_Fit$vertical_model$exp_estimates,
    vcov = TCT_Fit$vertical_model$vcov[c(1:n_points, n_points + 1:length(coef(TCT_Fit))),
                                       c(1:n_points, n_points + 1:length(coef(TCT_Fit)))],
    TCT_vcov = vcov,
    interpolation = TCT_Fit$interpolation,
    B = B,
    bs_fix_vcov = bs_fix_vcov,
    return_se = TRUE,
    gls_est = gls_est,
    select_coef = select_coef,
    constraints = constraints
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
    coefficients = gamma_common_estimate,
    inference = inference,
    vcov = gamma_common_vcov,
    bootstrap_estimates = bs_estimates,
    bootstrap_estimates_null = bs_estimates_null,
    interpolation = TCT_Fit$interpolation,
    lht_common = lht_common,
    type = type,
    weights = weights,
    vertical_model = TCT_Fit$vertical_model,
    select_coef = select_coef
  )
}

new_TCT_common = function(coefficients,
                          inference,
                          vcov,
                          bootstrap_estimates,
                          bootstrap_estimates_null,
                          interpolation,
                          lht_common,
                          type,
                          weights,
                          vertical_model,
                          select_coef
                          ) {
  structure(
    list(
      coefficients = coefficients,
      inference = inference,
      vcov = vcov,
      bootstrap_estimates = bootstrap_estimates,
      bootstrap_estimates_null = bootstrap_estimates_null,
      interpolation = interpolation,
      lht_common = lht_common,
      type = type,
      weights = weights,
      vertical_model = vertical_model,
      select_coef = select_coef
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
  # Wald-based inference
  if (x$inference == "wald") {
    if (delta_transformation == "identity") {
      gamma_common_se = sqrt(diag(x$vcov))
      z_value = (1 - coef(x)) / gamma_common_se
      gamma_ci_lower = coef(x) - stats::qnorm(1 - (alpha / 2)) * gamma_common_se
      gamma_ci_upper = coef(x) + stats::qnorm(1 - (alpha / 2)) * gamma_common_se
    }
    else if (delta_transformation == "log") {
      gamma_common_se = (1 / coef(x)) * sqrt(diag(x$vcov))
      z_value = (log(coef(x))) / gamma_common_se
      gamma_ci_lower = exp(log(coef(x)) - stats::qnorm(1 - (alpha / 2)) * gamma_common_se)
      gamma_ci_upper = exp(log(coef(x)) + stats::qnorm(1 - (alpha / 2)) * gamma_common_se)
    }
    else if (delta_transformation == "log10") {
      gamma_common_se = log10(exp(1)) * (1 / coef(x)) * sqrt(diag(x$vcov))
      z_value = log10(coef(x)) / gamma_common_se
      gamma_ci_lower = 10 ** (log10(coef(x)) - stats::qnorm(1 - (alpha / 2)) * gamma_common_se)
      gamma_ci_upper = 10 ** (log10(coef(x)) + stats::qnorm(1 - (alpha / 2)) * gamma_common_se)
    }
    gamma_common_ci = matrix(
      data = c(gamma_ci_lower, ci_delta_upper),
      ncol = 2,
      byrow = FALSE
    )
    p_value =  (1 - stats::pnorm(abs(z_value))) * 2
  }
  # Score based inference
  if (x$inference == "score") {
    # (Re)construct reference trajectory.
    ref_fun = ref_fun_constructor(
      x$vertical_model$time_points,
      x$vertical_model$ctrl_estimates,
      x$interpolation
    )
    # Compute confidence interval
    gamma_common_ci = matrix(
      data = score_conf_int_common(
        time_points = x$vertical_model$time_points,
        ctrl_estimates = x$vertical_model$ctrl_estimates,
        exp_estimates = x$vertical_model$exp_estimates,
        ref_fun = ref_fun,
        interpolation = x$interpolation,
        vcov = x$vertical_model$vcov,
        gamma_est = x$coefficients,
        type = x$type,
        j = x$select_coef,
        weights = x$weights,
        alpha = alpha
      ),
      ncol = 2
    )
    # Compute score test statistic and p-value
    temp = score_test_common(
      time_points = x$vertical_model$time_points,
      ctrl_estimates = x$vertical_model$ctrl_estimates,
      exp_estimates = x$vertical_model$exp_estimates,
      ref_fun = ref_fun,
      interpolation = x$interpolation,
      vcov = x$vertical_model$vcov,
      gamma_0 = 1,
      type = x$type,
      j = x$select_coef,
      weights = x$weights)
    z_value = temp[1]
    p_value = temp[2]
    # SE not yet implemented
    gamma_common_se = NA
  }


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
    gamma_common_se = gamma_common_se,
    z_value = z_value,
    gamma_common_ci = gamma_common_ci,
    p_value = p_value,
    vcov_bootstrap = vcov_bootstrap,
    se_bootstrap = se_bootstrap,
    ci_bootstrap = ci_bootstrap,
    p_bootstrap = p_bootstrap,
    alpha = alpha
  )
}

new_summary.TCT_common = function(
    x,
    gamma_common_se,
    z_value,
    gamma_common_ci,
    p_value,
    vcov_bootstrap,
    se_bootstrap,
    ci_bootstrap,
    p_bootstrap,
    alpha
) {
  structure(append(
    x,
    list(
      gamma_common_se = gamma_common_se,
      z_value = z_value,
      gamma_common_ci = gamma_common_ci,
      p_value = p_value,
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
      `Std. Error (delta)` = x$gamma_common_se,
      `z-value (delta)` = x$z_value,
      `p-value (delta)` = x$p_value,
      `CI (delta)` = paste0("(",
                            format(x$gamma_common_ci[1], digits = 5),
                            ", ",
                            format(x$gamma_common_ci[2], digits = 5),
                            ")"),
      check.names = FALSE
    )
  }
  else {
    coefficients_df = data.frame(
      Value = coef(x),
      `Std. Error (delta)` = x$gamma_common_se,
      `z-value (delta)` = x$z_value,
      `p-value (delta)` = x$p_value,
      `CI (delta)` = paste0("(",
                            format(x$gamma_common_ci[1], digits = 5),
                            ", ",
                            format(x$gamma_common_ci[2], digits = 5),
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





