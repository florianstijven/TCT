#' Transform vertical Treatment Effect Estimates to the Time Scale
#'
#' The [TCT_meta()] function transforms so-called vertical parameter estimates
#' to parameter estimates on the time scale. These treatment effect estimates on
#' the time scale are so-called acceleration factors. This is an application of
#' the Time-Component Test (TCT) methodology.
#'
#' @param time_points Ordered vector that contains the times corresponding to
#'   the estimated means in the `ctrl_estimates` vector. The first element
#'   should be zero which generally corresponds to the time of randomization.
#' @param ctrl_estimates Estimated mean outcome in the control group at the
#'   fixed occasions corresponding to the times in `time_points`.
#' @param exp_estimates Estimated mean outcomes in the experimental group at
#'   fixed occasions corresponding to the times in `time_points[-1]`. Note that
#'   first element in `exp_estimates` should correspond to some time after
#'   randomization.
#' @param vcov The variance-covariance matrix for the means. In order to map to
#'   the correct estimates, this matrix should be the variance-covariance matrix
#'   of `c(ctrl_means, exp_means)`. If an element of `c(ctrl_means, exp_means)`
#'   is known (e.g., mean at baseline is zero when using change from baseline as
#'   outcome), then the corresponding row and column in `vcov` should be set to
#'   zero.
#' @param inference Which approach is used for estimation and inference? Should
#'   be `"wald"`, `"score"`, or `"least-squares"`. The `"wald"`-approach is explained and
#'   implemented in [DeltaMethod()], the `"score"`-approach is explained and
#'   documented in [score_test()] and [score_test_common()].
#' @param interpolation Which interpolation method to use?
#'  * `"linear"`: linear interpolation.
#'  * `"spline"`: natural cubic spline interpolation. This interpolation method has been most
#'   thoroughly tested is most stable.
#'  * `"monoH.FC`: monotone Hermite spline according to the method of Fritsch
#'   and Carlson.
#' @param B Number of parametric bootstrap replications. If `B = 0`, no
#'   bootstrap is performed (default).
#' @param constraints Use the constrained generalized least squares estimator
#'   for the vertical treatment effects.
#' @return S3 object of class `"TCT_meta"`
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
#' \eqn{f_0(t; \boldsymbol{\alpha})} and \eqn{f_1(t; \boldsymbol{\beta})}.
#'
#' The treatment effects on the time scale are acceleration factors, analogous
#' to accelerated failure time models. These are defined as follows at
#' \eqn{t_j}, \deqn{f_1(t; \boldsymbol{\beta}) = f_0(\gamma_j \cdot t;
#' \boldsymbol{\alpha})}
#' where \eqn{\gamma_j} is the so-called acceleration factor at \eqn{t_j}, i.e.,
#' treatment causes an acceleration of \eqn{\gamma_j}. For example, if
#' \eqn{\gamma_j = 0.5}, patients in the active treatment group progress half as
#' slow as patients in the control group.
#'
#' # Estimation and Inference
#'
#' Following options for estimation and inference are available:
#' * Delta method. More information in [DeltaMethod()].
#' * Parametric bootstrap. More information in
#' [pm_bootstrap_vertical_to_horizontal()].
#' * Score test: More information in [score_test()] and [score_test_common()].
#'
#' Note that the estimators in the above three methods are equivalent. The
#' estimates for the acceleration factor at \eqn{t_j} will thus be identical.
#' The difference lies in the procedures to computing standard errors, p-values,
#' and confidence intervals.
#'
#' Inference based on a direct application of the delta method can be unstable
#' and should be avoided. The parametric bootstrap and score test generally
#' yield very similar inferences. The main advanatage of the score test is that
#' it reduces to "classical" hypothesis tests for \eqn{H_0: \alpha_j = \beta_j}.
#' The p-value of the score test will thus equal that of classical hypothesis
#' tests.
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
#' TCT_Fit = TCT_meta(
#'   time_points = 0:4,
#'   ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#'   exp_estimates = coef(mmrm_fit)[5:8],
#'   vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#'   interpolation = "spline",
#'   B = 1e3
#' )
#' # The summary() generic can be used to obtain the most useful quantities from
#' # the Meta-TCT.
#' summary(TCT_Fit)
TCT_meta = function(time_points,
               ctrl_estimates,
               exp_estimates,
               vcov,
               interpolation = "spline",
               inference = "wald",
               B = 0,
               constraints = FALSE) {
  # Apply constraint GLS estimator if asked.
  if (constraints) {
    # Known mean parameters has not been implemented for the constrained GLS
    # estimator. An error should be raised when the user tries to use the
    # constrained estimator with known means.
    if (any(rowSums(vcov()) == 0) || any(colSums(vcov()) == 0)) {
      stop("`vcov` cannot contain zero variances when `constraints = TRUE`. This feature has not been implemented yet.")
    }
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

  return(
    new_TCT_meta(
      coefficients = coefficients,
      vcov = vcov_matrix,
      inference = inference,
      bootstrap_estimates = bootstrap_estimates,
      interpolation = interpolation,
      type = "time-based treatment effects",
      vertical_model = list(
        time_points = time_points,
        ctrl_estimates = ctrl_estimates,
        exp_estimates = exp_estimates,
        vcov = vcov
      )
    )
  )
}

new_TCT_meta = function(coefficients,
                   vcov,
                   inference,
                   bootstrap_estimates,
                   interpolation,
                   type,
                   vertical_model) {
  # All options regarding inference are put into a single list.
  inference_options = list(
    interpolation = interpolation,
    type = type,
    inference = inference
  )

  structure(list(coefficients = coefficients,
                 vcov = vcov,
                 inference_options = inference_options,
                 bootstrap_estimates = bootstrap_estimates,
                 vertical_model = vertical_model),
            class = "TCT_meta")
}


#' Print Meta-Time Component Test object
#'
#' @param x Object returned by [TCT_meta()]
#' @inheritDotParams base::print
#'
#' @export
#' @return NULL
#' @importFrom stats coef
print.TCT_meta = function(x, ...) {
  cat(
    paste0(
      "Meta-Time Component Test:",
      "\n\n"
    )
  )
  cat("Estimated Acceleration Factors: \n")
  print(coef(x))
  cat("\n Interpolation Method: ")
  cat(x$inference_options$interpolation)
}


new_summary_TCT_meta = function(
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
  structure(unclass(append(
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
  )),
  class = "summary_TCT_meta")
}

#' Summarize fitted Meta-Time Component Test model
#'
#' @param object Object returned by [TCT_meta()].
#' @param alpha Two-sided confidence level for confidence intervals.
#' @param delta_transformation Transformation when applying the delta-method to
#'   obtain confidence intervals.
#' @inheritDotParams base::summary
#'
#' @return S3 object of class `"summary_TCT_meta"`
#' @export
#' @inherit TCT_meta examples
#' @inheritParams score_conf_int
#' @importFrom stats coef
summary.TCT_meta = function(object,
                            alpha = 0.05,
                            delta_transformation = "identity",
                            bounds = c(-5, 5),
                            ...) {
  # Extract information from the TCT_meta object that is used further on.
  inference = object$inference_options$inference
  ctrl_estimates = object$vertical_model$ctrl_estimates
  exp_estimates = object$vertical_model$exp_estimates
  time_points = object$vertical_model$time_points
  interpolation = object$inference_options$interpolation
  vcov_vertical = object$vertical_model$vcov

  # Wald-based inference if the inference option is equal to "wald".
  if (inference == "wald") {
    if (delta_transformation == "identity") {
      se_delta = sqrt(diag(object$vcov))
      z_values = (1 - coef(object)) / se_delta
      ci_delta_lower = coef(object) - stats::qnorm(1 - (alpha / 2)) * se_delta
      ci_delta_upper = coef(object) + stats::qnorm(1 - (alpha / 2)) * se_delta
    }
    else if (delta_transformation == "log") {
      se_delta = (1 / coef(object)) * sqrt(diag(object$vcov))
      z_values = (log(coef(object))) / se_delta
      ci_delta_lower = exp(log(coef(object)) - stats::qnorm(1 - (alpha / 2)) * se_delta)
      ci_delta_upper = exp(log(coef(object)) + stats::qnorm(1 - (alpha / 2)) * se_delta)
    }
    else if (delta_transformation == "log10") {
      se_delta = log10(exp(1)) * (1 / coef(object)) * sqrt(diag(object$vcov))
      z_values = (log10(coef(object))) / se_delta
      ci_delta_lower = 10**(log10(coef(object)) - stats::qnorm(1 - (alpha / 2)) * se_delta)
      ci_delta_upper = 10**(log10(coef(object)) + stats::qnorm(1 - (alpha / 2)) * se_delta)
    }
    ci_matrix = matrix(
      data = c(ci_delta_lower, ci_delta_upper),
      ncol = 2,
      byrow = FALSE
    )
  }
  # Score-based inference if the inference option is equal to "score".
  if (inference == "score") {
    se_delta = sqrt(diag(object$vcov))
    # (Re)construct reference trajectory.
    ref_fun = ref_fun_constructor(
      time_points,
      ctrl_estimates,
      interpolation
    )
    # Compute confidence intervals for the measurement occasion-specific
    # acceleration factors based on the score test.
    ci_list = lapply(
      X = 1:length(exp_estimates),
      FUN = function(j) {
        score_conf_int(
          time_points = time_points,
          ctrl_estimates = ctrl_estimates,
          exp_estimates = exp_estimates,
          ref_fun = ref_fun,
          interpolation = interpolation,
          vcov = vcov_vertical,
          j = j,
          alpha = alpha,
          bounds = bounds
        )
      }
    )
    # Convert list of confidence intervals to a matrix. The first column
    # contains the lower limits, the second column contains the upper limits.
    ci_matrix = matrix(unlist(ci_list), ncol = 2, byrow = TRUE)
    # Compute the z-statistics for the measurement occasion-specific
    # acceleration factors based on the score test.
    z_values = vapply(
      X = 1:length(exp_estimates),
      FUN = function(j) {
        score_test(
          time_points = time_points,
          ctrl_estimates = ctrl_estimates,
          exp_estimates = exp_estimates,
          ref_fun = ref_fun,
          interpolation = interpolation,
          vcov = vcov_vertical,
          j = j,
          gamma_0 = 1
        )[1]
      },
      FUN.VALUE = 1.1
    )
  }

  lht_delta = linearHypothesis.default(
    model = object,
    vcov. = object$vcov,
    coef. = coef(object),
    rhs = rep(1, length(coef(object))),
    hypothesis.matrix = diag(1, nrow = length(coef(object)), ncol = length(coef(object)))
  )
  p_values =  (1 - stats::pnorm(abs(z_values))) * 2

  # inference based on parametric bootstrap
  if (!(is.null(object$bootstrap_estimates))) {
    vcov_bootstrap = stats::var(object$bootstrap_estimates, na.rm = TRUE)
    se_bootstrap = sqrt(diag(vcov_bootstrap))
    ci_bootstrap = t(apply(
      X = object$bootstrap_estimates,
      MARGIN = 2,
      FUN = stats::quantile,
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE
    ))
    p_bootstrap = apply(
      X = object$bootstrap_estimates,
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
    new_summary_TCT_meta(
      x = object,
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

#' Print summary_TCT_meta object
#'
#' @param x Object returned by [summary.TCT_meta()].
#' @inheritDotParams base::print
#'
#' @export
#' @return NULL
#' @importFrom stats coef
print.summary_TCT_meta = function(x, ...) {
  cat(
    paste0(
      "Meta-Time Component Test: ",
      "\n\n"
    )
  )
  cat("Coefficients: \n")
  if (is.null(x$ci_bootstrap)) {
    coefficients_df = data.frame(
      "Estimate" = coef(x),
      "Std. Error" = x$se_delta,
      "z value" = x$z_values,
      "p value" = x$p_values,
      "CI" = paste0(
        "(",
        format(x$ci_matrix[, 1], digits = 5),
        ", ",
        format(x$ci_matrix[, 2], digits = 5),
        ")"
      ),
      check.names = FALSE
    )
  }
  else {
    coefficients_df = data.frame(
      Value = coef(x),
      `Std. Error` = x$se_delta,
      `z value` = x$z_values,
      `p value` = x$p_values,
      `CI` = paste0(
        "(",
        format(x$ci_matrix[, 1], digits = 5),
        ", ",
        format(x$ci_matrix[, 2], digits = 5),
        ")"
      ),
      `CI (bootstrap)` = paste0(
        "(",
        format(x$ci_bootstrap[, 1], digits = 5),
        ", ",
        format(x$ci_bootstrap[, 2], digits = 5),
        ")"
      ),
      check.names = FALSE
    )

  }

  print(coefficients_df, digits = 5)

  cat(paste0("alpha = ", x$alpha))
  cat("\n")
  cat("\nInterpolation Method: ")
  cat(x$inference_options$interpolation)
  cat("\nEstimation and Inference: ")
  cat(x$inference_options$inference)
  cat("\n")


}

#' Estimate Common Acceleration Factor
#'
#' The [TCT_meta_common()] function estimates the common acceleration factor.
#'
#' @param TCT_Fit Object returned by [TCT_meta()]
#' @param select_coef Estimates from the `exp_estimates` in [TCT_meta()] to use
#'  in estimating the common acceleration factor. If there is reason to believe
#'  that the proportional slowing assumption does not hold, e.g., for the first
#'  measurement after randomization, then the corresponding estimate should not
#'  be used in estimation the common acceleration factor.
#' @inheritParams TCT_meta
#' @inheritParams score_test_common
#' @inheritParams pm_bootstrap_vertical_to_common
#' @inheritParams nonlinear_gls_estimator
#'
#' @return S3 object of class `"TCT_meta_common"`
#' @export
#' @importFrom stats coef
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
#' TCT_Fit = TCT_meta(
#'   time_points = 0:4,
#'   ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
#'   exp_estimates = coef(mmrm_fit)[5:8],
#'   vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
#'   interpolation = "spline",
#'   B = 1e3
#' )
#' TCT_fit_common = TCT_meta_common(
#'   TCT_Fit = TCT_Fit,
#'   inference = "wald"
#' )
TCT_meta_common = function(TCT_Fit,
                           inference = "wald",
                           B = 0,
                           bs_fix_vcov = FALSE,
                           select_coef = 1:length(coef(TCT_Fit)),
                           constraints = FALSE,
                           type = NULL,
                           weights = NULL,
                           start_gamma = 0.75)
{
  # Extract information from the TCT_meta object that is used further on.
  ctrl_estimates = TCT_Fit$vertical_model$ctrl_estimates
  exp_estimates = TCT_Fit$vertical_model$exp_estimates
  time_points = TCT_Fit$vertical_model$time_points
  interpolation = TCT_Fit$inference_options$interpolation
  vcov_vertical = TCT_Fit$vertical_model$vcov
  # (Re)construct reference trajectory.
  ref_fun = ref_fun_constructor(time_points,
                                ctrl_estimates,
                                interpolation)

  # Use wald-based inference if this is asked by the user.
  n_points = length(TCT_Fit$vertical_model$time_points)
  if (inference == "wald") {
    estimates = coef(TCT_Fit)[select_coef]
    vcov = TCT_Fit$vcov[select_coef, select_coef]
    # delta method
    p = length(estimates)
    vec_1 = matrix(1, nrow = p, ncol = 1)
    gamma_common_estimate = (t(vec_1) %*% solve(vcov) %*% matrix(estimates, ncol = 1)) /
      (t(vec_1) %*% solve(vcov) %*% vec_1)
    gamma_common_vcov = (t(vec_1) %*% solve(vcov) %*% vec_1) ** (-1)
  }
  else if (inference == "score") {
    # Compute optimal weights if the user did not provide weights themselves.
    if ((type == "custom") & (is.null(weights))) {
      weights = optimize_weights(
        time_points = time_points,
        ctrl_estimates = ctrl_estimates,
        exp_estimates = exp_estimates,
        ref_fun = ref_fun,
        interpolation = interpolation,
        vcov = vcov_vertical,
        j = select_coef
      )
    }
    gamma_common_estimate = score_estimate_common(
      time_points = time_points,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      ref_fun = ref_fun,
      interpolation = interpolation,
      vcov = vcov_vertical,
      type = type,
      j = select_coef,
      weights = weights
    )

    gamma_common_vcov = score_estimate_common_se(
      gamma_est = gamma_common_estimate,
      time_points = time_points,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      interpolation = interpolation,
      vcov = vcov_vertical,
      type = type,
      j = select_coef,
      weights = weights
    )
  }
  else if (inference == "least-squares") {
    # Estimate parameter vector. The first length(ctrl_estimates) elements are
    # estimated for the mean parameters in the control group. The last element
    # is the estimate for the common acceleration factor.
    estimates_vec = nonlinear_gls_estimator(
      time_points = time_points,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      interpolation = interpolation,
      vcov = vcov_vertical,
      j = select_coef,
      start_gamma = start_gamma
    )$estimate
    # Extract estimate for the common acceleration factor.
    gamma_common_estimate = estimates_vec[length(ctrl_estimates) + 1]
    # Extract estimates for the mean parameters in the control group.
    alpha_est = estimates_vec[1:length(ctrl_estimates)]
    # Compute SE for the common acceleration factor.
    gamma_common_vcov = nonlinear_gls_estimator_se(
      time_points = time_points,
      interpolation = interpolation,
      vcov = vcov_vertical,
      j = select_coef,
      gamma_est = gamma_common_estimate,
      alpha_est = alpha_est
    )
  }

  bs_estimates = pm_bootstrap_vertical_to_common(
    TCT_Fit = TCT_Fit,
    inference = inference,
    B = B,
    bs_fix_vcov = bs_fix_vcov,
    return_se = TRUE,
    select_coef = select_coef,
    constraints = constraints,
    type = type,
    weights = weights,
    start_gamma = start_gamma
  )

  # Test for common slowing parameter (OLD)
  # lht_matrix = matrix(0,
  #                     nrow = length(estimates) - 1,
  #                     ncol = length(estimates) - 1)
  # diag(lht_matrix) = -1
  # lht_matrix = cbind(1, lht_matrix)

  # lht_common = linearHypothesis.default(
  #   model = TCT_Fit,
  #   vcov. = vcov,
  #   coef. = estimates,
  #   rhs = rep(0, length(estimates) - 1),
  #   hypothesis.matrix = lht_matrix
  # )
  # Test for common slowing factor based on the general least squares criterion.
  proportional_slowing_test = proportional_slowing_gls_test(
    time_points = time_points,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    interpolation = interpolation,
    vcov = vcov_vertical,
    j = select_coef,
    start_gamma = start_gamma
  )


  lht_common = NULL


  new_TCT_meta_common(
    coefficients = gamma_common_estimate,
    inference = inference,
    vcov = gamma_common_vcov,
    bootstrap_estimates = bs_estimates,
    interpolation = interpolation,
    lht_common = lht_common,
    type = type,
    weights = weights,
    vertical_model = TCT_Fit$vertical_model,
    select_coef = select_coef,
    proportional_slowing_test = proportional_slowing_test
  )
}

new_TCT_meta_common = function(coefficients,
                          inference,
                          vcov,
                          bootstrap_estimates,
                          interpolation,
                          lht_common,
                          type,
                          weights,
                          vertical_model,
                          select_coef,
                          proportional_slowing_test
                          ) {
  # All options regarding inference are put into a single list.
  inference_options = list(
    interpolation = interpolation,
    type = type,
    inference = inference,
    weights = weights,
    select_coef = select_coef
  )

  structure(
    list(
      coefficients = coefficients,
      vcov = vcov,
      inference_options = inference_options,
      bootstrap_estimates = bootstrap_estimates,
      lht_common = lht_common,
      vertical_model = vertical_model,
      proportional_slowing_test = proportional_slowing_test
    ),
    class = "TCT_meta_common"
  )
}


#' @export
vcov.TCT_meta_common = function(object, ...) {
  object$vcov
}

#' Print TCT_meta_common object
#'
#' @param x Object returned by [TCT_meta_common()].
#' @inheritDotParams base::print
#'
#' @export
#' @return NULL
#' @importFrom stats coef
print.TCT_meta_common = function(x, ...) {
  cat(
    paste0(
      "Meta-Time Component Test - ",
      "Common Acceleration Factor:",
      "\n\n"
    )
  )
  cat("Estimated Common Acceleration Factor: \n")
  print(coef(x))
  cat("\nInterpolation Method: ")
  cat(x$inference_options$interpolation)
  cat("\n")

}

#' Summarize fitted Meta-Time Component Test model with Proportional Slowing
#'
#' @param object Object returned by [TCT_meta_common()].
#' @inheritParams summary.TCT_meta
#' @inheritDotParams base::summary
#'
#' @return S3 object of class `"summary_TCT_meta_common"`
#' @export
#' @inherit TCT_meta_common examples
#' @importFrom stats coef
summary.TCT_meta_common = function(object,
                                   alpha = 0.05,
                                   delta_transformation = "identity",
                                   ...) {
  # Extract information from the TCT_meta object that is used further on.
  inference = object$inference_options$inference
  ctrl_estimates = object$vertical_model$ctrl_estimates
  exp_estimates = object$vertical_model$exp_estimates
  time_points = object$vertical_model$time_points
  interpolation = object$inference_options$interpolation
  vcov_vertical = object$vertical_model$vcov
  type = object$inference_options$type
  select_coef = object$inference_options$select_coef
  weights = object$inference_options$weights
  proportional_slowing_test = object$proportional_slowing_test
  # Wald-based inference
  if (inference == "wald") {
    if (delta_transformation == "identity") {
      gamma_common_se = sqrt(diag(object$vcov))
      z_value = (1 - coef(object)) / gamma_common_se
      gamma_ci_lower = coef(object) - stats::qnorm(1 - (alpha / 2)) * gamma_common_se
      gamma_ci_upper = coef(object) + stats::qnorm(1 - (alpha / 2)) * gamma_common_se
    }
    else if (delta_transformation == "log") {
      gamma_common_se = (1 / coef(object)) * sqrt(diag(object$vcov))
      z_value = (log(coef(object))) / gamma_common_se
      gamma_ci_lower = exp(log(coef(object)) - stats::qnorm(1 - (alpha / 2)) * gamma_common_se)
      gamma_ci_upper = exp(log(coef(object)) + stats::qnorm(1 - (alpha / 2)) * gamma_common_se)
    }
    else if (delta_transformation == "log10") {
      gamma_common_se = log10(exp(1)) * (1 / coef(object)) * sqrt(diag(object$vcov))
      z_value = log10(coef(object)) / gamma_common_se
      gamma_ci_lower = 10 ** (log10(coef(object)) - stats::qnorm(1 - (alpha / 2)) * gamma_common_se)
      gamma_ci_upper = 10 ** (log10(coef(object)) + stats::qnorm(1 - (alpha / 2)) * gamma_common_se)
    }
    gamma_common_ci = matrix(
      data = c(gamma_ci_lower, gamma_ci_upper),
      ncol = 2,
      byrow = FALSE
    )
    p_value =  (1 - stats::pnorm(abs(z_value))) * 2
  }
  # Score based inference
  if (inference == "score") {
    # (Re)construct reference trajectory.
    ref_fun = ref_fun_constructor(
      time_points,
      ctrl_estimates,
      interpolation
    )
    # Compute confidence interval
    gamma_common_ci = matrix(
      data = score_conf_int_common(
        time_points = time_points,
        ctrl_estimates = ctrl_estimates,
        exp_estimates = exp_estimates,
        ref_fun = ref_fun,
        interpolation = interpolation,
        vcov = vcov_vertical,
        gamma_est = coef(object),
        type = type,
        j = select_coef,
        weights = weights,
        alpha = alpha
      ),
      ncol = 2
    )
    # Compute score test statistic and p-value
    temp = score_test_common(
      time_points = time_points,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      ref_fun = ref_fun,
      interpolation = interpolation,
      vcov = vcov_vertical,
      gamma_0 = 1,
      type = type,
      j = select_coef,
      weights = weights
    )
    z_value = temp[1]
    p_value = temp[2]
    # Estimated standard error.
    gamma_common_se = object$vcov
  }
  # Nonlinear GLS based inference
  if (inference == "least-squares") {
    # (Re)construct reference trajectory.
    ref_fun = ref_fun_constructor(
      time_points,
      ctrl_estimates,
      interpolation
    )
    # Compute confidence interval
    gamma_common_ci = matrix(
      data = nonlinear_gls_conf_int_common(
        time_points = time_points,
        ctrl_estimates = ctrl_estimates,
        exp_estimates = exp_estimates,
        interpolation = interpolation,
        vcov = vcov_vertical,
        j = select_coef,
        alpha = alpha,
        start_gamma = coef(object)
      ),
      ncol = 2
    )
    # Compute test statistic and p-value
    temp = nonlinear_gls_test(
      time_points = time_points,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      interpolation = interpolation,
      vcov = vcov_vertical,
      gamma_0 = 1,
      j = select_coef
    )
    z_value = temp[1]
    p_value = temp[2]
    # Estimated standard error.
    gamma_common_se = object$vcov
  }


  # inference based on parametric bootstrap
  if (!(is.null(object$bootstrap_estimates))) {
    vcov_bootstrap = stats::var(object$bootstrap_estimates[[1]], na.rm = TRUE)
    se_bootstrap = sqrt(vcov_bootstrap)
    ci_bootstrap = stats::quantile(
      x = object$bootstrap_estimates[[1]],
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE
    )
    p_bootstrap = min(mean(object$bootstrap_estimates[[1]] > 1, na.rm = TRUE) * 2,
                      (1 - mean(object$bootstrap_estimates[[1]] > 1, na.rm = TRUE)) * 2)
  }
  else {
    vcov_bootstrap = NULL
    se_bootstrap = NULL
    ci_bootstrap = NULL
    p_bootstrap = NULL
  }


  new_summary_TCT_meta_common(
    x = object,
    gamma_common_se = gamma_common_se,
    z_value = z_value,
    gamma_common_ci = gamma_common_ci,
    p_value = p_value,
    vcov_bootstrap = vcov_bootstrap,
    se_bootstrap = se_bootstrap,
    ci_bootstrap = ci_bootstrap,
    p_bootstrap = p_bootstrap,
    alpha = alpha,
    proportional_slowing_test = proportional_slowing_test
  )
}

new_summary_TCT_meta_common = function(
    x,
    gamma_common_se,
    z_value,
    gamma_common_ci,
    p_value,
    vcov_bootstrap,
    se_bootstrap,
    ci_bootstrap,
    p_bootstrap,
    alpha,
    proportional_slowing_test
) {
  structure(unclass(append(
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
      alpha = alpha,
      proportional_slowing_test = proportional_slowing_test
    )
  )),
  class = "summary_TCT_meta_common")
}

#' Print summary_TCT_meta_common object
#'
#' @param x Object returned by [summary.TCT_meta_common()].
#' @inheritDotParams base::print
#'
#' @export
#' @return NULL
#' @importFrom stats coef
print.summary_TCT_meta_common = function(x, ...) {
  cat(
    paste0(
      "Meta-Time Component Test - ",
      "Common Acceleration Factor:",
      "\n\n"
    )
  )
  cat("Estimated Common Acceleration Factor: \n")
  if (is.null(x$vcov_bootstrap)) {
    coefficients_df = data.frame(
      Estimate = coef(x),
      `Std. Error` = x$gamma_common_se,
      `z value` = x$z_value,
      `p value` = x$p_value,
      `CI` = paste0("(",
                            format(x$gamma_common_ci[1], digits = 5),
                            ", ",
                            format(x$gamma_common_ci[2], digits = 5),
                            ")"),
      check.names = FALSE
    )
  }
  else {
    coefficients_df = data.frame(
      Estimate = coef(x),
      `Std. Error` = x$gamma_common_se,
      `z value` = x$z_value,
      `p value` = x$p_value,
      `CI` = paste0("(",
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
      `Std. Error (bootstrap)` = x$se_bootstrap,
      check.names = FALSE
    )
  }
  # If the omnibus score-test is used, the test statistic is not a z-value, but
  # rather chi-squared statistic.
  if (x$inference_options$inference == "score") {
    if (x$inference_options$type == "omnibus") {
      colnames(coefficients_df)[3] = "chi-squared"
    }
  }
  if (x$inference_options$inference == "least-squares") {
    colnames(coefficients_df)[3] = "chi-squared"
  }

  print(coefficients_df, digits = 5, row.names = FALSE)
  cat(paste0("alpha = ", x$alpha))
  cat("\n")
  cat("\nInterpolation Method: ")
  cat(x$inference_options$interpolation)
  cat("\nTime Points Used in Estimator: ")
  cat(x$vertical_model$time_points[x$inference_options$select_coef + 1])
  cat("\nEstimation and Inference: ")
  cat(x$inference_options$inference)
  if (x$inference_options$inference == "score") {
    cat("\n")
    cat("\t")
    cat(paste0("Type of Score Test: ", x$inference_options$type))
    if (x$inference_options$type == "custom") {
      cat("\n")
      cat("\t")
      cat("Weights: ")
      cat(round(x$inference_options$weights, 3))
    }
  }
  cat("\n\n")




  # Print test for proportional slowing (OLD)
  # cat("Test for proportional slowing factor:\n")
  # print(
  #   data.frame("Df" = x$lht_common$Df[2],
  #              "Chisq" = x$lht_common$Chisq[2],
  #              "p-value" = x$lht_common$`Pr(>Chisq)`[2]),
  #   digits = 5
  # )

  # Print test for proportional slowing based on minimized GLS criterion.
  cat("Test for proportional slowing factor (at the selected time points):\n")
  print(
    data.frame(
      "Df" = x$proportional_slowing_test["df"],
      "Chisq" = x$proportional_slowing_test["chi_squared"],
      "p-value" = x$proportional_slowing_test["p_value"]
    ),
    digits = 5,
    row.names = FALSE
  )
}





