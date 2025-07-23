#' Compute TCT contrast test z-value
#'
#' The [contrast_test()] function computes the z-value for the contrast test at
#' measurement `j`. This corresponds to `exp_estimates[j]`. For technical
#' details on the contrast test, see Details.
#'
#' @param j Measurement occasion to test acceleration factor for. This
#'   corresponds to `exp_estimates[j]`.
#' @param gamma_0 Value for the acceleration factor under the null hypothesis.
#' @inheritParams DeltaMethod
#'
#' @details
#'
#' In the following sections, we use notation introduced in the documentation of
#' [TCT_meta()].
#'
#' # Contrast Test
#'
#' For constructing a contrast test, we start from the null hypothesis,
#' \deqn{H_{\gamma_0}: \gamma_j = \gamma_0,} where \eqn{\gamma_j} is the
#' acceleration factor at time \eqn{t_j}. Under this null hypothesis, and taking
#' the control group as reference, we can
#' compute the experimental treatment group mean at \eqn{t_j} as \deqn{\beta_{0,
#' j} = f_0( \gamma_0 \cdot t_{j}; \boldsymbol{\alpha})} given the true control
#' group mean vector, \eqn{\boldsymbol{\alpha}}. Consequently, a test for
#' \eqn{H_0: \beta_j = \beta_{0, j}} is also a valid test for
#' \eqn{H_{\gamma_0}}. Well-established tests are available for the former
#' hypothesis. However, these cannot be applied directly since
#' \eqn{\boldsymbol{\alpha}} is not known, but only estimated. Indeed, \eqn{\beta_{0,
#' j}} is itself estimated by \eqn{\hat{\beta}_{0, j} =  f_0( \gamma_0 \cdot t_{j};
#' \hat{\boldsymbol{\alpha}})}. Still, this provides the starting point for a
#' contrast test that draws upon well-established statistical tests.
#'
#' We first define the contrast vector as follows, \deqn{\boldsymbol{\Phi}(\boldsymbol{\gamma}; \hat{\boldsymbol{\alpha}}_n, \hat{\boldsymbol{\beta}}_n) =
#' \hat{\boldsymbol{\beta}}_n - \left(f_0(\gamma_1 \cdot t_1;
#' \hat{\boldsymbol{\alpha}}_n), \dots, f_0(\gamma_K \cdot t_K;
#' \hat{\boldsymbol{\alpha}}_n)\right)^\top.}
#' Under \eqn{H_{\gamma_0}} for all time points and using the delta
#' method, we have that
#' \deqn{n^{1/2}\boldsymbol{\Phi}(\boldsymbol{\gamma_0}; \hat{\boldsymbol{\alpha}}_n,
#' \hat{\boldsymbol{\beta}}_n) \overset{d}{\to} \mathcal{N}(\boldsymbol{0}, \Omega)} where
#' \deqn{\Omega =
#' \dot{\boldsymbol{\Phi}}_{\boldsymbol{\alpha}, \boldsymbol{\beta}}(\boldsymbol{\gamma_0}; \boldsymbol{\alpha_0}, \boldsymbol{\beta_0}) \cdot \Sigma \cdot
#' \dot{\boldsymbol{\Phi}}_{\boldsymbol{\alpha}, \boldsymbol{\beta}}(\boldsymbol{\gamma_0}; \boldsymbol{\alpha_0}, \boldsymbol{\beta_0})^\top} with
#' \eqn{\dot{\boldsymbol{\Phi}}_{\boldsymbol{\alpha}, \boldsymbol{\beta}}(\boldsymbol{\gamma_0}; \boldsymbol{\alpha_0}, \boldsymbol{\beta_0})} the
#' Jacobian of \eqn{(\boldsymbol{\alpha}, \boldsymbol{\beta}) \mapsto \boldsymbol{\Phi}(\boldsymbol{\gamma_0}; \boldsymbol{\alpha}, \boldsymbol{\beta})} evaluated in the true parameter vector,
#' \eqn{(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})'}.
#' In practice, \eqn{\Omega} is unknown because the true parameter vector is unknown.
#' We, therefore, replace the unkown parameters wit the corresponding estimates in
#' the above expression for \eqn{\Omega}; the corresponding estimated matrix is
#' denoted by \eqn{\hat{\Omega}_n}. As long as \eqn{\hat{\Omega}_n} is consistent
#' for \eqn{\Omega}, we have the following distributional result, which forms the
#' basis for hypothesis tests and confidence intervals:
#' \deqn{n^{1/2}\boldsymbol{\Phi}(\boldsymbol{\gamma_0}; \hat{\boldsymbol{\alpha}}_n,
#' \hat{\boldsymbol{\beta}}_n) \overset{d}{\to} \mathcal{N}(\boldsymbol{0}, \hat{\Omega}_n).}
#' Specifically, for tests for time-specific acceleration factors, we use the
#' following distributional result:
#' \deqn{n^{1/2}\Phi_j(\gamma_0; \hat{\boldsymbol{\alpha}}_n,
#' \hat{\boldsymbol{\beta}}_n) \overset{d}{\to} \mathcal{N}(\boldsymbol{0}, \hat{\Omega}_n)}
#' where
#'
#'
#' @return Named (numeric) vector with two elements:
#'
#' 1. `"z"`. This test-statistic follows a standard normal distribution under
#' the null hypothesis.
#' 2. `"p-value"`: two-sided p-value.
contrast_test = function(time_points,
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
  deriv_f0_alpha_t_j  = -1 * deriv_f0_alpha(
    t_m = gamma_0 * time_points[j + 1],
    x_ref = time_points,
    y_ref = ctrl_estimates,
    method = interpolation
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
  # Return z-value and corresponding two-sided p-value.
  return(c("z" = as.numeric(z),
           "p-value" = 2 * (1 - stats::pnorm(abs(as.numeric(z))))
           )
         )
}

#' Compute Confidence Interval Based on Contrast Test
#'
#' The [contrast_conf_int()] function computes the confidence interval for the
#' contrast test at measurement `j`. This confidence interval is based on the
#' contrast test implemented in [contrast_test()].
#'
#' @param alpha `1 - alpha` represents the two-sided confidence level. Defaults
#'   to `0.05`.
#' @param bounds (numeric) vector with two elements that defines the bounds of
#'   the search interval for the (contrast-based) confidence limits or estimate.
#'   Defaults to `c(-5, 5)`.
#' @inheritParams contrast_test
#'
#' @details
#'
#' # Contrast Test Confidence Intervals
#'
#' For the construction of a confidence interval based on the contrasts, we make
#' use of the relationship between confidence intervals and hypothesis tests. A
#' \eqn{1 − α} confidence interval can be defined as follows, \deqn{\left\{
#' \gamma : p(\gamma) > \alpha \right\}} where \eqn{p(\gamma)} is the p-value
#' for the null that the acceleration factor is equal to \eqn{\gamma}. This
#' p-value can be based on the time-specific contrast test implemented in
#' [contrast_test()], or can be based on the joint contrast test implemented in
#' [contrast_test_common()]. The former gives a confidence interval for the
#' time-specific acceleration factor while the latter gives a confidence
#' interval for the common acceleration factor. The latter interval should be
#' used with some care. Indeed, the confidence interval for a common
#' acceleration factor is only sensible if there is no evidence against the
#' assumptions of constant slowing of the disease progression.
#'
#' @return (numeric) vector with two elements. The first element is the lower
#'   confidence limit, the second element is the upper confidence limit.
contrast_conf_int = function(time_points,
                          ctrl_estimates,
                          exp_estimates,
                          ref_fun,
                          interpolation,
                          vcov,
                          j,
                          alpha = 0.05,
                          bounds = c(-5, 5)) {
  # Construct function of gamma that return the z-value.
  z_value = function(gamma) {
    return(contrast_test(time_points,
                      ctrl_estimates,
                      exp_estimates,
                      ref_fun,
                      interpolation,
                      vcov,
                      j,
                      gamma_0 = gamma)[1])
  }
  # Below, a confidence region is computed based on the equivalence between
  # hypothesis testing and confidence regions. We assume that the confidence
  # region is a confidence interval. The limits can then be obtained as the
  # parameter values for which the corresponding p-value for the null that the
  # true value is equal to the parameter value is equal to alpha, or
  # equivalently, the critical z-value (typically 1.96 and 1.96). Depending on
  # the context, the lower confidence limit is the parameter with a
  # corresponding z-value equal to -1.96 or 1.96. Below, we assume that this is
  # -1.96 for the lower limit and 1.96 for the upper limit. If this assumption
  # was wrong, we just switch the computed limits.

  # Find upper limit
  z_critical = stats::qnorm(p = 1 - alpha / 2)

  # Test whether the z-value crosses z_critical going from the lower to the
  # upper bound of the search interval.
  if (sign(z_value(bounds[1]) + z_critical) != sign(z_value(bounds[2]) + z_critical)) {
    # If the sign changes, then the z-value crosses the critical value.
    # Therefore, we can compute the upper confidence limit.
    upper_limit = stats::uniroot(
      f = function(gamma)
        z_value(gamma) + z_critical,
      interval = bounds,
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }
  else {
    # If the sign does not change, then the z-value does not cross the critical
    # value. Therefore, we set the upper limit to infinity and raise a warning.
    upper_limit = +Inf
    warning("The z-value does not cross the critical value in the search interval for computing the contrast-based confidence interval. \\
            The corresponding confidence limit is set to Inf. Consider increasing the search interval bounds using the `bounds` argument.")
  }

  # Test whether the z-value crosses z_critical going from the lower to the
  # upper bound of the search interval.
  if (sign(z_value(bounds[1]) - z_critical) != sign(z_value(bounds[2]) - z_critical)) {
    # If the sign changes, then the z-value crosses the critical value.
    # Therefore, we can compute the upper confidence limit.
    lower_limit = stats::uniroot(
      f = function(gamma)
        z_value(gamma) - z_critical,
      interval = bounds,
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }
  else {
    # If the sign does not change, then the z-value does not cross the critical
    # value. Therefore, we set the upper limit to infinity and raise a warning.
    lower_limit = -Inf
    warning("The z-value does not cross the critical value in the search interval for computing the contrast-based confidence interval. \\
            The corresponding confidence limit is set to Inf. Consider increasing the search interval bounds using the `bounds` argument.")
  }

  # Return estimated confidence interval. As mentioned in the comment above, we
  # switch the confidence limits if required.
  return(
    c(min(lower_limit, upper_limit),
      max(lower_limit, upper_limit))
  )
}

#' contrast test for common acceleration factor
#'
#' The [contrast_test_common()] function implements the contrast test under the
#' assumption that there exists a common acceleration factor. Multiple variants
#' to this test exist and are implemented.
#'
#' @param type Which type of test statistic should be used. See Test Statistic
#'   Variants. Should be one of
#'   1. `type = "omnibus"`
#'   2. `type = "directional"`
#'   3. `type = "inverse variance"`
#'   4. `type = "custom"`
#' @param j (Integer) vector that indicates which elements in `exp_estimates`
#'   should be used for the contrast test. Defaults to `1:length(exp_estimates)`,
#'   i.e., all elements are used.
#' @param weights If `type == "custom"`, the user should specify a weight
#'   vector for weighting estimates at different time points differently.
#' @inheritParams contrast_test
#' @inheritSection contrast_test Contrast Test
#' @details
#'
#' # Test Statistic Variants
#'
#' Two types of test statistics are implemented in the [contrast_test_common()]
#' function:
#'
#' 1. `type = "omnibus"`
#' 2. `type = "custom"`
#'
#' These are discussed in more detail next.
#'
#' ## Omnibus
#'
#' The omnibus contrast test is based on the "classic" chi-squared statistic that is
#' defined as follows,
#' \deqn{t^2 = \boldsymbol{s}(\gamma_0 \cdot \boldsymbol{t}; \hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})^t \cdot \Sigma_s^{-1} \cdot \boldsymbol{s}(\gamma_0 \cdot \boldsymbol{t}; \hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}}).}
#' which follows a chi-squared distribution with K degrees of freedom under
#' \eqn{H_{\gamma_0}} for all time points. An important feature of this test
#' statistic is that it reduces to the chi-squared statistic for \eqn{H_0 :
#' \alpha_j = \beta_j} for all \eqn{j} in `j` when `gamma_0 = 1`. However, this
#' also means that there is no gain in power.
#'
#' ## Custom
#'
#' The "custom" contrast test allows for user-specified weights, \eqn{w}. The test
#' statistic is defined as follows,
#' \deqn{z = \frac{v(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}}; \gamma_0)}{\sqrt{\boldsymbol{w}^t \Sigma_s \boldsymbol{w}}} \; \dot\sim \; N(0, 1).}
#' As indicated above, this test statistic follows a standard normal
#' distribution under the null hypothesis.
#'
#' Note that [optimize_weights()] allows one to find optimal weights in the
#' sense of minimizing the estimated standard error for the corresponding
#' estimator of the common acceleration factor.
#'
#'
#' @return Named (numeric) vector with two elements:
#'
#' 1. `"z"` or `"chi-squared"`. The test statistic.
#' 2. `"p-value"`: two-sided p-value.
#'
#' @export
contrast_test_common = function(time_points,
                             ctrl_estimates,
                             exp_estimates,
                             ref_fun,
                             interpolation,
                             vcov,
                             gamma_0 = 1,
                             type = "omnibus",
                             j = 1:length(exp_estimates),
                             weights = NULL){
  K = length(j)
  J = contrast_vector_jacobian(
    time_points,
    ctrl_estimates,
    interpolation,
    gamma_0,
    j
  )
  # Compute the variance of the "test statistic "contrast vector", g_gamma_0 under
  # the null. The inverse of this matrix is also computed.
  Sigma_g = J %*% vcov[c(1:length(time_points), length(time_points) + j),
                          c(1:length(time_points), length(time_points) + j)] %*% t(J)
  Sigma_g_inv = solve(Sigma_g)
  # Compute the contrast vector.
  g = (exp_estimates[j] - ref_fun(gamma_0 * time_points[j + 1]))
  if (type == "omnibus") {
    # Compute test-statistic
    t_sq =  g %*% Sigma_g_inv %*% g
    return(
      c(
        "chi-squared" = as.numeric(t_sq),
        "p-value" = 1 - stats::pchisq(as.numeric(t_sq), K)
      )
    )
  }
  else if (type == "directional") {
    # Colmun vector of ones.
    ones = matrix(1, nrow = K)
    t_sq = (1 / t(ones) %*% Sigma_g_inv %*% ones) * (t(ones) %*% Sigma_g_inv %*% g)**2
    return(
      c(
        "chi-squared" = as.numeric(t_sq),
        "p-value" = 1 - stats::pchisq(as.numeric(t_sq), 1)
      )
    )
  }
  else if (type == "inverse variance") {
    # Colmun vector of ones.
    ones = matrix(1, nrow = K)
    if (nrow(Sigma_g) == 1) D = 1 / as.numeric(Sigma_g)
    else D = solve(diag(diag(Sigma_g)))
    z_value = (1 / sqrt(t(ones) %*% D %*% Sigma_g %*% t(D) %*% ones)) * t(ones) %*% D %*% g
    return(
      c(
        "z" = as.numeric(z_value),
        "p-value" = 2 * (1 - stats::pnorm(abs(as.numeric(z_value))))
      )
    )
  }
  else if (type == "custom") {
    weights = matrix(weights, ncol = 1)
    z_value = ( 1 / sqrt(t(weights) %*% Sigma_g %*% weights)) * t(weights) %*% g
    return(
      c(
        "z" = as.numeric(z_value),
        "p-value" = 2 * (1 - stats::pnorm(abs(as.numeric(z_value))))
      )
    )
  }
}



#' Confidence interval based on contrast test
#'
#' The [contrast_conf_int_common()] function computes the confidence interval for
#' common acceleration factor. This confidence interval is based on the contrast
#' test implemented in [contrast_test_common()].
#'
#' @inheritParams contrast_estimate_common
#' @inheritParams contrast_conf_int
#' @param gamma_est Estimate for the common acceleration factor.
#' @inheritSection contrast_conf_int Contrast Test Confidence Intervals
#' @inherit contrast_conf_int return
contrast_conf_int_common = function(time_points,
                                 ctrl_estimates,
                                 exp_estimates,
                                 ref_fun,
                                 interpolation,
                                 vcov,
                                 gamma_est,
                                 type = "omnibus",
                                 j = 1:length(exp_estimates),
                                 weights = NULL,
                                 alpha = 0.05,
                                 bounds = c(-5, 5)){
  # Force argument values. This is required because we're using these arguments
  # in a function factory.
  force(time_points); force(ctrl_estimates); force(exp_estimates)
  force(ref_fun); force(interpolation); force(vcov); force(gamma_est)
  force(type); force(j); force(weights); force(alpha)
  if (type %in% c("omnibus", "directional")) {
    # Construct function of gamma that return the z-value.
    t_sq_value = function(gamma) {
      t_sq = contrast_test_common(
        time_points,
        ctrl_estimates,
        exp_estimates,
        ref_fun,
        interpolation,
        vcov,
        gamma,
        type,
        j,
        weights
      )[1]
      return(t_sq)
    }
  }
  else {
    # Construct function of gamma that return the z-value.
    t_sq_value = function(gamma) {
      t_sq = contrast_test_common(
        time_points,
        ctrl_estimates,
        exp_estimates,
        ref_fun,
        interpolation,
        vcov,
        gamma,
        type,
        j,
        weights
      )[1] ** 2
      return(t_sq)
    }
  }
  # Degrees of freedom for chi-squared statistic.
  if (type == "omnibus") df = length(j)
  else df = 1
  # Compute critical value.
  t_sq_critical = stats::qchisq(p = 1 - alpha, df = df)
  # If the test statistic evaluated in the estimated value is larger than the
  # critical value, then we cannot compute the confidence interval. This means
  # that a common acceleration factor is not consistent with the data.
  # Nonetheless, the contrast-based estimator remains well-defined. A confidence
  # interval can still be computed from the bootstrap. However, this should be
  # done with care.
  if (t_sq_value(gamma_est) > t_sq_critical) return(NA)
  # Find upper limit.

  # If the right limit of the test statistic, as a function of gamma, does not
  # cross the critical value, the upper confidence limit is infinity. The same
  # principle applies to the lower limit.
  if (t_sq_value(bounds[2]) < t_sq_critical) {
    warning("The test statistic does not cross the critical value in the search interval for computing the contrast-based confidence interval. \\
            The corresponding confidence limit is set to Inf. Consider increasing the search interval bounds using the `bounds` argument.")
    upper_limit = +Inf
  }
  else {
    upper_limit = stats::uniroot(
      f = function(gamma)
        sqrt(t_sq_value(gamma)) - sqrt(t_sq_critical),
      interval = c(gamma_est,
                   bounds[2]),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }

  # Find lower limit
  if (t_sq_value(bounds[1]) < t_sq_critical) {
    warning("The test statistic does not cross the critical value in the search interval for computing the contrast-based confidence interval. \\
            The corresponding confidence limit is set to Inf. Consider increasing the search interval bounds using the `bounds` argument.")
    lower_limit = -Inf
  }
  else {
    lower_limit = stats::uniroot(
      f = function(gamma)
        sqrt(t_sq_value(gamma)) - sqrt(t_sq_critical),
      interval = c(bounds[1],
                   gamma_est),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }

  # Return estimated confidence interval.
  return(
    c(lower_limit, upper_limit)
  )
}

#' Estimate the common acceleration factor by minimizing the squared contrast
#' statistic
#'
#' The [contrast_estimate_common()] function estimates the common acceleration
#' factor. The estimate is the common acceleration factor which minimizes the
#' test statistic for the corresponding null hypothesis, i.e., the most likely
#' value of \eqn{\gamma} where "likelihood" is defined in terms of a test
#' statistic.
#'
#' @inheritParams contrast_test_common
#' @inheritParams contrast_conf_int
#' @param ... Tuning parameters that are passed to [stats::optim()].
#' @param penalty This a function that is added to the (squared) test
#'   statistics. Defaults to a constant function. This is mostly useful for
#'   small samples to put a penalty on values of gamma outside the unit
#'   interval.
#'
#' @return (numeric) estimated for the common acceleration factor.
contrast_estimate_common = function(time_points,
                                 ctrl_estimates,
                                 exp_estimates,
                                 ref_fun,
                                 interpolation,
                                 vcov,
                                 type = "omnibus",
                                 j = 1:length(exp_estimates),
                                 weights = NULL,
                                 penalty = function(x) 0,
                                 bounds = c(-5, 5),
                                 ...) {
  # Help function that computes the squared test statistics for a given gamma.
  # If the test statistic is a chi-squared value, then we do not have to square.
  # For a z-statistic, we compute the square.
  if (type %in% c("omnibus", "directional")) {
    objective_function = function(gamma) {
      test_statistic = contrast_test_common(
        time_points,
        ctrl_estimates,
        exp_estimates,
        ref_fun,
        interpolation,
        vcov,
        gamma,
        type,
        j,
        weights
      )[1]
      return(test_statistic + penalty(gamma))
    }
  }
  else {
    objective_function = function(gamma) {
      test_statistic = contrast_test_common(
        time_points,
        ctrl_estimates,
        exp_estimates,
        ref_fun,
        interpolation,
        vcov,
        gamma,
        type,
        j,
        weights
      )[1] ** 2
      return(test_statistic + penalty(gamma))
    }
  }

  # Find the gamma-value that minimizes the (squared) test-statistic. Starting
  # value is found by simple grid search.
  gammas = seq(from = -0.5, to = 1.5, length.out = 20)
  objectives = sapply(X = gammas, FUN = objective_function)
  # Select starting value
  gamma_start = gammas[which.min(objectives)]
  gamma_est = stats::optim(
    par = gamma_start,
    fn = objective_function,
    hessian = FALSE,
    method = "L-BFGS-B",
    lower = bounds[1],
    upper = bounds[2],
    ...
  )$par

  # Raise warning of the estimate equals the bounds approximately.


  if (abs(gamma_est - bounds[1]) < 1e-4 || abs(gamma_est - bounds[2]) < 1e-4) {
    warning("The estimated common acceleration factor is equal to the lower or upper bound of the search interval. Consider increasing the search interval bounds using the `bounds` argument.")
  }
  return(gamma_est)
}

#' Standard Error of contrast-based estimator of the common acceleration factor
#'
#'
#' @param gamma_est Estimated value for the common acceleration factor.
#' @inheritParams contrast_estimate_common
#'
#' @return (numeric) Estimated SE of the estimator.
contrast_estimate_common_se = function(gamma_est,
                                    time_points,
                                    ctrl_estimates,
                                    exp_estimates,
                                    interpolation,
                                    vcov,
                                    type = "omnibus",
                                    j = 1:length(exp_estimates),
                                    weights = NULL,
                                    ...) {
  # Only the SE for the weights-based estimator has been implemented. For
  # estimators of which the SE has not yet been implemented, NA is returned.
  if (type != "custom") {
    return(NA)
  }

  # Compute derivative of psi with respect to gamma.
  deriv_psi_gamma = v_deriv_gamma(time_points,
                                  ctrl_estimates,
                                  interpolation,
                                  gamma_est,
                                  j,
                                  weights)

  # Compute derivative of psi with respect to (alpha, beta).
  grad_psi_alpha_beta = v_deriv_alpha_beta(time_points,
                                           ctrl_estimates,
                                           interpolation,
                                           gamma_est,
                                           j,
                                           weights)

  se = (deriv_psi_gamma**-2) * grad_psi_alpha_beta %*%
    vcov %*%
      t(grad_psi_alpha_beta)
  return(as.numeric(se))
}











