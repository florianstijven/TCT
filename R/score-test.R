#' Compute TCT score test z-value
#'
#' The [score_test()] function computes the z-value for the score test at
#' measurement `j`. This corresponds to `exp_estimates[j]`. For technical
#' details on the score test, see Details.
#'
#' @param j Measurement occasion to test acceleration factor for. This
#'   corresponds to `exp_estimates[j]`.
#' @param gamma_0 Value for the acceleration factor under the null hypothesis.
#' @inheritParams DeltaMethod
#'
#' @details
#'
#' # Score Test
#'
#' For constructing a score test, we start from the null hypothesis,
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
#' score test that draws upon well-established statistical tests.
#'
#' We first define the score vector as follows, \deqn{\boldsymbol{s}(\gamma
#' \cdot \boldsymbol{t}; \hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}}) =
#' \hat{\boldsymbol{\beta}} - \boldsymbol{f_0}(\gamma \cdot \boldsymbol{t};
#' \hat{\boldsymbol{\alpha}})}
#' where \eqn{\boldsymbol{t} = (t_1, ..., t_K)'} and \eqn{\boldsymbol{f_0}(\gamma \cdot
#' \boldsymbol{t}; \hat{\boldsymbol{\alpha}}) = (f_0(\gamma \cdot t_1;
#' \hat{\boldsymbol{\alpha}}), ..., f_0(\gamma \cdot t_K;
#' \hat{\boldsymbol{\alpha}}))'}. Under \eqn{H_{\gamma_0}} for all time points and using the delta
#' method, we have that \deqn{\boldsymbol{s}(\gamma \cdot \boldsymbol{t};
#' \hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}}) \; \dot\sim \; N \left(
#' \boldsymbol{0}, \Sigma_s \right)} where \deqn{\Sigma_s =
#' J_s(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0}; \gamma_0) \cdot \Sigma \cdot
#' J_s(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0}; \gamma_0)^t} with
#' \eqn{J_s(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0}; \gamma_0)} the
#' Jacobian of \eqn{\boldsymbol{s}(\gamma \cdot \boldsymbol{t}; \boldsymbol{\alpha},
#' \boldsymbol{\beta})} evaluated in the true parameter vector,
#' \eqn{(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})'}, with \eqn{\gamma} held
#' constant at \eqn{\gamma_0}. Note that this matrix is a function of
#' \eqn{\boldsymbol{\alpha}}, \eqn{\boldsymbol{\beta}}, and \eqn{\gamma_0}. The
#' Jacobian matrix with the true parameter vector replaced by the corresponding
#' estimates, \eqn{(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})'}, is
#' denoted by \eqn{\hat{\Sigma}_s}.
#'
#' # Time-Specific Acceleration Factor
#'
#' From the distributional result in the previous section follows that \deqn{z =
#' \frac{s(\gamma_0 \cdot t_j; \hat{\boldsymbol{\alpha}},
#' \hat{\boldsymbol{\beta}})}{\sqrt{\hat{\Sigma}_{s,jj}}} \; \dot\sim \;N(0, 1)}
#' under \eqn{H_{\gamma_0}} where \eqn{\hat{\Sigma}_{s,jj}} is the \eqn{j}'th
#' diagonal element of \eqn{\hat{\Sigma}_s}. Note that this result relies on the
#' delta method and is thus only asymptotic. The finite sample properties of
#' this test will thus depend to a large extent on the accuracy of the delta
#' method, i.e., the degree of non-linearity of
#' \eqn{g_{\gamma_0}(\boldsymbol{\alpha}, \boldsymbol{\beta}; t_j)} around the
#' true parameter values. The [score_test()] function returns the above z-value.
#'
#' @return Named (numeric) vector with two elements:
#'
#' 1. `"z"`. This test-statistic follows a standard normal distribution under
#' the null hypothesis.
#' 2. `"p-value"`: two-sided p-value.
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
  # Return z-value and corresponding two-sided p-value.
  return(c("z" = as.numeric(z),
           "p-value" = 2 * (1 - stats::pnorm(abs(as.numeric(z))))
           )
         )
}

#' Compute Confidence Interval Based on Score Test
#'
#' The [score_conf_int()] function computes the confidence interval for the
#' score test at measurement `j`. This confidence interval is based on the
#' score test implemented in [score_test()].
#'
#' @param alpha `1 - alpha` represents the two-sided confidence level. Defaults
#'   to `0.05`.
#' @inheritParams score_test
#'
#' @details
#'
#' # Score Test Confidence Intervals
#'
#' For the construction of a confidence interval based on the scores, we make
#' use of the relationship between confidence intervals and hypothesis tests. A
#' \eqn{1 − α} confidence interval can be defined as follows, \deqn{\left\{
#' \gamma : p(\gamma) > \alpha \right\}} where \eqn{p(\gamma)} is the p-value
#' for the null that the acceleration factor is equal to \eqn{\gamma}. This
#' p-value can be based on the time-specific score test implemented in
#' [score_test()], or can be based on the joint score test implemented in
#' [score_test_common()]. The former gives a confidence interval for the
#' time-specific acceleration factor while the latter gives a confidence
#' interval for the common acceleration factor. The latter interval should be
#' used with some care. Indeed, the confidence interval for a common
#' acceleration factor is only sensible if there is no evidence against the
#' assumptions of constant slowing of the disease progression.
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
                      gamma_0 = gamma)[1])
  }
  # Find upper limit
  z_critical = stats::qnorm(p = 1 - alpha / 2)
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

#' Score test for common acceleration factor
#'
#' The [score_test_common()] function implements the score test under the
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
#'   should be used for the score test. Defaults to `1:length(exp_estimates)`,
#'   i.e., all elements are used.
#' @param weights If `type == "custom"`, the user should specify a weight
#'   vector for weighting estimates at different time points differently.
#' @inheritParams score_test
#' @inheritSection score_test Score Test
#' @details
#'
#' # Test Statistic Variants
#'
#' Two types of test statistics are implemented in the [score_test_common()]
#' function:
#'
#' 1. `type = "omnibus"`
#' 2. `type = "custom"`
#'
#' These are discussed in more detail next.
#'
#' ## Omnibus
#'
#' The omnibus score test is based on the "classic" chi-squared statistic that is
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
#' The "custom" score test allows for user-specified weights, \eqn{w}. The test
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
score_test_common = function(time_points,
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
  J = score_vector_jacobian(
    time_points,
    ctrl_estimates,
    interpolation,
    gamma_0,
    j
  )
  # Compute the variance of the "test statistic "score vector", g_gamma_0 under
  # the null. The inverse of this matrix is also computed.
  Sigma_g = J %*% vcov[c(1:length(time_points), length(time_points) + j),
                          c(1:length(time_points), length(time_points) + j)] %*% t(J)
  Sigma_g_inv = solve(Sigma_g)
  # Compute the score vector.
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



#' Confidence interval based on score test
#'
#' The [score_conf_int_common()] function computes the confidence interval for
#' common acceleration factor. This confidence interval is based on the score
#' test implemented in [score_test_common()].
#'
#' @inheritParams score_estimate_common
#' @inheritParams score_conf_int
#' @param gamma_est Estimate for the common acceleration factor.
#' @inheritSection score_conf_int Score Test Confidence Intervals
#' @inherit score_conf_int return
score_conf_int_common = function(time_points,
                                 ctrl_estimates,
                                 exp_estimates,
                                 ref_fun,
                                 interpolation,
                                 vcov,
                                 gamma_est,
                                 type = "omnibus",
                                 j = 1:length(exp_estimates),
                                 weights = NULL,
                                 alpha = 0.05){
  # Force argument values. This is required because we're using these arguments
  # in a function factory.
  force(time_points); force(ctrl_estimates); force(exp_estimates)
  force(ref_fun); force(interpolation); force(vcov); force(gamma_est)
  force(type); force(j); force(weights); force(alpha)
  if (type %in% c("omnibus", "directional")) {
    # Construct function of gamma that return the z-value.
    t_sq_value = function(gamma) {
      t_sq = score_test_common(
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
      t_sq = score_test_common(
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
  # Nonetheless, the score-based estimator remains well-defined. A confidence
  # interval can still be computed from the bootstrap. However, this should be
  # done with care.
  if (t_sq_value(gamma_est) > t_sq_critical) return(NA)
  # Find upper limit.

  # If the right limit of the test statistic, as a function of gamma, does not
  # cross the critical value, the upper confidence limit is infinity. The same
  # principle applies to the lower limit.
  if (t_sq_value(10) < t_sq_critical) {
    upper_limit = +Inf
  }
  else {
    upper_limit = stats::uniroot(
      f = function(gamma)
        sqrt(t_sq_value(gamma)) - sqrt(t_sq_critical),
      interval = c(gamma_est,
                   10),
      tol = .Machine$double.eps ^ 0.5,
      maxiter = 1e3
    )$root
  }

  # Find lower limit
  if (t_sq_value(-10) < t_sq_critical) {
    lower_limit = -Inf
  }
  else {
    lower_limit = stats::uniroot(
      f = function(gamma)
        sqrt(t_sq_value(gamma)) - sqrt(t_sq_critical),
      interval = c(-10,
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

#' Estimate the common acceleration factor by minimizing the squared score
#' statistic
#'
#' The [score_estimate_common()] function estimates the common acceleration
#' factor. The estimate is the common acceleration factor which minimizes the
#' test statistic for the corresponding null hypothesis, i.e., the most likely
#' value of \eqn{\gamma} where "likelihood" is defined in terms of a test
#' statistic.
#'
#' @inheritParams score_test_common
#' @param ... Tuning parameters that are passed to [stats::optim()].
#' @param penalty This a function that is added to the (squared) test
#'   statistics. Defaults to a constant function. This is mostly useful for
#'   small samples to put a penalty on values of gamma outside the unit
#'   interval.
#'
#' @return (numeric) estimated for the common acceleration factor.
score_estimate_common = function(time_points,
                                 ctrl_estimates,
                                 exp_estimates,
                                 ref_fun,
                                 interpolation,
                                 vcov,
                                 type = "omnibus",
                                 j = 1:length(exp_estimates),
                                 weights = NULL,
                                 penalty = function(x) 0,
                                 ...) {
  # Help function that computes the squared test statistics for a given gamma.
  # If the test statistic is a chi-squared value, then we do not have to square.
  # For a z-statistic, we compute the square.
  if (type %in% c("omnibus", "directional")) {
    objective_function = function(gamma) {
      test_statistic = score_test_common(
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
      test_statistic = score_test_common(
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
  stats::optim(
    par = gamma_start,
    fn = objective_function,
    hessian = FALSE,
    method = "L-BFGS-B",
    lower = -2,
    upper = 4,
    ...
  )$par
}

#' Standard Error of Score-based estimator of the common acceleration factor
#'
#'
#' @param gamma_est Estimated value for the common acceleration factor.
#' @inheritParams score_estimate_common
#'
#' @return (numeric) Estimated SE of the estimator.
score_estimate_common_se = function(gamma_est,
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
  gr_gamma_w = gradient_gamma_w_analytical(
    time_points = time_points,
    ctrl_estimates = ctrl_estimates,
    exp_estimates = exp_estimates,
    interpolation = interpolation,
    gamma_0 = gamma_est,
    j = j,
    weights = weights
  )
  se = t(gr_gamma_w) %*%
    vcov %*%
    gr_gamma_w
  return(as.numeric(se))
}











