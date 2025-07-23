#' Parametric Bootstrap for Time-Specific Acceleration Factors
#'
#' The [pm_bootstrap_vertical_to_horizontal()] function implements a parametric
#' bootstrap for the time-specific acceleration factors.
#'
#' @param B Number of bootstrap replications.
#' @param null (boolean): conduct the bootstrap under the null hypothesis of no
#'   treatment effect? Defaults to `FALSE`.
#' @inheritParams TCT_meta
#'
#' @details
#' We can write the vector of time-specific acceleration factors
#' as a vector-valued function of the mean parameters:
#' \eqn{
#' (\boldsymbol{\alpha}, \boldsymbol{\beta})
#' \mapsto
#' \boldsymbol{\gamma}(\boldsymbol{\alpha}, \boldsymbol{\beta})
#' },
#' where we assume asymptotic normality of the mean parameter estimator:
#'\deqn{
#' n^{1/2}
#' \left((\hat{\boldsymbol{\alpha}}_n, \hat{\boldsymbol{\beta}}_n)' -
#' (\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})' \right)
#' \overset{d}{\to}
#' \mathcal{N}
#' \left( \boldsymbol{0}, \Sigma \right),
#' }
#' where \eqn{(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})'} is the true
#' parameter vector. We let \eqn{\Sigma_n} be a consistent estimator of
#' \eqn{\Sigma}.
#'
#' The [pm_bootstrap_vertical_to_horizontal()] function implements a parametric
#' bootstrap in the following two steps:
#' 1. Sample \eqn{B} bootstrap replications of the vertical parameters from
#' \deqn{
#' (\hat{\boldsymbol{\alpha}}_b, \hat{\boldsymbol{\beta}}_b)'
#' \sim
#' \mathcal{N} \left( (\hat{\boldsymbol{\alpha}}_n, \hat{\boldsymbol{\beta}}_n)', \Sigma_n  \right)
#' }
#' where \eqn{b} refers to the \eqn{b}'th bootstrap replicate.
#' 2. Transform the bootstrap replicates to the time scale,
#' \deqn{
#' \hat{\boldsymbol{\gamma}}_b :=
#' \boldsymbol{\gamma}(\hat{\boldsymbol{\alpha}}_b, \hat{\boldsymbol{\beta}}_b).
#' }
#' This transformation is implemented in [g_Delta_bis()] for the time-specific
#' acceleration factors.
#'
#' For the common acceleration factor, we replace the vector-valued function
#' \eqn{
#' (\boldsymbol{\alpha}, \boldsymbol{\beta})
#' \mapsto
#' \boldsymbol{\gamma}(\boldsymbol{\alpha}, \boldsymbol{\beta})
#' }
#' with the real-valued function
#' \eqn{
#' (\boldsymbol{\alpha}, \boldsymbol{\beta})
#' \mapsto
#' \gamma(\boldsymbol{\alpha}, \boldsymbol{\beta})
#' }
#' and repeat the same steps as described above. The nature of this real-valued
#' function depends on the estimator for the common acceleration factor.
#'
#' @return A \eqn{(B \times K)} matrix where the \eqn{b}'th row corresponds to
#' \eqn{\hat{\boldsymbol{\gamma}}^b} and \eqn{K = } `length(exp_estimates)`.
#'
#' @importFrom mnormt pd.solve
#' @importFrom mvtnorm rmvnorm
pm_bootstrap_vertical_to_horizontal = function(time_points,
                                               ctrl_estimates,
                                               exp_estimates,
                                               vcov,
                                               interpolation = "spline",
                                               B = 100,
                                               null = FALSE,
                                               constraints) {
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

  # If constraints option is selected, first apply the constrained GLS
  # modification.
  if (constraints) {
    length_alpha = length(ctrl_estimates)
    par_sampled = apply(
      X = par_sampled,
      MARGIN = 1,
      FUN = function(x) {
        constrained_vertical_estimator(x[1:length_alpha],
                                       x[-1 * (1:length_alpha)],
                                       vcov)
      }
    )
    par_sampled = t(par_sampled)
  }

  estimates = matrix(0, nrow = B, ncol = length(exp_estimates))
  for (i in 1:B) {
    estimates[i, ] = g_Delta_bis(par = par_sampled[i, ],
                                 interpolation = interpolation,
                                 time_points = time_points)
  }
  return(estimates)
}


#' Parametric Bootstrap for a Common Acceleration Factor
#'
#' The [pm_bootstrap_vertical_to_common()] function implements a parametric
#' bootstrap for the common acceleration factor.
#'
#' @param return_se (boolean) Return the estimated standard error from each
#' bootstrap replication? This standard error is computed with the delta method.
#' @inheritParams pm_bootstrap_vertical_to_horizontal
#' @inheritParams TCT_meta_common
#' @return A list with two elements:
#'  * `estimates_bootstrap`: (numeric) vector of length `B` that contains the
#'  bootstrap replicates for the common acceleration factor.
#'  * `se_bootstrap`: (numeric) vector of length `B` that contains the bootstrap
#'  replicates of the estimated standard error. Is a vector of `NA`s when
#'  `return_se = FALSE`.
#'
#' @inherit pm_bootstrap_vertical_to_horizontal
#'

pm_bootstrap_vertical_to_common = function(TCT_Fit,
                                           inference,
                                           B = 100,
                                           return_se = TRUE,
                                           null = FALSE,
                                           select_coef,
                                           constraints = FALSE,
                                           weights,
                                           start_gamma) {
  if (B == 0)
    return(NULL)

  # Extract information from the TCT_meta object that is used further on.
  ctrl_estimates = TCT_Fit$vertical_model$ctrl_estimates
  exp_estimates = TCT_Fit$vertical_model$exp_estimates
  time_points = TCT_Fit$vertical_model$time_points
  interpolation = TCT_Fit$inference_options$interpolation
  vcov_vertical = TCT_Fit$vertical_model$vcov
  TCT_vcov = TCT_Fit$vcov

  n_points = length(time_points)
  p = length(exp_estimates[select_coef])
  vec_1 = matrix(1, nrow = p, ncol = 1)
  estimates_bootstrap = 1:B
  se_bootstrap = rep(NA, B)

  vcov_vertical = vcov_vertical[c(1:n_points, n_points + 1:length(coef(TCT_Fit))),
                                c(1:n_points, n_points + 1:length(coef(TCT_Fit)))]
  if (null) {
    # estimated covariance matrix of arm-specific vertical estimates
    vcov_gls = vcov_vertical[-1, -1]
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
      sigma = vcov_vertical
    )
  }
  else {
    par_sampled = mvtnorm::rmvnorm(
      n = B,
      mean = c(ctrl_estimates, exp_estimates),
      sigma = vcov_vertical
    )
  }

  # If constraints option is selected, first apply the constrained GLS
  # modification.
  if (constraints) {
    length_alpha = length(ctrl_estimates)
    par_sampled = apply(
      X = par_sampled,
      MARGIN = 1,
      FUN = function(x) {
        constrained_vertical_estimator(x[1:n_points],
                                       x[n_points + 1:length(coef(TCT_Fit))],
                                       vcov_vertical)
      }
    )
    par_sampled = t(par_sampled)
  }

  # Re-estimate common acceleration factor for Wald-based inference
  if (inference == "delta-method") {
    # Do not fix the variance-covariance matrix for the time-specific
    # acceleration factors across bootstrap replications.
    bs_fix_vcov = FALSE
    for (i in seq_along(estimates_bootstrap)) {
      if (bs_fix_vcov) {
        vcov_gls = TCT_vcov[select_coef, select_coef]
        coef_gls = g_Delta_bis(par = par_sampled[i, ],
                               time_points = time_points,
                               interpolation = interpolation)[select_coef]
      }
      else {
        tct_results = TCT_meta(
          time_points = time_points,
          ctrl_estimates = par_sampled[i, 1:n_points],
          exp_estimates = par_sampled[i, (n_points + 1):length(par_sampled[1, ])],
          vcov = vcov_vertical,
          interpolation = interpolation,
          B = 0,
          constraints = constraints
        )
        if (!is.na(tct_results$vcov[[1]])) {
          vcov_gls = tct_results$vcov[select_coef, select_coef]
        }
        coef_gls = stats::coef(tct_results)[select_coef]
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
  }

  else if (inference %in% c("contrast", "least-squares")){
    bootstrap_replications_model = apply(
      X = par_sampled,
      MARGIN = 1,
      FUN = function(par_sampled) {
        ctrl_estimates = par_sampled[1:n_points]
        exp_estimates = par_sampled[-1 * (1:n_points)]
        # Estimate TCT meta model.
        TCT_Fit_replicate = TCT_meta(
          time_points = time_points,
          ctrl_estimates = ctrl_estimates,
          exp_estimates = exp_estimates,
          vcov = vcov_vertical,
          interpolation = interpolation,
          # Inference in TCT_meta() does not matter when calling
          # TCT_meta_common() further on.
          inference = "delta-method",
          B = 0,
          # We can set the constraints to zero because this has already been
          # checked while drawing the vertical parameter estimates.
          constraints = FALSE
        )
        # Estimate meta TCT common model.
        TCT_common_fit = TCT_meta_common(
          TCT_Fit = TCT_Fit_replicate,
          inference = inference,
          B = 0,
          select_coef = select_coef,
          constraints = FALSE,
          weights = weights,
          start_gamma = start_gamma
          )
        return(TCT_common_fit)
      }
    )
    # Extract Estimates and corresponding standard errors.
    estimates_bootstrap = sapply(
      X = bootstrap_replications_model,
      FUN = function(x) stats::coef(x)
    )
    se_bootstrap = sapply(
      X = bootstrap_replications_model,
      FUN = function(x) sqrt(stats::vcov(x))
    )
  }

  return(list(estimates_bootstrap  = estimates_bootstrap,
              se_bootstrap = se_bootstrap))
}
