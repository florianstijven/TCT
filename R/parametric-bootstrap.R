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
#' We can write the slowing factor at \eqn{\boldsymbol{t} = (t_1, ..., t_K)'}
#' as a function of the the mean parameters,
#' \deqn{
#' \boldsymbol{\gamma} = \boldsymbol{\gamma}(\boldsymbol{t}; \boldsymbol{\alpha}, \boldsymbol{\beta}).
#' }
#' Further assume that \eqn{\hat{\boldsymbol{\alpha}}} and \eqn{\hat{\boldsymbol{\beta}}}
#' are normally distributed,
#'\deqn{
#' (\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})' \sim N\left( (\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})', D \right)
#' }
#' where \eqn{(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})'} is the true
#' parameter vector. We can obtain an approximate sampling distribution for
#' \eqn{\hat{\boldsymbol{\gamma}}} by a parametric bootstrap in two steps,
#' 1. Sample \eqn{B} bootstrap replications of the vertical parameters from
#' \deqn{
#' (\hat{\boldsymbol{\alpha}}^b, \hat{\boldsymbol{\beta}}^b)' \sim N \left( (\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})', \hat{D}  \right)
#' }
#'  where \eqn{b} refers to the \eqn{b}'th bootstrap replicate.
#' 2. Transform the bootstrap replicates to the time scale,
#' \deqn{
#' \hat{\boldsymbol{\gamma}}^b = \boldsymbol{\gamma}(\boldsymbol{t}; \hat{\boldsymbol{\alpha}}^b, \hat{\boldsymbol{\beta}}^b).
#' }
#' This transformation is implemented in [g_Delta_bis()].
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
                                               constraints,
                                               type) {
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
#' bootstrap for the time-specific acceleration factors.
#'
#' @param bs_fix_vcov (boolean) Fix the estimated variance-covariance matrix for
#' the estimated acceleration factors? This speeds up computations, but can have
#' a negative impact of the method's properties.
#' @param return_se (boolean) Return the estimated standard error from each
#' bootstrap replication? This standard error is computed with the delta method.
#' @inheritParams pm_bootstrap_vertical_to_horizontal
#' @inheritParams TCT_meta_common
#'
#' @return A list with two element:
#'  * `estimates_bootstrap`: (numeric) vector of length `B` that contains the
#'  bootstrap replicates for the common acceleration factor.
#'  * `se_bootstrap`: (numeric) vector of length `B` that contains the bootstrap
#'  replicates of the estimated standard error. Is a vector of `NA`s when
#'  `return_se = FALSE`.
pm_bootstrap_vertical_to_common = function(TCT_Fit,
                                           inference,
                                           B = 100,
                                           bs_fix_vcov = TRUE,
                                           return_se = TRUE,
                                           null = FALSE,
                                           select_coef,
                                           constraints = FALSE,
                                           type,
                                           weights) {
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
        constrained_vertical_estimator(x[1:length_alpha],
                                       x[-1 * (1:length_alpha)],
                                       vcov_vertical)
      }
    )
    par_sampled = t(par_sampled)
  }

  # Re-estimate common acceleration factor for Wald-based inference
  if (inference == "wald") {
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

  else if (inference %in% c("score", "least-squares")){
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
          inference = "wald",
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
          bs_fix_vcov = FALSE,
          select_coef = select_coef,
          constraints = FALSE,
          type = type,
          weight = weights
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
