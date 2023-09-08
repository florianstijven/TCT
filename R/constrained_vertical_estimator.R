#' Re-estimate Vertical Parameters Under Constraints for TCT
#'
#' The [constrained_vertical_estimator()] function re-estimates the vertical
#' parameters under the constraints that are implicitly assumed by TCT-methods.
#'
#' @param alpha_obs Estimates for vertical parameters in the control group.
#' @param beta_obs Estimates for the vertical parameters in the experimental
#'   group.
#' @param Sigma_obs Variance-covariance matrix for
#'   \eqn{(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})'}.
#'
#' @details
#'
#' # Implicit Assumptions Made in TCT
#'
#' There are two implicit assumptions that underlie the TCT methodology:
#'
#'  * \eqn{\boldsymbol{\alpha}} and \eqn{\boldsymbol{\beta}} are strictly
#'    monotone sequences, i.e.,
#'    \deqn{\alpha_0 < \alpha_1 < ... < \alpha_k \; \text{and} \; \beta_0 < \beta_1 < ... < \beta_k}
#'  * The values for the trajectory in the active treatment group lie in the
#'    range of values for the trajectory in the placebo group. Especially
#'    important is the minimum: \eqn{\min(\boldsymbol{\alpha}) < \min(\boldsymbol{\beta})}.
#'
#'  The available vertical parameter estimates may not satisfy the above
#'  assumptions. We can "re-estimate" the vertical parameters such that they
#'  satisfy the above assumptions.  This re-estimation procedure is explained
#'  in the next section.
#'
#'  # Re-Estimation Procedure
#'
#'  Let \eqn{N((\boldsymbol{\alpha}_0, \boldsymbol{\beta}_0)', \Sigma)} be the
#'  sampling distribution of the (unmodified) vertical parameter estimators. Let
#'  \eqn{(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})' = \boldsymbol{y}}
#'  be a sample from this sampling distribution. The generalized least squares
#'  estimator for \eqn{(\boldsymbol{\alpha}_0, \boldsymbol{\beta}_0)'} is simply
#'  \eqn{\boldsymbol{y}}. However, \eqn{\boldsymbol{y}} may not satisfy the
#'  above two assumptions, even though they hold for the true parameter vector
#'  \eqn{(\boldsymbol{\alpha}_0, \boldsymbol{\beta}_0)'}. Next, we modify the
#'  generalized least squares estimator to an estimator that satisfies those
#'  assumptions. This is termed the *constrained generalized least squares estimator*,
#'  denoted by \eqn{(\hat{\boldsymbol{\alpha}}_{c}, \hat{\boldsymbol{\beta}}_{c})'}.
#'
#'  The generalized least squares estimator minimizes the following objective function
#'  \deqn{O(\boldsymbol{\alpha}, \boldsymbol{\beta}) = (\boldsymbol{y} - (\boldsymbol{\alpha}, \boldsymbol{\beta})')' \cdot \Sigma^{-1} \cdot (\boldsymbol{y} - (\boldsymbol{\alpha}, \boldsymbol{\beta})').}
#'  This immediately leads us to the constrained generalized least-squares
#'  estimator. This estimator is found as the solution to the following
#'  minimiation problem,
#'  \deqn{(\hat{\boldsymbol{\alpha}}_{c}, \hat{\boldsymbol{\beta}}_{c})' = \arg \; \min_{(\boldsymbol{\alpha}, \boldsymbol{\beta})'}   O(\boldsymbol{\alpha}, \boldsymbol{\beta}),}
#'  under the following linear constraints,
#'  * \eqn{\alpha_0 < \alpha_1 < ... < \alpha_k}
#'  * \eqn{\beta_0 < \beta_1 < ... < \beta_k}
#'  * \eqn{\alpha_0 \le \beta_0}
#'
#'  If the original estimates already satisfy the above constraints, then the
#'  original estimates are of course returned.
#'
#' @return (numeric) Re-estimated parameter vector under the linear
#' constraints that are implicitly assumed by TCT.
#'   \eqn{(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})'} under the
#'   appropriate constraints.
#' @export
#'
#' @examples
#'
#' # transform example data set to desired format
#' library(dplyr)
#' data = simulated_test_trial %>%
#'   mutate(time_int = (Week %/% 25)) %>%
#'   arrange(trial_number, SubjId, time_int) %>%
#'   mutate(time_int = as.integer(time_int) + 1L) %>%
#'   mutate(arm_time = ifelse(time_int == 1L,
#'                            "baseline",
#'                            paste0(arm, ":", time_int)))
#' # fit e.g. MMRM model to obtain estimates of profiles
#' mmrm_fit = analyze_mmrm(data)
#' # Re-estimate parameter vector under the constrainst that are implictly
#' # assumed by TCT.
#' re_estimated_vector = constrained_vertical_estimator(
#'   alpha_obs = coef(mmrm_fit)[c(9, 1:4)],
#'   beta_obs = coef(mmrm_fit)[5:8],
#'   Sigma_obs = vcov(mmrm_fit)
#' )
constrained_vertical_estimator = function(alpha_obs, beta_obs, Sigma_obs) {
  # Compute the inverse covariance matrix. This leads to efficiency gains since
  # the inverse would otherwise be computed multiple time.
  Sigma_inv = solve(Sigma_obs)
  length_alpha = length(alpha_obs)
  obj_function = function(theta, alpha_beta_obs) {
    alpha = theta[1:length_alpha]
    beta = theta[(length_alpha + 1):length(theta)]
    mahalanobis(
      x = alpha_beta_obs,
      center = c(alpha, beta),
      cov = Sigma_inv,
      inverted = TRUE
    )
  }
  # Constraints matrix for first assumption: strict monotonically increasing
  # trajectory.
  A1 = constraint_matrix_A1(length(alpha_obs), length(beta_obs))
  # Constrain matrix for second assumption: min(alpha) < min(beta).
  A2 = constraint_matrix_A2(length(alpha_obs), length(beta_obs))
  # Combine both constraint matrices.
  constraint_matrix = rbind(A1, A2)
  # Define constraint vector. Since all constraint are inequalities without
  # constants, this is just the zero vector.
  constraint_vector = matrix(0, nrow = nrow(constraint_matrix), ncol = 1)

  # Custom starting values that automatically satisfy the linear inequalities.
  starting_values = starting_values_constrained_optimization(alpha_obs, beta_obs)

  stats::constrOptim(
    theta = starting_values,
    f = obj_function,
    grad = gradient_mahalanobis,
    ui = constraint_matrix,
    ci = constraint_vector,
    alpha_beta_obs = c(alpha_obs, beta_obs)
  )$par
}

# Return constraint matrix for Assumption 1
constraint_matrix_A1 = function(length_alpha, length_beta) {
  A1_alpha = -1 * diag(1, nrow = length_alpha)+
    matrix(c(0, rep(c(1, rep(0, length_alpha)), length_alpha - 1)),
           ncol = length_alpha, byrow = TRUE)
  A1_alpha = A1_alpha[-length_alpha, ]
  A1_beta = -1 * diag(1, nrow = length_beta) +
    matrix(c(0, rep(c(1, rep(0, length_beta)), length_beta - 1)),
           ncol = length_beta, byrow = TRUE)
  A1_beta = A1_beta[-length_beta, ]
  # Define empty A1 matrix with the correct dimensions.
  A1_matrix = matrix(0,
                     ncol = length_alpha + length_beta,
                     nrow = length_alpha + length_beta - 2)
  # Plug in the previously defined block matrices.
  A1_matrix[1:(length_alpha - 1), 1:length_alpha] = A1_alpha
  A1_matrix[(length_alpha):(length_alpha + length_beta - 2),
            (length_alpha + 1):(length_alpha + length_beta)] = A1_beta
  return(A1_matrix)
}

# Return Constraint matrix for Assumption 2
constraint_matrix_A2 = function(length_alpha, length_beta) {
  A2_matrix = matrix(
    c(-1, rep(0, length_alpha - 1), 1, rep(0, length_beta - 1)),
    nrow = 1
  )
  return(A2_matrix)
}

gradient_mahalanobis = function(alpha_beta, alpha_beta_obs) {
  y = alpha_beta_obs
  return(-2 * y + 2 * alpha_beta)
}

starting_values_constrained_optimization = function(alpha_obs, beta_obs) {
  # order both vectors to satisify assumption 1
  alpha_obs = sort(alpha_obs)
  beta_obs = sort(beta_obs)
  for (i in seq_along(beta_obs)) {
    if (beta_obs[i] < alpha_obs[1]) {
      beta_obs[i] = alpha_obs[1] + i * 1e-5
    }
  }
  return(c(alpha_obs, beta_obs))
}
