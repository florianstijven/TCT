#' Compute TCT score test z-value
#'
#' The [score_test()] function computes the z-value for the score test at
#' measurement `j`. This corresponds to `exp_estimates[j]`.
#'
#' @param j Measurement occasion to test acceleration factor for.
#' @param gamma_0 Value under the null hypothesis for the acceleration factor.
#' @inheritParams DeltaMethod
#'
#'
#' @return (numeric) z-value
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
  # Return z-value
  return(as.numeric(z))
}

score_conf_int = function(){

}

score_test_common = function(){

}

score_conf_int_common = function(){

}

#' Compute Gradient for Score Test
#'
#' The [gradient_g()] function computes
#' \eqn{\nabla g_{\gamma}(\boldsymbol{\alpha}, \boldsymbol{\beta}; t_j)^t} where
#' \eqn{g_{\gamma}(\boldsymbol{\alpha}, \boldsymbol{\beta}; t_j) = \beta_j - f_0( \gamma \cdot t_{j}; \boldsymbol{\alpha})}.
#'
#' @return
#' @export
#'
#' @examples
gradient_g = function(t_m, x_ref, y_ref, method = "spline" ){

}
