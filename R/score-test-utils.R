gradient_gamma_w = function(time_points,
                            ctrl_estimates,
                            exp_estimates,
                            interpolation,
                            vcov,
                            gamma_0,
                            j = 1:length(exp_estimates),
                            weights,
                            analytical = FALSE){
  K = length(time_points)
  force(time_points); force(ctrl_estimates); force(exp_estimates)
  force(interpolation); force(vcov);
  force(j); force(weights); force(K); force(ref_fun); force(gamma_0)
  v = function(theta) {
    v_function(
      time_points,
      theta[1:K],
      theta[K + j],
      interpolation,
      vcov,
      gamma_0,
      j,
      weights)
  }

  f = function(theta) {
    v_gamma = function(gamma) {
      v_function(
        time_points,
        theta[1:K],
        theta[K + j],
        interpolation,
        vcov,
        gamma,
        j,
        weights)
    }
    v_deriv_gamma = numDeriv::grad(
      func = v_gamma,
      x = gamma_0
    )
    v(theta) * (1 / v_deriv_gamma)
  }
  numDeriv::grad(
    func = f,
    x = c(ctrl_estimates, exp_estimates[j])
  )
}

gradient_gamma_w_analytical = function(time_points,
                                       ctrl_estimates,
                                       exp_estimates,
                                       interpolation,
                                       vcov,
                                       gamma_0,
                                       j = 1:length(exp_estimates),
                                       weights) {

}

#' Compute Jacobian of the score vector
#'
#' The [score_vector_jacobian()] function computes the Jacobian of the score
#' vector with respect to \eqn{(\boldsymbol{\alpha},\boldsymbol{\beta})'},
#' \deqn{\frac{\partial \boldsymbol{s}(\gamma \cdot \boldsymbol{t}; \boldsymbol{\alpha}, \boldsymbol{\beta})}{\partial (\boldsymbol{\alpha}, \boldsymbol{\beta})}}.
#'
#' @inheritParams score_test_common
#' @return
score_vector_jacobian = function(time_points,
                                 ctrl_estimates,
                                 interpolation,
                                 gamma_0,
                                 j) {
  # Number of measurement occasions in the experimental group.
  K = length(j)
  # The Jacobian matrix is computed in two parts. First, the (K x K + 1) part is
  # computed.
  J_f0_alpha_t  = -1 * attr(
    deriv_f0_alpha(
      t_m = gamma_0 * time_points[j + 1],
      x_ref = time_points,
      y_ref = ctrl_estimates,
      method = interpolation
    ),
    "gradient"
  )
  # Next, the second part of the Jacobian. This corresponds to the identity
  # matrix. We join both parts to get the Jacobian.
  J = cbind(J_f0_alpha_t, diag(x = 1, nrow = K))
  return(J)
}


#' Derivative of v-function with respect to \eqn{(\boldsymbol{\alpha}, \boldsymbol{\beta})'}
#'
#' This function computes
#' \eqn{\frac{\partial v(\boldsymbol{\alpha}, \boldsymbol{\beta}; \gamma_0)}{\partial (\boldsymbol{\alpha}, \boldsymbol{\beta})}}.
#'
#' @inheritParams score_test_common
#'
#' @return
v_deriv_alpha_beta = function(time_points,
                              ctrl_estimates,
                              exp_estimates,
                              interpolation,
                              vcov,
                              gamma_0,
                              j,
                              weights) {
  ref_fun = ref_fun_constructor(time_points,
                                ctrl_estimates,
                                interpolation)
  # Compute the jacobian of the score vector.
  J = score_vector_jacobian(time_points,
                            ctrl_estimates,
                            exp_estimates,
                            interpolation,
                            vcov,
                            gamma_0,
                            j,
                            weights)
  v_deriv = t(weights) %*% J
  return(v_deriv)
}

v_function = function(time_points,
                              ctrl_estimates,
                              exp_estimates,
                              interpolation,
                              vcov,
                              gamma_0,
                              j,
                              weights) {
  ref_fun = ref_fun_constructor(time_points,
                                ctrl_estimates,
                                interpolation)
  # Compute the jacobian of the score vector.
  J = score_vector_jacobian(time_points,
                            ctrl_estimates,
                            exp_estimates,
                            interpolation,
                            vcov,
                            gamma_0,
                            j,
                            weights)
  # Compute the variance of the "test statistic "score vector", g_gamma_0 under
  # the null. The inverse of this matrix is also computed.
  Sigma_g = J %*% vcov[c(1:length(time_points), length(time_points) + j),
                       c(1:length(time_points), length(time_points) + j)] %*% t(J)
  Sigma_g_inv = solve(Sigma_g)
  # Compute the score vector.
  g = (exp_estimates[j] - ref_fun(gamma_0 * time_points[j + 1]))
  weights = matrix(weights, ncol = 1)
  v = t(weights) %*% g
  return(as.numeric(v))
}
