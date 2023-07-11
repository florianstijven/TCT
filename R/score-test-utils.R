gradient_gamma_w_numerical = function(time_points,
                            ctrl_estimates,
                            exp_estimates,
                            interpolation,
                            gamma_0,
                            j = 1:length(exp_estimates),
                            weights){
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
        gamma,
        j,
        weights)
    }
    v_deriv_gamma = numDeriv::grad(
      func = v_gamma,
      x = gamma_0
    )
    -1 * v(theta) * (1 / v_deriv_gamma)
  }
  gradient = numDeriv::grad(
    func = f,
    x = c(ctrl_estimates, exp_estimates[j])
  )
  return(matrix(gradient, ncol = 1))
}

# v_prime_deriv = function(time_points,
#                             ctrl_estimates,
#                             exp_estimates,
#                             interpolation,
#                             gamma_0,
#                             j = 1:length(exp_estimates),
#                             weights,
#                             analytical = FALSE){
#   K = length(time_points)
#   force(time_points); force(ctrl_estimates); force(exp_estimates)
#   force(interpolation); force(vcov);
#   force(j); force(weights); force(K); force(ref_fun); force(gamma_0)
#   v = function(theta) {
#     v_function(
#       time_points,
#       theta[1:K],
#       theta[K + j],
#       interpolation,
#       gamma_0,
#       j,
#       weights)
#   }
#
#   f = function(theta) {
#     v_gamma = function(gamma) {
#       v_function(
#         time_points,
#         theta[1:K],
#         theta[K + j],
#         interpolation,
#         gamma,
#         j,
#         weights)
#     }
#     v_deriv_gamma = numDeriv::grad(
#       func = v_gamma,
#       x = gamma_0
#     )
#     (1 / v_deriv_gamma)
#   }
#   numDeriv::grad(
#     func = f,
#     x = c(ctrl_estimates, exp_estimates[j])
#   )
# }
#
# v_gradient = function(time_points,
#                          ctrl_estimates,
#                          exp_estimates,
#                          interpolation,
#                          gamma_0,
#                          j = 1:length(exp_estimates),
#                          weights,
#                          analytical = FALSE){
#   K = length(time_points)
#   force(time_points); force(ctrl_estimates); force(exp_estimates)
#   force(interpolation); force(vcov);
#   force(j); force(weights); force(K); force(ref_fun); force(gamma_0)
#   v = function(theta) {
#     v_function(
#       time_points,
#       theta[1:K],
#       theta[K + j],
#       interpolation,
#       gamma_0,
#       j,
#       weights)
#   }
#
#   numDeriv::grad(
#     func = v,
#     x = c(ctrl_estimates, exp_estimates[j])
#   )
# }

gradient_gamma_w_analytical = function(time_points,
                                       ctrl_estimates,
                                       exp_estimates,
                                       interpolation,
                                       gamma_0,
                                       j = 1:length(exp_estimates),
                                       weights) {
  # Diagonal matrix with time_points vector as diagonal.
  D_t = diag(time_points[j + 1])
  # Construct reference function trajectory.
  ref_fun = ref_fun_constructor(x_ref = time_points,
                                y_ref = ctrl_estimates,
                                method = interpolation)
  # Derivative of the reference trajectory at the time_points multiplied with
  # gamma_0. These are essentially the "corresponding" time points for the
  # active treatment group.
  f0_gradient_t_gamma0 = f0_gradient_t(ref_fun, time_points[j + 1] * gamma_0)
  # Jacobian of the above function with respect to the vertical parameters.
  f0_gradient_t_gamma0_jacobian = f0_gradient_t_jacobian(time_points[j + 1] * gamma_0,
                                                         time_points,
                                                         ctrl_estimates,
                                                         exp_estimates,
                                                         interpolation)
  # Row vector with weights.
  weights_row = matrix(weights, nrow = 1)
  # Gradient of v-function with respect to vertical parameters.
  v_gradient = v_deriv_alpha_beta(time_points,
                                  ctrl_estimates,
                                  interpolation,
                                  gamma_0,
                                  j,
                                  weights)
  # v-function
  v = v_function(time_points,
                 ctrl_estimates,
                 exp_estimates,
                 interpolation,
                 gamma_0,
                 j,
                 weights)
  # Compute expression that is used twice in the gradient.
  denominator = as.numeric(weights_row %*% D_t %*% f0_gradient_t_gamma0)

  # Return the analytical gradient.
  return(
    (t(v_gradient) / denominator) -
      v * t((weights_row %*% D_t %*% f0_gradient_t_gamma0_jacobian) / (denominator ** 2))
  )

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
                              interpolation,
                              gamma_0,
                              j,
                              weights) {
  # Compute the jacobian of the score vector.
  J = score_vector_jacobian(time_points,
                            ctrl_estimates,
                            interpolation,
                            gamma_0,
                            j)
  v_deriv = matrix(weights, nrow = 1) %*% J
  return(v_deriv)
}

#' Gradient of the reference trajectory at fixed time point
#'
#' The [f0_gradient_t()] function is a simple wrapper that returns
#' \eqn{\frac{\partial \boldsymbol{f_0} (\boldsymbol{t^*}; \boldsymbol{\alpha}, \boldsymbol{\beta})}{\partial \boldsymbol{t^*}}}.
#'
#' @param ref_fun
#' @param t_vec
#'
#' @return
#' @export
#'
#' @examples
f0_gradient_t = function(ref_fun, t_vec) {
  return(matrix(ref_fun(t_vec, deriv = 1), ncol = 1))
}

f0_gradient_t_jacobian = function(t_vec,
                                  time_points,
                                  ctrl_estimates,
                                  exp_estimates,
                                  interpolation) {
  myenv <- new.env()
  myenv$x_ref <- time_points
  myenv$y_ref <- ctrl_estimates
  myenv$method <- interpolation
  myenv$t_vec <- t_vec
  myenv$ref_fun_constructor <- ref_fun_constructor
  myenv$f0_gradient_t <- f0_gradient_t
  J_object = stats::numericDeriv(expr = quote({
    ref_fun = ref_fun_constructor(x_ref, y_ref, method)
    f0_gradient_t(ref_fun, t_vec)
  }),
  theta = c("y_ref"),
  rho = myenv)
  # Next, the second part of the Jacobian. This corresponds to the identity
  # matrix. We join both parts to get the Jacobian.
  J = cbind(attr(J_object, "gradient"), matrix(0, nrow = length(t_vec), ncol = length(exp_estimates)))
  return(J)
}



v_function = function(time_points,
                      ctrl_estimates,
                      exp_estimates,
                      interpolation,
                      gamma_0,
                      j,
                      weights
) {
  ref_fun = ref_fun_constructor(time_points,
                                ctrl_estimates,
                                interpolation)
  # Compute the jacobian of the score vector.
  J = score_vector_jacobian(time_points,
                            ctrl_estimates,
                            interpolation,
                            gamma_0,
                            j)
  # Compute the score vector.
  g = (exp_estimates[j] - ref_fun(gamma_0 * time_points[j + 1]))
  weights = matrix(weights, ncol = 1)
  v = t(weights) %*% g
  return(as.numeric(v))
}
