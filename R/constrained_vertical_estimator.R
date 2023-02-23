#' Constrained vertical Estimator
#'
#' @return
#' @export
#'
#' @examples
constrained_vertical_estimator = function(alpha_obs, beta_obs, Sigma_obs) {
  # Compute the inverse covariance matrix. This leads to efficiency gains since
  # the inverse would otherwise be computed multiple time.
  Sigma_inv = solve(Sigma_obs)
  length_alpha = length(alpha_obs)
  obj_function = function(theta) {
    alpha = theta[1:length_alpha]
    beta = theta[(length_alpha + 1):length(theta)]
    mahalanobis(x = c(alpha_obs, beta_obs),
                center = c(alpha, beta),
                cov = Sigma_inv, inverted = TRUE)
  }
  # Constraint matrix for first assumption: strict monotonically increasing
  # trajectory.
  A1 = constraint_matrix_A1(length(alpha_obs), length(beta_obs))
  # Constrain matrix for second assumption: min(alpha) < min(beta).
  A2 = constraint_matrix_A2(length(alpha_obs), length(beta_obs))
  # Combine both constraint matrices.
  constraint_matrix = rbind(A1, A2)
  # Define constraint vector. Since all constraint are inequalities without
  # constants, this is just the zero vector.
  constraint_vector = matrix(0, nrow = nrow(constraint_matrix), ncol = 1)

  starting_values = c(alpha_obs, beta_obs)
  constrOptim(
    theta = starting_values,
    f = obj_function,
    grad = gradient_mahalanobis,
    ui = constraint_matrix,
    ci = constraint_vector
  )
}

#' Return constrain matrix for Assumption 1
#'
#' @param length_alpha
#' @param length_beta
#'
#' @return
#' @export
#'
#' @examples
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

#' Return Constraint matrix for Assumption 2
#'
#' @param length_alpha
#' @param length_beta
#'
#' @return
#' @export
#'
#' @examples
constraint_matrix_A2 = function(length_alpha, length_beta) {
  A2_matrix = matrix(
    c(-1, rep(0, length_alpha - 1), 1, rep(0, length_beta - 1)),
    nrow = 1
  )
  return(A2_matrix)
}

gradient_mahalanobis = function(alpha_obs, beta_obs, alpha, beta) {
  y = c(alpha_obs, beta_obs)
  alpha_beta = c(alpha, beta)
  return(-2 * y + 2 * alpha_beta)
}
