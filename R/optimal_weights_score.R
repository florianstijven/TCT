update_gamma = function(time_points,
                        ctrl_estimates,
                        exp_estimates,
                        ref_fun,
                        interpolation,
                        vcov,
                        type,
                        j,
                        weights) {
  score_estimate_common(
    time_points,
    ctrl_estimates,
    exp_estimates,
    ref_fun,
    interpolation,
    vcov,
    type,
    j,
    weights
  )
}


#' Find the "optimal" weights for the score estimator
#'
#' The [optimize_weights()] function finds the weights for the score estimator
#' that minimize the estimated variance of this estimator. The estimated
#' variance is based in the delta method.
#'
#'
#' @inheritParams score_test_common
#' @param weights Staring values for the weights vector.
#' @param epsilon Tolerance to stop optimization algorithm. If the l2-norm of
#'   the difference in weights from two subsequent iterations is smaller than
#'   `epsilon`, the optimization algorithm stops.
#'
#' @return (numeric) vector with optimal weights.
optimize_weights = function(time_points,
                            ctrl_estimates,
                            exp_estimates,
                            ref_fun,
                            interpolation,
                            vcov,
                            j = 1:length(exp_estimates),
                            weights = rep(1, length(exp_estimates)),
                            epsilon = 1e-6) {
  # The weights multiplied by a constant are equivalent. So, we let the weights
  # sum to one.
  w_new = weights / sum(weights)
  # Because the weights sum to one, the last weight parameter is redundant.
  w_new = w_new[-length(weights)]
  # Log-transform the weight "odds" with the last category is reference. This
  # ensure that the weights are always positive and sum to one. In principle,
  # negative weights still lead to valid estimators and inferences. However,
  # they make no sense in practice. We therefore do not allow for negative
  # weights in this procedure.
  log_odds_w_new = log(w_new / (1 - sum(w_new)))
  gamma_new = update_gamma(
    time_points,
    ctrl_estimates,
    exp_estimates,
    ref_fun,
    interpolation,
    vcov,
    "custom",
    j,
    weights
  )
  sigma_squared = function(log_odds_w) {
    w_K = 1 / (1 + sum(exp(log_odds_w)))
    w = c(w_K * exp(log_odds_w),
          w_K)
    gr_gamma_w = gradient_gamma_w_analytical(
      time_points = time_points,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      interpolation = interpolation,
      gamma_0 = gamma_old,
      j = j,
      weights = w
    )
    se = t(gr_gamma_w) %*%
      vcov %*%
      gr_gamma_w
    return(as.numeric(se))
  }
  stopping_criterion = FALSE

  gamma_old = gamma_new
  optim_value_new = sigma_squared(log_odds_w_new)
  while (!(stopping_criterion)) {
    log_odds_w_old = log_odds_w_new
    optim_value_old = optim_value_new
    w_K_old = 1 / (1 + sum(exp(log_odds_w_old)))
    w_old = w_K_old * exp(log_odds_w_old)

    optim_object = stats::optim(
      f = sigma_squared,
      par = log_odds_w_old,
      control = list(maxit = 30)
    )
    log_odds_w_new = optim_object$par
    optim_value_new = optim_object$value
    w_K_new = 1 / (1 + sum(exp(log_odds_w_new)))
    w_new = w_K_new * exp(log_odds_w_new)
    gamma_new = update_gamma(
      time_points,
      ctrl_estimates,
      exp_estimates,
      ref_fun,
      interpolation,
      vcov,
      "custom",
      j,
      c(w_new, w_K_new)
    )
    stopping_criterion = (sum((w_old - w_new) ** 2) < epsilon) |
      ((optim_value_old - optim_value_new) < epsilon)
  }
  return(c(w_new, w_K_new))
}

