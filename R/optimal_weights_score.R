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

update_weights = function() {

}


optimize_weights = function(time_points,
                            ctrl_estimates,
                            exp_estimates,
                            ref_fun,
                            interpolation,
                            vcov,
                            j = 1:length(exp_estimates),
                            weights = NULL,
                            penalty = function(x) 0,
                            ...) {
  epsilon = 0.001
  # Let the weights sum to one and leave out the last element. The last element
  # is known if the other elements are known and the vector sums to one.
  w_new = weights / sum(weights)
  w_new = w_new[-length(weights)]
  # Transform weights to logit scale. This ensure that the weights are always in
  # the unit interval.
  logit_w_new = log(w_new / ( 1 - w_new))
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
  sigma_squared = function(logit_w) {
    w = 1 / (1 + exp(-logit_w))
    # ANALYTICAL GRADIENT SHOULD STILL BE IMPLEMENTED!!!
    gr_gamma_w = gradient_gamma_w(
      time_points = time_points,
      ctrl_estimates = ctrl_estimates,
      exp_estimates = exp_estimates,
      vcov = vcov,
      interpolation = interpolation,
      gamma_0 = gamma_old,
      j = j,
      weights = c(w, 1 - sum(w))
    )
    se = matrix(gr_gamma_w, nrow = 1) %*%
      vcov %*%
      matrix(gr_gamma_w, ncol = 1)
    return(as.numeric(se))
  }
  stoppping_criterion = FALSE

  while (!(stoppping_criterion)) {
    logit_w_old = logit_w_new
    w_old = 1 / (1 + exp(-logit_w_old))
    gamma_old = gamma_new
    logit_w_new = optim(
      f = sigma_squared,
      par = logit_w_old,
      control = list(maxit = 30)
    )$par
    w_new = 1 / (1 + exp(-logit_w_new))
    gamma_new = update_gamma(
      time_points,
      ctrl_estimates,
      exp_estimates,
      ref_fun,
      interpolation,
      vcov,
      "custom",
      j,
      c(w_new, 1 - sum(w_new))
    )
    stoppping_criterion = sum((w_old - w_new) ** 2) < epsilon
  }
  return(w_new)
}

