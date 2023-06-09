analyze_qr = function(data_trial, tau = 0.5, R = 100) {
  rq_fit = quantreg::rq(formula = ADAScog_integer ~ arm_time + 0,
                        data = data_trial,
                        tau = tau)
  summary(rq_fit,
          se = "boot",
          cluster = data_sim$SubjId,
          covariance = TRUE,
          R = 100)
}


