data = simulated_test_trial %>%
  dplyr::mutate(time_int = (Week %/% 25)) %>%
  dplyr::arrange(trial_number, SubjId, time_int) %>%
  dplyr::mutate(time_int = as.integer(time_int) + 1L) %>%
  dplyr::mutate(arm_time = ifelse(time_int == 1L,
                                  "baseline",
                                  paste0(arm, ":", time_int)))
library(tidyverse)
data %>%
  group_by(Week, arm) %>%
  summarize(mean = mean(ADAScog_integer)) %>%
  ggplot(aes(x = Week, y = mean, color = as.factor(arm))) +
  geom_line()

mmrm_fit = analyze_mmrm(data)

constrained_vertical_estimator(
  coef(mmrm_fit)[c(9, 1:4)],
  coef(mmrm_fit)[5:8],
  vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)]
)
TCT_Fit = TCT(
  time_points = 0:4,
  ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)] + change,
  exp_estimates = coef(mmrm_fit)[5:8],
  vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
  interpolation = "spline",
  B = 5e3,
  constraints = FALSE
)
summary.TCT(TCT_Fit)

change = c(0, -1, 0, 0, 0)
TCT_Fit_c = TCT(
  time_points = 0:4,
  ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)] + change,
  exp_estimates = coef(mmrm_fit)[5:8],
  vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
  interpolation = "spline",
  B = 5e3,
  constraints = TRUE
)
summary.TCT(TCT_Fit_c)

set.seed(1)


TCT_Fit = TCT(
  time_points = 0:4,
  ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
  exp_estimates = coef(mmrm_fit)[5:8],
  vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
  interpolation = "fourPL",
  B = 0
)
summary(TCT_Fit)

library(quantreg)
library(dplyr)
# control_medians = data %>%
#   filter(arm == 0) %>%
#   group_by(time_int) %>%
#   summarize(median = median(ADAScog_integer))
trajectory_5 = function(t,
                        z,
                        gamma_log,
                        alpha1,
                        alpha2,
                        alpha3,
                        alpha4,
                        alpha5,
                        x_ref,
                        method) {
  alpha = c(alpha1, alpha2, alpha3, alpha4, alpha5)
  trajectory(t, z, exp(gamma_log), alpha, x_ref, method)
}
# nlrq_fit = nlrq(
#   formula = ADAScog_integer ~ trajectory_5(Week,
#                                            arm,
#                                            1 - gamma,
#                                            alpha1,
#                                            alpha2,
#                                            alpha3,
#                                            alpha4,
#                                            alpha5,
#                                            x_ref = c(0, 25, 50, 75, 100),
#                                            "spline"),
#   data = data,
#   start = list(gamma = 0.75,
#                alpha1 = control_medians$median[1],
#                alpha2 = control_medians$median[2],
#                alpha3 = control_medians$median[3],
#                alpha4 = control_medians$median[4],
#                alpha5 = control_medians$median[5]),
#   tau = 0.3
# )
#
# summary(nlrq_fit)
#
# nls_fit = nls(
#   formula = ADAScog_integer ~ trajectory_5(Week,
#                                            arm,
#                                            1 - gamma,
#                                            alpha1,
#                                            alpha2,
#                                            alpha3,
#                                            alpha4,
#                                            alpha5,
#                                            x_ref = c(0, 25, 50, 75, 100),
#                                            "spline"),
#   data = data,
#   start = list(gamma = 1,
#                alpha1 = control_medians$median[1],
#                alpha2 = control_medians$median[2],
#                alpha3 = control_medians$median[3],
#                alpha4 = control_medians$median[4],
#                alpha5 = control_medians$median[5]),
# )
# summary(nls_fit)


library(quantreg)
library(tidyverse)
# Simulate "easy" data
set.seed(1)
vcov_ref = matrix(c(45.1, 40.0, 45.1, 54.9, 53.6,
                    40.0, 57.8, 54.4, 66.3, 64.1,
                    45.1, 54.4, 72.0, 80.0, 77.6,
                    54.9, 66.3, 80.0, 109.8, 99.3,
                    53.6, 64.1, 77.6, 99.3, 111.4), nrow = 5, byrow = TRUE)

time_means_list = list(
  short = c(19.6, 20.5, 20.9, 22.7, 23.8),
  long = c(18, 19.7, 20.9, 22.7, 24.7)
)
n_sim = 2e3
duration = c("long")
gamma_slowing = c(1)
n = c(2e2)
N_trials = 1
gamma_nlrq = c()
gamma_tct = c()
gamma_tct_c = c()
gamma_nls = c()
gamma_nlrq_0.25 = c()
gamma_nlrq_0.75 = c()
gamma_nlrq_se <- gamma_nlrq_se_0.25 <- gamma_nlrq_se_0.75 <- rep(1e6, n_sim)

for (i in 1:n_sim) {
  time_means = time_means_list[[1]]
  time_means_trt = spline(x = 0:4, y = time_means, xout = gamma_slowing * 0:4)$y

  data_trial = rbind(
    mvtnorm::rmvnorm(n = (n / 2),
                     mean = time_means,
                     sigma = vcov_ref),
    mvtnorm::rmvnorm(n = (n / 2),
                     mean = time_means_trt,
                     sigma = vcov_ref)
  )
  colnames(data_trial) = c("w0", "w25", "w50", "w75", "w100")
  data_trial = data.frame(data_trial)
  data_trial$arm = rep(0:1, each = (n / 2) * N_trials)
  data_trial$trial_number = rep(1:N_trials, times = n)
  data_trial$SubjId = 1:(n * N_trials)
  data_trial = data_trial %>%
    pivot_longer(cols = c("w0", "w25", "w50", "w75", "w100"),
                 values_to = "ADAScog_integer")
  data_trial$time_int = rep(1:5, N_trials * n)
  data_trial = data_trial %>%
    mutate(arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)),
           Week = (time_int - 1) * 25)

  # quantiles = data_trial %>%
  #   filter(arm == 0) %>%
  #   group_by(Week) %>%
  #   summarise(quantile = median(ADAScog_integer))
  #
  # nlrq_fit = nlrq(
  #   formula = ADAScog_integer ~ trajectory_5(Week,
  #                                            arm,
  #                                            gamma,
  #                                            alpha1,
  #                                            alpha2,
  #                                            alpha3,
  #                                            alpha4,
  #                                            alpha5,
  #                                            x_ref = c(0, 25, 50, 75, 100),
  #                                            "spline"),
  #   data = data_trial,
  #   start = list(gamma = log(gamma_slowing),
  #                alpha1 = quantiles$quantile[1],
  #                alpha2 = quantiles$quantile[2],
  #                alpha3 = quantiles$quantile[3],
  #                alpha4 = quantiles$quantile[4],
  #                alpha5 = quantiles$quantile[5]),
  #   tau = 0.5
  # )
  #
  # gamma_nlrq[i] = coef(nlrq_fit)[1]
  # gamma_nlrq_se[i] = try(sqrt(summary(nlrq_fit)$cov[1, 1]))
  #
  # quantiles = data_trial %>%
  #   filter(arm == 0) %>%
  #   group_by(Week) %>%
  #   summarise(quantile = quantile(ADAScog_integer, 0.25))
  #
  # nlrq_fit_0.25 = nlrq(
  #   formula = ADAScog_integer ~ trajectory_5(Week,
  #                                            arm,
  #                                            gamma,
  #                                            alpha1,
  #                                            alpha2,
  #                                            alpha3,
  #                                            alpha4,
  #                                            alpha5,
  #                                            x_ref = c(0, 25, 50, 75, 100),
  #                                            "spline"),
  #   data = data_trial,
  #   start = list(gamma = log(gamma_slowing),
  #                alpha1 = quantiles$quantile[1],
  #                alpha2 = quantiles$quantile[2],
  #                alpha3 = quantiles$quantile[3],
  #                alpha4 = quantiles$quantile[4],
  #                alpha5 = quantiles$quantile[5]),
  #   tau = 0.25
  # )
  #
  # gamma_nlrq_0.25[i] = coef(nlrq_fit_0.25)[1]
  # gamma_nlrq_se_0.25[i] = try(sqrt(summary(nlrq_fit_0.25)$cov[1, 1]))
  #
  # quantiles = data_trial %>%
  #   filter(arm == 0) %>%
  #   group_by(Week) %>%
  #   summarise(quantile = quantile(ADAScog_integer, 0.75))
  #
  # nlrq_fit_0.75 = nlrq(
  #   formula = ADAScog_integer ~ trajectory_5(Week,
  #                                            arm,
  #                                            gamma,
  #                                            alpha1,
  #                                            alpha2,
  #                                            alpha3,
  #                                            alpha4,
  #                                            alpha5,
  #                                            x_ref = c(0, 25, 50, 75, 100),
  #                                            "spline"),
  #   data = data_trial,
  #   start = list(gamma = log(gamma_slowing),
  #                alpha1 = time_means[1],
  #                alpha2 = time_means[2],
  #                alpha3 = time_means[3],
  #                alpha4 = time_means[4],
  #                alpha5 = time_means[5]),
  #   tau = 0.75
  # )
  #
  # gamma_nlrq_0.75[i] = coef(nlrq_fit_0.75)[1]
  # gamma_nlrq_se_0.75[i] = try(sqrt(summary(nlrq_fit_0.75)$cov[1, 1]))
  #
  # gamma_nls[i] = NA
  # try(
  #   {
  #     nls_fit = nls(
  #       formula = ADAScog_integer ~ trajectory_5(Week,
  #                                                arm,
  #                                                gamma,
  #                                                alpha1,
  #                                                alpha2,
  #                                                alpha3,
  #                                                alpha4,
  #                                                alpha5,
  #                                                x_ref = c(0, 25, 50, 75, 100),
  #                                                "spline"),
  #       data = data_trial,
  #       start = list(gamma = log(gamma_slowing),
  #                    alpha1 = quantiles$quantile[1],
  #                    alpha2 = quantiles$quantile[2],
  #                    alpha3 = quantiles$quantile[3],
  #                    alpha4 = quantiles$quantile[4],
  #                    alpha5 = quantiles$quantile[5]),
  #     )
  #
  #     gamma_nls[i] = coef(nls_fit)[1]
  #   }
  # )
  #


  mmrm_fit = analyze_mmrm(data_trial)
  # set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0,
    constraints = FALSE
  )
  # summary.TCT(TCT_Fit)
  TCT_common_fit = TCT_common(TCT_Fit, B = 0, constraints = FALSE)
  # median(TCT_common_fit$bootstrap_estimates[[1]])
  # mean(TCT_common_fit$bootstrap_estimates[[1]])
  # hist(TCT_common_fit$bootstrap_estimates[[1]])
  # summary.TCT_common(TCT_common_fit)
  gamma_tct[i] = coef(TCT_common_fit)[1]

  # set.seed(1)
  TCT_Fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0,
    constraints = TRUE
  )
  # summary.TCT(TCT_Fit)
  TCT_common_fit = TCT_common(TCT_Fit, B = 0, constraints = TRUE)
  # median(TCT_common_fit$bootstrap_estimates[[1]])
  # mean(TCT_common_fit$bootstrap_estimates[[1]])
  # hist(TCT_common_fit$bootstrap_estimates[[1]])
  # summary.TCT_common(TCT_common_fit)
  gamma_tct_c[i] = coef(TCT_common_fit)[1]
}

gamma_nlrq_se_0.25 = as.numeric(gamma_nlrq_se_0.25)
gamma_nlrq_se_0.25[is.na(gamma_nlrq_se_0.25)] = 1e10
results_df = data.frame(gamma_tct = gamma_tct,
                        gamma_nlrq,
                        gamma_nls,
                        gamma_nlrq_0.25,
                        gamma_nlrq_0.75) %>%
  mutate(
    mean_gamma_nlrq = (gamma_nlrq + gamma_nlrq_0.25 + gamma_nlrq_0.75) / 3,
    wmean_gamma_nlrq =
      (gamma_nlrq_se^-1 * gamma_nlrq + gamma_nlrq_se_0.25^-1 * gamma_nlrq_0.25 + gamma_nlrq_se_0.75^-1 * gamma_nlrq_0.75) /
      (gamma_nlrq_se^-1 + gamma_nlrq_se_0.25^-1 + gamma_nlrq_se_0.75^-1)
  ) %>%
  pivot_longer(cols = 1:7, values_to = "estimate", names_to = "type") %>%
  mutate(estimate = ifelse(type == "gamma_tct",
                           estimate,
                           exp(estimate)))

results_df %>%
ggplot(aes(x = estimate, fill = type, color = type)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  xlim(-2, 3)

results_df %>%
  ggplot(aes(x = estimate, color = type)) +
  geom_density() +
  xlim(-2, 3)

pairs({data.frame(gamma_tct, gamma_nlrq, gamma_nls, gamma_nlrq_0.25, gamma_nlrq_0.75) %>%
        mutate(mean_gamma_nlrq = (gamma_nlrq + gamma_nlrq_0.25 + gamma_nlrq_0.75) / 3)})

results_df %>%
  group_by(type) %>%
  summarise(mean = mean(estimate, na.rm = TRUE),
            MSE = sqrt(mean((estimate - gamma_slowing)**2, na.rm = TRUE)),
            mean_log = mean(estimate, na.rm = TRUE),
            MSE_log = sqrt(mean(abs((estimate - gamma_slowing)))))
summary({results_df %>%
          mutate(id = rep(1:n_sim, each = 7)) %>%
          pivot_wider(id_cols = "id", names_from = "type", values_from = "estimate")})
results_df_wide = results_df %>%
  mutate(id = rep(1:n_sim, each = 7)) %>%
  pivot_wider(id_cols = "id", names_from = "type", values_from = "estimate")
