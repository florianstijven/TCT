#' Trajectory function
#'
#' This function returns the value of the trajectory at `t`.
#'
#' @param t Time at which the trajectory's value is computed.
#' @param z Group indicater:
#'  * `z = 0`: indicates the reference group
#'  * `z = 1`: indicates the treatment group
#' @param gamma Acceleration factor parameter.
#' @param alpha (Estimated) points on the reference profile.
#' @inheritParams get_new_time
#'
#' @return (numeric) value of the trajctory at `t`.
trajectory = function(t, z, gamma, alpha, x_ref, method) {
  # Reference profile
  f0 = ref_fun_constructor(x_ref, alpha, method = method)
  f0(t - z * (1 - gamma) * t)
}
