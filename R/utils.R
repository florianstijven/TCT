#' Maps observed outcomes to a reference profile using interpolation
#'
#' The [get_new_time()] function maps observed (mean) outcomes to the time-axis
#' of a reference profile. The reference profile is defined by a set of
#' time-outcome pairs in the reference group. These pairs are usually the
#' estimated mean outcomes at fixed time-points in the control group.
#'
#' @details Let \eqn{f_{0}(t)} denote the reference profile. This profiles is
#'   estimated as the interpolation curve between the points specified in
#'   `x_ref` and `y_ref`. The interpolating curve is constructed with natural
#'   cubic splines or a piecewise linear function.
#'
#'   Let \eqn{y^{\ast}} denote the observed (mean) outcome. The time-mapping is
#'   then defined as \eqn{t^{\ast} = f^{-1}_{0}(y^{\ast})}. This is equivalent
#'   to solving \eqn{f_{0}(t^{\ast}) = y^{\ast}} for \eqn{t^{\ast}}. The latter
#'   approach is implemented numerically with [stats::uniroot()]. Although
#'   \eqn{f^{-1}_{0}(\cdot)} can be obtained analytically in some cases, the
#'   numerical method is fully general.
#'
#' @param y_ref Vector of reference (mean) outcome values.
#' @param x_ref Vector of reference time values.
#' @param y_obs Vector of (mean) outcome values to be mapped to the reference
#'   profile.
#' @param method Method for interpolation:
#' * `"linear"`: piecewise linear interpolation
#' * `"spline"`: natural cubic spline interpolation
#'
#' @return A vector of the time-mapped values.
#'
get_new_time = function(y_ref, x_ref, y_obs, method = "linear") {
  # amount of extrapolation
  extrapol = 1e4
  # Only the first first value in y_obs is used. This enables us to use
  # straightforward data handling methods from dplyr.
  y_obs = y_obs
  if (method == "linear") {
    ref_fun = function(x,
                       .x_ref = x_ref,
                       .y_ref = y_ref) {
      if ((x < min(.x_ref))) {
        # standard function does not do extrapolation, this extrapolation is
        # implemented here.
        y2 = .y_ref[2]; y1 = .y_ref[1]
        x2 = .x_ref[2]; x1 = .x_ref[1]
        # slope
        a = (y2 - y1) / (x2 - x1)
        # intercept
        b = y2 - a * x2
        y_interpolated = a * x + b
      }
      else if (x > max(.x_ref)) {
        p = length(.x_ref)
        y2 = .y_ref[p]; y1 = .y_ref[p - 1]
        x2 = .x_ref[p]; x1 = .x_ref[p - 1]
        # slope
        a = (y2 - y1) / (x2 - x1)
        # intercept
        b = y2 - a * x2
        y_interpolated = a * x + b
      }
      else {
        y_interpolated = approx(.x_ref,
                         .y_ref,
                         method = "linear",
                         rule = 1,
                         xout = x)$y
      }
      return(y_interpolated)
    }
  }
  else if (method == "spline") {
    ref_fun = splinefun(x_ref,
                        y_ref,
                        method = "natural")
  }
  # x-value that maps to y_obs with the reference function. The reference
  # function represents the reference mean trajectory. To find this x-value, we
  # solve the following equation for x_mapped: f_ref (x_mapped) = y_obs. This is
  # found by finding the root of g(x) = f_ref (x) - y_obs.

  # If y_obs is larger than the reference function everywhere, or smaller than
  # the reference function everywhere. Then the root of the g(x) above does not
  # exist within the time range considered here. If that is the case, NA is
  # returned.
  x_mapped = rep(NA, length(y_obs))
  for (i in seq_along(x_mapped)) {
    if (sign(ref_fun(min(x_ref) - extrapol) - y_obs[i]) ==
        sign(ref_fun(max(x_ref) + extrapol) - y_obs[i])) {
      x_mapped[i] = NA
    }
    else {
      # The uniroot() function finds the root of g(x).
      x_mapped[i] = uniroot(
        f = function(x)
          ref_fun(x) - y_obs[i],
        interval = c(min(x_ref) - extrapol,
                     max(x_ref) + extrapol),
        tol = .Machine$double.eps^0.5,
        maxiter = 10e4
      )$root
    }
  }

  return(x_mapped)
}
