#' Bootstrap confidence intervals
#'
#' The [bootstrap_CI()] function computes bootstrap confidence intervals given
#' the bootstrap replicates.
#'
#' @param boot_replicates A numeric vector with bootstrap replicates of the
#' estimate of interest.
#' @param estimate The point estimate of the target parameter.
#' @param alpha Significance level for the confidence interval. Default is `0.05`.
#' @param type Type of bootstrap confidence interval to compute. Currently only
#' two types are implemented: `"BC percentile"` and `"percentile"`.
#'
#' @returns (numeric) vector with lower and upper limits
bootstrap_CI = function(boot_replicates, estimate, alpha = 0.05, type = "BCa") {
  # Compute the required type of bootstrap CI.
  if (type == "BC percentile") {
    # Compute BCa interval
    return(
      BC_percentile_CI(
        boot_replicates = boot_replicates,
        estimate = estimate,
        alpha = alpha
      )
    )
  } else if (type == "percentile") {
    # Compute percentile interval
    return(percentile_CI(boot_replicates, alpha))
  } else {
    stop("Invalid type. Must be 'BCa' or 'percentile'.")
  }
}

#' Percentile confidence intervals
#'
#' The [percentile_CI()] function computes the percentile bootstrap confidence
#' interval given the bootstrap replicates.
#'
#' @inherit bootstrap_CI
percentile_CI = function(boot_replicates, alpha = 0.05) {
  # Check for missing values or NaN in the bootstrap replicates. A warning is
  # raised if there are any missing values. The warning also details the number
  # of problematic values such that the user can decide whether to ignore this.
  if (any(is.na(boot_replicates)) | any(is.nan(boot_replicates))) {
    warning_message = paste(
      "Some bootstrap replicates have missing values. These are removed for computing percentile confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates)))
    warning(warning_message)
    # Remove missing values from the bootstrap replicates and standard errors.
    boot_replicates = boot_replicates[!is.na(boot_replicates) & !is.nan(boot_replicates)]
  }
  ci_lower = stats::quantile(x = boot_replicates, probs = alpha / 2, na.rm = TRUE)
  ci_upper = stats::quantile(x = boot_replicates, probs = 1 - (alpha / 2), na.rm = TRUE)

  return(c(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

#' Bias-corrected percentile confidence intervals
#'
#' The [BC_percentile_CI()] function computes the bias-corrected percentile
#' bootstrap confidence interval given the bootstrap replicates.
#'
#' @inherit bootstrap_CI
BC_percentile_CI = function(estimate, boot_replicates, alpha = 0.05) {
  # Check for missing values or NaN in the bootstrap replicates. A warning is
  # raised if there are any missing values. The warning also details the number
  # of problematic values such that the user can decide whether to ignore this.
  if (any(is.na(boot_replicates)) | any(is.nan(boot_replicates))) {
    warning_message = paste(
      "Some bootstrap replicates have missing values. These are removed for computing percentile confidence intervals.",
      "\nNumber of missing values in bootstrap replicates: ",
      sum(is.na(boot_replicates) | is.nan(boot_replicates)))
    warning(warning_message)
    # Remove missing values from the bootstrap replicates and standard errors.
    boot_replicates = boot_replicates[!is.na(boot_replicates) & !is.nan(boot_replicates)]
  }
  # Compute bias-correction value.
  p0 = mean(boot_replicates < estimate) + 0.5 * mean(boot_replicates == estimate)
  z0 = stats::qnorm(p0)

  alpha_lower = stats::pnorm(2 * z0 + stats::qnorm(alpha / 2))
  alpha_upper = stats::pnorm(2 * z0 + stats::qnorm(1 - alpha / 2))

  ci_lower = stats::quantile(x = boot_replicates, probs = alpha_lower, na.rm = TRUE)
  ci_upper = stats::quantile(x = boot_replicates, probs = alpha_upper, na.rm = TRUE)

  return(c(ci_lower = ci_lower, ci_upper = ci_upper))
}
