#' Maps observed outcomes to a reference profile using interpolation
#'
#' The [get_new_time()] function maps observed (mean) outcomes to the time-axis
#' of a reference profile. The reference profile is defined by a set of
#' time-outcome pairs in the reference group. These pairs are usually the
#' estimated mean outcomes at fixed time-points in the control group.
#'
#' @details
#' # Reference Profile
#'  Let \eqn{f_{0}(t)} denote the reference profile. This
#' profile is estimated as the interpolating curve between the points specified
#' in `x_ref` and `y_ref`. The interpolation method is specified in the `method`
#' argument.
#'
#' # Time Mapping
#'  Let \eqn{y^{\ast}} denote the observed (mean) outcome. The time-mapping is
#'  then defined as \eqn{t^{\ast} = f^{-1}_{0}(y^{\ast})}. This is equivalent to
#'  solving \eqn{f_{0}(t^{\ast}) = y^{\ast}} for \eqn{t^{\ast}}. The latter
#'  approach is implemented numerically with [stats::uniroot()] for all
#'  `method`s, expect for `"fourPL"`. For the latter method,
#'  \eqn{f^{-1}_{0}(\cdot)} is obtained analytically.
#'
#' @param y_ref Vector of reference (mean) outcome values.
#' @param x_ref Vector of reference time values.
#' @param y_obs Vector of (mean) outcome values to be mapped to the reference
#'   profile.
#' @param method Method for interpolation:
#' * `"linear"`: piecewise linear interpolation
#' * `"spline"`: natural cubic spline interpolation
#' * `"monoH.FC`: monotone Hermite spline according to the method of Fritsch and
#'      Carlson.
#' * `"fourPL`: four parameter logistic regression. This can only be used when
#'      the reference profile contains four or more means, i.e.,
#'      `length(y_ref) >= 4`. When there are more than four means, this is
#'      strictly speaking no longer interpolation. The corresponding parameters
#'      are found by minimizing the sum of squared deviations from the four PL
#'      curve.
#'
#' @return A vector of the time-mapped values.
get_new_time = function(y_ref, x_ref, y_obs, method = "linear") {
  # For fourPL, analytic solution exists.
  if (method == "fourPL") {
    theta = fit_4PL(x_ref, y_ref)
    inv_extraplation = inv_linear(x_ref, y_ref)
    x_mapped = ifelse(
      test = (y_obs < y_ref[1]) | (y_obs > y_ref[length(y_ref)]),
      yes = inv_extraplation(y_obs),
      no = inv_4PL(y_obs, theta)
    )
    return(x_mapped)
  }

  # amount of extrapolation
  extrapol = 1e5
  # Only the first first value in y_obs is used. This enables us to use
  # straightforward data handling methods from dplyr.
  p = length(x_ref)
  ref_fun = ref_fun_constructor(x_ref, y_ref, method)

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
      # The uniroot() function finds the root of g(x). If error is returned by
      # uniroot, then the NA value remains at position i. Error are most likely
      # due to mulitple possible time mappings.
      try({
        x_mapped[i] = stats::uniroot(
          f = function(x)
            ref_fun(x) - y_obs[i],
          interval = c(min(x_ref) - extrapol,
                       max(x_ref) + extrapol),
          tol = .Machine$double.eps ^ 0.5,
          maxiter = 1e3
        )$root
      }

      )
    }
  }
  return(x_mapped)
}

#' Function factory for the extrapolating part of the reference profile
#'
#' This functions returns a function that extrapolates beyond the range of the
#' observed time points in `.x_ref`. This function is the linear function
#' through the first and last point of the reference profile. This ensures some
#' stability in the extrapolation which is not present when some interpolation
#' methods are used for extrapolation.
#'
#' @inheritParams get_new_time
#'
#' @return A (linear) R function.
extrapol_fun_factory = function(x_ref, y_ref) {
  # standard function does not do extrapolation, this extrapolation is
  # implemented here.
  p = length(y_ref)
  y2 = y_ref[p]
  y1 = y_ref[1]
  x2 = x_ref[p]
  x1 = x_ref[1]
  # slope
  a = (y2 - y1) / (x2 - x1)
  # intercept
  b = y2 - a * x2
  return(
    function(x, deriv = 0) {
      a * (x ** (1 - deriv))  + b * (1 - deriv)
    }
  )
}

#' Function factory for the interpolating part of the reference profile
#'
#' This function returns a function that interpolates within the range of the
#' observed time points in `x_ref`. This is done using the interpolation method
#' in the `method` argument.
#'
#'
#' @return An interpolating R function.
#'
#' @inheritParams get_new_time
interpol_fun_factory = function(x_ref, y_ref, method) {
  switch(
    method,
    linear = linear_interpolation_f_factory(x_ref,
                                            y_ref),
    spline = stats::splinefun(x_ref,
                              y_ref,
                              method = "natural"),
    monoH.FC = stats::splinefun(x_ref,
                                y_ref,
                                method = "monoH.FC"),
    fourPL = four_pl_f_factory(x_ref,
                               y_ref)# Methods can be added here
  )
}

#' Function factory for linear interpolation function
#'
#' @inheritParams get_new_time
#'
#' @return A linear interpolation R function.
linear_interpolation_f_factory = function(x_ref, y_ref) {
  linear_spline_fun = stats::approxfun(x_ref,
                                       y_ref,
                                       method = "linear",
                                       rule = 1)
  linear_spline_deriv_fun = function(x) {
    epsilon = 1e-6
    (linear_spline_fun(x + epsilon) - linear_spline_fun(x)) / epsilon
  }
  function(x, deriv = 0) {
    if (deriv == 0) {
      linear_spline_fun(x)
    }
    else {
     linear_spline_deriv_fun(x)
    }
  }
}

#' Function factory for 4 PL function
#'
#' @inheritParams get_new_time
#'
#' @return A 4 PL function.
four_pl_f_factory = function(x_ref, y_ref) {
  # fit 4PL model
  theta = fit_4PL(x_ref, y_ref)
  four_pl_fun = function(x) {
    predict_4PL(x, theta)
  }
  # compute derivative of 4PL model
  four_pl_deriv_fun = function(x) {
    deriv_4PL(x, theta)
  }
  # combine both into a single function
  function(x, deriv = 0) {
    ifelse(
      deriv == 0,
      four_pl_fun(x),
      four_pl_deriv_fun(x)
    )
  }
}

#' Fit 4 PL function to observed reference profile.
#'
#' The [fit_4PL()] function fits the four parameter logistic function to the
#' data points in `x_ref` and `y_ref`. This is done by minimizing the sum of
#' deviations between the 4PL function and the values in `y_ref`.
#'
#' @inheritParams get_new_time
#'
#' @return A vector containing the four estimates for the 4PL function
#'   parameters.
fit_4PL = function(x_ref, y_ref) {
  # Fit the 4PL function by minimizing the squared loss.

  # Data-based starting values
  p = length(y_ref)
  a = (y_ref[p] - y_ref[1]) / (x_ref[p] - x_ref[1])
  inits = c(min(y_ref), mean(x_ref), (4 * a) / (y_ref[1] - y_ref[p]), max(y_ref))
  optim_return = stats::optim(
    par = inits,
    fn = function(par) {
      sum((predict_4PL(x_ref, par) - y_ref) ** 2)
    },
    method = "BFGS",
    control = list(reltol = 1e-5),
    gr = function(par) {
      grad_4PL(x_ref, y_ref, par)
    }
  )
  # Return estimated 4PL model parameters
  return(optim_return$par)
}


grad_4PL = function(x_ref, y_ref, theta) {
  # save 2 * (f(x_i|theta) - y_i)
  A = 2 * (predict_4PL(x_ref, theta) - y_ref)
  # save {1 - (exp(x_i)/theta_2)^theta_3}^-1
  B = (1 + (exp(x_ref) / theta[2])**theta[3])**(-1)

  grad1 = A %*% (1 - B)
  grad2 = sum(A *(theta[4] - theta[1]) * (B**2) * exp(x_ref * theta[3]) * theta[3] *
    theta[2]**(-theta[3] - 1))
  grad3 = sum(A * (theta[1] - theta[4]) * (B**2) * (exp(x_ref) / theta[2])**theta[3] *
                (x_ref - log(theta[2])))
  grad4 = sum(A * B)
  return(c(grad1, grad2, grad3, grad4))
}

predict_4PL = function(x, theta) {
  # The 4PL function is evaluated in x, given the four parameters in theta.
  z = exp(x)
  return(
    theta[1] + (theta[4] - theta[1]) / (1 + (z / theta[2])**theta[3])
  )
}

inv_4PL = function(y, theta) {
  (1 / theta[3]) * log(((theta[4] - theta[1]) / (y - theta[1])) - 1) + log(theta[2])
}

inv_linear = function(x_ref, y_ref) {
  # standard function does not do extrapolation, this extrapolation is
  # implemented here. In addition, we compute the inverse function here.
  p = length(y_ref)
  y2 = y_ref[p]
  y1 = y_ref[1]
  x2 = x_ref[p]
  x1 = x_ref[1]
  # slope
  a = (y2 - y1) / (x2 - x1)
  # intercept
  b = y2 - a * x2
  return(
    function(y) {
      (y - b) / a
    }
  )
}

deriv_4PL = function(x, theta) {
  # Return the derivative of the 4PL function parameterized by theta at the
  # points in x.
  z = exp(x)
  # Save intermediate computation in derivative formula.
  a = (z / theta[2])**theta[3]
  return(
    ((theta[1] - theta[4]) / (1 + a)**2) * theta[3] * a
  )
}


ref_fun_factory = function(.x_ref, extrapol_fun, interpol_fun) {
  ref_fun = function(x, deriv = 0) {
    ifelse(test = (x < min(.x_ref)) | (x > max(.x_ref)),
           yes = extrapol_fun(x, deriv),
           no = interpol_fun(x, deriv))
  }
  return(ref_fun)
}

deriv_f0_t = function(ref_fun, t_m) {
  ref_fun(t_m, 1)
}

deriv_tm_beta = function(t_m, beta, ref_fun) {
  ((ref_fun(t_m, 1)) ** -1) * 1
}

deriv_gamma_beta = function(t_m, t_j, ref_fun) {
  diag_J = ((ref_fun(t_m, 1)) ** -1) / t_j
  return(diag(x = diag_J))
}

deriv_gamma_alpha = function(t_m, t_j, x_ref, y_ref, method = "spline") {
  ref_fun = ref_fun_constructor(x_ref, y_ref, method)
  -1 *diag(1 / (deriv_f0_t(ref_fun, t_m) * t_j)) %*%
    attr(deriv_f0_alpha(t_m, x_ref, y_ref, method), "gradient")
}

deriv_f0_alpha = function(t_m, x_ref, y_ref, method = "spline") {
  myenv <- new.env()
  myenv$x_ref <- x_ref
  myenv$y_ref <- y_ref
  myenv$method <- method
  myenv$t_m <- t_m
  myenv$ref_fun_constructor <- ref_fun_constructor
  stats::numericDeriv(expr = quote({
    ref_fun = ref_fun_constructor(x_ref, y_ref, method)
    ref_fun(t_m)
  }),
  theta = c("y_ref"),
  rho = myenv)
}

# deriv_f0_alpha_bis = function(t_m, x_ref, y_ref, finite_diff = 1e-6, method = "spline") {
#   numDeriv::jacobian(
#     func = function(y_ref) {
#       ref_fun = ref_fun_constructor(x_ref, y_ref, method)
#       ref_fun(t_m)
#     },
#     x = y_ref
#   )
# }

ref_fun_constructor = function(x_ref, y_ref, method) {
  # if (method == "linear") {
  #   interpol_fun = stats::approxfun(x_ref,
  #                                   y_ref,
  #                                   method = "linear",
  #                                   rule = 1)
  # }
  # else if (method == "spline") {
  #   interpol_fun = stats::splinefun(x_ref,
  #                                   y_ref,
  #                                   method = "natural")
  # }
  # else if (method == "monoH.FC") {
  #   interpol_fun = stats::splinefun(x_ref,
  #                                   y_ref,
  #                                   method = "monoH.FC")
  # }
  interpol_fun = interpol_fun_factory(x_ref, y_ref, method)
  extrapol_fun = extrapol_fun_factory(x_ref, y_ref)
  ref_fun = ref_fun_factory(x_ref, extrapol_fun, interpol_fun)
  return(ref_fun)
}

jacobian_tct = function(t_m, t_j, x_ref, y_ref, ref_fun, method = "spline") {
 return(
   cbind(
     deriv_gamma_alpha(t_m, t_j, x_ref, y_ref, method),
     deriv_gamma_beta(t_m, t_j, ref_fun)
   )
 )
}


# copy from linearHypothesis() from car package. That package has many
# dependencies that can cause issue with installation. It is therefore better
# not to depend on that package.
linearHypothesis.default = function (model,
                                     hypothesis.matrix,
                                     rhs = NULL,
                                     test = c("Chisq",
                                              "F"),
                                     vcov. = NULL,
                                     singular.ok = FALSE,
                                     verbose = FALSE,
                                     coef. = coef(model),
                                     suppress.vcov.msg = FALSE,
                                     error.df,
                                     ...)

{
  if (missing(error.df)) {
    df <- stats::df.residual(model)
    test <- match.arg(test)
    if (test == "F" && (is.null(df) || is.na(df))) {
      test <- "Chisq"
      message("residual df unavailable, test set to 'Chisq'")
    }
  }
  else {
    df <- error.df
  }
  if (is.null(df))
    df <- Inf
  if (df == 0)
    stop("residual df = 0")
  V <- if (is.null(vcov.))
    stas::vcov(model, complete = FALSE)
  else if (is.function(vcov.))
    vcov.(model)
  else vcov.
  b <- coef.
  if (any(aliased <- is.na(b)) && !singular.ok)
    stop("there are aliased coefficients in the model")
  b <- b[!aliased]
  if (is.null(b))
    stop(paste("there is no coef() method for models of class",
               paste(class(model), collapse = ", ")))

  L <- if (is.null(dim(hypothesis.matrix)))
    t(hypothesis.matrix)
  else
    hypothesis.matrix
  if (is.null(rhs))
    rhs <- rep(0, nrow(L))
  q <- NROW(L)
  value.hyp <- L %*% b - rhs
  vcov.hyp <- L %*% V %*% t(L)
  if (verbose) {
    cat("\nHypothesis matrix:\n")
    print(L)
    cat("\nRight-hand-side vector:\n")
    print(rhs)
    cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs)\n")
    print(drop(value.hyp))
    cat("\n")
    if (length(vcov.hyp) == 1)
      cat("\nEstimated variance of linear function\n")
    else cat("\nEstimated variance/covariance matrix for linear function\n")
    print(drop(vcov.hyp))
    cat("\n")
  }
  SSH <- as.vector(t(value.hyp) %*% solve(vcov.hyp) %*% value.hyp)
  test <- match.arg(test)
  if (!(is.finite(df) && df > 0))
    test <- "Chisq"
  name <- try(stats::formula(model), silent = TRUE)
  if (inherits(name, "try-error"))
    name <- substitute(model)
  title <- "Linear hypothesis test\n\nHypothesis:"
  topnote <- paste("Model 1: restricted model", "\n", "Model 2: ",
                   paste(deparse(name), collapse = "\n"), sep = "")
  note <- if (is.null(vcov.) || suppress.vcov.msg)
    ""
  else "\nNote: Coefficient covariance matrix supplied.\n"
  rval <- matrix(rep(NA, 8), ncol = 4)
  colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test,
                                                  ")", sep = ""))
  rownames(rval) <- 1:2
  rval[, 1] <- c(df + q, df)
  if (test == "F") {
    f <- SSH/q
    p <- stats::pf(f, q, df, lower.tail = FALSE)
    rval[2, 2:4] <- c(q, f, p)
  }
  else {
    p <- stats::pchisq(SSH, q, lower.tail = FALSE)
    rval[2, 2:4] <- c(q, SSH, p)
  }
  if (!(is.finite(df) && df > 0))
    rval <- rval[, -1]
  result <- structure(as.data.frame(rval), heading = c(title,
                                                       printHypothesis(L, rhs, names(b)), "", topnote, note),
                      class = c("anova", "data.frame"))
  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp
  result
}

printHypothesis = function (L, rhs, cnames)
{
  hyp <- rep("", nrow(L))
  for (i in 1:nrow(L)) {
    sel <- L[i, ] != 0
    h <- L[i, sel]
    h <- ifelse(h < 0, as.character(h), paste("+", h, sep = ""))
    nms <- cnames[sel]
    h <- paste(h, nms)
    h <- gsub("-", " - ", h)
    h <- gsub("+", "  + ", h, fixed = TRUE)
    h <- paste(h, collapse = "")
    h <- gsub("  ", " ", h, fixed = TRUE)
    h <- sub("^\\ \\+", "", h)
    h <- sub("^\\ ", "", h)
    h <- sub("^-\\ ", "-", h)
    h <- paste(" ", h, sep = "")
    h <- paste(h, "=", rhs[i])
    h <- gsub(" 1([^[:alnum:]_.]+)[ *]*", "", gsub("-1([^[:alnum:]_.]+)[ *]*",
                                                   "-", gsub("- +1 +", "-1 ", h)))
    h <- sub("Intercept)", "(Intercept)", h)
    h <- gsub("-", " - ", h)
    h <- gsub("+", "  + ", h, fixed = TRUE)
    h <- gsub("  ", " ", h, fixed = TRUE)
    h <- sub("^ *", "", h)
    hyp[i] <- h
  }
  hyp
}
