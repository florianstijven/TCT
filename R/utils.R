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
#' * `"monoH.FC`: monotone Hermite spline according to the method of Fritsch and
#'      Carlson.
#'
#' @return A vector of the time-mapped values.
#'
get_new_time = function(y_ref, x_ref, y_obs, method = "linear") {
  # amount of extrapolation
  extrapol = 1e3
  # Only the first first value in y_obs is used. This enables us to use
  # straightforward data handling methods from dplyr.
  if (method == "linear") {
    interpol_fun = stats::approxfun(x_ref,
                                    y_ref,
                                    method = "linear",
                                    rule = 1)
  }
  else if (method == "spline") {
    interpol_fun = stats::splinefun(x_ref,
                                    y_ref,
                                    method = "natural")
  }
  else if (method == "monoH.FC") {
    interpol_fun = stats::splinefun(x_ref,
                                    y_ref,
                                    method = "monoH.FC")
  }

  p = length(x_ref)
  extrapol_fun = extrapol_fun_factory(x_ref, y_ref)
  ref_fun = ref_fun_factory(x_ref, extrapol_fun, interpol_fun)


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
      x_mapped[i] = stats::uniroot(
        f = function(x)
          ref_fun(x) - y_obs[i],
        interval = c(min(x_ref) - extrapol,
                     max(x_ref) + extrapol),
        tol = .Machine$double.eps^0.5,
        maxiter = 1e3
      )$root
    }
  }
  return(x_mapped)
}

extrapol_fun_factory = function(.x_ref, .y_ref) {
  # standard function does not do extrapolation, this extrapolation is
  # implemented here.
  p = length(.y_ref)
  y2 = .y_ref[p]
  y1 = .y_ref[1]
  x2 = .x_ref[p]
  x1 = .x_ref[1]
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

ref_fun_factory = function(.x_ref, extrapol_fun, interpol_fun, deriv = 0) {
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
  numericDeriv(expr = quote({
    ref_fun = ref_fun_constructor(x_ref, y_ref, method)
    ref_fun(t_m)
  }),
  theta = c("y_ref"),
  rho = myenv)
}

deriv_f0_alpha_bis = function(t_m, x_ref, y_ref, finite_diff = 1e-6, method = "spline") {
  numDeriv::jacobian(
    func = function(y_ref) {
      ref_fun = ref_fun_constructor(x_ref, y_ref, method)
      ref_fun(t_m)
    },
    x = y_ref
  )
}

ref_fun_constructor = function(x_ref, y_ref, method) {
  if (method == "linear") {
    interpol_fun = stats::approxfun(x_ref,
                                    y_ref,
                                    method = "linear",
                                    rule = 1)
  }
  else if (method == "spline") {
    interpol_fun = stats::splinefun(x_ref,
                                    y_ref,
                                    method = "natural")
  }
  else if (method == "monoH.FC") {
    interpol_fun = stats::splinefun(x_ref,
                                    y_ref,
                                    method = "monoH.FC")
  }

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
    df <- df.residual(model)
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
    vcov(model, complete = FALSE)
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
  if (is.character(hypothesis.matrix)) {
    L <- makeHypothesis(names(b), hypothesis.matrix, rhs)
    if (is.null(dim(L)))
      L <- t(L)
    rhs <- L[, NCOL(L)]
    L <- L[, -NCOL(L), drop = FALSE]
    rownames(L) <- hypothesis.matrix
  }
  else {
    L <- if (is.null(dim(hypothesis.matrix)))
      t(hypothesis.matrix)
    else hypothesis.matrix
    if (is.null(rhs))
      rhs <- rep(0, nrow(L))
  }
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
  name <- try(formula(model), silent = TRUE)
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
    p <- pf(f, q, df, lower.tail = FALSE)
    rval[2, 2:4] <- c(q, f, p)
  }
  else {
    p <- pchisq(SSH, q, lower.tail = FALSE)
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
