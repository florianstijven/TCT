% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parametric-bootstrap.R
\name{pm_bootstrap_vertical_to_common}
\alias{pm_bootstrap_vertical_to_common}
\title{Parametric Bootstrap for a Common Acceleration Factor}
\usage{
pm_bootstrap_vertical_to_common(
  TCT_Fit,
  inference,
  B = 100,
  bs_fix_vcov = TRUE,
  return_se = TRUE,
  null = FALSE,
  select_coef,
  constraints = FALSE,
  type,
  weights,
  start_gamma
)
}
\arguments{
\item{TCT_Fit}{Object returned by \code{\link[=TCT_meta]{TCT_meta()}}}

\item{inference}{Which approach is used for estimation and inference? Should
be \code{"wald"}, \code{"score"}, or \code{"least-squares"}. The \code{"wald"}-approach is explained and
implemented in \code{\link[=DeltaMethod]{DeltaMethod()}}, the \code{"score"}-approach is explained and
documented in \code{\link[=score_test]{score_test()}} and \code{\link[=score_test_common]{score_test_common()}}.}

\item{B}{Number of bootstrap replications.}

\item{bs_fix_vcov}{(boolean) Fix the estimated variance-covariance matrix for
the estimated acceleration factors? This speeds up computations, but can have
a negative impact of the method's properties.}

\item{return_se}{(boolean) Return the estimated standard error from each
bootstrap replication? This standard error is computed with the delta method.}

\item{null}{(boolean): conduct the bootstrap under the null hypothesis of no
treatment effect? Defaults to \code{FALSE}.}

\item{select_coef}{Estimates from the \code{exp_estimates} in \code{\link[=TCT_meta]{TCT_meta()}} to use
in estimating the common acceleration factor. If there is reason to believe
that the proportional slowing assumption does not hold, e.g., for the first
measurement after randomization, then the corresponding estimate should not
be used in estimation the common acceleration factor.}

\item{constraints}{Use the constrained generalized least squares estimator
for the vertical treatment effects.}

\item{type}{Which type of test statistic should be used. See Test Statistic
Variants. Should be one of
\enumerate{
\item \code{type = "omnibus"}
\item \code{type = "directional"}
\item \code{type = "inverse variance"}
\item \code{type = "custom"}
}}

\item{weights}{If \code{type == "custom"}, the user should specify a weight
vector for weighting estimates at different time points differently.}

\item{start_gamma}{Starting value for the acceleration factor used for for the
generalized least squares estimator.}
}
\value{
A list with two element:
\itemize{
\item \code{estimates_bootstrap}: (numeric) vector of length \code{B} that contains the
bootstrap replicates for the common acceleration factor.
\item \code{se_bootstrap}: (numeric) vector of length \code{B} that contains the bootstrap
replicates of the estimated standard error. Is a vector of \code{NA}s when
\code{return_se = FALSE}.
}
}
\description{
The \code{\link[=pm_bootstrap_vertical_to_common]{pm_bootstrap_vertical_to_common()}} function implements a parametric
bootstrap for the time-specific acceleration factors.
}
