% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score-test.R
\name{score_conf_int}
\alias{score_conf_int}
\title{Compute Confidence Interval Based on Score Test}
\usage{
score_conf_int(
  time_points,
  ctrl_estimates,
  exp_estimates,
  ref_fun,
  interpolation,
  vcov,
  j,
  alpha = 0.05,
  bounds = c(-5, 5)
)
}
\arguments{
\item{time_points}{Ordered vector that contains the times corresponding to
the estimated means in the \code{ctrl_estimates} vector. The first element
should be zero which generally corresponds to the time of randomization.}

\item{ctrl_estimates}{Estimated mean outcome in the control group at the
fixed occasions corresponding to the times in \code{time_points}.}

\item{exp_estimates}{Estimated mean outcomes in the experimental group at
fixed occasions corresponding to the times in \code{time_points[-1]}. Note that
first element in \code{exp_estimates} should correspond to some time after
randomization.}

\item{ref_fun}{Intra- and extrapolation function that is returned by
\code{ref_fun_constructor()}.}

\item{interpolation}{Which interpolation method to use?
\itemize{
\item \code{"linear"}: linear interpolation.
\item \code{"spline"}: natural cubic spline interpolation. This interpolation method has been most
thoroughly tested is most stable.
\item \verb{"monoH.FC}: monotone Hermite spline according to the method of Fritsch
and Carlson.
}}

\item{vcov}{The variance-covariance matrix for the means. In order to map to
the correct estimates, this matrix should be the variance-covariance matrix
of \code{c(ctrl_means, exp_means)}. If an element of \code{c(ctrl_means, exp_means)}
is known (e.g., mean at baseline is zero when using change from baseline as
outcome), then the corresponding row and column in \code{vcov} should be set to
zero.}

\item{j}{Measurement occasion to test acceleration factor for. This
corresponds to \code{exp_estimates[j]}.}

\item{alpha}{\code{1 - alpha} represents the two-sided confidence level. Defaults
to \code{0.05}.}

\item{bounds}{(numeric) vector with two elements that defines the bounds of
the search interval for the (score-based) confidence limits or estimate.
Defaults to \code{c(-5, 5)}.}
}
\value{
(numeric) vector with two elements. The first element is the lower
confidence limit, the second element is the upper confidence limit.
}
\description{
The \code{\link[=score_conf_int]{score_conf_int()}} function computes the confidence interval for the
score test at measurement \code{j}. This confidence interval is based on the
score test implemented in \code{\link[=score_test]{score_test()}}.
}
\section{Score Test Confidence Intervals}{
For the construction of a confidence interval based on the scores, we make
use of the relationship between confidence intervals and hypothesis tests. A
\eqn{1 − α} confidence interval can be defined as follows, \deqn{\left\{
\gamma : p(\gamma) > \alpha \right\}} where \eqn{p(\gamma)} is the p-value
for the null that the acceleration factor is equal to \eqn{\gamma}. This
p-value can be based on the time-specific score test implemented in
\code{\link[=score_test]{score_test()}}, or can be based on the joint score test implemented in
\code{\link[=score_test_common]{score_test_common()}}. The former gives a confidence interval for the
time-specific acceleration factor while the latter gives a confidence
interval for the common acceleration factor. The latter interval should be
used with some care. Indeed, the confidence interval for a common
acceleration factor is only sensible if there is no evidence against the
assumptions of constant slowing of the disease progression.
}

