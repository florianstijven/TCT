% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time_scale_effects.R
\name{summary.TCT_meta}
\alias{summary.TCT_meta}
\title{Summarize fitted Meta-Time Component Test model}
\usage{
\method{summary}{TCT_meta}(
  object,
  alpha = 0.05,
  delta_transformation = "identity",
  bounds = c(-5, 5),
  ...
)
}
\arguments{
\item{object}{Object returned by \code{\link[=TCT_meta]{TCT_meta()}}.}

\item{alpha}{Two-sided confidence level for confidence intervals.}

\item{delta_transformation}{Transformation when applying the delta-method to
obtain confidence intervals.}

\item{bounds}{(numeric) vector with two elements that defines the bounds of
the search interval for the (score-based) confidence limits or estimate.
Defaults to \code{c(-5, 5)}.}

\item{...}{
  Arguments passed on to \code{\link[base:summary]{base::summary}}
  \describe{
    \item{\code{}}{}
  }}
}
\value{
S3 object of class \code{"summary_TCT_meta"}
}
\description{
Summarize fitted Meta-Time Component Test model
}
\examples{
# transform example data set to desired format
library(dplyr)
data = simulated_test_trial \%>\%
mutate(time_int = (Week \%/\% 25)) \%>\%
  arrange(trial_number, SubjId, time_int) \%>\%
  mutate(time_int = as.integer(time_int) + 1L) \%>\%
  mutate(arm_time = ifelse(time_int == 1L,
                           "baseline",
                           paste0(arm, ":", time_int)))
# fit e.g. MMRM model to obtain estimates of profiles
mmrm_fit = analyze_mmrm(data)
set.seed(1)
TCT_Fit = TCT_meta(
  time_points = 0:4,
  ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
  exp_estimates = coef(mmrm_fit)[5:8],
  vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
  interpolation = "spline",
  B = 1e3
)
# The summary() generic can be used to obtain the most useful quantities from
# the Meta-TCT.
summary(TCT_Fit)
}
