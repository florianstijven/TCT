% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmrm.R
\name{analyze_mmrm}
\alias{analyze_mmrm}
\title{Analyze trial data with MMRM}
\usage{
analyze_mmrm(data_trial, method = "ML", type = "full", package = "mmrm")
}
\arguments{
\item{data_trial}{A data frame from a single clinical trial. The following
variables should be present:
\itemize{
\item \code{ADAScog_integer}: integer-valued ADAScog score
\item \code{arm_time}: factor that characterizes the treatment arm-time combination.
\item \code{time_int}: integer that indicates the measurement occasion. Should start
with \code{1} for the first occasion, values for subsequent occasions should be
consecutive integers.
}}

\item{method}{Which method should be used for fitting the MMRM? Should be one
of the following strings:
\itemize{
\item \code{"ML"} (default): maximum likelihood
\item \code{"REML"}: restricted maximum likelihood
}}

\item{type}{Which MMRM model should be use?
\itemize{
\item \code{"full"} (default): full model with treatment-time interaction.
\item \code{"null"}: null model with no treatment effect. The only use case for
fitting this model is to perform likelihood ratio tests.
}}

\item{package}{Which package to use for fitting the mmrm? Shoud be one of:
\itemize{
\item \code{"mmrm"}
\item \code{"nlme"}
}}
}
\value{
fitted-model object from \code{\link[nlme:gls]{nlme::gls()}} or \code{\link[mmrm:mmrm]{mmrm::mmrm()}}
}
\description{
The \code{\link[=analyze_mmrm]{analyze_mmrm()}} function analyzes data from a single trial with the
mixed model for repeated measures (See Details).
}
\examples{

# load example data and add required variables (time_int and arm_time)
library(dplyr)
data = simulated_test_trial \%>\%
  mutate(time_int = (Week \%/\% 25)) \%>\%
  arrange(trial_number, SubjId, time_int) \%>\%
  mutate(time_int = as.integer(time_int) + 1L) \%>\%
  mutate(arm_time = ifelse(time_int == 1L,
                           "baseline",
                           paste0(arm, ":", time_int)))
mmrm_fit = analyze_mmrm(data)


}
