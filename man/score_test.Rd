% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score-test.R
\name{score_test}
\alias{score_test}
\title{Compute TCT score test z-value}
\usage{
score_test(
  time_points,
  ctrl_estimates,
  exp_estimates,
  ref_fun,
  interpolation,
  vcov,
  j,
  gamma_0 = 1
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
of \code{c(ctrl_means, exp_means)}.}

\item{j}{Measurement occasion to test acceleration factor for. This
corresponds to \code{exp_estimates[j]}.}

\item{gamma_0}{Value for the acceleration factor under the null hypothesis.}
}
\value{
Named (numeric) vector with two elements:
\enumerate{
\item \code{"z"}. This test-statistic follows a standard normal distribution under
the null hypothesis.
\item \code{"p-value"}: two-sided p-value.
}
}
\description{
The \code{\link[=score_test]{score_test()}} function computes the z-value for the score test at
measurement \code{j}. This corresponds to \code{exp_estimates[j]}. For technical
details on the score test, see Details.
}
\section{Score Test}{
For constructing a score test, we start from the null hypothesis,
\deqn{H_{\gamma_0}: \gamma_j = \gamma_0,} where \eqn{\gamma_j} is the
acceleration factor at time \eqn{t_j}. Under this null hypothesis, and taking
the control group as reference, we can
compute the experimental treatment group mean at \eqn{t_j} as \deqn{\beta_{0,
j} = f_0( \gamma_0 \cdot t_{j}; \boldsymbol{\alpha})} given the true control
group mean vector, \eqn{\boldsymbol{\alpha}}. Consequently, a test for
\eqn{H_0: \beta_j = \beta_{0, j}} is also a valid test for
\eqn{H_{\gamma_0}}. Well-established tests are available for the former
hypothesis. However, these cannot be applied directly since
\eqn{\boldsymbol{\alpha}} is not known, but only estimated. Indeed, \eqn{\beta_{0,
j}} is itself estimated by \eqn{\hat{\beta}_{0, j} =  f_0( \gamma_0 \cdot t_{j};
\hat{\boldsymbol{\alpha}})}. Still, this provides the starting point for a
score test that draws upon well-established statistical tests.

We first define the score vector as follows, \deqn{\boldsymbol{s}(\gamma
\cdot \boldsymbol{t}; \hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}}) =
\hat{\boldsymbol{\beta}} - \boldsymbol{f_0}(\gamma \cdot \boldsymbol{t};
\hat{\boldsymbol{\alpha}})}
where \eqn{\boldsymbol{t} = (t_1, ..., t_K)'} and \eqn{\boldsymbol{f_0}(\gamma \cdot
\boldsymbol{t}; \hat{\boldsymbol{\alpha}}) = (f_0(\gamma \cdot t_1;
\hat{\boldsymbol{\alpha}}), ..., f_0(\gamma \cdot t_K;
\hat{\boldsymbol{\alpha}}))'}. Under \eqn{H_{\gamma_0}} for all time points and using the delta
method, we have that \deqn{\boldsymbol{s}(\gamma \cdot \boldsymbol{t};
\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}}) \; \dot\sim \; N \left(
\boldsymbol{0}, \Sigma_s \right)} where \deqn{\Sigma_s =
J_s(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0}; \gamma_0) \cdot \Sigma \cdot
J_s(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0}; \gamma_0)^t} with
\eqn{J_s(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0}; \gamma_0)} the
Jacobian of \eqn{\boldsymbol{s}(\gamma \cdot \boldsymbol{t}; \boldsymbol{\alpha},
\boldsymbol{\beta})} evaluated in the true parameter vector,
\eqn{(\boldsymbol{\alpha_0}, \boldsymbol{\beta_0})'}, with \eqn{\gamma} held
constant at \eqn{\gamma_0}. Note that this matrix is a function of
\eqn{\boldsymbol{\alpha}}, \eqn{\boldsymbol{\beta}}, and \eqn{\gamma_0}. The
Jacobian matrix with the true parameter vector replaced by the corresponding
estimates, \eqn{(\hat{\boldsymbol{\alpha}}, \hat{\boldsymbol{\beta}})'}, is
denoted by \eqn{\hat{\Sigma}_s}.
}

\section{Time-Specific Acceleration Factor}{
From the distributional result in the previous section follows that \deqn{z =
\frac{s(\gamma_0 \cdot t_j; \hat{\boldsymbol{\alpha}},
\hat{\boldsymbol{\beta}})}{\sqrt{\hat{\Sigma}_{s,jj}}} \; \dot\sim \;N(0, 1)}
under \eqn{H_{\gamma_0}} where \eqn{\hat{\Sigma}_{s,jj}} is the \eqn{j}'th
diagonal element of \eqn{\hat{\Sigma}_s}. Note that this result relies on the
delta method and is thus only asymptotic. The finite sample properties of
this test will thus depend to a large extent on the accuracy of the delta
method, i.e., the degree of non-linearity of
\eqn{g_{\gamma_0}(\boldsymbol{\alpha}, \boldsymbol{\beta}; t_j)} around the
true parameter values. The \code{\link[=score_test]{score_test()}} function returns the above z-value.
}

