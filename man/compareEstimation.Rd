% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leptokurticMixture.R
\name{compareEstimation}
\alias{compareEstimation}
\title{Compare the two methods of estimation}
\usage{
compareEstimation(
  mod = NULL,
  data = NULL,
  G = NULL,
  n = 10^4,
  tol = 1e-06,
  wt = NULL,
  n0 = 25,
  lab = NULL
)
}
\arguments{
\item{mod}{A character of length 4 such as "VVVV", indicating the model; the covariance and beta parameters.}

\item{data}{A n x p matrix of observations.}

\item{G}{The number of components to fit.}

\item{n}{The maximum number of EM iterations.}

\item{tol}{The tolerance for the stopping rule; lack of progress. The default is 1e-6 but it depends on the dataset.}

\item{wt}{a (n x d) matrix of weights for initialization if NULL, then a random weight matrix is generated.}

\item{n0}{Given wt, the number of iterations used to obtain the initial parameters}

\item{lab}{Using given labels (lab) as starting values.}
}
\value{
A vector of times, number of iterations and log-likelihood values.
}
\description{
Compare the two methods of estimation for fitting a finite mixture of multivariate elliptical leptokurtic-normal distributions; fixed point iterations and MM algorithm.
}
