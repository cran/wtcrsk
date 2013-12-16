\name{predict.crrwt}
\alias{predict.crrwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute predictive CIFs for given set of covariates.
}
\description{
This is a function to calculate prediction of the cumulative incidence function (CIF) as well as its variance at observed failure times given in original data.
}
\usage{
\method{predict}{crrwt}(object, z, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a 'crrwt' class object obtained from crr.wt function.}
  \item{z}{sets of covariates used for prediction, each row represents a new set of covariates.}
  \item{...}{additional arguments affecting the predictions produced.}
}
\details{
More derivations are given in the reference.
}
\value{
\item{z}{given sets of covariates.}
\item{time}{observed failure times.}
\item{F1}{predicted cumulative incidence probabilities at observed failure times.}
\item{F1sd}{standard errors of predicted cumulative incidence probabilities.}
}
\references{
He P, Scheike TH and Zhang MJ, A proportional hazards regression model for the subdistribution with covariates adjusted censoring weight for competing risks data, Technical Report #61, Division of Biostatistics, Medical College of Wisconsin, November 2013.

}
\author{
Peng He
}
