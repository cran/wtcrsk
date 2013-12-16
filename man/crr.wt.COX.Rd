\name{crr.wt.COX}
\alias{crr.wt.COX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression modeling of subdistribution function using the Cox weight.
}
\description{
Fits a subdistribution hazards regression model using the Cox model adjusted weight. Use predicted censoring probability based on the Cox model for each individual.
}
\usage{
crr.wt.COX(t, ic, z, zc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{a vector of observed times.}
  \item{ic}{a vector of cause indicators.}
  \item{z}{a matrix of covariates for fitting the model.}
  \item{zc}{a matrix of associated covariates for fitting the censoring distribution.}
}
\details{
More details are given in the reference.
}
\value{a "crrwt" class object.}
\references{
He P, Scheike TH and Zhang MJ, A proportional hazards regression model for the subdistribution with covariates adjusted censoring weight for competing risks data, Technical Report #61, Division of Biostatistics, Medical College of Wisconsin, November 2013.
}
\author{Peng He}

