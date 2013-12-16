\name{crr.wt.COX.str}
\alias{crr.wt.COX.str}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression modeling of subdistribution functions using the stratified Cox weight.
}
\description{
Fits a subdistribution hazards regression model using the stratified Cox model adjusted weight. Use predicted censoring probability based on the stratified Cox model for each individual.
}
\usage{
crr.wt.COX.str(t, ic, z, zc, strata.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{a vector of observed times.}
  \item{ic}{a vector of cause indicators.}
  \item{z}{a matrix of covariates for fitting the model.}
  \item{zc}{a matrix of covariates for fitting the censoring distribution.}
  \item{strata.var}{a vector contains information about strata membership.}
}
\details{
Same estimating procedure with a different weight.
}
\value{a "crrwt" class object.}

\author{Peng He}
