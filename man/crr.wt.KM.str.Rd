\name{crr.wt.KM.str}
\alias{crr.wt.KM.str}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression modeling of subdistribution functions using the stratified Kaplan-Meier weight.
}
\description{
Fits a subdistribution hazards regression model using the stratified Kaplan-Meier adjusted weight. Use predicted censoring probability based on the stratified Kaplan-Meier estimate for each individual.
}
\usage{
crr.wt.KM.str(t, ic, z, strata.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{a vector of observed times.}
  \item{ic}{a vector of cause indicators.}
  \item{z}{a matrix of covariates for fitting the model.}
  \item{strata.var}{a vector contains information about strata membership.}
}
\details{
Same estimating procedure with a different weight.
}
\value{a "crrwt" class object.}

\author{Peng He}