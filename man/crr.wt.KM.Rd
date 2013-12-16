\name{crr.wt.KM}
\alias{crr.wt.KM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression modeling of subdistribution functions using Kaplan-Meier weight.
}
\description{
Fits a subdistribution hazards regression model using the Kaplan-Meier weight. This is the same model as Fine and Gray (1999).
}
\usage{
crr.wt.KM(t, ic, z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{a vector of observed times.}
  \item{ic}{a vector of cause indicators.}
  \item{z}{a matrix of covariates for fitting the model.}
}
\details{
More details are given in the reference.
}
\value{a "crrwt" class object.}
\references{
He P, Scheike TH and Zhang MJ, A proportional hazards regression model for the subdistribution with covariates adjusted censoring weight for competing risks data, Technical Report #61, Division of Biostatistics, Medical College of Wisconsin, November 2013.

Fine JP and Gray RJ (1999), A proportional hazards model for the subdistribution of a competing risk. Journal of The American Statistical Association 94, 496-509.
}
\author{
Peng He
}