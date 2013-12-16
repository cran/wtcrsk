\name{plot.crrwt.pred}
\alias{plot.crrwt.pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot function for a 'crrwt.pred' class object.
}
\description{
Making plots of cumulative incidence function over time.
}
\usage{
\method{plot}{crrwt.pred}(x, multiple = 0, se = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a 'crrwt.pred' class object.}
  \item{multiple}{if multiple = 0, plots are overlayed. If multiple = 1, plots are separated.}
  \item{se}{if se = 0, no 95\% point-wise confidence interval. If se = 1, a 95\% point-wise confidence interval is included. "se" option only works when multiple = 1 (separate plots).}
  \item{...}{same as the generic 'plot' function}
}
\author{Peng He}
