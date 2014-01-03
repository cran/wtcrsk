\name{crr.wt}
\alias{crr.wt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Competing risks regression using different weights.
}
\description{
Regression modeling of subdistribution functions in competing risks with adjustments for covariate-dependent censoring.
}
\usage{
crr.wt(formula, data, weight = c("KM", "COX", "KM.Strata", "COX.Strata"), cens.formula, 
cause = 1, strata.var, variance = TRUE, var.conservative = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a competing risks object as returned by the Crsk function.}
  \item{data}{a data.frame in which to interpret the variables named in the formula.}
  \item{weight}{type of weight for modeling. "KM" - Kaplan-Meier weight. "COX" - Cox regression adjusted weight. "KM.strata" - stratified Kaplan-Meier adjusted weight. "COX.strata" - stratified COX regression adjusted weight.}
  \item{cens.formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a censoring object in the format of \code{Cens(t,ic)}. This option is needed only for weight="COX" or weight="COX.strata". See the example for more details}
  \item{cause}{code of cause that is of interest.}
  \item{strata.var}{a categorical stratification variable, needed only for weight="KM.strata" or 
  weight="COX.strata".}
  \item{variance}{Includes variance calculation if TRUE, Skips variance calculation if FALSE.}
  \item{var.conservative}{Includes minor term of variance of beta if FALSE, no minor term for variance of beta if TRUE.}
}
\details{
Fits the subdistribution hazards regression model with adjustments to covariate-dependent censoring as described in He et al. (2013). In addition, this function includes two new weights - "KM.strata" and "COX.strata" for fitting stratified models for the censoring distribution. \code{strata.var} is an option for users to specify the stratification variable, which has to be categorical with levels. When fitting the regression model using the "COX.strata" weight, the stratification variable cannot be included in your censoring-dependent covariates set. If \code{cens.formula} was not specified but either "COX" or "COX.strata" weight was used, by default the function uses "z" as the set of covariates for \code{cens.formula} with necessary removal of the specified stratification variable. 

Return values W.lambda and W.beta are used for prediction calculations. Incomplete observations are removed before the model fitting.
}
\value{
  \item{weight}{type of weight used for modeling the censoring distribution.}
  \item{varname}{covariate names.}
  \item{converge}{1 if the iterative algorithm coverged. 0 if not.}
  \item{cens.det}{for weight="COX" and weight="COX.strata" only, 0 if there is a problem with the determinant of the information matrix when fitting the censoring distribution, 1 if no problem with determinant.}
  \item{beta}{estimated regression coefficients.}
  \item{beta_se}{estimated standard errors for regression coefficients.}
  \item{time}{a vector of observed failure times}
  \item{a10}{a vector of estimated cumulative baseline subdistribution hazards at observed failure times.}
  \item{a10se}{a vector of estimated standard errors for estimated cumulative baseline subdistribution hazards.}
  \item{W_lambda}{used for prediction function.}
  \item{W_beta}{used for prediction function.}
  \item{n}{total number of observations.}
  \item{n.missing}{number of observations with missing values.}
}
\references{
He P, Scheike TH and Zhang MJ, A proportional hazards regression model for the subdistribution with covariates adjusted censoring weight for competing risks data, Technical Report #61, Division of Biostatistics, Medical College of Wisconsin, November 2013.
}
\author{
Peng He
}

\seealso{
\code{\link{predict.crrwt}},\code{\link{print.crrwt}},\code{\link{summary.crrwt}},\code{\link{crr.wt.KM}},\code{\link{crr.wt.COX}},\code{\link{crr.wt.KM.str}},\code{\link{crr.wt.COX.str}}
}
\examples{
# Simulated data
set.seed(123321)
t <- rexp(200)
ic <- sample(0:2,200,replace=TRUE)
z <- matrix(runif(600),nrow=200)
colnames(z) <- c("z1","z2","z3")
dat <- data.frame(t,ic,z)

# Model fitting with a Kaplan-Meier weight
print(out.KM <- crr.wt(Crsk(t,ic)~z1+z2+z3,data=dat,weight="KM",cause=1))

# Model fitting with a Cox weight
print(out.COX <- crr.wt(Crsk(t,ic)~z1+z2+z3,data=dat,weight="COX",cause=1,Cens(t,ic)~z1+z2))

# Summary information
summary(out.COX)

# Prediction of cumulative incidence proability
newdata <- data.frame(c(.2,.3,.5),c(.9,.1,.6),c(.3,.9,.2))
colnames(newdata) <- c("z1","z2","z3")
pred.COX <- predict(out.COX,z=newdata)

# Plots of prediction
plot(pred.COX)
plot(pred.COX,multiple=1)
plot(pred.COX,multiple=1,se=1)
}
