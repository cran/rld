\name{rld}
\alias{rld}
\alias{print.rld}
\alias{summary.rld}
\alias{print.summary.rld}

\title{
Fit a Discrete-Time Survival Model
}

\description{
This function fits a discrete-time survival model with and without random effects.
}

\usage{
rld(formula, data, na.action, initial=NULL, lower=NULL, upper=NULL, frailty=TRUE)
}

\arguments{
  \item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by \code{Surv} function. The terms is a series of terms which specify linear predictors for response.}

  \item{data}{a data.frame in which to interpret the variables named in the formula. This augmented data frame can be returned by function \code{transdata}.}

  \item{na.action}{a function which indicates what should happen when the data contain NAs.}

  \item{initial}{a vector of initial values for the paramters to be optimized over. If NULL, the default initial values will be used.}

  \item{lower}{a vector of lower bound values for the paramters. If NULL, the default lower bound will be used.}

  \item{upper}{a vector of upper bound values for the parameters. If NULL, the default upper bound will be used.}

  \item{frailty}{logic value: if TRUE, the discrete-time survival model with random effects will be run. Otherwise it is assumed that there is no random effect. The default is TRUE.}
}

\details{
Kang et al. (2015) proposed to use a discrete-time survival model with gamma-distributed random effects and a complementary log-log link function to model data from repeated low-dose challenge studies, assuming an animal's risks of infection across challenges are independent of each other conditional on random effects. Please refer to Kang et al.(2015) for more details.
}

\value{
\code{rld} returns an object of class ``rld''. The functions \code{summary} is used to obtain and print a summary of the results.

  \item{coefficients}{a vector of parameter estimates.}

  \item{hessian}{the hessian matrix returned from the function \code{optim}.}

  \item{X}{the design matrix created based on the input formula.}

  \item{VEexpr}{the formula expression on the right of ~ operator.}

  \item{loglikvalue}{the log-likelihood value.}

  \item{call}{the matched call.}

  \item{frailty}{the chosen model.}

  \item{augdata}{the augmented data set.}
}

\references{
Kang, C., Huang, Y., and Miller, C. (2015). A discrete-time survival model with random effects for designing and analyzing repeated low-dose challenge experiments. \emph{Biostatistics}, 16(2): 295-310.
}

\author{
Bin Yao, Ying Huang and Chaeryon Kang
}

\seealso{
\code{\link{rld.fit}}
}

\examples{
data(SampleData)
newdata <- transdata(data = SampleData, ndlevel = 3, nexposure = c(10, 10, 2))

#interaction between the hightest dose level and treatment assignment
#under the discrete-time survival model with random effects

fitout1 <- rld(formula = survival::Surv(time, delta)~factor(dose)+trt+I(I(dose==3)*trt),
              data = newdata, frailty = TRUE)
#summary(fitout1)

\dontrun{
#main effects model without random effectss
ini <- rep(0.5, 4)
lwr <- rep(-Inf, 4)
upr <- rep(Inf, 4)

fitout2 <- rld(formula = survival::Surv(time, delta)~factor(dose)+trt,
               initial = ini, lower = lwr, upper = upr, data = newdata,
               frailty = FALSE)
#summary(fitout2)
}
}

\keyword{rld}
\keyword{Discrete-time survival model}
