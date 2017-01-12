\name{rld.fit}
\alias{rld.fit}

\title{
Discrete-Time Survival Model Fitting Function
}

\description{
This is the function called by \code{rld} that do the actual computation.
}

\usage{
rld.fit(X, C, delta, initial, lower, upper, frailty)
}

\arguments{
  \item{X}{a design matrix created based on the input formula.}

  \item{C}{a vector containing censoring times.}

  \item{delta}{a vector contatining censoring status in which"1" denotes failure, "0" denotes right-censoring.}

  \item{initial}{a vector containing initial values for the paramters to be optimized over.}

  \item{lower}{a vector containing lower bound values for the paramters to be optimized over.}

  \item{upper}{a vector containing upper bound values for the paramters to be optimized over.}

  \item{frailty}{logical value: If TRUE, a discrete-time survival model with random effects will be used.}
}

\details{
\code{optim} is used to maximize the log-likelihood function. Method "L-BFGS-B" is that of Byrd et. al (1995) which allows box constraints, that is each variable can be given a lower and/or upper bound.
}

\value{
  \item{coefficients}{parameter estimates.}

  \item{hessian}{hessian matrix.}

  \item{LikFunValue}{log-likelihood value.}
}

\references{
Kang, C., Huang, Y., and Miller, C. (2015). A discrete-time survival model with random effects for designing and analyzing repeated low-dose challenge experiments. \emph{Biostatistics}, 16(2): 295-310.

Byrd, R.H., Lu, P. Nocedal, J. and Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. \emph{SIAMJ. Scientific Computing}, 16, 1190-1208.

}

\author{
Bin Yao, Ying Huang and Chaeryon Kang
}

\note{
Don't use this function when you are not familiar with the whole computational procedure.
}

\seealso{
 \code{\link{rld}}, \code{\link{optim}}
}

\examples{
data(SampleData)
newdata <- transdata(data = SampleData, ndlevel = 3, nexposure = c(10, 10, 2))
formulaexp <- survival::Surv(time, delta)~factor(dose)+trt+I(I(dose==3)*trt)
designmat <- model.matrix(formulaexp, data = newdata)
time <- SampleData$time
status <- SampleData$delta
inival <- rep(0.5, 6)
lwrval <- c(rep(-Inf, 5), 0.01)
uprval <- rep(Inf, 6)
frailty <- TRUE

rld.fit(X = designmat, C = time, delta = status, initial = inival,
        lower = lwrval, upper = uprval, frailty = frailty)
}

\keyword{fitting function}
