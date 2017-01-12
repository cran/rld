\name{calcpk}
\alias{calcpk}
\alias{print.calcpk}
\alias{summary.calcpk}
\alias{print.summary.calcpk}

\title{
Calculate Per-Challenge Probability of Infection
}

\description{
This function estimates the per-challenge probability of infection.
}

\usage{
calcpk(object, predlevel, CIlevel = 0.95)
}

\arguments{
  \item{object}{a fitted object of class inheriting from "rld".}

  \item{predlevel}{a value or a vector indicating covariate for prediction of per-challenge probability of infection.}

  \item{CIlevel}{confidence level. The default is 0.95.}
}

\details{
Calculate the per-challenge risk of infection for the contrast group and reference group.
}

\value{
  \item{pk}{risk of infection for the contrast group.}

  \item{pk_se}{standard error of the estimated pk.}

  \item{lwr}{lower bound value of the confidence interval of the estimated pk.}

  \item{upr}{upper bound value of the confidence interval of the estimated pk.}

}

\references{
Kang, C., Huang, Y., and Miller, C. (2015). A discrete-time survival model with random effects for designing and analyzing repeated low-dose challenge experiments. \emph{Biostatistics}, 16(2): 295-310.
}

\author{
Bin Yao, Ying Huang and Chaeryon Kang
}

\seealso{
\code{\link{calcVEk}}, \code{\link{calcVEt}}
}

\note{The variable names in the predlevel should match the names in the data frame.}

\examples{
data(SampleData)
newdata <- transdata(data = SampleData, ndlevel = 3, nexposure = c(10, 10, 2))
fitout <- rld(formula = survival::Surv(time, delta)~factor(dose)+trt+I(I(dose==3)*trt),
              data = newdata, frailty = TRUE)

predictdata <- 1
names(predictdata) <- c("trt")

pkout <- calcpk(object = fitout, predlevel = predictdata, CIlevel = 0.95)
summary(pkout)
}

\keyword{per-challenge probability of infection}
