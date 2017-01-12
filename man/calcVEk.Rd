\name{calcVEk}
\alias{calcVEk}
\alias{print.calcVEk}
\alias{summary.calcVEk}
\alias{print.summary.calcVEk}

\title{
Calculate Per-Challenge Vaccine Efficacy
}

\description{
This function calculate the estimated per-challenge vaccine efficacy.
}

\usage{
calcVEk(object, newdata, CIlevel=0.95)
}

\arguments{
  \item{object}{a fitted object of class inheriting from "rld".}

  \item{newdata}{a data list for predicting vaccine efficacy where "contrgroup" and "refgroup" list names must be included.}

  \item{CIlevel}{a confidence level. The default is 0.95.}
}

\details{
Per-challenge vaccine efficacy is defined as the relative reduction in the risk of infection caused by vaccination at a particular challenge, conditional on non-infection before the challenge. Please refer to Kang et al.(2015) for more details about the formula.
}

\value{
  \item{VE}{the vaccine efficacy estimates for contrast group and reference group.}

  \item{se}{standard deviations of per-challenge vaccine efficacy estimates.}

  \item{lwr}{a vector containing the lower bound values of confidence interval for VE.}

  \item{upr}{a vector containing the upper bound values of confidence interval for VE.}
}

\references{
Kang, C., Huang, Y., and Miller, C. (2015). A discrete-time survival model with random effects for designing and analyzing repeated low-dose challenge experiments. \emph{Biostatistics}, 16(2): 295-310.
}

\author{
Bin Yao, Ying Huang and Chaeryon Kang
}

\seealso{
\code{\link{calcVEt}}, \code{\link{calcpk}}
}

\note{In the newdata list for vaccine efficacy prediction, users must assign variable names in the contrast group and reference group.}

\examples{
data(SampleData)
augdata <- transdata(data = SampleData, ndlevel = 3, nexposure = c(10, 10, 2))
fitout <- rld(formula = survival::Surv(time, delta)~factor(dose)+trt+I(I(dose==3)*trt),
              data = augdata, frailty = TRUE)

contrgroup <- 1
refgroup <- 0
predata <- list(contrgroup, refgroup)
names(predata) <- c("contrgroup", "refgroup")
names(predata$contrgroup) <- c("trt")
names(predata$refgroup) <- c("trt")

VEkout <- calcVEk(object = fitout, newdata = predata, CIlevel = 0.95)
summary(VEkout)
}

\keyword{per-challenge vaccine efficacy }
