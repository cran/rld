\name{calcVEt}
\alias{calcVEt}
\alias{print.calcVEt}
\alias{summary.calcVEt}
\alias{print.summary.calcVEt}

\title{
Calculate Vaccine Efficacy Before or At the Time of Challenge t
}

\description{
This function estimates the vaccine efficacy before or at the time of challenge t. VE(t)>0 indicates that the vaccine is effective in reducing the risk of infection before or at time t, whereas VE(t)<=0 indicate that the vaccine is not effective or has a negative effect.
}

\usage{
calcVEt(object, nexposure, newdata, CIlevel = 0.95)
}

\arguments{
  \item{object}{a fitted object of class inheriting from "rld".}

  \item{nexposure}{a vector of challenges or exposures for all dose levels for predicting VE(t).}

 \item{newdata}{a data list for predicting vaccine efficacy where "contrgroup" and "refgroup" list names must be included.}

  \item{CIlevel}{the confidence level. The default is 0.95.}
}

\details{
Vaccine efficacy for preventing infection before or at the time of challenge t, VE(t), is defined as the relative reduction in the risk of infection before or at time t for the vaccine group compared to the placebo group. Please refer to Kang et al.(2015) for more details.
}

\value{
 \item{VE}{a vector containing vaccine efficacy estimates for contrast group and reference group.}

  \item{se}{a vector containing standard deviations of per-challenge vaccine efficacy estimates.}

  \item{lwr}{a vector containing lower bound of confidence interval for VE(t).}

  \item{upr}{a vector containing upper bound value of confidence interval for VE(t).}

  \item{time}{a vector containing challenge times.}
}

\references{
Kang, C., Huang, Y., and Miller, C. (2015). A discrete-time survival model with random effects for designing and analyzing repeated low-dose challenge experiments. \emph{Biostatistics}, 16(2): 295-310.
}

\author{
Bin Yao, Ying Huang and Chaeryon Kang
}

\seealso{
\code{\link{calcVEk}}, \code{\link{calcpk}}
}

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

VEtout <- calcVEt(object = fitout, nexposure = c(10, 10, 2), newdata = predata,
                  CIlevel = 0.95)
summary(VEtout)
}

\keyword{vaccine efficacy}
