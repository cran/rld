\name{lrtest}
\alias{lrtest}
\alias{print.lrtest}

\title{
Likelihood Ratio Test
}

\description{
This function performs likelihood ratio test (LRT) to test regression paramters in the mean model and frailty variance paramter nu.
}

\usage{
lrtest(model1, model2, TestNu=TRUE, Siglevel=0.05)
}

\arguments{
  \item{model1}{a result of the nested model which is returned by \code{rld}.}

  \item{model2}{a result of the full model which is returned by \code{rld}.}

  \item{TestNu}{logic value: If TRUE, the frailty variance paramter will be tested. Otherwise, the regression
  parameters will be tested. The default is TRUE.}

  \item{Siglevel}{a value: significance level. The default is 0.05.}
}

\details{
Under the null hypothesis of frailty variance parameter \code{nu=0}, the test statistic converges to a mixture of chi-squared distribution. For regression parameters, under the null hypothesis \code{beta=0}, the test statistic converges to chi-squared distribution. Please refer to Kang et al.(2015) for more details.
}

\value{
  \item{statistic}{a value: test statistic.}

  \item{df}{a value: degree of freedom.}

  \item{pvalue}{the p-value.}
}

\references{
Kang, C., Huang, Y., and Miller, C. (2015). A discrete-time survival model with random effects for designing and analyzing repeated low-dose challenge experiments. \emph{Biostatistics}, 16(2): 295-310.
}

\author{
Bin Yao, Ying Huang and Chaeryon Kang
}

\seealso{
\code{\link{rld}}
}

\examples{
#test frailty variance parameter
data(SampleData)
newdata <- transdata(data = SampleData, ndlevel = 3, nexposure = c(10, 10, 2))
fitout1 <- rld(formula = survival::Surv(time, delta)~factor(trt)+factor(dose), data = newdata,
               frailty = FALSE)
fitout2 <- rld(formula = survival::Surv(time, delta)~factor(trt)+factor(dose), data = newdata,
               frailty = TRUE)
testnu <- lrtest(model1 = fitout1, model2 = fitout2, TestNu=TRUE)

\dontrun{
#test regression parameters
fitout3 <- rld(formula = survival::Surv(time, delta)~factor(dose), data = newdata,
               frailty = TRUE)
fitout4 <- rld(formula = survival::Surv(time, delta)~factor(trt)+factor(dose), data = newdata,
               frailty = TRUE)
testbeta <- lrtest(model1 = fitout3, model2 = fitout4, TestNu=FALSE)
}
}

\keyword{LRT}
