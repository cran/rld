\name{rld.design}
\alias{rld.design}
\alias{print.rld.design}

\title{
Power analysis for the Repeated-Low Dose Experiments
}

\description{
This function performs power analyis to design a repeated low-dose challenge experiment with a vaccine and a placebo
arm.
}

\usage{
rld.design(nsim, nv, np, ndlevel, nexposure, rho, p0, RR,
           method=c("LRT", "log-rank"), Siglevel)
}

\arguments{
  \item{nsim}{a value indicating the number of simulations to run.}

  \item{nv}{a value indicating the number of subjects in vaccine group.}

  \item{np}{a value indicating the number of subjects in placebo group.}

  \item{ndlevel}{a value indicating the number of dose levels.}

  \item{nexposure}{a vector of challenges or exposures for all dose levels.}

  \item{rho}{a value: intracluster correlation between underlying continuous responses.}

  \item{p0}{probability of infection in placebo group at each challenge dose level.}

  \item{RR}{a value: relative risk of vaccine relative to placebo at each challenge dose level.}

  \item{method}{"LRT": likelihood ratio test; "log-rank": log-rank test.}

  \item{Siglevel}{a value indicating significance level.}
}

\details{
Users need to specify the paramters of the experiment. The function will generate the data from the discrete-time survival model with random effects. The power is defined as the proportion of rejecting the null hypothesis that treatment has no effect. There are only two groups in the study, i.e. vaccine group and placebo group. There are two types of test availabe for use, likelihood ratio test and log-rank test. Note that likelihood ratio test takes more simulation time than log-rank test because of model fitting.
}

\value{
  \item{method}{a character which is either a "LRT" or "log-rank".}

  \item{power}{a value: statistcal power.}
}

\references{
Yao,B and Huang, Y. (2016+). rld: An R package for designing and analyzing repeated low-dose experiments.
}

\author{
Bin Yao, Ying Huang and Chaeryon Kang
}

\examples{
designout <- rld.design(nsim = 50, nv = 25, np = 25, ndlevel = 3, nexposure = c(10, 10, 2),
           rho = 0.2, p0 = c(0.16, 0.22, 0.27), RR = c(0.3, 0.45, 0.55),
           method = "log-rank", Siglevel = 0.05)
}

\keyword{design}
