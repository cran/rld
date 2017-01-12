\name{rld-package}
\alias{rld-package}
\alias{spBayesSurv}
\docType{package}

\title{
Designing and Analyzing Repeated Low-Dose Challenge Experiments
}

\description{
This package can analyze the data from repeated low-dose (RLD) challenge experiments for evaluating effect of treatment(vaccine). Based on the fit of a discrete-time survival model with gamma distributed random effect and a complementary log-log link function, this package can provide estimates of regression parameters for effects of treatment and challenge dose, as well as estimates of treatment(vaccine) efficacy. In addition, this package can test regression parameters and frailty variance parameter using likelihood ratio test. Also it can help users design the RLD study by performing power analysis.
}

\details{
\tabular{ll}{
Package: \tab rld\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-01-11\cr
License: \tab GPL (>= 2)\cr
}
This package is created based on a paper written by Kang et al. (2015). In this paper, the authors proposd to use a discrete-time survival model with random effects to model the data from repeated low-dose challenge experiments.
}

\author{
Bin Yao \email{byao@fredhutch.org}, Ying Huang \email{yhuang@fhcrc.org}
and Chaeryon Kang \email{crkang@pitt.edu}
}

\references{
Kang, C., Huang, Y., and Miller, C. (2015). A discrete-time survival model with random effects for designing and analyzing repeated low-dose challenge experiments. \emph{Biostatistics}, 16(2): 295-310.

Yao, B. and Huang, Y. (2016). rld: An R package for designing and analyzing repeated low-dose experiments.
}

\keyword{package, discrete-time survival model}
