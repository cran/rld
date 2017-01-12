\name{SampleData}
\alias{SampleData}
\docType{data}

\title{
Simulated Sample Data for Illustration
}

\description{
This is a simulated data set containing 3 variables:
}

\usage{data(SampleData)}

\format{
  A data frame with 50 observations on the following 3 variables.
  \describe{
    \item{\code{time}}{a numeric vector: failure or censoring times, whichever comes first}
    \item{\code{delta}}{a numeric vector: censoring status}
    \item{\code{trt}}{a numeric vector: vaccination status}
  }
}

\details{
The dataset is artificial and used to test out some of the features of rld.
}

\examples{
data(SampleData)
}

\keyword{datasets}
