\name{transdata}
\alias{transdata}

\title{
Transfer Input Data Frame to Model Fitting Data Frame
}

\description{
This function transfers the original input data frame to required data frame for model fitting.
}

\usage{
transdata(data, ndlevel, nexposure)
}

\arguments{
  \item{data}{a data frame which must include variables named "time" and "delta".}

  \item{ndlevel}{a value: number of dose levels.}

  \item{nexposure}{a vector: number of challenges or exposures for each dose level.}
}

\details{
The original data frame from users include the time points when the subjects are right-censored or failed, the censoring status, vaccination status and baseline information. In addition, users must provide dose information including dose level and number of challenges or exposures for each dose level. However, the original data frame cannot be used directly to create design matrix because dose levels are time-dependent. Therefore, this function lists all time points before or at the right-censoring time or failure time for each subject. Then for each subject the function replicates each dose level multiple times according to the number of challenges the subject receives at each dose level. Finally, the function replicates the vaccination status and baseline information and returns a whole data frame.
}

\value{
a data frame containing variables "id", "time", "delta", "dose", vaccination status and baseline information.
}

\author{
Bin Yao, Ying Huang and Chaeryon Kang
}

\note{
The original data frame must have "time" and "delta" variables.
}

\examples{
data(SampleData)
newdata <- transdata(data = SampleData, ndlevel = 3, nexposure = c(10, 10, 2))
}

\keyword{transfer data}
