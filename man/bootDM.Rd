% Generated by roxygen2 (4.0.1): do not edit by hand
\name{bootDM}
\alias{bootDM}
\title{Take a bootstrap sample from a square distance matrix.}
\usage{
bootDM(dm, subject_group_vector)
}
\arguments{
\item{dm}{a square distance matrix}

\item{subject_group_vector}{number of subjects in each group to sample, as a vector.}
}
\value{
A square distance matrix.
}
\description{
Random sample with replacement from a distance matrix, with resulting matrix specified by number of subjects per group.
}
\examples{
bootDM(calcUJstudy(simStudy()),c(3,4,5))
}
\seealso{
\code{\link{calcUJstudy}}, \code{\link{calcWJstudy}}
}

