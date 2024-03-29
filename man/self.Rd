% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossing.R
\name{self}
\alias{self}
\title{Self individuals}
\usage{
self(pop, nProgeny = 1, parents = NULL, keepParents = TRUE, simParam = NULL)
}
\arguments{
\item{pop}{an object of \code{\link{Pop-class}}}

\item{nProgeny}{total number of selfed progeny per individual}

\item{parents}{an optional vector of indices for allowable parents}

\item{keepParents}{should previous parents be used for mother and
father.}

\item{simParam}{an object of \code{\link{SimParam}}}
}
\value{
Returns an object of \code{\link{Pop-class}}
}
\description{
Creates selfed progeny from each individual in a
population. Only works when sexes is "no".
}
\examples{
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
\dontshow{SP$nThreads = 1L}

#Create population
pop = newPop(founderPop, simParam=SP)

#Self pollinate each individual
pop2 = self(pop, simParam=SP)

}
