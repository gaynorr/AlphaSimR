% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popSummary.R
\name{genicVarAA}
\alias{genicVarAA}
\title{Additive-by-additive genic variance}
\usage{
genicVarAA(pop, simParam = NULL)
}
\arguments{
\item{pop}{an object of \code{\link{Pop-class}}}

\item{simParam}{an object of \code{\link{SimParam}}}
}
\description{
Returns additive-by-additive epistatic
genic variance for all traits
}
\examples{
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)
SP$setVarE(h2=0.5)
\dontshow{SP$nThreads = 1L}

#Create population
pop = newPop(founderPop, simParam=SP)
genicVarAA(pop, simParam=SP)

}
