% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popSummary.R
\name{genParam}
\alias{genParam}
\title{Sumarize genetic parameters}
\usage{
genParam(pop, indValues = FALSE, simParam = NULL)
}
\arguments{
\item{pop}{an object of \code{\link{Pop-class}}}

\item{indValues}{should breeding values, dominance deviations 
and allele subsitution effects be returned}

\item{simParam}{an object of \code{\link{SimParam}}}
}
\value{
\describe{
\item{varA}{an nTrait by nTrait matrix of additive genetic variances}
\item{varD}{an nTrait by nTrait matrix of dominance genetic variances}
\item{varG}{an nTrait by nTrait matrix of total genetic variances}
\item{genicVarA}{an nTrait vector of additive genic variances}
\item{genicVarD}{an nTrait vector of dominance genic variances}
\item{genicVarG}{an nTrait vector of total genic variances}
\item{bv}{an nInd by nTrait matrix of breeding values with dimensions nInd by nTraits}
\item{dd}{an nInd by nTrait matrix of dominance deviations with dimensions nInd by nTraits}
\item{alpha}{an nTrait list of allele subsitution effects}
}
}
\description{
Calculates genetic and genic additive and dominance variances 
for an object of \code{\link{Pop-class}}
}