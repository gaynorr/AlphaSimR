% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hybrids.R
\name{setPhenoGCA}
\alias{setPhenoGCA}
\title{Set GCA as phenotype}
\usage{
setPhenoGCA(
  pop,
  testers,
  use = "pheno",
  varE = NULL,
  reps = 1,
  fixEff = 1L,
  p = NULL,
  inbred = FALSE,
  onlyPheno = FALSE,
  simParam = NULL
)
}
\arguments{
\item{pop}{an object of \code{\link{Pop-class}}}

\item{testers}{an object of \code{\link{Pop-class}}}

\item{use}{true genetic value (\code{gv}) or phenotypes (\code{pheno}, default)}

\item{varE}{error variances for phenotype if \code{use="pheno"}. A vector
of length nTraits for independent error or a square matrix of 
dimensions nTraits for correlated errors.}

\item{reps}{number of replications for phenotype. See details.}

\item{fixEff}{fixed effect to assign to the population. Used 
by genomic selection models only.}

\item{p}{the p-value for the environmental covariate 
used by GxE traits. If NULL, a value is
sampled at random.}

\item{inbred}{are both pop and testers fully inbred. They are only 
fully inbred if created by \code{\link{newPop}} using inbred founders 
or by the \code{\link{makeDH}} function}

\item{onlyPheno}{should only the phenotype be returned, see return}

\item{simParam}{an object of \code{\link{SimParam}}}
}
\value{
Returns an object of \code{\link{Pop-class}} or 
a matrix if onlyPheno=TRUE
}
\description{
Calculates general combining ability from a set of testers and 
returns these values as phenotypes for a population.
}
\details{
The reps parameter is for convenient representation of replicated data. 
It was intended for representation of replicated yield trials in plant 
breeding programs. In this case, varE is set to the plot error and 
reps is set to the number plots per entry. The resulting phenotype 
would reflect the mean of all replications.
}
\examples{
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10, inbred=TRUE)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)

#Set phenotype to average per
pop2 = setPhenoGCA(pop, pop, use="gv", inbred=TRUE, simParam=SP)

}
