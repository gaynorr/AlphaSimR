% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/founderPop.R
\name{newMapPop}
\alias{newMapPop}
\title{New MapPop}
\usage{
newMapPop(genMap, haplotypes, inbred = FALSE, ploidy = 2L)
}
\arguments{
\item{genMap}{a list of genetic maps}

\item{haplotypes}{a list of matrices or data.frames that
can be coerced to matrices. See details.}

\item{inbred}{are individuals fully inbred}

\item{ploidy}{ploidy level of the organism}
}
\value{
an object of \code{\link{MapPop-class}}
}
\description{
Creates a new \code{\link{MapPop-class}} from user supplied
genetic maps and haplotypes.
}
\details{
Each item of genMap must be a vector of ordered genetic lengths in
Morgans. The first value must be zero. The length of the vector
determines the number of segregating sites on the chromosome.

Each item of haplotypes must be coercible to a matrix. The columns
of this matrix correspond to segregating sites. The number of rows
must match the number of individuals times the ploidy if using
inbred=FALSE. If using inbred=TRUE, the number of rows must equal
the number of individuals. The haplotypes can be stored as numeric,
integer or raw. The underlying C++ function will use raw.
}
\examples{
# Create genetic map for two chromosomes, each 1 Morgan long
# Each chromosome contains 11 equally spaced segregating sites
genMap = list(seq(0,1,length.out=11),
               seq(0,1,length.out=11))

# Create haplotypes for 10 outbred individuals
chr1 = sample(x=0:1,size=20*11,replace=TRUE)
chr1 = matrix(chr1,nrow=20,ncol=11)
chr2 = sample(x=0:1,size=20*11,replace=TRUE)
chr2 = matrix(chr2,nrow=20,ncol=11)
haplotypes = list(chr1,chr2)

founderPop = newMapPop(genMap=genMap, haplotypes=haplotypes)

}
