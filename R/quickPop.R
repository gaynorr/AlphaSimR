#' @title Quick founder genotypes
#'
#' @description 
#' Creates founder genotypes with evely spaced segregating 
#' sites and allele frequencies of 0.5.
#' 
#' @param nInd number of individuals
#' @param nChr number of chromosomes
#' @param segSites number of segregating sites per chromosome
#' @param genLen genetic length of chromosomes in Morgans
#' @param inbred should founder individuals be inbred
#' 
#' @return an object of \code{\link{MapPop-class}}
#' 
#' @export
quickPop = function(nInd,nChr,segSites,genLen=1,inbred=TRUE){
  geno = vector("list",nChr)
  genMaps = vector("list",nChr)
  for(i in 1:nChr){
    genMaps[[i]] = seq(0,genLen,length.out=segSites)
    geno[[i]] = array(sample(as.raw(c(0,1)),nInd*segSites*2,
                                    replace=TRUE),
                      dim=c(segSites,2,nInd))
    if(inbred){
      geno[[i]][,2,] = geno[[i]][,1,]
    }
  }
  output = new("MapPop",nInd=as.integer(nInd),nChr=as.integer(nChr),
               ploidy=2L,nLoci=as.integer(rep(segSites,nChr)),
               gender=rep("H",nInd),geno=as.matrix(geno),genMaps=as.matrix(genMaps))
  return(output)
}