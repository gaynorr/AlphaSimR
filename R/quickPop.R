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
               geno=as.matrix(geno),genMaps=as.matrix(genMaps))
  return(output)
}

#' @title Haplotype tracking population
#' 
#' @description
#' Creates a population for tracking haplotypes.
#'
#' @param genMaps a list of genetic maps
#' @param nInd number of individuals
#' @param inbred should individuals be fully inbred
#' 
#' @details
#' The number of chromosomes is determined by the length of genMaps. 
#' Each item of genMaps must be a vector of ordered numeric values and 
#' the first value must be zero. The length of the vector determines the 
#' number of segregating sites on the chromosome.
#' 
#' If inbred=FALSE, the value of nInd must be less than or equal to 
#' 128. Otherwise, it must be less than or equal to 256.
#' 
#' @examples
#' genMaps = list(seq(0,1,length.out=11))
#' FOUNDERPOP = trackHaploPop(genMaps=genMaps,nInd=10)
#' SIMPARAM = createSimulation(founderPop=FOUNDERPOP)
#' pop = newPop(FOUNDERPOP,simParam=SIMPARAM)
#' segSites = pullSegSiteHaplo(pop,simParam=SIMPARAM)
#' 
#' @export
trackHaploPop = function(genMaps,nInd,inbred=FALSE){
  stopifnot(is.list(genMaps))
  if(inbred){
    stopifnot(nInd<=128)
  }else{
    stopifnot(nInd<=256)
  }
  nInd = as.integer(nInd)
  nChr = length(genMaps)
  nLoci = unlist(lapply(genMaps,length))
  geno = vector("list",nChr)
  for(i in 1:nChr){
    tmpGeno = as.raw(0:(2*nInd-1))
    tmpGeno = array(raw(),dim=c(nLoci[i],2L,nInd))
    tmp=-1
    for(j in 1:nInd){
      if(inbred){
        tmp=tmp+1
        tmpGeno[,1:2,j] = as.raw(tmp)
      }else{
        for(k in 1:2){
          tmp=tmp+1
          tmpGeno[,k,j] = as.raw(tmp)
        } 
      }
    }
    geno[[i]] = tmpGeno
  }
  output = new("MapPop",nInd=nInd,nChr=nChr,ploidy=2L,
               nLoci=nLoci,geno=as.matrix(geno),
               genMaps=as.matrix(genMaps))
  return(output)
}