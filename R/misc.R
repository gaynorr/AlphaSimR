#' @title Selection Intensity
#' 
#' @description 
#' Calculates the standardized selection intensity
#' 
#' @param p the proportion of individuals selected
#' 
#' @export
selInt = function(p){
  return(dnorm(qnorm(1-p))/p)
}

#' @title Calculate Smith-Hazel Weights
#' 
#' @description
#' Calculates weights for Smith-Hazel index given economice weights 
#' and phenotypic and genotypic variance-covariance matrices.
#' 
#' @param econWt vector of economic weights
#' @param varG the genetic variance-covariance matrix
#' @param varP the phenotypic variance-covariance matrix
#' 
#' @return a vector of weight for calculating index values
#' 
#' @export
smithHazel = function(econWt,varG,varP){
  return(solve(varP)%*%varG%*%econWt)
}

#' @title Selection Index
#' 
#' @description
#' Calculates values of a selection index given trait values and 
#' weights. This function is intended to be used in combination with 
#' selection functions working on populations such as 
#' \code{\link{selectInd}}.
#' 
#' @param Y a matrix of trait values
#' @param b a vector of weights
#' 
#' @export
selIndex = function(Y,b){
  return(Y%*%b)
}

#' @title Edit Genome
#' 
#' @description
#' Edits selected loci of selected individuals to a homozygous 
#' state for either the 1 or 0 allele. The gv slot is recalculated to 
#' reflect the any changes due to editing, but other slots remain the same.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param ind a vector of individuals to edit
#' @param chr a vector of chromosomes to edit. Length must match 
#' length of segSites.
#' @param segSites a vector of segregating sites to edit. Length must 
#' match length of chr.
#' @param allele either 0 or 1 for desired allele
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
editGenome = function(pop,ind,chr,segSites,allele,simParam=SIMPARAM){
  ind = unique(as.integer(ind))
  stopifnot(all(ind%in%(1:pop@nInd)))
  chr = as.integer(chr)
  segSites = as.integer(segSites)
  stopifnot(length(chr)==length(segSites))
  allele = as.integer(allele)
  stopifnot(allele==0L | allele==1L)
  allele = as.raw(allele)
  for(selChr in unique(chr)){
    selSegSites = segSites[chr==selChr]
    pop@geno[[selChr]][selSegSites,,ind] = allele
  }
  gv = lapply(simParam@traits,getGv,pop=pop,w=0.5)
  gv = do.call("cbind",gv)
  pop@gv = gv
  validObject(pop)
  return(pop)
}


