#' @title Write AlphaSuite style genotypes
#' 
#' @description
#' To Do
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param file name of output file
#' @param chr which chromosome
#' @param chips ??
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
writeAlphaGenotypes = function(pop,file,chr=1,chips=rep(0,pop@nInd),
                               simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  allSnps = numeric(0)
  uniqueChips = unique(chips)
  if (0 %in% uniqueChips){
    for (c in 1:simParam@nSnpChips) {
      allSnps = sort(union(allSnps,simParam@snpChips[[c]]@lociLoc))
    }
  }
  else {
    for (i in uniqueChips){
      allSnps = sort(union(allSnps,simParam@snpChips[[i]]@lociLoc))
    }
  }
  positions = list()
  for (i in 1:simParam@nSnpChips) {
    start = 1
    end = simParam@snpChips[[i]]@lociPerChr[1]
    if (chr > 1) {
      for (j in 2:chr) {
        start = start + simParam@snpChips[[i]]@lociPerChr[j-1]
        end = end + simParam@snpChips[[i]]@lociPerChr[j]
      }
    }
    positions[[i]] = simParam@snpChips[[i]]@lociLoc[start:end]
  }

  writeASGenotypes(pop@geno[[chr]],positions,allSnps,chips,pop@id,'9',normalizePath(file, mustWork=FALSE))
}

#' @title Write AlphaSuite style haplotypes
#' 
#' @description
#' To Do
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param file name of output file
#' @param chr which chromosome
#' @param chips ??
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
writeAlphaHaplotypes = function(pop,file,chr=1,chips=rep(0,pop@nInd),simParam) {
  allSnps = numeric(0)
  uniqueChips = unique(chips)
  if (0 %in% uniqueChips){
    for (c in 1:simParam@nSnpChips) {
      allSnps = sort(union(allSnps,simParam@snpChips[[c]]@lociLoc))
    }
  }
  else {
    for (i in uniqueChips){
      allSnps = sort(union(allSnps,simParam@snpChips[[i]]@lociLoc))
    }
  }

  positions = list()
  for (i in 1:simParam@nSnpChips) {
    start = 1
    end = simParam@snpChips[[i]]@lociPerChr[1]
    if (chr > 1) {
      for (j in 2:chr) {
        start = start + simParam@snpChips[[i]]@lociPerChr[j-1]
        end = end + simParam@snpChips[[i]]@lociPerChr[j]
      }
    }
    positions[[i]] = simParam@snpChips[[i]]@lociLoc[start:end]
  }

  writeASHaplotypes(pop@geno[[chr]],positions,allSnps,chips,sprintf("%s",pop@id),'9',normalizePath(file, mustWork=FALSE))
}
