#' @export
writeAlphaGenotypes = function(pop,file,chr=1,chips=rep(0,pop@nInd),simParam=SIMPARAM) {
  allSnps = numeric(0)
  uniqueChips = unique(chips)
  if (0 %in% uniqueChips){
    for (c in 1:simParam@nSnpChips) {
      allSnps = sort(union(allSnps,simParam@snpChips[[c]]@lociLoc))
    }
  }
  else {
    for (c in uniqueChips){
      allSnps = sort(union(allSnps,simParam@snpChips[[c]]@lociLoc))
    }
  }
  
  positions = list()
  for (i in 1:simParam@nSnpChips) {
    start = 1
    end = params@snpChips[[i]]@lociPerChr[1]
    if (chr > 1) {
      for (j in 2:chr) {
        start = start + params@snpChips[[i]]@lociPerChr[j-1]
        end = end + params@snpChips[[i]]@lociPerChr[j]
      }
    }
    positions[[i]] = params@snpChips[[i]]@lociLoc[start:end]
  }

  writeASGenotypes(pop@geno[[chr]],positions,allSnps,chips,sprintf("%s",newpop@id),'9',normalizePath(file, mustWork=FALSE))
}

#' @export
writeAlphaHaplotypes = function(pop,file,chr=1,chips=rep(0,pop@nInd),simParam=SIMPARAM) {
  allSnps = numeric(0)
  uniqueChips = unique(chips)
  if (0 %in% uniqueChips){
    for (c in 1:simParam@nSnpChips) {
      allSnps = sort(union(allSnps,simParam@snpChips[[c]]@lociLoc))
    }
  }
  else {
    for (c in uniqueChips){
      allSnps = sort(union(allSnps,simParam@snpChips[[c]]@lociLoc))
    }
  }

  positions = list()
  for (i in 1:simParam@nSnpChips) {
    start = 1
    end = params@snpChips[[i]]@lociPerChr[1]
    if (chr > 1) {
      for (j in 2:chr) {
        start = start + params@snpChips[[i]]@lociPerChr[j-1]
        end = end + params@snpChips[[i]]@lociPerChr[j]
      }
    }
    positions[[i]] = params@snpChips[[i]]@lociLoc[start:end]
  }

  writeASHaplotypes(pop@geno[[chr]],positions,allSnps,chips,sprintf("%s",newpop@id),'9',normalizePath(file, mustWork=FALSE))
}