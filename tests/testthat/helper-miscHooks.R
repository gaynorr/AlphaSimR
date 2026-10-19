# Shared populations with miscellaneous-data hooks for integration tests.
miscFunctionFixture = function(){
  founder = quickHaplo(nInd=4, nChr=1, segSites=10)
  sp = SimParam$new(founder)
  sp$nThreads = 1L
  sp$subsetMiscPop = function(pop, i, simParam){
    list(total=sum(pop@misc$value[i]), ids=pop@id[i])
  }
  sp$mergeMiscPop = function(popList, simParam){
    list(total=sum(vapply(popList, function(x) sum(x@misc[[1]]), numeric(1))),
         ids=unlist(lapply(popList, function(x) x@id), use.names=FALSE))
  }
  pop = newPop(founder, simParam=sp)
  pop@misc = list(value=c(2, 4, 6, 8), mat=matrix(1:8, nrow=4))
  pop@miscPop = list(total=20, ids=pop@id)
  list(pop=pop, sp=sp, founder=founder)
}
