#' @title Merge list of populations
#'
#' @description Rapidly merges a list of populations into a
#' single population
#'
#' @param popList a list containing \code{\link{Pop-class}} elements
#' or a \code{\link{MultiPop-class}}
#'
#' @return Returns a \code{\link{Pop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create a list of populations and merge list
#' pop = newPop(founderPop, simParam=SP)
#' popList = list(pop, pop)
#' pop2 = mergePops(popList)
#'
#' @export
mergePops = function(popList){
  if(is(popList,"MultiPop")){
    for(i in 1:length(popList@pops)){
      if(is(popList@pops[i],"MultiPop")){
        popList@pops[i] = mergePops(popList@pops[i])
      }
    }
    popList = popList@pops
  }

  classes = do.call("c",lapply(popList,
                               function(x) class(x)))
  if(any(classes=="NULL")){
    remove = which(classes=="NULL")
    popList = popList[-remove]
    classes = classes[-remove]
  }
  stopifnot(all(classes=="Pop"))
  #nChr
  nChr = do.call("c",lapply(popList,
                            function(x) x@nChr))
  stopifnot(all(nChr==nChr[1]))
  nChr = nChr[1]
  #ploidy
  ploidy = do.call("c",lapply(popList,
                              function(x) x@ploidy))
  stopifnot(all(ploidy==ploidy[1]))
  ploidy = ploidy[1]
  #nLoci
  nLoci = do.call("c",lapply(popList,
                             function(x){
                               all(x@nLoci==popList[[1]]@nLoci)
                             }))
  stopifnot(all(nLoci))
  nLoci = popList[[1]]@nLoci
  #id
  id = do.call("c",
               lapply(popList,
                      function(x) x@id))
  #iid
  iid = do.call("c",
                lapply(popList,
                       function(x) x@iid))
  #mother
  mother = do.call("c",
                   lapply(popList,
                          function(x) x@mother))
  #father
  father= do.call("c",
                  lapply(popList,
                         function(x) x@father))
  #fixEff
  fixEff= do.call("c",
                  lapply(popList,
                         function(x) x@fixEff))

  #misc
  misc = do.call("c",
                 lapply(popList,
                        function(x) x@misc))

  #sex
  sex = do.call("c",
                   lapply(popList,
                          function(x) x@sex))
  #nTraits
  nTraits = do.call("c",lapply(popList,
                               function(x) x@nTraits))
  stopifnot(all(nTraits==nTraits[1]))
  nTraits = nTraits[1]
  #nInd
  nInd = do.call("c",lapply(popList,
                            function(x) x@nInd))
  #gv
  gv = do.call("rbind",lapply(popList,
                              function(x) x@gv))
  #pheno
  pheno = do.call("rbind",lapply(popList,
                                 function(x) x@pheno))
  #ebv
  ebv = do.call("c",lapply(popList,
                           function(x) ncol(x@ebv)))
  if(all(ebv==ebv[1])){
    ebv = do.call("rbind",lapply(popList,
                                 function(x) x@ebv))
  }else{
    ebv = matrix(NA_real_,nrow=sum(nInd),ncol=0)
  }
  #gxe
  if(nTraits>=1){
    gxe = vector("list",length=nTraits)
    for(trait in 1:nTraits){
      if(!is.null(popList[[1]]@gxe[[trait]])){
        tmp = lapply(popList,function(x) x@gxe[[trait]])
        tmp = do.call("c",tmp)
        gxe[[trait]] = tmp
      }
    }
  }else{
    gxe = list()
  }
  #geno
  nBin = as.integer(nLoci%/%8L + (nLoci%%8L > 0L))
  geno = mergeMultGeno(popList,nInd=nInd,nBin=nBin,ploidy=ploidy)
  dim(geno) = NULL # Account for matrix bug in RcppArmadillo
  nInd = sum(nInd)
  return(new("Pop",
             nInd=nInd,
             nChr=nChr,
             ploidy=ploidy,
             nLoci=nLoci,
             sex=sex,
             geno=geno,
             id=id,
             iid=iid,
             mother=mother,
             father=father,
             fixEff=fixEff,
             nTraits=nTraits,
             gv=gv,
             gxe=gxe,
             pheno=pheno,
             ebv=ebv,
             misc=misc,
             miscPop=list()))
}
