# fmt: skip file

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
#' pop@misc$tmp = rnorm(n=10)
#' pop@misc$tmp2 = rnorm(n=10)
#'
#' popList = list(pop, pop)
#' pop2 = mergePops(popList)
#'
#' @export
mergePops = function(popList){
  if(is(popList,"MultiPop")){
    for(i in seq_len(length(popList@pops))){
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
  tmp = sapply(popList, function(x) length(x@misc))
  if(!all(tmp == tmp[1])) {
    warning("number of misc elements differs - setting misc to an empty list!")
    misc = list()
  } else {
    if(tmp[1]>0) {
      tmp = lapply(popList, function(x) names(x@misc))
      allMatch = TRUE
      if(length(tmp)>1){
        for(i in 2:length(tmp)){
          if(!all(tmp[[1]]==tmp[[i]])){
            allMatch = FALSE
            break
          }
        }
      }
      if(allMatch){
        misc = vector("list", length=length(tmp[[1]]))
        for(i in seq_len(length(tmp[[1]]))){
          miscTmp = lapply(popList, function(x) x@misc[[i]])
          if (is.matrix(miscTmp[[1]])) {
            misc[[i]] = do.call("rbind", miscTmp)
          } else {
            misc[[i]] = do.call("c", miscTmp)
          }
        }
        names(misc) = tmp[[1]]
      }else{
        warning("misc element names do not match - setting misc to an empty list!")
        misc = list()
      }
    } else {
      misc = list()
    }
  }

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
    for(trait in seq_len(nTraits)){
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

  #wrap it all up into a Pop
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
             misc=misc,
             miscPop=list(),
             nTraits=nTraits,
             gv=gv,
             gxe=gxe,
             pheno=pheno,
             ebv=ebv))
}
