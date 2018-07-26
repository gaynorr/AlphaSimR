#' @title Merge list of populations
#' 
#' @description Rapidly merges a list of populations into a
#' single population
#' 
#' @param popList a list containing \code{\link{Pop-class}} elements
#' 
#' @return Returns a \code{\link{Pop-class}}
#' 
#' @export
mergePops = function(popList){
  classes = do.call("c",lapply(popList,
                               function(x) class(x)))
  stopifnot(all(classes=="Pop"))
  #nChr
  nChr = do.call("c",lapply(popList,
                            function(x) x@nChr))
  stopifnot(all(nChr==nChr[1]))
  nChr = nChr[1]
  #ploidy
  ploidy = do.call("c",lapply(popList,
                              function(x) x@ploidy))
  stopifnot(all(ploidy==2L))
  ploidy = 2L
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
  #reps
  reps= do.call("c",
                lapply(popList,
                       function(x) x@reps))
  #gender
  gender = do.call("c",
                   lapply(popList,
                          function(x) x@gender))
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
                                 function(x) x@pheno))
  }else{
    ebv = matrix(NA_real_,nrow=nInd,ncol=0)
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
  geno = mergeMultGeno(popList,nInd=nInd,nLoci=nLoci,ploidy=ploidy)
  nInd = sum(nInd)
  return(new("Pop",
             nInd=nInd,
             nChr=nChr,
             ploidy=ploidy,
             nLoci=nLoci,
             gender=gender,
             geno=geno,
             id=id,
             mother=mother,
             father=father,
             fixEff=fixEff,
             reps=reps,
             nTraits=nTraits,
             gv=gv,
             gxe=gxe,
             pheno=pheno,
             ebv=ebv))
}
