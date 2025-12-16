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
    warning("Populations have different numbers of EBV columns; EBVs removed!")
    ebv = matrix(NA_real_,nrow=sum(nInd),ncol=0,
                 dimnames=list(NULL, NULL))
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

#' @title Flatten a MultiPop object to a specified depth
#'
#' @description
#' Recursively flatten a \code{\link{MultiPop-class}} object into a shallower
#' \code{MultiPop} containing only \code{\link{Pop-class}} objects
#' down to a requested level (see Details).
#'
#' @param x A \code{\link{Pop-class}} or \code{\link{MultiPop-class}} object.
#' @param level Integer scalar >= 1. Number of \code{MultiPop} levels to
#'   preserve.
#'
#' @details
#' The \code{level} argument controls how many levels of nesting are preserved.
#' \code{level = 1}: flatten structure of \code{MultiPop-class} so that
#'   all items are \code{Pop-class} objects (that is, flatten everything
#'   below the top level).
#' \code{level > 1}: preserve the top \code{level} levels of nesting; any
#'   deeper \code{MultiPop-class} objects are flattened.
#'
#' If \code{level} is greater than or equal to the nesting depth, the original
#' object is returned unchanged.
#'
#' @return If \code{x} is a \code{\link{Pop-class}}, the same \code{x}
#' object is returned. Otherwise a \code{\link{MultiPop-class}} is returned
#' whose \code{x@pops} slot contains \code{\link{Pop-class}} (and possibly
#' \code{\link{MultiPop-class}}) objects flattened according to \code{level}.
#'
#' @seealso \code{\link{mergeMultiPops}} and \code{\link{mergePops}}
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=12, nChr=1, segSites=10)
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' # Create a multi-population with down to level 3 nesting
#' mp_nested = newMultiPop(pop[1:2],
#'                         newMultiPop(pop[3:4],
#'                                     newMultiPop(pop[5:7], pop[8:12])))
#' mp_nested
#'
#' # Completely flatten to a single top-level MultiPop
#' flattenMultiPop(mp_nested)
#'
#' # Preserve two levels of nesting
#' flattenMultiPop(mp_nested, level=2)
#'
#' @export
flattenMultiPop = function(x, level=1) {
  if (isPop(x)) return(x)
  stopifnot(isMultiPop(x))
  multi = which(sapply(x@pops, isMultiPop))
  while (level > 1) {
    level = level - 1
    for (i in multi) {
      x@pops[[i]] = flattenMultiPop(x@pops[[i]], level = level)
    }
    return(do.call(newMultiPop, x@pops))
  }
  flatPopList = .flattenMultiPop(x)
  return(do.call(newMultiPop, flatPopList))
}

#' Helper function to recursively extract Pop objects from a MultiPop
#'
#' @param mp \code{\link{MultiPop-class}} object
#'
#' @keywords internal
.flattenMultiPop = function(mp) {
  popList = list()
  for (item in mp@pops) {
    if (isPop(item)) {
      popList = c(popList, list(item))
    } else if (isMultiPop(item)) {
      popList = c(popList, .flattenMultiPop(item))
    }
  }
  return(popList)
}

#' @title Merge Pop and MultiPop objects
#'
#' @description
#' Merge one or more \code{\link{Pop-class}} and \code{\link{MultiPop-class}}
#' objects. Because a \code{MultiPop} can have a nested structure
#' merging is controlled by the \code{level} argument (see Details).
#'
#' @param ... \code{\link{Pop-class}} or \code{\link{MultiPop-class}} objects;
#'   \code{NULL} values are ignored.
#' @param level Integer scalar >= 0 to merge at a sepecific level of nesting;
#'   see Details.
#'
#' @details
#' The function accepts multiple inputs and merges them according to the
#' \code{level} argument.
#' \code{level = 0}: merge all \code{Pop} objects in the inputs into a single
#'   \code{Pop} object (the inputs are first completely flattened then merged).
#' \code{level = 1}: merge inputs into a \code{MultiPop} so that each item is a
#'   \code{Pop} (level 1). Each input is first flattened to level 1 and its
#'   \code{Pop} objects merged into one \code{Pop}, and then all these \code{Pop}
#'   objects are merged into a single \code{MultiPop}
#' \code{level > 1}: merge inputs into a \code{MultiPop} while preserving top
#'   \code{level} structure. Each input is first flattened to the requested
#'   \code{level} and its items merged into a single \code{MultiPop}, and then
#'   inputs are merged.
#'
#' Important behavioral notes: If a single \code{Pop} is provided, it
#' is returned unchanged. If a single \code{MultiPop} is provided, it
#' is processed according to \code{level}. If multiple objects are provided,
#' they are combined into a \code{MultiPop} (preserving the order of inputs)
#' and then processed according to \code{level}. If neither of the multiple inputs
#' is a \code{MultiPop}, an error is raised. An error is also raised if inputs
#' are not a \code{Pop} or a \code{MultiPop}. \code{NULL} inputs are ignored.
#'
#' @seealso \code{\link{MultiPop-class}}, \code{\link{flattenMultiPop}}, and
#'   \code{\link{mergePops}}
#'
#' @return If \code{level == 0}, or when merging yields a single \code{Pop},
#' a \code{\link{Pop-class}} object is returned. Otherwise a
#' \code{\link{MultiPop-class}} object is returned with the requested level of
#' nesting preserved.
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=11, nChr=1, segSites=10)
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' # Create two multi-population with different levels of nesting
#' mp1 = newMultiPop(pop[1:2], pop[3:5])
#' mp1
#' mp2 = newMultiPop(pop[6:7],
#'                   newMultiPop(pop[8:10], pop[11]))
#' mp2
#'
#' # Fully merge all Pops in inputs into one Pop
#' mergeMultiPops(mp1[[1]]) # nothing happens with a single Pop
#' mergeMultiPops(mp1)
#' mergeMultiPops(mp1, mp2)
#'
#' # Merge into a MultiPop where each level 1 item is a Pop
#' mergeMultiPops(mp1[[1]], level=1) # nothing happens with a single Pop
#' mergeMultiPops(mp1, level=1)
#' mergeMultiPops(mp1, mp2[[1]], level=1)
#' mergeMultiPops(mp2, level=1)
#' mergeMultiPops(mp1, mp2, level=1)
#'
#' # Merge into a MultiPop and preserve level 1 and 2 structure
#' mergeMultiPops(mp2, level=2)
#' mergeMultiPops(mp2, mp1, level=2)
#' mergeMultiPops(mp2, mp1[[1]], level=2)
#'
#' @export
mergeMultiPops = function(..., level=0){

  popList = list(...)
  classes = do.call("c", lapply(popList, class))

  if(any(classes == "NULL")){
    remove = which(classes == "NULL")
    popList = popList[-remove]
    classes = classes[-remove]
  }

  # If popList contains a single object
  if (length(classes) == 1) {
    if (classes == "Pop") {
      # If the object is Pop, return it without the list wrapping
      return(popList[[1]])
    } else if (classes == "MultiPop") {
      # Else, remove the list wrapping from the single multiPop object
      multiPop = popList[[1]]
      popList = multiPop@pops
    } else {
      stop("One or more objects are not of Pop or Multi-Pop class!")
    }
  } else if (all(classes == "Pop")) {
    stop("Use mergePops() to merge multiple Pop objects!")
  } else {
    # Combine all arguments into a single MultiPop
    multiPop = do.call('c', popList)
    popList = multiPop@pops
  }

  multi = which(sapply(popList, isMultiPop))
  while (level > 0) {
    level = level - 1
    for (i in multi) {
      popList[[i]] = mergeMultiPops(popList[[i]], level = level)
    }
    return(do.call(newMultiPop, popList))
  }

  flatMultiPop = flattenMultiPop(multiPop)
  return(mergePops(flatMultiPop))
}
