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
      if(is(popList@pops[[i]],"MultiPop")){
        popList@pops[[i]] = mergePops(popList@pops[[i]])
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
#' @param preserveNames Character scalar controlling how to handle names when
#'   flattening. One of:  \cr
#'   - \code{"auto"} (default): Keep names only if all are present and
#'       unique after flattening; otherwise drop them.  \cr
#'   - \code{"concatenate"}: Prefix child names with parent names using
#'       underscore separator. Keep names only if resulting names are unique.  \cr
#'   - \code{"force"}: Always keep concatenated names, making them unique
#'       via \code{\link{make.unique}} if duplicates are found (with a warning).  \cr
#'   - \code{"none"}: Drop all names.
#'
#' @details
#' The \code{level} argument controls how many levels of nesting are preserved.  \cr
#' - \code{level = 1}: flatten structure of \code{MultiPop-class} so that
#'   all items are \code{Pop-class} objects (that is, flatten everything
#'   below the top level).  \cr
#' - \code{level > 1}: preserve the top \code{level} levels of nesting; any
#'   deeper \code{MultiPop-class} objects are flattened.
#'
#' If \code{level} is greater than or equal to the nesting depth, the original
#' object is returned unchanged.
#' 
#' The \code{preserveNames} argument allows control over name preservation during
#' flattening. This is useful when you want to track the source of flattened
#' populations through their hierarchical names. The default \code{"auto"} mode
#' is conservative: names are kept only when they are meaningful (all present
#' and unique). Use \code{"concatenate"} to always build hierarchical names
#' (when possible), \code{"force"} to guarantee names with uniqueness enforcement,
#' or \code{"none"} to explicitly discard all names.
#'
#' @return If \code{x} is a \code{\link{Pop-class}}, the same \code{x}
#' object is returned. Otherwise a \code{\link{MultiPop-class}} is returned
#' whose \code{x@pops} slot contains \code{\link{Pop-class}} (and possibly
#' \code{\link{MultiPop-class}}) objects flattened according to \code{level}
#' and \code{preserveNames}.
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
#' mp_nested = newMultiPop(pop1 = pop[1:2],
#'                         mp1 = newMultiPop(pop2 = pop[3:4],
#'                                           mp2 = newMultiPop(pop3 = pop[5:7], pop4 = pop[8:12])))
#' mp_nested
#'
#' # Completely flatten to a single top-level MultiPop
#' # With default "auto" mode: names are kept if unique
#' flattenMultiPop(mp_nested)
#'
#' # With "concatenate": build hierarchical names
#' flattenMultiPop(mp_nested, preserveNames = "concatenate")
#'
#' # With "force": guarantee unique names with suffixes
#' # (not needed here but useful if duplicates were present)
#' flattenMultiPop(mp_nested, preserveNames = "force")
#'
#' # With "none": drop all names
#' flattenMultiPop(mp_nested, preserveNames = "none")
#'
#' # Preserve two levels of nesting (level 2 names are kept)
#' flattenMultiPop(mp_nested, level = 2)
#'
#' @export
flattenMultiPop = function(x, level = 1,
                           preserveNames = c("auto", "concatenate", "force", "none")) {
  preserveNames = match.arg(preserveNames)
  if (isPop(x)) {
    return(x)
  }
  stopifnot(isMultiPop(x))

  if (level > 1L) {
    for (i in seq_along(x@pops)) {
      if (isMultiPop(x@pops[[i]])) {
        x@pops[[i]] = flattenMultiPop(
          x@pops[[i]],
          level = level - 1L,
          preserveNames = preserveNames
        )
      }
    }
    validObject(x)
    return(x)
  }
  flatPopList = .flattenMultiPop(x, preserveNames = preserveNames)
  multiPop = do.call(newMultiPop, flatPopList)
  validObject(multiPop)
  return(multiPop)
}

#' Helper function to recursively extract Pop objects from a MultiPop
#'
#' @param mp \code{\link{MultiPop-class}} object
#'
#' @keywords internal
.flattenMultiPop = function(mp, preserveNames = c("auto", "concatenate", "force", "none")) {
  preserveNames = match.arg(preserveNames)

  if (.depthMultiPop(mp) == 1L) {
    res = mp@pops
    if (preserveNames == "none") {
      names(res) = NULL
    }
    return(res)
  }

  nm = names(mp@pops)
  popList = list()
  nameVec = character()

  for (i in seq_along(mp@pops)) {
    item = mp@pops[[i]]
    label = if (!is.null(nm) && nzchar(nm[i])) nm[i] else as.character(i)

    if (isPop(item)) {
      popList = c(popList, list(item))
      nameVec = c(nameVec, label)
    } else if (isMultiPop(item)) {
      child = .flattenMultiPop(item, preserveNames = preserveNames)
      childNames = names(child)

      if (preserveNames %in% c("concatenate", "force")) {
        if (is.null(childNames)) {
          childNames = as.character(seq_along(child))
        }
        childNames = paste(label, childNames, sep = "_")
        names(child) = childNames
      }

      popList = c(popList, child)
      nameVec = c(
        nameVec,
        if (is.null(names(child))) {
          rep(NA_character_, length(child))
        } else {
          names(child)
        }
      )
    }
  }

  if (preserveNames == "none") {
    names(popList) = NULL
    return(popList)
  }

  # decide whether to keep/transform names
  if (preserveNames == "force") {
    # nm_final = nameVec
    nameVec[is.na(nameVec) | !nzchar(nameVec)] = as.character(which(
      is.na(nameVec) | !nzchar(nameVec)
    ))
    if (any(duplicated(nameVec))) {
      warning(
        "Duplicate names found in 'force' mode. Making names unique by appending suffixes.",
        call. = FALSE
      )
    }
    names(popList) = make.unique(nameVec)
    return(popList)
  }

  # auto or concatenate: keep only if no NA and no duplicates
  if (length(nameVec) > 0L && 
       !any(is.na(nameVec)) && 
       !any(duplicated(nameVec))
  ) {
    names(popList) = nameVec
  } else {
    names(popList) = NULL
  }

  popList
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
#' - \code{level = 0}: merge all \code{Pop} or \code{MultiPop} objects in the 
#'   inputs into a single \code{Pop} object (the inputs are first completely 
#'   flattened then merged).  \cr
#' - \code{level = 1}: merge inputs into a \code{MultiPop} so that each item is a
#'   \code{Pop} (level 1). Each input is first flattened to level 1 and its
#'   \code{Pop} objects merged into one \code{Pop}, and then all these \code{Pop}
#'   objects are merged into a single \code{MultiPop}.  \cr
#' - \code{level > 1}: merge inputs into a \code{MultiPop} while preserving top
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
  if (length(classes) == 1L) {
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
    multiPop = do.call(newMultiPop, popList)
    validObject(multiPop)
    return(multiPop)
  }

  flatMultiPop = flattenMultiPop(multiPop)
  return(mergePops(flatMultiPop))
}

#' @title Split Pop or MultiPop
#'
#' @description
#' Split a \code{\link{Pop-class}} or \code{\link{MultiPop-class}} object
#' into a \code{MultiPop} using one or more grouping specifications passed
#' through \code{by}.
#'
#' @param x a \code{\link{Pop-class}} or \code{\link{MultiPop-class}} object
#' @param by a vector, function, list, or data frame defining groupings
#' @param level a positive integer, a vector of positive integers, or
#'   \code{Inf}. Only relevant when \code{x} is a \code{\link{MultiPop-class}}
#'   object. See Details.
#'
#' @details
#' The \code{by} argument can take several forms:  \cr
#' - Atomic vectors: used directly as grouping labels. The
#'   vector is passed to \code{\link[base]{split}} and may be recycled as
#'   needed. If the vector length is not a multiple of the population size,
#'   \code{\link[base]{split}} issues a warning.  \cr
#' - Functions: called on each \code{Pop} object and must return
#'   an atomic vector of grouping labels.  \cr
#' - Lists: each element must be a vector or function. Elements
#'   are applied recursively to create nested \code{\link{MultiPop-class}}
#'   objects.  \cr
#' - Data frames: each column defines one nesting level. In this
#'   case, \code{x} must be a \code{\link{Pop-class}} object. Rows of the data
#'   frame must correspond to individuals in \code{x}. If row names are present
#'   and match \code{x@id}, they are used to align rows to individuals;
#'   otherwise rows are assumed to be in population order and a warning is
#'   issued. \code{NA} values are allowed in grouping columns and can be used
#'   to represent uneven nesting structures. If all grouping values for a
#'   subpopulation are \code{NA} at a given level, splitting stops for that
#'   branch.
#'
#' The \code{level} argument is only relevant when \code{x} is a
#' \code{\link{MultiPop-class}} object. If \code{level = Inf}, all
#' \code{\link{Pop-class}} objects contained in \code{x} are split according
#' to \code{by}. If \code{level = c(a, b)}, only \code{Pop-class} objects at
#' nesting levels \code{a} and \code{b} (with the top level being \code{1}) are
#' split. An error is raised if any requested level exceeds the maximum nesting
#' depth of \code{x}.
#'
#' @return Returns a \code{\link{MultiPop-class}} object.
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd = 10, nChr = 1, segSites = 10)
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#'
#' # Create population
#' pop = newPop(founderPop, simParam = SP)
#'
#' # Split pop into groups A/B deterministically
#' mp1 = splitPop(pop, by = rep(c("A", "B"), length.out = nInd(pop)))
#' mp1
#'
#' # Split using a data frame, with rows aligned by individual ID
#' by_df = data.frame(
#'   level1 = sample(LETTERS[1:3], nInd(pop), replace = TRUE),
#'   level2 = sample(1:2, nInd(pop), replace = TRUE),
#'   row.names = pop@id
#' )
#' splitPop(pop, by = by_df)
#'
#' # Uneven nested structure using NA values in lower levels
#' by_df2 = data.frame(
#'   level1 = c(rep("A", 4), rep("B", 6)),
#'   level2 = c(rep(NA_character_, 4), rep(c("C", "D"), each = 3)),
#'   row.names = pop@id
#' )
#' splitPop(pop, by = by_df2)
#'
#' # Nested split: first by a random grouping vector, then by family
#' mp2 = splitPop(
#'   pop,
#'   by = list(
#'     sample(LETTERS[1:3], nInd(pop), replace = TRUE),
#'     function(x) paste(x@mother, x@father, sep = "_")
#'   )
#' )
#' mp2
#'
#' # When x is a MultiPop, control which nesting levels are split
#' mp_nested = newMultiPop(pop[1:3], newMultiPop(pop[4:6], pop[7:9]))
#' splitPop(
#'   mp_nested,
#'   by = function(p) rep(c("A", "B"), length.out = nInd(p)),
#'   level = 1
#' )
#'
#' @export
splitPop = function(x, by, level = Inf) {
  if (!isPop(x) && !isMultiPop(x)) {
    stop("`x` must be a Pop or MultiPop object")
  }
  
  # Validate by argument
  if (!is.list(by)) {
    by = list(by)
  }
  if (length(by) == 0) {
    stop("`by` must have at least one grouping spec.")
  }
  if (is.data.frame(by) && isPop(x)){
    if (.row_names_info(by) > 0) {
      if (setequal(rownames(by), x@id)){
        by = by[x@id, , drop = FALSE]
      } else {
        warning("Row names of data frame `by` don't match `x@id`. Mapping rows by order instead of names.")
      }
    } else {
      warning("Mapping rows of data frame (`by`) to individuals' identifiers (`x@id`) by order.\n",
              "Consider setting row names of `by` to match `x@id` for clarity.")
    }
    res = .splitPop(x, by)
    return(res)
  }
  
  # Get max depth of nesting in MultiPop
  md = ifelse(isMultiPop(x), .depthMultiPop(x), 1L)

  # Validate level argument
  if (length(level) == 1L) {
    if (is.infinite(level)) {
      levels = Inf
    } else {
      if (!is.numeric(level) || is.na(level) ||
          level < 1 || level != as.integer(level)) {
        stop("`level` must be a positive integer or Inf")
      }
      if (level > md) {
        stop(sprintf("requested `level` exceeds max depth of `x` (%d)", md))
      }
      levels = as.integer(level)
    }
  } else {
    if (!is.numeric(level) || any(is.na(level)) || any(level <= 0)) {
      stop("`level` must be a numeric vector of positive integers (no NA)")
    }
    if (any(is.infinite(level))) {
      stop("cannot mix Inf with integer levels")
    }
    if (any(level != trunc(level))) {
      stop("`level` must be a numeric vector of positive integers (no NA)")
    }
    if (any(level > md)) {
      stop(sprintf("requested level(s) exceed max depth of `x` (%d)", md))
    }
    levels = as.integer(unique(level))
  }

  # Recursively split at requested levels
  res = .splitAtLevels(x, currentLevel = 1L, levels = levels, by = by)
  return(res)
}

#' Helper to apply splitting only at requested nesting levels
#'
#' @param obj A \code{\link{Pop-class}} or \code{\link{MultiPop-class}} object.
#' @param currentLevel Integer; current depth during recursion (top-level = 1).
#' @param levels Integer vector or \code{Inf}; levels at which to apply splits.
#' @param by A list of grouping specs (vectors or functions).
#'
#' @return The input object with splits applied at the requested levels.
#'
#' @keywords internal
.splitAtLevels = function(obj, currentLevel, levels, by) {
  if (isPop(obj)) {
    # A Pop at the currentLevel: split only if currentLevel requested (or Inf)
    if (any(is.infinite(levels)) || currentLevel %in% levels) {
      return(.splitPop(obj, by))
    } else {
      return(obj)
    }
  }

  if (isMultiPop(obj)) {
    # For each child: if child is MultiPop, its children are one level deeper;
    # if child is Pop, it sits at currentLevel.
    obj@pops = lapply(obj@pops, function(child) {
      if (isMultiPop(child)) {
        .splitAtLevels(child, currentLevel = currentLevel + 1L, levels = levels, by = by)
      } else {
        .splitAtLevels(child, currentLevel = currentLevel, levels = levels, by = by)
      }
    })
    validObject(obj)
    return(obj)
  }
}

#' Helper function to recursively split a Pop object
#'
#' @param pop A \code{\link{Pop-class}} object.
#' @param by A grouping specification. This may be:  \cr
#' - an atomic vector.  \cr
#' - a function returning an atomic vector.  \cr
#' - a list of vectors/functions for recursive nested splitting.  \cr
#' - a data frame whose columns define nested grouping levels.
#' 
#' For data frames, rows are aligned to \code{pop@id} when possible, otherwise
#' they are assumed to already be in population order. During recursion, the
#' data frame is subset to each child population so that row alignment is
#' preserved across levels.
#'
#' @return A \code{\link{MultiPop-class}} produced by recursively applying \code{by}.
#'
#' @keywords internal
.splitPop = function(pop, by) {
  if (length(by) == 0) {
    return(pop)
  }

  f = by[[1]]
  groups = if (is.function(f)) f(pop) else f

  if (!is.atomic(groups)) {
    stop("Grouping spec must be an atomic vector or a function returning one.")
  }
  if (all(is.na(groups))) {
    return(pop)
  }
  if (any(is.na(groups))) {
    stop("Grouping vector contains NA values.")
  }

  popList = split(pop, groups)
  mp = newEmptyMultiPop()
  if (is.data.frame(by)){
    groups = lapply(split(by, groups), `[`, -1)
    mp@pops = mapply(FUN = .splitPop, pop = popList, by = groups)
  } else {
    mp@pops = lapply(popList, .splitPop, by = by[-1])
  }
  return(mp)
}