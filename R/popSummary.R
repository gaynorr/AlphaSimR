# fmt: skip file

#' @title Mean genetic values
#'
#' @description
#' Returns mean genetic values for all traits in a population.
#' Supports \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, and
#' \code{\link{MultiPop-class}} inputs. For \code{MultiPop} objects, output can
#' optionally be simplified to a requested nesting level.
#'
#' @param pop A \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, or
#'   \code{\link{MultiPop-class}} object.
#' @param simplify Logical. Only used when \code{pop} is a
#'   \code{\link{MultiPop-class}}. If \code{TRUE}, flatten \code{pop} to the
#'   requested \code{level} with \code{\link{flattenMultiPop}}, combine results
#'   across populations with \code{\link{rbind}}, and attach a \code{"source"}
#'   attribute describing row origins.
#' @param level Integer scalar \eqn{\ge 1}. Only used when \code{pop} is a
#'   \code{\link{MultiPop-class}} and \code{simplify=TRUE}. Number of
#'   \code{MultiPop} levels to preserve. Passed to
#'   \code{\link{flattenMultiPop}}. Ignored when \code{simplify=FALSE}.
#'
#' @details
#' When \code{simplify=FALSE}, \code{MultiPop} structure is preserved and output
#' follows the same nesting as \code{pop}.
#'
#' When \code{simplify=TRUE}, the result is simplified to the requested
#' \code{level}. If \code{level} exceeds the nesting depth of \code{pop}, the
#' structure is returned unchanged. If \code{level} \eqn{\lt 1}, the function 
#' sets \code{level=1} issuing a warning.
#'
#' Simplified output includes a \code{"source"} attribute. This is a data frame
#' with columns \code{level1}, \code{level2}, etc., indicating the name (or
#' index) of the source population for each row at each retained nesting level.
#' 
#' @seealso
#' \code{\link{meanGPop}}, \code{\link{gv}}
#'
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or \code{\link{HybridPop-class}},
#' returns a numeric vector of mean genetic values (one value per trait).
#'
#' If \code{pop} is a \code{\link{MultiPop-class}} and \code{simplify=FALSE},
#' returns a list matching the \code{MultiPop} nesting structure, with one
#' numeric vector per \code{Pop} object. If \code{simplify=TRUE}, the nested
#' list is simplified down to a requested \code{level} by binding individual
#' \code{Pop} results with \code{\link{rbind}}: rows correspond to populations
#' and columns to traits and a \code{"source"} attribute records row origins.
#' If \code{level=1} the whole structure is simplified to a numeric matrix.
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#' meanG(pop)
#'
#' pop2 = randCross(pop, nCrosses = 3, nProgeny = 4)
#'
#' # Create a nested MultiPop
#' mp1 = splitPop(
#'   pop2,
#'   by = list(
#'     function(x) rep(LETTERS[1:2], length.out = length(x)),
#'     function(x) paste(x@mother, x@father, sep = "_")
#'   )
#' )
#'
#' # Calculate mean genetic values and simplify at different levels
#' meanG(mp1, simplify = FALSE)
#' meanG(mp1, simplify = TRUE, level = 2)
#' meanG(mp1, simplify = TRUE, level = 1)
#'
#' @export
meanG = function(pop, simplify = FALSE, level = 1L){
  if (simplify && level < 1L) {
    warning("`level` should be >= 1. Setting default `level=1`")
    level = 1L
  }
  calcPopValue(
    pop,
    FUN = function(x) colMeans(x@gv),
    simplify = simplify,
    level = level
  )
}

#' @title Mean phenotype values
#'
#' @description
#' Returns mean phenotype values for all traits in a population.
#' Supports \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, and
#' \code{\link{MultiPop-class}} inputs. For \code{MultiPop} objects, output can
#' optionally be simplified to a requested nesting level.
#'
#' @param pop A \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, or
#'   \code{\link{MultiPop-class}} object.
#' @param simplify Logical. Only used when \code{pop} is a
#'   \code{\link{MultiPop-class}}. If \code{TRUE}, flatten \code{pop} to the
#'   requested \code{level} with \code{\link{flattenMultiPop}}, combine results
#'   across populations with \code{\link{rbind}}, and attach a \code{"source"}
#'   attribute describing row origins.
#' @param level Integer scalar \eqn{\ge 1}. Only used when \code{pop} is a
#'   \code{\link{MultiPop-class}} and \code{simplify=TRUE}. Number of
#'   \code{MultiPop} levels to preserve. Passed to
#'   \code{\link{flattenMultiPop}}. Ignored when \code{simplify=FALSE}.
#'
#' @details
#' When \code{simplify=FALSE}, \code{MultiPop} structure is preserved and output
#' follows the same nesting as \code{pop}.
#'
#' When \code{simplify=TRUE}, the result is simplified to the requested
#' \code{level}. If \code{level} exceeds the nesting depth of \code{pop}, the
#' structure is returned unchanged. If \code{level} \eqn{\lt 1}, the function 
#' sets \code{level=1} issuing a warning.
#'
#' Simplified output includes a \code{"source"} attribute. This is a data frame
#' with columns \code{level1}, \code{level2}, etc., indicating the name (or
#' index) of the source population for each row at each retained nesting level.
#' 
#' @seealso
#' \code{\link{meanPPop}}, \code{\link{pheno}}, \code{\link{setPheno}}
#'
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or \code{\link{HybridPop-class}},
#' returns a numeric vector of mean phenotype values (one value per trait).
#'
#' If \code{pop} is a \code{\link{MultiPop-class}} and \code{simplify=FALSE},
#' returns a list matching the \code{MultiPop} nesting structure, with one
#' numeric vector per \code{Pop} object. If \code{simplify=TRUE}, the nested
#' list is simplified down to a requested \code{level} by binding individual
#' \code{Pop} results with \code{\link{rbind}}: rows correspond to populations
#' and columns to traits and a \code{"source"} attribute records row origins.
#' If \code{level=1} the whole structure is simplified to a numeric matrix.
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#' meanP(pop)
#'
#' pop2 = randCross(pop, nCrosses = 3, nProgeny = 4)
#'
#' # Create a nested MultiPop
#' mp1 = splitPop(
#'   pop2,
#'   by = list(
#'     function(x) rep(LETTERS[1:2], length.out = length(x)),
#'     function(x) paste(x@mother, x@father, sep = "_")
#'   )
#' )
#'
#' # Calculate mean phenotype values and simplify at different levels
#' meanP(mp1, simplify = FALSE)
#' meanP(mp1, simplify = TRUE, level = 2)
#' meanP(mp1, simplify = TRUE, level = 1)
#'
#' @export
meanP = function(pop, simplify = FALSE, level = 1L){
  if (simplify && level < 1L) {
    warning("`level` should be >= 1. Setting default `level=1`")
    level = 1L
  }
  calcPopValue(
    pop,
    FUN = function(x) colMeans(x@pheno),
    simplify = simplify,
    level = level
  )
}

#' @title Mean estimated breeding values
#'
#' @description
#' Returns mean estimated breeding values for all traits in a population.
#' Supports \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, and
#' \code{\link{MultiPop-class}} inputs. For \code{MultiPop} objects, output can
#' optionally be simplified to a requested nesting level.
#'
#' @param pop A \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, or
#'   \code{\link{MultiPop-class}} object.
#' @param simplify Logical. Only used when \code{pop} is a
#'   \code{\link{MultiPop-class}}. If \code{TRUE}, flatten \code{pop} to the
#'   requested \code{level} with \code{\link{flattenMultiPop}}, combine results
#'   across populations with \code{\link{rbind}}, and attach a \code{"source"}
#'   attribute describing row origins.
#' @param level Integer scalar \eqn{\ge 1}. Only used when \code{pop} is a
#'   \code{\link{MultiPop-class}} and \code{simplify=TRUE}. Number of
#'   \code{MultiPop} levels to preserve. Passed to
#'   \code{\link{flattenMultiPop}}. Ignored when \code{simplify=FALSE}.
#'
#' @details
#' When \code{simplify=FALSE}, \code{MultiPop} structure is preserved and output
#' follows the same nesting as \code{pop}.
#'
#' When \code{simplify=TRUE}, output is simplified to the requested
#' \code{level}. If \code{level} exceeds the nesting depth of \code{pop}, the
#' structure is returned unchanged. If \code{level} \eqn{< 1}, the function
#' sets \code{level=1} and issues a warning.
#'
#' Simplified output requires compatible \code{@ebv} columns 
#' (see \code{\link{ebv}}) across terminal populations so that results can be
#' combined with \code{\link{rbind}}.
#'
#' Simplified output includes a \code{"source"} attribute. This is a data frame
#' with columns \code{level1}, \code{level2}, etc., indicating the name (or
#' index) of the source population for each row at each retained nesting level.
#' 
#' @seealso
#' \code{\link{ebv}}, \code{\link{setEBV}}
#'
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or \code{\link{HybridPop-class}},
#' returns a numeric vector of mean estimated breeding values (one value per
#' column in \code{ebv(pop)}).
#'
#' If \code{pop} is a \code{\link{MultiPop-class}} and \code{simplify=FALSE},
#' returns a list matching the \code{MultiPop} nesting structure, with one
#' numeric vector per \code{Pop} object. If \code{simplify=TRUE}, the nested
#' list is simplified down to a requested \code{level} by binding individual
#' \code{Pop} results with \code{\link{rbind}}: rows correspond to populations
#' and columns to traits in \code{ebv(pop)} and a \code{"source"} attribute 
#' records row origins. If \code{level=1} the whole structure is simplified 
#' to a numeric matrix.
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' trtH2 = 0.5
#' SP$setVarE(h2=trtH2)
#' \dontshow{SP$nThreads = 1L}
#'
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop2 = randCross(pop, nCrosses = 3, nProgeny = 4)
#'
#' # Individual performance based EBV
#' pop2@ebv = trtH2 * (pheno(pop2) - meanP(pop2))
#' meanEBV(pop2)
#'
#' # Create a nested MultiPop
#' mp1 = splitPop(
#'   pop2,
#'   by = list(
#'     function(x) rep(LETTERS[1:2], length.out = length(x)),
#'     function(x) paste(x@mother, x@father, sep = "_")
#'   )
#' )
#'
#' # Calculate mean estimated breeding values and simplify at different levels
#' meanEBV(mp1, simplify = FALSE)
#' meanEBV(mp1, simplify = TRUE, level = 2)
#' meanEBV(mp1, simplify = TRUE, level = 1)
#'
#' @export
meanEBV = function(pop, simplify = FALSE, level = 1L){
  if (simplify && level < 1L) {
    warning("`level` should be >= 1. Setting default `level=1`")
    level = 1L
  }
  calcPopValue(
    pop,
    FUN = function(x) colMeans(x@ebv),
    simplify = simplify,
    level = level
  )
}

#' @title Mean genetic values between \code{Pops} in a \code{MultiPop}
#'
#' @description Computes mean genetic values for individual
#' \code{\link{Pop-class}} objects in a \code{\link{MultiPop-class}} object.
#' Means are then recursively aggregated upward through the nested
#' \code{MultiPop} structure by taking the mean of child-population means
#' at each level. Aggregation is intentionally unweighted across child 
#' populations.
#'
#' @param x A \code{\link{Pop-class}} or \code{\link{MultiPop-class}} object.
#' @param level Integer scalar \code{>= 0} indicating the requested aggregation
#'   level (see Details).
#' @param .req_level Internal argument used during recursion. Do not set
#'   manually.
#'
#' @details
#' The \code{level} argument controls the nesting level in a
#' \code{MultiPop} structure at which mean aggregation stops:  \cr
#' - When \code{level>0}, the function returns a numeric matrix with one
#'   row per \code{Pop} or \code{MultiPop} unit at the requested \code{level}.
#'   The returned matrix includes a \code{"source"} attribute with columns
#'   \code{level1}, \code{level2}, etc., indicating the origin of each row,
#'   and values indicating the name or index of the population at each
#'   nesting level.  \cr
#' - When \code{level=0}, the function returns a numeric vector giving the
#'   overall mean-of-means across all branches for each trait. Any
#'   \code{level<0} is treated as \code{level=0}.
#' 
#' @seealso
#' \code{\link{meanG}}, \code{\link{gv}}
#'
#' @return
#' If \code{x} is a \code{\link{Pop-class}} object, or if \code{level=0}, a
#' numeric vector of mean genetic values for each trait.
#'
#' If \code{x} is a \code{\link{MultiPop-class}} object and \code{level>0}, a
#' numeric matrix of (possibly aggregated) mean genetic values for each trait.
#' The origin of each row is described by the \code{"source"} attribute.
#'
#' @examples
#' founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' \dontshow{SP$nThreads = 1L}
#' pop = newPop(founderPop, simParam = SP)
#'
#' mp = splitPop(
#'   pop,
#'   by = list(
#'     sample(rep(LETTERS[1:2], length.out = pop@nInd)),
#'     function(x) sample(letters[5:6], length(x), replace = TRUE)
#'   )
#' )
#'
#' meanGPop(mp, level = 1)
#' meanGPop(mp, level = 0)
#'
#' @export
meanGPop = function(x, level = 0, .req_level) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level)) {
    stop("`level` must be a single non-NA integer value.")
  }
  level = as.integer(level)

  if (isPop(x)) {
    if (nrow(x@gv) == 0L) {
      stop("One of the populations in `x` is empty")
    }
    return(colMeans(x@gv))
  }
  stopifnot(isMultiPop(x))
  if (length(x) == 0L) {
    stop("`x` contains no populations.")
  }

  if (missing(.req_level)) {
    .req_level = level
  }

  popValueList = lapply(
    x@pops,
    meanGPop,
    level = level - 1,
    .req_level = .req_level
  )
  popValues = do.call('rbind', popValueList)

  if (level < 1) {
    return(colMeans(popValues))
  } else {
    if (level == .req_level) {
      md = .depthMultiPop(x)
      cols = ifelse(level > md, md, level)
      paths = .collectLeafPaths(x)
      src = .formatPopSource(
        paths = paths,
        nRows = rep(1, length(paths)),
        level_offset = 0L
      )
      src = src[, seq_len(cols), drop = FALSE]
      src = unique(src)
      rownames(src) = NULL
      attr(popValues, "source") = src
      rownames(popValues) = NULL
    }
    return(popValues)
  }
}

#' @title Mean phenotype values between \code{Pops} in a \code{MultiPop}
#'
#' @description Computes mean phenotype values for individual
#' \code{\link{Pop-class}} objects in a \code{\link{MultiPop-class}} object.
#' Means are then recursively aggregated upward through the nested
#' \code{MultiPop} structure by taking the mean of child-population means
#' at each level. Aggregation is intentionally unweighted across child 
#' populations.
#'
#' @param x A \code{\link{Pop-class}} or \code{\link{MultiPop-class}} object.
#' @param level Integer scalar \code{>= 0} indicating the requested aggregation
#'   level (see Details).
#' @param .req_level Internal argument used during recursion. Do not set
#'   manually.
#'
#' @details
#' The \code{level} argument controls the nesting level in a
#' \code{MultiPop} structure at which mean aggregation stops:  \cr
#' - When \code{level>0}, the function returns a numeric matrix with one
#'   row per \code{Pop} or \code{MultiPop} unit at the requested \code{level}.
#'   The returned matrix includes a \code{"source"} attribute with columns
#'   \code{level1}, \code{level2}, etc., indicating the origin of each row,
#'   and values indicating the name or index of the population at each
#'   nesting level.  \cr
#' - When \code{level=0}, the function returns a numeric vector giving the
#'   overall mean-of-means across all branches for each trait. Any
#'   \code{level<0} is treated as \code{level=0}.
#' 
#' @seealso
#' \code{\link{meanP}}, \code{\link{pheno}}, \code{\link{setPheno}}
#'
#' @return
#' If \code{x} is a \code{\link{Pop-class}} object, or if \code{level=0}, a
#' numeric vector of mean phenotype values for each trait.
#'
#' If \code{x} is a \code{\link{MultiPop-class}} object and \code{level>0}, a
#' numeric matrix of (possibly aggregated) mean phenotype values for each trait.
#' The origin of each row is described by the \code{"source"} attribute.
#'
#' @examples
#' founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' \dontshow{SP$nThreads = 1L}
#' pop = newPop(founderPop, simParam = SP)
#'
#' mp = splitPop(
#'   pop,
#'   by = list(
#'     sample(rep(LETTERS[1:2], length.out = pop@nInd)),
#'     function(x) sample(letters[5:6], length(x), replace = TRUE)
#'   )
#' )
#'
#' meanPPop(mp, level = 1)
#' meanPPop(mp, level = 0)
#'
#' @export
meanPPop = function(x, level = 0, .req_level) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level)) {
    stop("`level` must be a single non-NA integer value.")
  }
  level = as.integer(level)

  if (isPop(x)) {
    if (nrow(x@pheno) == 0L) {
      stop("One of the populations in `x` is empty")
    }
    return(colMeans(x@pheno))
  }
  stopifnot(isMultiPop(x))
  if (length(x) == 0L) {
    stop("`x` contains no populations.")
  }

  if (missing(.req_level)) {
    .req_level = level
  }

  popValueList = lapply(
    x@pops,
    meanPPop,
    level = level - 1,
    .req_level = .req_level
  )
  popValues = do.call('rbind', popValueList)

  if (level < 1) {
    return(colMeans(popValues))
  } else {
    if (level == .req_level) {
      md = .depthMultiPop(x)
      cols = ifelse(level > md, md, level)
      paths = .collectLeafPaths(x)
      src = .formatPopSource(
        paths = paths,
        nRows = rep(1, length(paths)),
        level_offset = 0L
      )
      src = src[, seq_len(cols), drop = FALSE]
      src = unique(src)
      rownames(src) = NULL
      attr(popValues, "source") = src
      rownames(popValues) = NULL
    }
    return(popValues)
  }
}

#' @title Total genetic variance
#'
#' @description Returns total genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop)
#'
#' @export
varG = function(pop){
  G = popVar(pop@gv)
  rownames(G) = colnames(G) = colnames(pop@gv)
  return(G)
}

#' @title Phenotypic variance
#'
#' @description Returns phenotypic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varP(pop)
#'
#' @export
varP = function(pop){
  P = popVar(pop@pheno)
  rownames(P) = colnames(P) = colnames(pop@pheno)
  return(P)
}

#' @title Variance of estimated breeding values
#'
#' @description Returns variance of estimated breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' trtH2 = 0.5
#' SP$setVarE(h2=trtH2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@ebv = trtH2 * (pop@pheno - meanP(pop)) #ind performance based EBV
#' varA(pop)
#' varEBV(pop)
#'
#' @export
varEBV = function(pop){
  ebv = popVar(pop@ebv)
  rownames(ebv) = colnames(ebv) = colnames(pop@ebv)
  return(ebv)
}

#' @title Sumarize genetic parameters
#'
#' @description
#' Calculates genetic and genic additive and dominance variances
#' for an object of \code{\link{Pop-class}}
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return
#' \describe{
#' \item{varA}{an nTrait by nTrait matrix of additive genetic variances}
#' \item{varD}{an nTrait by nTrait matrix of dominance genetic variances}
#' \item{varAA}{an nTrait by nTrait matrix of additive-by-additive genetic variances}
#' \item{varG}{an nTrait by nTrait matrix of total genetic variances}
#' \item{genicVarA}{an nTrait vector of additive genic variances}
#' \item{genicVarD}{an nTrait vector of dominance genic variances}
#' \item{genicVarAA}{an nTrait vector of additive-by-additive genic variances}
#' \item{genicVarG}{an nTrait vector of total genic variances}
#' \item{covA_HW}{an nTrait vector of additive covariances due to non-random mating}
#' \item{covD_HW}{an nTrait vector of dominance covariances due to non-random mating}
#' \item{covAA_HW}{an nTrait vector of additive-by-additive covariances due to non-random mating}
#' \item{covG_HW}{an nTrait vector of total genic covariances due to non-random mating}
#' \item{covA_L}{an nTrait vector of additive covariances due to linkage disequilibrium}
#' \item{covD_L}{an nTrait vector of dominance covariances due to linkage disequilibrium}
#' \item{covAA_L}{an nTrait vector of additive-by-additive covariances due to linkage disequilibrium}
#' \item{covAD_L}{an nTrait vector of additive by dominance covariances due to linkage disequilibrium}
#' \item{covAAA_L}{an nTrait vector of additive by additive-by-additive covariances due to linkage disequilibrium}
#' \item{covDAA_L}{an nTrait vector of dominance by additive-by-additive covariances due to linkage disequilibrium}
#' \item{covG_L}{an nTrait vector of total genic covariances due to linkage disequilibrium}
#' \item{mu}{an nTrait vector of trait means}
#' \item{mu_HW}{an nTrait vector of expected trait means under random mating}
#' \item{gv}{a matrix of genetic values with dimensions nInd by nTraits}
#' \item{bv}{a matrix of breeding values with dimensions nInd by nTraits}
#' \item{dd}{a matrix of dominance deviations with dimensions nInd by nTraits}
#' \item{aa}{a matrix of additive-by-additive epistatic deviations with dimensions nInd by nTraits}
#' \item{gv_mu}{an nTrait vector of intercepts with dimensions nInd by nTraits}
#' \item{gv_a}{a matrix of additive genetic values with dimensions nInd by nTraits}
#' \item{gv_d}{a matrix of dominance genetic values with dimensions nInd by nTraits}
#' \item{gv_aa}{a matrix of additive-by-additive genetic values with dimensions nInd by nTraits}
#' \item{alpha}{a list of average allele subsitution effects with length nTraits}
#' \item{alpha_HW}{a list of average allele subsitution effects at Hardy-Weinberg equilibrium with length nTraits}
#' }
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' ans = genParam(pop, simParam=SP)
#'
#' @export
genParam = function(pop,simParam=NULL,nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  nInd = nInd(pop)
  nTraits = simParam$nTraits
  traitNames = simParam$traitNames

  # Blank nInd x nTrait matrices
  gv = matrix(NA_real_, nrow=nInd, ncol=nTraits)
  colnames(gv) = traitNames
  bv = dd = aa = gv_a = gv_d = gv_aa = gv

  # Blank nTrait vectors
  genicVarA = rep(NA_real_, nTraits)
  names(genicVarA) = traitNames
  genicVarD = genicVarAA = covA_HW = covD_HW = covAA_HW =
    covG_HW = mu = mu_HW = gv_mu = covAAA_L = covDAA_L =
    covAD_L = genicVarA

  # Average effect of an allele substitution
  alpha = vector("list", length=nTraits)
  names(alpha) = traitNames
  alpha_HW = alpha

  #Loop through trait calculations
  for(i in seq_len(nTraits)){
    trait = simParam$traits[[i]]
    tmp = calcGenParam(trait,pop,nThreads)
    genicVarA[i] = tmp$genicVarA2
    covA_HW[i] = tmp$genicVarA-tmp$genicVarA2
    gv[,i] = tmp$gv
    bv[,i] = tmp$bv
    mu[i] = tmp$mu
    mu_HW[i] = tmp$mu_HWE
    gv_a[,i] = tmp$gv_a
    gv_mu[i] = tmp$gv_mu
    if(.hasSlot(trait,"domEff")){
      genicVarD[i] = tmp$genicVarD2
      covD_HW[i] = tmp$genicVarD-tmp$genicVarD2
      dd[,i] = tmp$dd
      gv_d[,i] = tmp$gv_d
    }else{
      genicVarD[i] = 0
      covD_HW[i] = 0
      dd[,i] = rep(0,pop@nInd)
      gv_d[,i] = rep(0,pop@nInd)
    }
    if(.hasSlot(trait,"epiEff")){
      genicVarAA[i] = tmp$genicVarAA2
      covAA_HW[i] = tmp$genicVarAA-tmp$genicVarAA2
      aa[,i] = tmp$aa
      gv_aa[,i] = tmp$gv_aa
    }else{
      genicVarAA[i] = 0
      covAA_HW[i] = 0
      aa[,i] = rep(0,pop@nInd)
      gv_aa[,i] = rep(0,pop@nInd)
    }
    if(nInd==1){
      covAD_L[i] = 0
      covAAA_L[i] = 0
      covDAA_L[i] = 0
    } else {
      covAD_L[i] = popVar(cbind(bv[,i],dd[,i]))[1,2]
      covAAA_L[i] = popVar(cbind(bv[,i],aa[,i]))[1,2]
      covDAA_L[i] = popVar(cbind(dd[,i],aa[,i]))[1,2]
    }
    alpha[[i]] = tmp$alpha
    alpha_HW[[i]] = tmp$alpha_HW
  }

  varA = popVar(bv)
  rownames(varA) = colnames(varA) = traitNames

  varD = popVar(dd)
  rownames(varD) = colnames(varD) = traitNames

  varAA = popVar(aa)
  rownames(varAA) = colnames(varAA) = traitNames

  varG = popVar(gv)
  rownames(varG) = colnames(varG) = traitNames

  genicVarG = genicVarA + genicVarD + genicVarAA
  covG_HW = covA_HW + covD_HW + covAA_HW

  output = list(varA=varA,
                varD=varD,
                varAA=varAA,
                varG=varG,
                genicVarA=genicVarA,
                genicVarD=genicVarD,
                genicVarAA=genicVarAA,
                genicVarG=genicVarG,
                covA_HW=covA_HW,
                covD_HW=covD_HW,
                covAA_HW=covAA_HW,
                covG_HW=covG_HW,
                covA_L=diag(varA)-genicVarA-covA_HW,
                covD_L=diag(varD)-genicVarD-covD_HW,
                covAA_L=diag(varAA)-genicVarAA-covAA_HW,
                covAD_L=covAD_L,
                covAAA_L=covAAA_L,
                covDAA_L=covDAA_L,
                covG_L=diag(varG)-genicVarG-covG_HW,
                mu=mu,
                mu_HW=mu_HW,
                gv=gv,
                bv=bv,
                dd=dd,
                aa=aa,
                gv_mu=gv_mu,
                gv_a=gv_a,
                gv_d=gv_d,
                gv_aa=gv_aa,
                alpha=alpha,
                alpha_HW=alpha_HW)
  return(output)
}

#' @title Additive variance
#'
#' @description Returns additive variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varA(pop, simParam=SP)
#'
#' @export
varA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varA
}

#' @title Dominance variance
#'
#' @description Returns dominance variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varD(pop, simParam=SP)
#'
#' @export
varD = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varD
}

#' @title Additive-by-additive epistatic variance
#'
#' @description Returns additive-by-additive epistatic
#' variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varAA(pop, simParam=SP)
#'
#' @export
varAA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varAA
}

#' @title Breeding value
#'
#' @description Returns breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' bv(pop, simParam=SP)
#'
#' @export
bv = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$bv
}

#' @title Dominance deviations
#'
#' @description Returns dominance deviations for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' dd(pop, simParam=SP)
#'
#' @export
dd = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$dd
}

#' @title Additive-by-additive epistatic deviations
#'
#' @description Returns additive-by-additive epistatic
#' deviations for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' aa(pop, simParam=SP)
#'
#' @export
aa = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$aa
}

#' @title Additive genic variance
#'
#' @description Returns additive genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#'
#' @export
genicVarA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarA
}

#' @title Dominance genic variance
#'
#' @description Returns dominance genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarD(pop, simParam=SP)
#'
#' @export
genicVarD = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarD
}

#' @title Additive-by-additive genic variance
#'
#' @description Returns additive-by-additive epistatic
#' genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarAA(pop, simParam=SP)
#'
#' @export
genicVarAA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarAA
}

#' @title Total genic variance
#'
#' @description Returns total genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#'
#' @export
genicVarG = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarG
}

#' @title Genetic value
#'
#' @description A wrapper for accessing the gv slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' gv(pop)
#'
#' @export
gv = function(pop){
  pop@gv
}

#' @title Phenotype
#'
#' @description A wrapper for accessing the pheno slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pheno(pop)
#'
#' @export
pheno = function(pop){
  pop@pheno
}

#' @title Estimated breeding value
#'
#' @description A wrapper for accessing the ebv slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@ebv = matrix(rnorm(pop@nInd), nrow=pop@nInd, ncol=1)
#' ebv(pop)
#'
#' @export
ebv = function(pop){
  pop@ebv
}

#' @title Calculate parent average
#'
#' @param pop \code{\link{Pop-class}} with individuals whose parent average
#'   will be calculated
#' @param parents \code{\link{Pop-class}} with mothers and fathers of individuals
#'   in \code{pop}; if \code{NULL} must provide \code{mothers} and \code{fathers}
#' @param mothers \code{\link{Pop-class}} with mothers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param fathers \code{\link{Pop-class}} with fathers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param use character, calculate using \code{"\link{gv}"}, \code{"\link{bv}"},
#'   \code{"\link{ebv}"}, or \code{"\link{pheno}"}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return a matrix of parent averages with dimensions nInd by nTraits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop2 = randCross(pop, nCrosses=10, nProgeny=2)
#' parentAverage(pop2, parents = pop)
#' parentAverage(pop2, mothers = pop, fathers = pop)
#'
#' @export
parentAverage = function(pop, parents = NULL, mothers = NULL, fathers = NULL,
                         use = "gv", simParam = NULL, nThreads=NULL) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if (!is.null(parents)) {
    matchMothers = match(x = pop@mother, table = parents@id)
    matchFathers = match(x = pop@father, table = parents@id)
  } else {
    if (is.null(mothers) | is.null(fathers)) {
      stop("must provide either 'parents' or both 'mothers' and 'fathers'!")
    }
    matchMothers = match(x = pop@mother, table = mothers@id)
    matchFathers = match(x = pop@father, table = fathers@id)
  }
  if (anyNA(matchMothers)) {
    stop("some parents/mothers not found!")
  }
  if (anyNA(matchFathers)) {
    stop("some parents/fathers not found!")
  }
  if (use %in% c("gv", "ebv", "pheno")) {
    if (!is.null(parents)) {
      ret = 0.5 * (slot(object = parents, name = use)[matchMothers, , drop = FALSE] +
                   slot(object = parents, name = use)[matchFathers, , drop = FALSE])
    } else {
      ret = 0.5 * (slot(object = mothers, name = use)[matchMothers, , drop = FALSE] +
                   slot(object = fathers, name = use)[matchFathers, , drop = FALSE])
    }
  } else if (use == "bv") {
    if (!is.null(parents)) {
      ret = 0.5 * (bv(parents, simParam = simParam,
                      nThreads=nThreads)[matchMothers, , drop = FALSE] +
                   bv(parents, simParam = simParam,
                      nThreads=nThreads)[matchFathers, , drop = FALSE])
    } else {
      ret = 0.5 * (bv(mothers, simParam = simParam,
                      nThreads=nThreads)[matchMothers, , drop = FALSE] +
                   bv(fathers, simParam = simParam,
                      nThreads=nThreads)[matchFathers, , drop = FALSE])
    }
  } else {
    stop("use must be one of 'gv', 'bv', 'ebv', or 'pheno'!")
  }
  return(ret)
}

#' @title Calculate Mendelian sampling
#'
#' @param pop \code{\link{Pop-class}} with individuals whose parent average
#'   will be calculated
#' @param parents \code{\link{Pop-class}} with mothers and fathers of individuals
#'   in \code{pop}; if \code{NULL} must provide \code{mothers} and \code{fathers}
#' @param mothers \code{\link{Pop-class}} with mothers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param fathers \code{\link{Pop-class}} with fathers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param use character, calculate using \code{"\link{gv}"}, \code{"\link{bv}"},
#'   \code{"\link{ebv}"}, or \code{"\link{pheno}"}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return a matrix of Mendelian samplings with dimensions nInd by nTraits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop2 = randCross(pop, nCrosses=10, nProgeny=2)
#' mendelianSampling(pop2, parents = pop)
#' mendelianSampling(pop2, mothers = pop, fathers = pop)
#'
#' @export
mendelianSampling = function(pop, parents = NULL, mothers = NULL, fathers = NULL,
                             use = "gv", simParam = NULL, nThreads=NULL) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  pa = parentAverage(pop = pop, parents = parents, mothers = mothers, fathers = fathers,
                     use = use, simParam = simParam, nThreads=nThreads)
  if (use %in% c("gv", "ebv", "pheno")) {
    ret = slot(object = pop, name = use) - pa
  } else if (use == "bv") {
    ret = bv(pop, simParam = simParam, nThreads=nThreads) - pa
  } else {
    stop("use must be one of 'gv', 'bv', 'ebv', or 'pheno'!")
  }
  return(ret)
}

#' @title Number of individuals
#'
#' @description A wrapper for accessing the nInd slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' nInd(pop)
#'
#' @export
nInd = function(pop){
  pop@nInd
}
