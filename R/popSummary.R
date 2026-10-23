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
#' structure is returned unchanged. If \code{level} \eqn{< 1}, the function
#' sets \code{level=1} issuing a warning.
#'
#' Simplified output includes a \code{"source"} attribute. This is a data frame
#' with columns \code{level1}, \code{level2}, etc., indicating the name (or
#' index) of the source population for each row at each retained nesting level.
#' 
#' The mean is calculated from the genetic values in \code{pop}.
#' It includes the trait intercept and can differ from both the intercept
#' and the requested founder mean.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
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
#' structure is returned unchanged. If \code{level} \eqn{< 1}, the function
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
#' This function calculated the mean of values stored in the \code{ebv} slot;
#' supplied by the user or \code{\link{setEBV}}.
#' Depending on how it was populated, it can contain estimated breeding
#' values, estimated genetic values, or other estimated/predicted values.
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
#' @param level Integer scalar \eqn{\ge 0} indicating the requested aggregation
#'   level (see Details).
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
#'   overall mean-of-means across all branches for each trait.
#'
#' @seealso
#' \code{\link{meanG}}, \code{\link{gv}}, \code{\link{calcPopSummary}}
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
meanGPop = function(x, level = 0) {
  calcPopSummary(
    x,
    withinPopFUN = function(x) colMeans(x@gv),
    acrossPopFUN = colMeans,
    level = level
  )
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
#' @param level Integer scalar \eqn{\ge 0} indicating the requested aggregation
#'   level (see Details).
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
#'   overall mean-of-means across all branches for each trait.
#' 
#' @seealso
#' \code{\link{meanP}}, \code{\link{pheno}}, \code{\link{setPheno}},
#'   \code{\link{calcPopSummary}}
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
#' SP$setVarE(h2 = 0.5)
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
meanPPop = function(x, level = 0) {
  calcPopSummary(
    x,
    withinPopFUN = function(x) colMeans(x@pheno),
    acrossPopFUN = colMeans,
    level = level
  )
}

#' @title Total genetic variance within a \code{Pop}
#'
#' @description
#' Returns total genetic variance for all traits within a \code{Pop}.
#' 
#' In contrast to \code{\link{varGPop}}, which computes variance between
#' population means, \code{varG} computes variance within populations from
#' individual genetic values.
#'
#' @param pop A \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, or
#'   \code{\link{MultiPop-class}} object.
#' 
#' @details
#' A variance-covariance matrix of genetic values is computed for each
#' \code{Pop} across traits.
#' 
#' This is the variance of genetic values in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @seealso \code{\link{varGPop}}
#' 
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or \code{\link{HybridPop-class}},
#' returns a variance-covariance matrix of genetic values across traits.
#'
#' If \code{pop} is a \code{\link{MultiPop-class}}, returns a nested list
#' matching the \code{MultiPop} structure, with one variance-covariance matrix
#' per \code{Pop} object.
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' \dontshow{SP$nThreads = 1L}
#'
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' # Return a variance-covariance matrix for the population
#' varG(pop)
#' 
#' # Randomly split population into two sub-populations
#' mp = splitPop(pop, sample(c('A','B'), size = 10, replace = TRUE))
#' 
#' # Return a list of variance-covariance matrices for the MultiPop
#' varG(mp)
#' 
#' # Compare with between-population variance
#' varGPop(mp, level = 1)
#'
#' @export
varG = function(pop) {
  calcPopValue(pop, FUN = function(x) {
    G = popVar(x@gv)
    rownames(G) = colnames(G) = colnames(x@gv)
    return(G)
  })
}

#' @title Phenotypic variance within a \code{Pop}
#'
#' @description
#' Returns phenotypic variance for all traits within a \code{Pop}.
#' 
#' In contrast to \code{\link{varPPop}}, which computes variance between
#' population means, \code{varP} computes variance within populations from
#' individual phenotype values.
#'
#' @param pop A \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, or
#'   \code{\link{MultiPop-class}} object.
#' 
#' @details
#' A variance-covariance matrix of phenotype values is computed for each 
#' \code{Pop} across traits.
#' 
#' @seealso \code{\link{varPPop}}
#' 
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or \code{\link{HybridPop-class}},
#' returns a variance-covariance matrix of phenotype values across traits.
#'
#' If \code{pop} is a \code{\link{MultiPop-class}}, returns a nested list
#' matching the \code{MultiPop} structure, with one variance-covariance matrix
#' per \code{Pop} object.
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
#' 
#' # Return a variance-covariance matrix for the population
#' varP(pop)
#' 
#' # Randomly split population into two sub-populations
#' mp = splitPop(pop, sample(c('A','B'), size = 10, replace = TRUE))
#' 
#' # Return a list of variance-covariance matrices for the MultiPop
#' varP(mp)
#' 
#' # Compare with between-population variance
#' varPPop(mp, level = 1)
#'
#' @export
varP = function(pop) {
  calcPopValue(pop, FUN = function(x) {
    G = popVar(x@pheno)
    rownames(G) = colnames(G) = colnames(x@pheno)
    return(G)
  })
}

#' @title Variance of estimated breeding values within a \code{Pop}
#'
#' @description
#' Returns variance of estimated breeding values for all traits within a
#' \code{Pop}.
#'
#' @param pop A \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, or
#'   \code{\link{MultiPop-class}} object.
#' 
#' @details
#' A variance-covariance matrix of estimated breeding values is computed for
#' each \code{Pop} across traits in \code{ebv(pop)}.
#' 
#' This function calculated the variance of values stored in the \code{ebv} slot;
#' supplied by the user or \code{\link{setEBV}}.
#' Depending on how it was populated, it can contain estimated breeding
#' values, estimated genetic values, or other estimated/predicted values.
#'
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or \code{\link{HybridPop-class}},
#' returns a variance-covariance matrix of estimated breeding values across
#' traits in \code{ebv(pop)}.
#'
#' If \code{pop} is a \code{\link{MultiPop-class}}, returns a nested list
#' matching the \code{MultiPop} structure, with one variance-covariance matrix
#' per \code{Pop} object.
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
#' pop@ebv = trtH2 * (pheno(pop) - meanP(pop))
#' 
#' # Return a variance-covariance matrix for the population
#' varEBV(pop)
#' 
#' # Randomly split population into two sub-populations
#' mp = splitPop(pop, sample(c('A','B'), size = 10, replace = TRUE))
#' 
#' # Return a list of variance-covariance matrices for the MultiPop
#' varEBV(mp)
#'
#' @export
varEBV = function(pop) {
  calcPopValue(pop, FUN = function(x) {
    G = popVar(x@ebv)
    rownames(G) = colnames(G) = colnames(x@ebv)
    return(G)
  })
}

#' @title Total genetic variance between populations
#'
#' @description
#' Returns total genetic variance between populations in a
#' \code{\link{MultiPop-class}} object.
#'
#' In contrast to \code{\link{varG}}, which computes variance within each
#' terminal population from individual genetic values, \code{varGPop} computes
#' variance among population means.
#'
#' @param x A \code{\link{MultiPop-class}} object.
#' @param level Integer scalar \eqn{\ge 1} indicating the aggregation level for
#'   population means (see Details).
#'
#' @details
#' This function first computes mean genetic values using
#' \code{\link{meanGPop}} at the requested \code{level}. The resulting matrix
#' (rows = populations, columns = traits) is then passed to
#' \code{\link{popVar}} to calculate a variance-covariance matrix of genetic
#' values across traits, representing variance between populations at the
#' requested \code{level}.
#'
#' If the requested \code{level} yields only one aggregated population, there is
#' no between-population variation to estimate. In this case, the function
#' returns a variance-covariance matrix of zeros and issues a warning.
#' 
#' @seealso \code{\link{varG}}
#'
#' @return
#' A variance-covariance matrix of genetic values across traits, representing
#' variance between populations at the requested aggregation \code{level}.
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
#' # Between-population variance at different aggregation levels
#' varGPop(mp, level = 2)
#' varGPop(mp, level = 1)
#'
#' # Compare with within-population variance
#' varG(mp)
#'
#' @export
varGPop = function(x, level = 1L) {
  m = meanGPop(x, level = level)
  if (is.vector(m)) {
    m = t(m)
  }
  if (nrow(m) == 1L) {
    warning(
      "Returning a variance-covariance matrix of zeros. You may want to increase the value of `level`."
    )
  }
  G = popVar(m)
  rownames(G) = colnames(G) = colnames(m)
  return(G)
}

#' @title Phenotypic variance between populations
#'
#' @description
#' Returns phenotypic variance between populations in a
#' \code{\link{MultiPop-class}} object.
#'
#' In contrast to \code{\link{varP}}, which computes variance within each
#' terminal population from individual phenotype values, \code{varPPop} computes
#' variance among population means.
#'
#' @param x A \code{\link{MultiPop-class}} object.
#' @param level Integer scalar \eqn{\ge 1} indicating the aggregation level for
#'   population means (see Details).
#'
#' @details
#' This function first computes mean phenotype values using
#' \code{\link{meanPPop}} at the requested \code{level}. The resulting matrix
#' (rows = populations, columns = traits) is then passed to
#' \code{\link{popVar}} to calculate a variance-covariance matrix of phenotype
#' values across traits, representing variance between populations at the
#' requested \code{level}.
#'
#' If the requested \code{level} yields only one aggregated population, there is
#' no between-population variation to estimate. In this case, the function
#' returns a variance-covariance matrix of zeros and issues a warning.
#' 
#' @seealso \code{\link{varP}}
#'
#' @return
#' A variance-covariance matrix of phenotype values across traits, representing
#' variance between populations at the requested aggregation \code{level}.
#'
#' @examples
#' founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 10)
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2 = 0.5)
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
#' # Between-population variance at different aggregation levels
#' varPPop(mp, level = 2)
#' varPPop(mp, level = 1)
#'
#' # Compare with within-population variance
#' varP(mp)
#'
#' @export
varPPop = function(x, level = 1L) {
  m = meanPPop(x, level = level)
  if (is.vector(m)) {
    m = t(m)
  }
  if (nrow(m) == 1L) {
    warning(
      "Returning a variance-covariance matrix of zeros. You may want to increase the value of `level`."
    )
  }
  G = popVar(m)
  rownames(G) = colnames(G) = colnames(m)
  return(G)
}

#' @title Genic covariance matrices between traits
#'
#' @description
#' Forms the genic covariances between traits, and their adjustments for
#' departures from Hardy-Weinberg genotype frequencies, from the per-locus
#' tables returned by \code{calcGenParamMulti} with
#' \code{locusTables=TRUE}. Used by \code{\link{genParam}}.
#'
#' @param tabs a list with one element per trait, the \code{tables}
#' returned for that trait
#' @param traits the list of traits, in the same order
#' @param pop the population the tables were calculated in
#'
#' @details
#' A genic covariance is a sum over loci of the product of two traits'
#' statistical values at each genotype, weighted by the genotype
#' frequencies, which belong to the locus rather than to either trait. A
#' trait contributes nothing at a locus where it has no QTL, so traits on
#' disjoint QTL have no genic covariance. Additive-by-additive deviations
#' are defined for pairs of loci, and two traits' deviations only covary
#' where both traits pair the same two loci, because a deviation averages
#' to zero over either locus of its pair. The sums are taken in R, in a
#' fixed order, so that they do not depend on BLAS threading.
#'
#' @return a list of nTraits by nTraits matrices: \code{A}, \code{D} and
#' \code{AA}, the genic covariances under Hardy-Weinberg and linkage
#' equilibrium, and \code{A_HW}, \code{D_HW} and \code{AA_HW}, the genic
#' covariances under the observed marginal genotype frequencies less those
#' under Hardy-Weinberg frequencies
#'
#' @keywords internal
.genicCovMatrices = function(tabs, traits, pop){
  nTraits = length(traits)
  nGeno = pop@ploidy + 1L

  # Genome-wide index of each trait's loci, in the order of its tables
  offset = c(0L, cumsum(pop@nLoci))[seq_len(pop@nChr)]
  locKey = lapply(traits, function(trait){
    chr = rep(seq_along(trait@lociPerChr), trait@lociPerChr)
    offset[chr] + trait@lociLoc
  })
  keys = sort(unique(unlist(locKey)))
  nKey = length(keys)
  take = lapply(locKey, match, table=keys)

  # Genotype frequencies are the same in every trait that has a locus
  freq = freqE = matrix(0, nGeno, nKey)
  for(t in seq_len(nTraits)){
    freq[,take[[t]]] = tabs[[t]]$freq
    freqE[,take[[t]]] = tabs[[t]]$freqE
  }

  # A trait's values spread over every locus, zero where it has no QTL
  spread = function(t, name){
    out = matrix(0, nGeno, nKey)
    if(!is.null(tabs[[t]][[name]])){
      out[,take[[t]]] = tabs[[t]][[name]]
    }
    c(out)
  }

  # Weighted sums of products for every pair of traits
  crossAll = function(cols, w){
    out = matrix(0, nTraits, nTraits)
    for(t in seq_len(nTraits)){
      wt = w*cols[[t]]
      for(s in t:nTraits){
        out[t,s] = out[s,t] = sum(wt*cols[[s]])
      }
    }
    out
  }

  covPair = function(name, nameE){
    obs = crossAll(lapply(seq_len(nTraits), spread, name=name), c(freq))
    hwe = crossAll(lapply(seq_len(nTraits), spread, name=nameE), c(freqE))
    list(genic=hwe, hw=obs-hwe)
  }
  covA = covPair("bv", "bvE")
  covD = covPair("dd", "ddE")

  # Epistatic deviations are keyed by their pair of loci, with the lower
  # index first. A table whose pair runs the other way is transposed, so
  # that the dosage at the lower locus varies fastest in every table.
  pairs = lapply(seq_len(nTraits), function(t){
    if(is.null(tabs[[t]]$pair1)) return(NULL)
    k1 = locKey[[t]][tabs[[t]]$pair1]
    k2 = locKey[[t]][tabs[[t]]$pair2]
    flip = k1 > k2
    aa = tabs[[t]]$aa
    aaE = tabs[[t]]$aaE
    if(any(flip)){
      # Row order of a transposed table
      tr = c(t(matrix(seq_len(nGeno*nGeno), nGeno, nGeno)))
      aa[,flip] = aa[tr,flip,drop=FALSE]
      aaE[,flip] = aaE[tr,flip,drop=FALSE]
    }
    list(lo=pmin(k1,k2), hi=pmax(k1,k2), aa=aa, aaE=aaE)
  })
  hasAA = !vapply(pairs, is.null, logical(1))
  if(any(hasAA)){
    pairId = function(p) paste(p$lo, p$hi)
    allPairs = unique(do.call(rbind, lapply(pairs[hasAA], function(p){
      cbind(p$lo, p$hi)
    })))
    allPairs = allPairs[order(allPairs[,1], allPairs[,2]),,drop=FALSE]
    allId = paste(allPairs[,1], allPairs[,2])
    nCell = nGeno*nGeno
    lo = match(allPairs[,1], keys)
    hi = match(allPairs[,2], keys)
    # Frequencies of the joint genotypes of each pair under linkage
    # equilibrium, from the marginal frequencies of its two loci, with the
    # dosage at the lower locus varying fastest
    rowLo = rep(seq_len(nGeno), times=nGeno)
    rowHi = rep(seq_len(nGeno), each=nGeno)
    wNoLD = freq[rowLo,lo,drop=FALSE]*freq[rowHi,hi,drop=FALSE]
    wNoLDE = freqE[rowLo,lo,drop=FALSE]*freqE[rowHi,hi,drop=FALSE]
    spreadAA = function(t, name){
      out = matrix(0, nCell, nrow(allPairs))
      if(hasAA[t]){
        out[,match(pairId(pairs[[t]]), allId)] = pairs[[t]][[name]]
      }
      c(out)
    }
    obs = crossAll(lapply(seq_len(nTraits), spreadAA, name="aa"), c(wNoLD))
    hwe = crossAll(lapply(seq_len(nTraits), spreadAA, name="aaE"), c(wNoLDE))
    covAA = list(genic=hwe, hw=obs-hwe)
  }else{
    zero = matrix(0, nTraits, nTraits)
    covAA = list(genic=zero, hw=zero)
  }

  return(list(A=covA$genic, A_HW=covA$hw,
              D=covD$genic, D_HW=covD$hw,
              AA=covAA$genic, AA_HW=covAA$hw))
}

#' @title Calculate quantitative genetic quantities
#'
#' @description
#' Calculates quantitative genetic quantities and their variances
#' for an object of \code{\link{Pop-class}}
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See \code{vignette("traits", package = "AlphaSimR")} and
#' the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory on trait and quantitative genetics implementation in AlphaSimR.
#' The text below is a very short summary of the theory.
#'
#' In quantitative genetics and AlphaSimR we use
#' different parameterizations of genetic values.
#' While these are equivalent
#' in the sense that they all parameterize the variation of genetic values,
#' the terms within each parameterization can have different meanings.
#'
#' To generate trait genetic values, AlphaSimR uses a trait model with
#' the genotypic parameterization of genetic values into
#' the trait intercept (\code{gv_mu}; fixed for a given trait definition) and
#' contributions from additive (\code{gv_a}), dominance (\code{gv_d}), and
#' additive-by-additive epistatic (\code{gv_aa}) genotypic effects.
#' For trait \code{t},
#' \code{gv[, t] = gv_mu[t] + gv_a[, t] + gv_d[, t] + gv_aa[, t] =
#' gv_mu[t] + gv_a[, t] + gv_n[, t]}.
#' Here, \code{gv_n = gv_d + gv_aa}. These contributions need not have
#' the population mean of zero. The intercept is calibrated to achieve the
#' requested founder mean and need not equal the population mean genetic value.
#' This genotypic parameterization is also called the "functional" or "biological" model,
#' though it is a simplified representation of the underlying trait biology.
#'
#' To generate classic quantitative genetic quantities, AlphaSimR uses
#' the breeding value parameterization of genetic values into
#' the population mean genetic value (\code{mu}; variable),
#' breeding value (\code{bv}), dominance deviation (\code{dd}), and
#' additive-by-additive epistatic deviation (\code{aa}).
#' For trait \code{t},
#' \code{gv[, t] = mu[t] + bv[, t] + dd[, t] + aa[, t] =
#' mu[t] + bv[, t] + nd[, t]}.
#' Here, \code{nd = dd + aa}.
#' The breeding values and all deviations have mean zero in \code{pop}.
#' The reference quantities are calculated from \code{pop}, so the same
#' individual can have different breeding values and deviations when
#' evaluated in different populations.
#' This breeding value parameterization is also called the "statistical" model
#' because it depends on statistical decomposition of genetic values.
#'
#' Note, the output \code{gv_mu} is the trait intercept \eqn{\mu_0},
#' whereas \code{mu} is the population mean of genetic value \eqn{\bar{G}}.
#' Despite its name, \code{gv_mu} is not the mean of \code{gv}.
#' For the returned quantities, the following identity holds
#' \code{mu = gv_mu + colMeans(gv_a + gv_d + gv_aa)}.
#' After scaling the genotypic effects,
#' the trait intercept is chosen as
#' the requested founder mean minus
#' the mean contribution of the genotypic effects in the founder population.
#' With the trait definition held fixed,
#' the intercept remains constant while the population mean can change.
#'
#' Each \code{covX_HW} is the difference between the corresponding genic
#' variance calculated using observed and Hardy-Weinberg genotype frequencies.
#' Each \code{covX_L}, for \code{X} in \code{G}, \code{A}, \code{D},
#' \code{AA}, and \code{N}, is calculated as \code{diag(varX) - genicVarX - covX_HW}.
#' These are within-trait variance contributions, which can be negative,
#' rather than covariances between traits.
#' The cross-component terms \code{covAD_L}, \code{covAAA_L}, \code{covDAA_L}, and
#' \code{covAN_L} are covariances.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than sample divisor \code{nInd - 1}.
#'
#' The same decomposition is also given between traits, as nTraits by
#' nTraits matrices.
#' \code{genicCovX} is the genic covariance under Hardy-Weinberg and
#' linkage equilibrium at the observed allele frequencies,
#' \code{covX_HWMat} is the adjustment for departures from Hardy-Weinberg
#' genotype frequencies, and \code{covX_LMat} is
#' \code{varX - genicCovX - covX_HWMat}, the contribution of linkage
#' disequilibrium, so the three add up to \code{varX}.
#' Their diagonals are the vectors \code{genicVarX}, \code{covX_HW} and
#' \code{covX_L}.
#' Genic covariances between traits come only from loci that both traits
#' have as QTL, so traits on disjoint QTL have no genic covariance, and any
#' observed covariance between them is due to linkage disequilibrium.
#' Additive-by-additive genic covariances come only from pairs of loci that
#' both traits treat as interacting.
#' These matrices are covariances between traits within one component; the
#' within-trait covariances between components, such as \code{covAD_L},
#' are only given as vectors.
#'
#' @return
#' \describe{
#' \item{varG}{an nTraits by nTraits matrix of total genetic variances and covariances between traits}
#' \item{varA}{an nTraits by nTraits matrix of additive genetic variances and covariances between traits}
#' \item{varD}{an nTraits by nTraits matrix of dominance genetic variances and covariances between traits}
#' \item{varAA}{an nTraits by nTraits matrix of additive-by-additive epistatic genetic variances and covariances between traits}
#' \item{varN}{an nTraits by nTraits matrix of non-additive genetic variances and covariances between traits}
#' \item{genicVarG}{an nTraits vector of total genic variances}
#' \item{genicVarA}{an nTraits vector of additive genic variances}
#' \item{genicVarD}{an nTraits vector of dominance genic variances}
#' \item{genicVarAA}{an nTraits vector of additive-by-additive epistatic genic variances}
#' \item{genicVarN}{an nTraits vector of non-additive genic variances}
#' \item{covG_HW}{an nTraits vector of adjustments to total genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covA_HW}{an nTraits vector of adjustments to additive genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covD_HW}{an nTraits vector of adjustments to dominance genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covAA_HW}{an nTraits vector of adjustments to additive-by-additive epistatic genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covN_HW}{an nTraits vector of adjustments to non-additive genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covG_L}{an nTraits vector of contributions to total genetic variances
#'   due to linkage disequilibrium}
#' \item{covA_L}{an nTraits vector of contributions to additive genetic variances
#'   due to linkage disequilibrium}
#' \item{covD_L}{an nTraits vector of contributions to dominance genetic variances
#'   due to linkage disequilibrium}
#' \item{covAA_L}{an nTraits vector of contributions to additive-by-additive epistatic genetic variances
#'   due to linkage disequilibrium}
#' \item{covN_L}{an nTraits vector of contributions to non-additive genetic variances
#'   due to linkage disequilibrium}
#' \item{covAD_L}{an nTraits vector of within-trait covariances
#'   between breeding values and dominance deviations}
#' \item{covAAA_L}{an nTraits vector of within-trait covariances
#'   between breeding values and additive-by-additive epistatic deviations}
#' \item{covDAA_L}{an nTraits vector of within-trait covariances
#'   between dominance deviations and additive-by-additive epistatic deviations}
#' \item{covAN_L}{an nTraits vector of within-trait covariances
#'   between breeding values and non-additive deviations}
#' \item{genicCovG}{an nTraits by nTraits matrix of total genic variances
#'   and covariances between traits}
#' \item{genicCovA}{an nTraits by nTraits matrix of additive genic variances
#'   and covariances between traits}
#' \item{genicCovD}{an nTraits by nTraits matrix of dominance genic variances
#'   and covariances between traits}
#' \item{genicCovAA}{an nTraits by nTraits matrix of additive-by-additive
#'   epistatic genic variances and covariances between traits}
#' \item{genicCovN}{an nTraits by nTraits matrix of non-additive genic
#'   variances and covariances between traits}
#' \item{covG_HWMat}{an nTraits by nTraits matrix of adjustments to total genic
#'   variances and covariances due to departures from Hardy-Weinberg
#'   genotype frequencies}
#' \item{covA_HWMat}{as \code{covG_HWMat}, for additive genic variances and covariances}
#' \item{covD_HWMat}{as \code{covG_HWMat}, for dominance genic variances and covariances}
#' \item{covAA_HWMat}{as \code{covG_HWMat}, for additive-by-additive epistatic
#'   genic variances and covariances}
#' \item{covN_HWMat}{as \code{covG_HWMat}, for non-additive genic variances and covariances}
#' \item{covG_LMat}{an nTraits by nTraits matrix of contributions to total genetic
#'   variances and covariances due to linkage disequilibrium}
#' \item{covA_LMat}{as \code{covG_LMat}, for additive genetic variances and covariances}
#' \item{covD_LMat}{as \code{covG_LMat}, for dominance genetic variances and covariances}
#' \item{covAA_LMat}{as \code{covG_LMat}, for additive-by-additive epistatic
#'   genetic variances and covariances}
#' \item{covN_LMat}{as \code{covG_LMat}, for non-additive genetic variances and covariances}
#' \item{gv}{a matrix of genetic values with dimensions nInd by nTraits}
#' \item{mu}{an nTraits vector of population mean of genetic values \eqn{\bar{G}} (breeding value parameterization)}
#' \item{mu_HW}{an nTraits vector of genetic value means expected under
#'   Hardy-Weinberg and linkage equilibrium at the observed allele frequencies}
#' \item{bv}{an nInd by nTraits matrix of breeding values (breeding value parameterization)}
#' \item{dd}{an nInd by nTraits matrix of dominance deviations (breeding value parameterization)}
#' \item{aa}{an nInd by nTraits matrix of additive-by-additive epistatic deviations (breeding value parameterization)}
#' \item{nd}{an nInd by nTraits matrix of non-additive deviations (breeding value parameterization)}
#' \item{gv_mu}{an nTraits vector of trait intercepts \eqn{\mu_0} (genotypic parameterization)}
#' \item{gv_a}{an nInd by nTraits matrix of additive genotypic effect
#'   contributions to genetic values (genotypic parameterization)}
#' \item{gv_d}{an nInd by nTraits matrix of dominance genotypic effect
#'   contributions to genetic values (genotypic parameterization)}
#' \item{gv_aa}{an nInd by nTraits matrix of additive-by-additive epistatic genotypic effect
#'   contributions to genetic values (genotypic parameterization)}
#' \item{gv_n}{an nInd by nTraits matrix of non-additive genotypic effects
#'   contributions to genetic values (genotypic parameterization)}
#' \item{alpha}{a list of average effects of allele substitution
#'   calculated using observed marginal genotype frequencies, with length nTraits}
#' \item{alpha_HW}{a list of average effects of allele substitution
#'   calculated using Hardy-Weinberg genotype frequencies at the observed
#'   allele frequencies, with length nTraits}
#' }
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genParam(pop, simParam=SP)
#'
#' @export
genParam = function(pop,simParam=NULL,nThreads=NULL){
  simParam = getSimParam(simParam)
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
  bv = dd = aa = nd = gv_a = gv_d = gv_aa = gv_n = gv

  # Blank nTrait vectors
  genicVarA = rep(NA_real_, nTraits)
  names(genicVarA) = traitNames
  genicVarD = genicVarAA = genicVarN =
    covG_HW = covA_HW = covD_HW = covAA_HW = covN_HW =
    covAD_L = covAAA_L = covDAA_L = covAN_L =
    mu = mu_HW = gv_mu = genicVarA

  # Average effect of an allele substitution
  alpha = vector("list", length=nTraits)
  names(alpha) = traitNames
  alpha_HW = alpha

  # Traits sharing loci, such as correlated traits, are worked out together
  # so their genotypes are only counted and read once. The per-locus tables
  # for genic covariances are only needed when there is more than one trait.
  allTmp = calcGenParamMulti(simParam$traits,pop,nThreads,
                             locusTables=nTraits>1)

  #Loop through trait calculations
  for(i in seq_len(nTraits)){
    trait = simParam$traits[[i]]
    tmp = allTmp[[i]]
    genicVarA[i] = tmp$genicVarA2
    genicVarN[i] = 0
    covA_HW[i] = tmp$genicVarA-tmp$genicVarA2
    covN_HW[i] = 0
    gv[,i] = tmp$gv
    mu[i] = tmp$mu # mean gv for the trait in the pop
    mu_HW[i] = tmp$mu_HWE
    bv[,i] = tmp$bv
    nd[,i] = rep(0,pop@nInd)
    gv_mu[i] = tmp$gv_mu # trait gv intercept
    gv_a[,i] = tmp$gv_a
    gv_n[,i] = rep(0,pop@nInd)
    if(.hasSlot(trait,"domEff")){
      genicVarD[i] = tmp$genicVarD2
      genicVarN[i] = genicVarN[i] + genicVarD[i]
      covD_HW[i] = tmp$genicVarD-tmp$genicVarD2
      covN_HW[i] = covN_HW[i] + covD_HW[i]
      dd[,i] = tmp$dd
      nd[,i] = nd[,i] + dd[,i]
      gv_d[,i] = tmp$gv_d
      gv_n[,i] = gv_n[,i] + gv_d[,i]
    }else{
      genicVarD[i] = 0
      covD_HW[i] = 0
      dd[,i] = rep(0,pop@nInd)
      gv_d[,i] = rep(0,pop@nInd)
    }
    if(.hasSlot(trait,"epiEff")){
      genicVarAA[i] = tmp$genicVarAA2
      genicVarN[i] = genicVarN[i] + genicVarAA[i]
      covAA_HW[i] = tmp$genicVarAA-tmp$genicVarAA2
      covN_HW[i] = covN_HW[i] + covAA_HW[i]
      aa[,i] = tmp$aa
      nd[,i] = nd[,i] + aa[,i]
      gv_aa[,i] = tmp$gv_aa
      gv_n[,i] = gv_n[,i] + gv_aa[,i]
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
      covAN_L[i] = 0
    } else {
      covAD_L[i] = popVar(cbind(bv[,i],dd[,i]))[1,2]
      covAAA_L[i] = popVar(cbind(bv[,i],aa[,i]))[1,2]
      covDAA_L[i] = popVar(cbind(dd[,i],aa[,i]))[1,2]
      covAN_L[i] = popVar(cbind(bv[,i],nd[,i]))[1,2]
    }
    alpha[[i]] = tmp$alpha
    alpha_HW[[i]] = tmp$alpha_HW
  }

  varG = popVar(gv)
  rownames(varG) = colnames(varG) = traitNames

  varA = popVar(bv)
  rownames(varA) = colnames(varA) = traitNames

  varD = popVar(dd)
  rownames(varD) = colnames(varD) = traitNames

  varAA = popVar(aa)
  rownames(varAA) = colnames(varAA) = traitNames

  varN = popVar(nd)
  rownames(varN) = colnames(varN) = traitNames

  genicVarG = genicVarA + genicVarD + genicVarAA
  covG_HW = covA_HW + covD_HW + covAA_HW

  # Genic covariance matrices between traits. Their diagonals are set to the
  # vectors above, which they equal up to rounding, so that the two forms of
  # the output agree exactly.
  if(nTraits>1){
    cm = .genicCovMatrices(lapply(allTmp, function(x) x$tables),
                           simParam$traits, pop)
  }else{
    cm = NULL
  }
  traitMat = function(m, v){
    if(is.null(m)){
      m = matrix(0, nTraits, nTraits)
    }
    diag(m) = v
    rownames(m) = colnames(m) = traitNames
    m
  }
  genicCovA = traitMat(cm$A, genicVarA)
  genicCovD = traitMat(cm$D, genicVarD)
  genicCovAA = traitMat(cm$AA, genicVarAA)
  # Summed in the same order as the vectors, so the diagonals match exactly
  genicCovN = genicCovD + genicCovAA
  genicCovG = genicCovA + genicCovD + genicCovAA
  covA_HWMat = traitMat(cm$A_HW, covA_HW)
  covD_HWMat = traitMat(cm$D_HW, covD_HW)
  covAA_HWMat = traitMat(cm$AA_HW, covAA_HW)
  covN_HWMat = covD_HWMat + covAA_HWMat
  covG_HWMat = covA_HWMat + covD_HWMat + covAA_HWMat

  output = list(varG=varG,
                varA=varA,
                varD=varD,
                varAA=varAA,
                varN=varN,
                genicVarG=genicVarG,
                genicVarA=genicVarA,
                genicVarD=genicVarD,
                genicVarAA=genicVarAA,
                genicVarN=genicVarN,
                covG_HW=covG_HW,
                covA_HW=covA_HW,
                covD_HW=covD_HW,
                covAA_HW=covAA_HW,
                covN_HW=covN_HW,
                covG_L=diag(varG)-genicVarG-covG_HW,
                covA_L=diag(varA)-genicVarA-covA_HW,
                covD_L=diag(varD)-genicVarD-covD_HW,
                covAA_L=diag(varAA)-genicVarAA-covAA_HW,
                covN_L=diag(varN)-genicVarN-covN_HW,
                covAD_L=covAD_L,
                covAAA_L=covAAA_L,
                covDAA_L=covDAA_L,
                covAN_L=covAN_L,
                genicCovG=genicCovG,
                genicCovA=genicCovA,
                genicCovD=genicCovD,
                genicCovAA=genicCovAA,
                genicCovN=genicCovN,
                covG_HWMat=covG_HWMat,
                covA_HWMat=covA_HWMat,
                covD_HWMat=covD_HWMat,
                covAA_HWMat=covAA_HWMat,
                covN_HWMat=covN_HWMat,
                covG_LMat=varG-genicCovG-covG_HWMat,
                covA_LMat=varA-genicCovA-covA_HWMat,
                covD_LMat=varD-genicCovD-covD_HWMat,
                covAA_LMat=varAA-genicCovAA-covAA_HWMat,
                covN_LMat=varN-genicCovN-covN_HWMat,
                gv=gv,
                mu=mu, # mean gv for the trait in the pop
                mu_HW=mu_HW,
                bv=bv,
                dd=dd,
                aa=aa,
                nd=nd,
                gv_mu=gv_mu, # trait gv intercept
                gv_a=gv_a,
                gv_d=gv_d,
                gv_aa=gv_aa,
                gv_n=gv_n,
                alpha=alpha,
                alpha_HW=alpha_HW)
  return(output)
}

#' @title Additive genetic variance
#'
#' @description Returns additive genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' This is the variance of additive genetic (breeding) values in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
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

#' @title Dominance genetic variance
#'
#' @description Returns dominance genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' This is the variance of dominance deviations in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
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

#' @title Additive-by-additive epistatic genetic variance
#'
#' @description Returns additive-by-additive epistatic genetic
#' variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' This is the variance of additive-by-additive epistatic deviations in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
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

#' @title Non-additive genetic variance
#'
#' @description Returns non-additive genetic variance for all traits
#'   (includes dominance and epistatic genetic variance)
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' This is the variance of total non-additive deviations in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#' Since \code{nd = dd + aa}, the variance includes twice the covariance
#' between dominance and additive-by-additive epistatic deviations.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop)
#' varA(pop, simParam=SP)
#' varN(pop, simParam=SP)
#'
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop)
#' varA(pop, simParam=SP)
#' varD(pop, simParam=SP)
#' varAA(pop, simParam=SP)
#' varN(pop, simParam=SP)
#'
#' @export
varN = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varN
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
#' @details
#' Breeding values derive from the breeding value parameterization of genetic values.
#' They have mean zero in \code{pop},
#' which supplies the reference quantities used by \code{\link{genParam}}.
#' Evaluating the same individual in a different population can therefore
#' change its returned values.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return An nInd by nTraits matrix.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
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
#' @details
#' Dominance deviations derive from the breeding value parameterization of genetic values.
#' They have mean zero in \code{pop},
#' which supplies the reference quantities used by \code{\link{genParam}}.
#' Evaluating the same individual in a different population can therefore
#' change returned values.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return An nInd by nTraits matrix.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
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
#' @details
#' Additive-by-additive epistatic deviations derive from the breeding value parameterization of genetic values.
#' They have mean zero in \code{pop},
#' which supplies the reference quantities used by \code{\link{genParam}}.
#' Evaluating the same individual in a different population can therefore
#' change its returned values.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return An nInd by nTraits matrix.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
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

#' @title Non-additive deviations
#'
#' @description Returns non-additive deviations for all traits
#'   (includes dominance and epistatic deviations)
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Non-additive deviations derive from the breeding value parameterization of genetic values.
#' They have mean zero in \code{pop},
#' which supplies the reference quantities used by \code{\link{genParam}}.
#' Evaluating the same individual in a different population can therefore
#' change its returned values.
#' The total non-additive deviation is \code{nd = dd + aa}.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return An nInd by nTraits matrix.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' cbind(dd(pop, simParam=SP),
#'       aa(pop, simParam=SP),
#'       nd(pop, simParam=SP))
#'
#' @export
nd = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$nd
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
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
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
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
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
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
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

#' @title Additive-by-additive epistatic genic variance
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
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
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

#' @title Non-additive genic variance
#'
#' @description Returns non-additive genic variance for all traits
#'   (includes dominance and epistatic genic variance)
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#' genicVarN(pop, simParam=SP)
#'
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#' genicVarD(pop, simParam=SP)
#' genicVarAA(pop, simParam=SP)
#' genicVarN(pop, simParam=SP)
#'
#' @export
genicVarN = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarN
}

#' @title Get genetic values from a population
#'
#' @description
#' Returns genetic values for all traits in a population by extracting the
#' \code{@gv} slot. Supports \code{\link{Pop-class}}, 
#' \code{\link{HybridPop-class}}, and \code{\link{MultiPop-class}} inputs.
#' For \code{MultiPop} objects, output can optionally be simplified to a
#' requested nesting level.
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
#' Simplified output includes a \code{"source"} attribute. This is a data frame
#' with columns \code{level1}, \code{level2}, etc., indicating the name (or
#' index) of the source population for each row at each retained nesting level.
#' 
#' Genetic values include the trait intercept and the contributions of
#' additive, dominance, and additive-by-additive epistatic genotypic effects
#' present in the trait model.
#' For a fixed trait definition, these values do not depend on the population,
#' unlike the breeding values and non-additive deviations.
#' For traits with genotype-by-environment effects, these genetic values
#' refer to the target environment, where the environmental covariate is zero.
#'
#' See \code{vignette("traits", package = "AlphaSimR")} for the trait model.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html}{Quantitative Genetics article} for
#' background theory and its demonstration.
#'
#' See the \href{https://gaynorr.github.io/AlphaSimR/articles/GxE.html}{Genotype-by-Environment article} for
#' the genotype-by-environment model
#'
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or
#' \code{\link{HybridPop-class}}, returns a numeric matrix of genetic
#' values (rows = individuals, columns = traits).
#'
#' If \code{pop} is a \code{\link{MultiPop-class}}:    \cr
#' - with \code{simplify=FALSE}, returns a nested list matching the
#'   \code{MultiPop} structure, with one numeric matrix per terminal
#'   \code{Pop};    \cr
#' - with \code{simplify=TRUE}, returns a simplified object where terminal
#'   \code{Pop} \code{@gv} matrices are combined (within each retained branch)
#'   using \code{\link{rbind}}. A \code{"source"} attribute records row origins.
#'   When \code{level=1}, output is a single numeric matrix.
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
#' gv(pop)
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
#' # Extract genetic values at different simplification levels
#' gv(mp1, simplify = FALSE)
#' gv(mp1, simplify = TRUE, level = 2)
#' gv(mp1, simplify = TRUE, level = 1)
#'
#' @export
gv = function(pop, simplify = FALSE, level = 1L){
  if (simplify && level < 1L) {
    warning("`level` should be >= 1. Setting default `level=1`")
    level = 1L
  }
  calcPopValue(
    pop,
    FUN = function(x) x@gv,
    simplify = simplify,
    level = level
  )
}

#' @title Get phenotype values from a population
#'
#' @description
#' Returns phenotype values for all traits in a population by extracting the
#' \code{@pheno} slot. Supports \code{\link{Pop-class}}, 
#' \code{\link{HybridPop-class}}, and \code{\link{MultiPop-class}} inputs.
#' For \code{MultiPop} objects, output can optionally be simplified to a
#' requested nesting level.
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
#' Simplified output includes a \code{"source"} attribute. This is a data frame
#' with columns \code{level1}, \code{level2}, etc., indicating the name (or
#' index) of the source population for each row at each retained nesting level.
#' 
#' @seealso \code{\link{setPheno}}
#'
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or
#' \code{\link{HybridPop-class}}, returns a numeric matrix of phenotype
#' values (rows = individuals, columns = traits).
#'
#' If \code{pop} is a \code{\link{MultiPop-class}}:    \cr
#' - with \code{simplify=FALSE}, returns a nested list matching the
#'   \code{MultiPop} structure, with one numeric matrix per terminal
#'   \code{Pop};    \cr
#' - with \code{simplify=TRUE}, returns a simplified object where terminal
#'   \code{Pop} \code{@pheno} matrices are combined (within each retained branch)
#'   using \code{\link{rbind}}. A \code{"source"} attribute records row origins.
#'   When \code{level=1}, output is a single numeric matrix.
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
#' pheno(pop)
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
#' # Extract phenotype values at different simplification levels
#' pheno(mp1, simplify = FALSE)
#' pheno(mp1, simplify = TRUE, level = 2)
#' pheno(mp1, simplify = TRUE, level = 1)
#'
#' @export
pheno = function(pop, simplify = FALSE, level = 1L){
  if (simplify && level < 1L) {
    warning("`level` should be >= 1. Setting default `level=1`")
    level = 1L
  }
  calcPopValue(
    pop,
    FUN = function(x) x@pheno,
    simplify = simplify,
    level = level
  )
}

#' @title Get estimated breeding values from a population
#'
#' @description
#' Returns estimated breeding values for all traits in a population by
#' extracting the \code{@ebv} slot.
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
#' Simplified output requires compatible \code{@ebv} columns across terminal
#' populations so that results can be combined with \code{\link{rbind}}.
#'
#' Simplified output includes a \code{"source"} attribute. This is a data frame
#' with columns \code{level1}, \code{level2}, etc., indicating the name (or
#' index) of the source population for each row at each retained nesting level.
#' 
#' The \code{ebv} slot stores predictions supplied by the user or
#' \code{\link{setEBV}}.
#' Depending on how it was populated, it can contain estimated breeding
#' values, estimated genetic values, or other estimated/predicted values.
#'
#' @seealso \code{\link{setEBV}}
#'
#' @return
#' If \code{pop} is a \code{\link{Pop-class}} or
#' \code{\link{HybridPop-class}}, returns a numeric matrix of estimated
#' breeding values (rows = individuals, columns = traits).
#'
#' If \code{pop} is a \code{\link{MultiPop-class}}:    \cr
#' - with \code{simplify=FALSE}, returns a nested list matching the
#'   \code{MultiPop} structure, with one numeric matrix per terminal
#'   \code{Pop};    \cr
#' - with \code{simplify=TRUE}, returns a simplified object where terminal
#'   \code{Pop} \code{@ebv} matrices are combined (within each retained branch)
#'   using \code{\link{rbind}}. A \code{"source"} attribute records row origins.
#'   When \code{level=1}, output is a single numeric matrix.
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
#' pop@ebv = trtH2 * (pheno(pop) - meanP(pop))
#' ebv(pop)
#'
#' # Create a nested MultiPop
#' pop2 = randCross(pop, nCrosses = 3, nProgeny = 4)
#' pop2@ebv = trtH2 * (pheno(pop2) - meanP(pop2))
#' mp1 = splitPop(
#'   pop2,
#'   by = list(
#'     function(x) rep(LETTERS[1:2], length.out = length(x)),
#'     function(x) paste(x@mother, x@father, sep = "_")
#'   )
#' )
#'
#' # Extract estimated breeding values at different simplification levels
#' ebv(mp1, simplify = FALSE)
#' ebv(mp1, simplify = TRUE, level = 2)
#' ebv(mp1, simplify = TRUE, level = 1)
#'
#' @export
ebv = function(pop, simplify = FALSE, level = 1L){
  if (simplify && level < 1L) {
    warning("`level` should be >= 1. Setting default `level=1`")
    level = 1L
  }
  calcPopValue(
    pop,
    FUN = function(x) x@ebv,
    simplify = simplify,
    level = level
  )
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
#' @param use character, calculate using \code{"\link{gv}"},
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
  simParam = getSimParam(simParam)
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
  # Commented out because bv() uses different reference populations for parents and offspring.
  # TODO: See https://github.com/gaynorr/AlphaSimR/issues/291
  # } else if (use == "bv") {
  #   if (!is.null(parents)) {
  #     ret = 0.5 * (bv(parents, simParam = simParam,
  #                     nThreads=nThreads)[matchMothers, , drop = FALSE] +
  #                  bv(parents, simParam = simParam,
  #                     nThreads=nThreads)[matchFathers, , drop = FALSE])
  #   } else {
  #     ret = 0.5 * (bv(mothers, simParam = simParam,
  #                     nThreads=nThreads)[matchMothers, , drop = FALSE] +
  #                  bv(fathers, simParam = simParam,
  #                     nThreads=nThreads)[matchFathers, , drop = FALSE])
  #   }
  } else {
    stop("use must be one of 'gv', 'ebv', or 'pheno'!")
  }
  return(ret)
}

#' @title Calculate Mendelian sampling
#'
#' @param pop \code{\link{Pop-class}} with individuals whose Mendelian samplings
#'   will be calculated
#' @param parents \code{\link{Pop-class}} with mothers and fathers of individuals
#'   in \code{pop}; if \code{NULL} must provide \code{mothers} and \code{fathers}
#' @param mothers \code{\link{Pop-class}} with mothers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param fathers \code{\link{Pop-class}} with fathers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param use character, calculate using \code{"\link{gv}"},
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
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  pa = parentAverage(pop = pop, parents = parents, mothers = mothers, fathers = fathers,
                     use = use, simParam = simParam, nThreads=nThreads)
  if (use %in% c("gv", "ebv", "pheno")) {
    ret = slot(object = pop, name = use) - pa
  # Commented out because bv() uses different reference populations for parents and offspring.
  # TODO: See https://github.com/gaynorr/AlphaSimR/issues/291
  # } else if (use == "bv") {
  #   ret = bv(pop, simParam = simParam, nThreads=nThreads) - pa
  } else {
    stop("use must be one of 'gv', 'ebv', or 'pheno'!")
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

#' @title Number of populations
#'
#' @description
#' Counts the number of \code{\link{Pop-class}} objects in a 
#' \code{\link{MultiPop-class}} object.
#'
#' @param x A \code{\link{Pop-class}} or \code{\link{MultiPop-class}} object.
#'
#' @return
#' Integer scalar giving the number of terminal \code{Pop} objects contained in
#' \code{x}.
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
#' nPop(pop)
#'
#' mp = splitPop(
#'   pop,
#'   by = list(
#'     function(x) rep(LETTERS[1:2], length.out = length(x)),
#'     function(x) sample(c("g1", "g2"), size = length(x), replace = TRUE)
#'   )
#' )
#' nPop(mp)
#'
#' @export
nPop = function (x) {
  if (isPop(x)) return(1L)
  stopifnot(isMultiPop(x))

  count = sum(vapply(x@pops, nPop, integer(1L)))
  return(count)
}