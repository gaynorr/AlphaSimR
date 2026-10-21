#' Returns a vector response from a population
#'
#' @param pop an object of class Pop or HybridPop
#' @param trait a vector or custom function
#' @param use a character ("rand", "gv", "ebv", "pheno", or "bv"; 
#' note that "bv" doesn't work on class HybridPop)
#' @param simParam simulation parameters. Only used when \code{use="bv"}.
#' If \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... are additional arguments passed to trait when trait is a function
#'
#' @keywords internal
getResponse = function(pop,trait,use,simParam=NULL,nThreads=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if(is(trait,"function")){
    if(is.character(use)){
      use = tolower(use)
      if(use=="rand"){
        # A one column matrix, to match what every other option
        # returns. Callers that read ncol() need the dimension.
        return(matrix(rnorm(pop@nInd), ncol=1L))
      }else if(use=="gv"){
        response = trait(pop@gv,...)
      }else if(use=="ebv"){
        response = trait(pop@ebv,...)
      }else if(use=="pheno"){
        response = trait(pop@pheno,...)
      }else if(use=="bv"){
        if(is(pop,"HybridPop")){
          stop("Use='bv' is not a valid option for HybridPop")
        }
        response = genParam(pop,simParam=simParam,nThreads=nThreads)$bv
        response = trait(response,...)
      }else{
        stop(paste0("Use=",use," is not an option"))
      }
    }else if(is(use,"function")){
      response = trait(use(pop, ...), ...)
    }else{
      stop("use must be a character or a function")
    }
  }else{ # trait is not a function, so must be numeric or character
    if(is.character(trait)){ # Suspect trait is a name
      take = match(trait, simParam$traitNames)
      if(any(is.na(take))){
        stop("'",trait[is.na(take)],"' did not match any trait names")
      }
      trait = take
    }
    if(is.character(use)){
      use = tolower(use)
      if(use=="rand"){
        # A one column matrix, to match what every other option
        # returns. Callers that read ncol() need the dimension.
        return(matrix(rnorm(pop@nInd), ncol=1L))
      }else if(use == "gv"){
        response = pop@gv[,trait,drop=FALSE]
      }else if(use=="ebv"){
        response = pop@ebv[,trait,drop=FALSE]
      }else if(use=="pheno"){
        response = pop@pheno[,trait,drop=FALSE]
      }else if(use=="bv"){
        if(is(pop,"HybridPop")){
          stop("Use='bv' is not a valid option for HybridPop")
        }
        response = genParam(pop,simParam=simParam,
                            nThreads=nThreads)$bv[,trait,drop=FALSE]
      }else{
        stop(paste0("Use=",use," is not an option"))
      }
    }else if(is(use,"function")){
      response = use(pop, trait=trait, ...)
    }else{
      stop("use must be a character or a function")
    }
  }
  if(any(is.na(response))){
    stop("selection trait has missing values, phenotype may need to be set")
  }
  return(response)
}

#' @title Calculate values from Pop or MultiPop objects
#'
#' @description
#' Apply a summary function to each \code{\link{Pop-class}} in a
#' \code{\link{MultiPop-class}} and optionally simplify the output to a
#' requested level of nesting.
#'
#' @param x A \code{\link{Pop-class}}, \code{\link{HybridPop-class}} or
#'   \code{\link{MultiPop-class}} object. A \code{MultiPop} can only hold
#'   \code{Pop} and \code{MultiPop} objects, so a \code{HybridPop} is only
#'   accepted on its own.
#' @param FUN Function to apply to each \code{Pop}. Must accept a \code{Pop}
#'   (or \code{HybridPop}, when \code{x} is one) as its first argument.
#' @param simplify Logical. If \code{TRUE}, simplify the output by flattening
#'   the \code{MultiPop} in \code{x} to the requested \code{level} using
#'   \code{\link{flattenMultiPop}}. The corresponding output matrices from
#'   each \code{Pop} are combined with \code{\link{rbind}}, and a
#'   \code{"source"} attribute is added to indicate the origin of each row
#'   (see Details).
#' @param level Integer scalar >= 1. Number of \code{MultiPop} levels to
#'   preserve when \code{simplify=TRUE}. Passed to \code{\link{flattenMultiPop}}.
#'   Ignored if \code{simplify=FALSE}.
#' @param ... Additional arguments passed to \code{FUN}.
#'
#' @details
#' The \code{level} argument controls the depth of nesting retained when
#' \code{simplify=TRUE}. If \code{level} exceeds the depth of \code{x}, the
#' output structure is returned unchanged.
#'
#' When the output is simplified (\code{simplify=TRUE}), a
#' \code{"source"} attribute is attached to the resulting matrix. This
#' attribute is a data frame with columns \code{level1}, \code{level2}, etc.,
#' that record the origin of each row. The values correspond to the name or
#' index of the population from which the row was derived at each nesting
#' level.
#'
#' @return If \code{x} is a \code{Pop}, returns the value of \code{FUN(x, ...)}.
#' Otherwise returns a list (or a simplified matrix) of results, optionally with
#' a \code{"source"} attribute.
#'
#' @seealso \code{\link{MultiPop-class}}, \code{\link{flattenMultiPop}}
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=6, nChr=1, segSites=10)
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10, mean = c(0, 0), var = c(1, 1))
#' SP$setVarE(h2 = c(0.6, 0.4))
#'
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
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
#' # Calculate phenotypes at different levels
#' calcPopValue(mp1, FUN = pheno, simplify = TRUE, level = 2)
#' calcPopValue(mp1, FUN = pheno, simplify = TRUE, level = 1)
#'
#' # Custom function returning a summary matrix
#' calcPopValue(mp1, FUN = function(x) colMeans(pheno(x)),
#'              simplify = TRUE, level = 1)
#'
#' @export
calcPopValue = function(
  x,
  FUN,
  simplify = FALSE,
  level = 1,
  ...
) {
  dots = list(...)
  .paths = dots[[".paths"]]
  .level_offset = dots[[".level_offset"]]
  if (is.null(.level_offset)) {
    .level_offset = 0L
  }
  dots[[".level_offset"]] = NULL
  dots[[".paths"]] = NULL

  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level != trunc(level) ||
      abs(level) > .Machine$integer.max) {
    stop("`level` must be a single non-NA integer value.")
  }
  if (!is.function(FUN)) {
    stop("`FUN` must be a function.")
  }

  if (isPop(x) || isHybridPop(x)) {
    FUN.ARGS = c(list(x), dots)
    return(do.call(FUN, FUN.ARGS))
  }

  if (!isMultiPop(x)) {
    stop("`x` must be a Pop, HybridPop or MultiPop object.")
  }

  if (is.null(.paths) && simplify) {
    .paths = .collectLeafPaths(x)
  }

  src = NULL
  if (simplify && level == 1L) {
    src = .paths
  }

  if (simplify) {
    x = flattenMultiPop(x, level = level)
  }

  popValueList = lapply(
    seq_along(x@pops),
    function(i) {
      child = x@pops[[i]]

      # Only a nested MultiPop uses the leaf paths, so filtering them for a
      # terminal Pop would scan every path once per child for nothing.
      childPaths = NULL
      if (isMultiPop(child)) {
        nm = names(x@pops)
        valid_name = !is.null(nm) && !is.na(nm[i]) && nzchar(nm[i])
        unique_name = valid_name && sum(nm == nm[i], na.rm = TRUE) == 1L
        child_label = if (unique_name) nm[i] else i
        keep = vapply(
          .paths,
          function(p) {
            length(p) > 0 && identical(p[[1]], child_label)
          },
          logical(1)
        )
        childPaths = lapply(.paths[keep], function(p) p[-1])
      }

      FUN.ARGS = c(
        list(
          child,
          FUN = FUN,
          simplify = simplify,
          level = level - 1L,
          .level_offset = .level_offset + 1L,
          .paths = childPaths
        ),
        dots
      )
      do.call(calcPopValue, FUN.ARGS)
    }
  )

  if (simplify && level == 1L) {
    nRows = vapply(
      popValueList,
      function(v) {
        if (is.null(v)) {
          return(0L)
        }
        if (is.data.frame(v) || is.matrix(v) || is.array(v)) {
          return(dim(v)[1]) # number of rows
        }
        if (is.atomic(v) && is.null(dim(v))) {
          return(1L) # vectors are one row
        }
        return(NA_integer_)
      },
      integer(1L)
    )
    if (any(is.na(nRows))) {
      warning(
        "Some values returned by FUN have unsupported types for simplification. Returning list output.",
        call. = FALSE
      )
      if (!is.null(names(x@pops))) {
        names(popValueList) = names(x@pops)
      }
      return(popValueList)
    }

    # Check consistent column names across matrices before combining
    same_columns = all(vapply(
      popValueList,
      function(x) {
        reference = popValueList[[1]]
        x_names = if (is.atomic(x) && is.null(dim(x))) names(x) else colnames(x)
        reference_names = if (is.atomic(reference) && is.null(dim(reference))) names(reference) else colnames(reference)
        identical(x_names, reference_names)
      },
      logical(1)
    ))
    if (!same_columns) {
      warning("Some values returned by FUN do not have consistent column names. Returning list output.")
      if (!is.null(names(x@pops))) {
        names(popValueList) = names(x@pops)
      }
      return(popValueList)
    }

    popValue = do.call(rbind, popValueList)
    attr(popValue, "source") = .formatPopSource(
      paths = src,
      nRows = nRows,
      level_offset = .level_offset
    )
    return(popValue)
  }

  if (!is.null(names(x@pops))) {
    names(popValueList) = names(x@pops)
  }

  return(popValueList)
}

#' @title Calculate summaries across a population hierarchy
#'
#' @description
#' Calculates summaries for \code{\link{Pop-class}} objects and recursively
#' aggregates them across a nested \code{\link{MultiPop-class}} object.
#'
#' @param x A \code{\link{Pop-class}}, \code{\link{HybridPop-class}} or
#'   \code{\link{MultiPop-class}} object. A \code{MultiPop} can only hold
#'   \code{Pop} and \code{MultiPop} objects, so a \code{HybridPop} is only
#'   accepted on its own.
#' @param withinPopFUN a function used to calculate a summary for each
#'   \code{Pop} (or \code{HybridPop}, when \code{x} is one) object. The
#'   function must accept a \code{Pop} object as its first argument and 
#'   return a numeric vector or one-row matrix.
#' @param withinPopArgs a list of additional arguments passed to
#'   \code{withinPopFUN}
#' @param acrossPopFUN a function used to aggregate summaries across the
#'   immediate children of each \code{MultiPop} object. The function must
#'   accept a numeric matrix as its first argument and return a numeric vector
#'   or one-row matrix with the same columns.
#' @param acrossPopArgs a list of additional arguments passed to
#'   \code{acrossPopFUN}
#' @param level an integer indicating the nesting level at which aggregation
#'   stops. See Details.
#'
#' @details
#' This function performs two operations. First, \code{withinPopFUN} is applied
#' to every terminal \code{Pop} object. Second, the resulting summaries are
#' recursively aggregated with \code{acrossPopFUN} through the nested
#' \code{MultiPop} structure.
#'
#' The rows of the matrix supplied to \code{acrossPopFUN} correspond to the
#' immediate children of a \code{MultiPop}. Aggregation is therefore unweighted
#' with respect to the number of individuals in each child population, unless
#' weighting is explicitly implemented by the user-supplied function.
#'
#' The \code{level} argument determines where recursive aggregation stops:   \cr
#'   - \code{level = 0} aggregates summaries across the complete object.    \cr
#'   - \code{level = 1} returns one summary for each immediate child of
#'     \code{x}.    \cr
#'   - Larger values return summaries from progressively deeper levels of
#'     the hierarchy.
#'
#' If a branch ends in a \code{Pop} object before the requested level is
#' reached, the summary for that terminal \code{Pop} is returned. Consequently,
#' levels greater than the maximum depth return one summary for every terminal
#' \code{Pop} object.
#'
#' The outputs of \code{withinPopFUN} must have consistent lengths and names.
#' The output of \code{acrossPopFUN} must preserve the number and names of
#' these summary values. Missing values are passed through unchanged, so
#' how they propagate is decided by the supplied functions. For example,
#' \code{colMeans} returns \code{NA} for a column containing a missing
#' value and leaves the other columns unaffected.
#'
#' @return If \code{x} is a \code{Pop} object, a numeric vector containing its
#'   summary.
#'
#' If \code{x} is a \code{MultiPop} object and \code{level = 0}, a numeric
#' vector containing the summary aggregated across the complete object.
#'
#' If \code{x} is a \code{MultiPop} object and \code{level > 0}, a numeric
#' matrix with one row per population unit at the requested level. The matrix
#' has a \code{"source"} attribute describing the origin of each row in the
#' population hierarchy.
#'
#' @examples
#' founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 20)
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(nQtlPerChr = 10)
#' SP$setVarE(h2 = 0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' pop = newPop(founderPop, simParam = SP)
#'
#' # Create a nested MultiPop
#' multiPop = splitPop(
#'   pop,
#'   by = list(
#'     rep(c("A", "B"), each = 8),
#'     function(x) {
#'       rep(c("one", "two"), length.out = nInd(x))
#'     }
#'   )
#' )
#'
#' # Calculate mean phenotypic values for each terminal Pop
#' popMeans = calcPopSummary(
#'   multiPop,
#'   withinPopFUN = function(x) colMeans(pheno(x)),
#'   acrossPopFUN = colMeans,
#'   level = 2
#' )
#' popMeans
#' attr(popMeans, "source")
#'
#' # Aggregate terminal Pop means within each level-1 MultiPop
#' groupMeans = calcPopSummary(
#'   multiPop,
#'   withinPopFUN = function(x) colMeans(pheno(x)),
#'   acrossPopFUN = colMeans,
#'   level = 1
#' )
#' groupMeans
#'
#' # Aggregate across the complete hierarchy
#' overallMean = calcPopSummary(
#'   multiPop,
#'   withinPopFUN = function(x) colMeans(pheno(x)),
#'   acrossPopFUN = colMeans,
#'   level = 0
#' )
#' overallMean
#'
#' @seealso \code{\link{selectMultiPop}}, \code{\link{meanGPop}},
#'   \code{\link{meanPPop}}
#'
#' @export
calcPopSummary = function(
  x,
  withinPopFUN,
  withinPopArgs = list(),
  acrossPopFUN,
  acrossPopArgs = list(),
  level = 0
) {
  # Validate the public arguments once, before starting recursion
  if (!isPop(x) && !isHybridPop(x) && !isMultiPop(x)) {
    stop("`x` must be a Pop, HybridPop or MultiPop object.")
  }
  if (length(x) == 0L) {
    stop("`x` contains no populations.")
  }
  if (!is.function(withinPopFUN)) {
    stop("`withinPopFUN` must be a function.")
  }
  if (!is.function(acrossPopFUN)) {
    stop("`acrossPopFUN` must be a function.")
  }
  if (!is.list(withinPopArgs)) {
    stop("`withinPopArgs` must be a list.")
  }
  if (!is.list(acrossPopArgs)) {
    stop("`acrossPopArgs` must be a list.")
  }
  if (
    !is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      !is.finite(level) ||
      level != trunc(level) ||
      level < 0 ||
      level > .Machine$integer.max
  ) {
    stop("`level` must be a single non-negative integer.")
  }

  level = as.integer(level)
  result = .calcPopSummary(
    x,
    levelsRemaining = level,
    withinPopFUN = withinPopFUN,
    withinPopArgs = withinPopArgs,
    acrossPopFUN = acrossPopFUN,
    acrossPopArgs = acrossPopArgs
  )

  # At level 0, recursion returns one summary vector for the complete object.
  if (level == 0L || isPop(x) || isHybridPop(x)) {
    return(result)
  }

  # At positive levels, recursion returns one row for each Pop or MultiPop
  # unit at the requested level, or for each terminal Pop when a branch ends
  # before that level.
  maxDepth = .depthMultiPop(x)
  sourceDepth = min(level, maxDepth)

  leafPaths = .collectLeafPaths(x)
  summaryPaths = unique(lapply(leafPaths, function(path) {
    path[seq_len(min(length(path), sourceDepth))]
  }))
  src = .formatPopSource(
    paths = summaryPaths,
    nRows = rep(1L, length(summaryPaths)),
    level_offset = 0L
  )
  rownames(src) = NULL

  if (nrow(src) != nrow(result)) {
    stop(
      "Internal error: the number of population summaries does not ",
      "match the number of source paths."
    )
  }
  attr(result, "source") = src
  rownames(result) = NULL

  return(result)
}

#' Identify candidate individuals
#' 
#' This function handles indexing by id and negative value indexing
#'
#' @param pop a population object
#' @param candidates a vector of candidates (numeric or character)
#'
#' @keywords internal
getCandidates = function(pop, candidates){
  if(is.character(candidates)){
    candidates = match(candidates, pop@id)
    if(any(is.na(candidates))){
      stop("Trying to select invalid individuals")
    }
    if(any(is.null(candidates))){
      stop("Not valid ids")
    }
  }else{
    if(any(abs(candidates)>pop@nInd)){
      stop("Trying to select invalid individuals")
    }
    candidates = (1:pop@nInd)[candidates]
  }
  return(candidates)
}

#' Find individuals of desired sex
#' 
#' This function ignores sex if sexes is "no" in SimParam.
#'
#' @param pop a population
#' @param sex the desired sex (M or F)
#' @param simParam simulation parameters. If \code{NULL}, the function uses
#' the object named \code{SP} from the global environment.
#' @param ...  captures use of old gender argument
#'
#' @keywords internal
checkSexes = function(pop,sex,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  sex = toupper(sex)
  eligible = 1:pop@nInd
  if(simParam$sexes=="no"){
    return(eligible)
  }else{
    # Check in gender is incorrectly being used
    args = list(...)
    if(any(names(args)=="gender")){
      stop("The discontinued 'gender' argument appears to be in use. This argument was renamed as 'sex' in AlphaSimR version 0.13.0.")
    }
    if(sex=="B"){
      return(eligible)
    }else if(sex=="F" | sex=="M"){
      return(eligible[pop@sex%in%sex])
    }else{
      stop("Invalid option supplied for 'sex'")
    }
  }
}

#' Check the onShortfall argument
#'
#' @param onShortfall what to do when fewer candidates are available than
#' requested: "warning", "error" or "none"
#'
#' @return the validated value
#'
#' @keywords internal
checkShortfall = function(onShortfall){
  match.arg(onShortfall, c("warning","error","none"))
}

#' Report a selection shortfall
#'
#' Raises a warning or an error, or does nothing, according to
#' onShortfall. The condition carries the call of the selection function
#' rather than of this helper, so the message names the function the user
#' called.
#'
#' @param msg the message for a warning, which says what is returned
#' @param onShortfall "warning", "error" or "none"
#' @param errorMsg the message for an error, which must not say that
#' anything is returned
#' @param call the call to report
#'
#' @keywords internal
reportShortfall = function(msg, onShortfall, errorMsg=msg,
                           call=sys.call(-1)){
  if(onShortfall=="warning"){
    warning(simpleWarning(msg, call))
  }else if(onShortfall=="error"){
    stop(simpleError(errorMsg, call))
  }
  invisible(NULL)
}

#' Stop if a population has no individuals
#'
#' The selection functions index candidates as 1:nInd, which is c(1, 0)
#' for an empty population, so without this check they fail later with
#' messages that say nothing about the cause.
#'
#' @param pop a population
#' @param call the call to report
#'
#' @keywords internal
checkNotEmpty = function(pop, call=sys.call(-1)){
  if(pop@nInd==0L){
    stop(simpleError("pop contains no individuals to select from", call))
  }
  invisible(NULL)
}

#' Is a population an empty member of a MultiPop
#'
#' A MultiPop may hold empty populations, which the selection functions
#' pass through unchanged rather than treating as an error.
#'
#' @param pop an element of a MultiPop
#'
#' @keywords internal
isEmptyMember = function(pop){
  !is(pop,"MultiPop") && pop@nInd==0L
}

#' Determine families
#'
#' @param pop a population object
#' @param famType type of family (fullsib "B", maternal "M", or paternal "F")
#'
#' @keywords internal
getFam = function(pop,famType){
  famType = toupper(famType)
  if(famType=="B"){
    return(paste(pop@mother,pop@father,sep="_"))
  }else if(famType=="F"){
    return(pop@mother)
  }else if(famType=="M"){
    return(pop@father)
  }else{
    stop(paste0("famType=",famType," is not a valid option"))
  }
}

#' @title Select individuals
#'
#' @description Selects a subset of nInd individuals from a
#' population.
#'
#' @param pop an object of \code{\link{Pop-class}},
#'   \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#' @param nInd the number of individuals to select
#' @param trait the trait for selection. Either a number indicating
#'   a single trait or a function returning a vector of length nInd.
#'   The function must work on a vector or matrix of \code{use} values as
#'   \code{trait(pop@use, ...)} - depending on what \code{use} is.
#'   See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values "gv", estimated breeding values "ebv", breeding values "bv",
#'   phenotypes "pheno", or randomly "rand") or
#'   a function returning a vector of length nInd.
#'   The function must work on \code{pop} as \code{use(pop, trait, ...)} or
#'   as \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param sex which sex to select. Use "B" for both, "F" for
#'   females and "M" for males. If the simulation is not using sexes,
#'   the argument is ignored.
#' @param selectTop selects highest values if true.
#'   Selects lowest values if false.
#' @param returnPop should results be returned as a
#'   \code{\link{Pop-class}}. If FALSE, only the index of selected
#'   individuals is returned.
#' @param candidates an optional vector of eligible selection candidates.
#' @param onShortfall what to do when fewer suitable candidates are
#'   available than requested. "warning" (the default) warns and returns
#'   all of them, "error" stops, and "none" returns them silently.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#'   \code{trait} or \code{use}
#'
#' @return Returns an object of \code{\link{Pop-class}},
#' \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#'
#' @family selection functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Select top 5 (directional selection)
#' pop2 = selectInd(pop, 5, simParam=SP)
#' hist(pop@pheno); abline(v=pop@pheno, lwd=2)
#' abline(v=pop2@pheno, col="red", lwd=2)
#'
#' #Select 5 most deviating from an optima (disruptive selection)
#' squaredDeviation = function(x, optima=0) (x - optima)^2
#' pop3 = selectInd(pop, 5, trait=squaredDeviation, selectTop=TRUE, simParam=SP)
#' hist(pop@pheno); abline(v=pop@pheno, lwd=2)
#' abline(v=pop3@pheno, col="red", lwd=2)
#'
#' #Select 5 least deviating from an optima (stabilizing selection)
#' pop4 = selectInd(pop, 5, trait=squaredDeviation, selectTop=FALSE, simParam=SP)
#' hist(pop@pheno); abline(v=pop@pheno, lwd=2)
#' abline(v=pop4@pheno, col="red", lwd=2)
#'
#' #Select 5 individuals based on miscelaneous information with use function
#' pop@misc = list(smth=rnorm(10), smth2=rnorm(10))
#' useFunc = function(pop, trait=NULL) pop@misc$smth + pop@misc$smth2
#' pop5 = selectInd(pop, 5, use=useFunc, simParam=SP)
#' pop5@id
#'
#' #... equivalent result with the use & trait function
#' useFunc2 = function(pop, trait=NULL) cbind(pop@misc$smth, pop@misc$smth2)
#' trtFunc = function(x) rowSums(x)
#' pop6 = selectInd(pop, 5, trait=trtFunc, use=useFunc2, simParam=SP)
#' pop6@id
#'
#' @export
selectInd = function(pop,nInd,trait=1,use="pheno",sex="B",
                     selectTop=TRUE,returnPop=TRUE,
                     candidates=NULL,simParam=NULL,nThreads=NULL,
                     onShortfall="warning",...){
  if(nInd<0) stop("nInd must be >= 0")
  onShortfall = checkShortfall(onShortfall)
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if(is(pop,"MultiPop")){
    if(!returnPop) stop("returnPop must be TRUE for a MultiPop")
    if(!is.null(candidates)) stop("candidates must be NULL for a MultiPop")
    pop@pops = lapply(pop@pops, function(x, ...){
      if(isEmptyMember(x)) return(x)
      selectInd(x, nInd=nInd, trait=trait,
                use=use, sex=sex, selectTop=selectTop,
                returnPop=TRUE, candidates=NULL,
                simParam=simParam, nThreads=nThreads,
                onShortfall=onShortfall, ...)
    }, ...)
    return(pop)
  }
  checkNotEmpty(pop)
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
    candidates = getCandidates(pop=pop,candidates=candidates)
    eligible = eligible[eligible%in%candidates]
  }
  if(length(eligible)<nInd){
    reportShortfall(paste0("Suitable candidates smaller than nInd, returning ",
                           length(eligible)," individuals"),
                    onShortfall,
                    paste0("Only ", length(eligible),
                           " suitable candidates for nInd = ", nInd))
    nInd = length(eligible)
  }
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,nThreads=nThreads,...)
  if(is.matrix(response)){
    if(ncol(response)!=1) stop("response must have a single column")
  }
  take = order(response,decreasing=selectTop)
  take = take[take%in%eligible]
  if(returnPop){
    return(pop[take[0:nInd]])
  }else{
    return(take[0:nInd])
  }
}

#' @title Select families
#'
#' @description Selects a subset of full-sib families from a
#' population.
#'
#' @param pop an object of \code{\link{Pop-class}},
#'   \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#' @param nFam the number of families to select
#' @param trait the trait for selection. Either a number indicating
#'   a single trait or a function returning a vector of length nInd.
#'   The function must work on a vector or matrix of \code{use} values as
#'   \code{trait(pop@use, ...)} - depending on what \code{use} is.
#'   See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values "gv", estimated breeding values "ebv", breeding values "bv",
#'   phenotypes "pheno", or randomly "rand") or
#'   a function returning a vector of length nInd.
#'   The function must work on \code{pop} as \code{use(pop, trait, ...)} or
#'   as \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param sex which sex to select. Use "B" for both, "F" for
#'   females and "M" for males. If the simulation is not using sexes,
#'   the argument is ignored.
#' @param famType which type of family to select. Use "B" for
#'   full-sib families, "F" for half-sib families on female side and "M"
#'   for half-sib families on the male side.
#' @param selectTop selects highest values if true.
#'   Selects lowest values if false.
#' @param returnPop should results be returned as a
#'   \code{\link{Pop-class}}. If FALSE, only the index of selected
#'   individuals is returned.
#' @param candidates an optional vector of eligible selection candidates.
#' @param onShortfall what to do when fewer suitable candidates are
#'   available than requested. "warning" (the default) warns and returns
#'   all of them, "error" stops, and "none" returns them silently.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#'   \code{trait} and \code{use}
#'
#' @return Returns an object of \code{\link{Pop-class}},
#' \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#'
#' @family selection functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Create 3 biparental families with 10 progeny
#' pop2 = randCross(pop, nCrosses=3, nProgeny=10, simParam=SP)
#'
#' #Select best 2 families
#' pop3 = selectFam(pop2, 2, simParam=SP)
#'
#' @export
selectFam = function(pop,nFam,trait=1,use="pheno",sex="B",
                     famType="B",selectTop=TRUE,returnPop=TRUE,
                     candidates=NULL,simParam=NULL,nThreads=NULL,
                     onShortfall="warning",...){
  if(nFam<0) stop("nFam must be >= 0")
  onShortfall = checkShortfall(onShortfall)
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if(is(pop,"MultiPop")){
    if(!returnPop) stop("returnPop must be TRUE for a MultiPop")
    if(!is.null(candidates)) stop("candidates must be NULL for a MultiPop")
    pop@pops = lapply(pop@pops, function(x, ...){
      if(isEmptyMember(x)) return(x)
      selectFam(x, nFam=nFam, trait=trait,
                use=use, sex=sex, famType=famType,
                selectTop=selectTop, returnPop=TRUE,
                candidates=NULL, simParam=simParam,
                nThreads=nThreads, onShortfall=onShortfall, ...)
    }, ...)
    return(pop)
  }
  checkNotEmpty(pop)
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
    candidates = getCandidates(pop=pop,candidates=candidates)
    eligible = eligible[eligible%in%candidates]
  }
  allFam = getFam(pop=pop,famType=famType)
  availFam = allFam[eligible]
  if(nFam>length(unique(availFam))){
    reportShortfall(paste0("Suitable families smaller than nFam, returning ",
                           length(unique(availFam)), " families"),
                    onShortfall,
                    paste0("Only ", length(unique(availFam)),
                           " suitable families for nFam = ", nFam))
    nFam = length(unique(availFam))
  }
  # aggregate fails on zero rows, and there is nothing to rank
  if(length(eligible)==0L){
    if(returnPop){
      return(pop[integer(0)])
    }else{
      return(integer(0))
    }
  }
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,nThreads=nThreads,...)
  if(is.matrix(response)){
    if(ncol(response)!=1) stop("response must have a single column")
  }
  response = response[eligible]
  #Calculate family means
  famMeans = aggregate(response,list(families=availFam),mean)
  response = famMeans$x
  #Select families
  bestFam = order(response,decreasing=selectTop)[0:nFam]
  bestFam = famMeans$families[bestFam]
  take = which(allFam%in%bestFam)
  take = take[take%in%eligible]
  if(returnPop){
    return(pop[take])
  }else{
    return(take)
  }
}

#' @title Select individuals within families
#'
#' @description Selects a subset of nInd individuals from each
#' full-sib family within a population. Will return all individuals
#' from a full-sib family if it has less than or equal to nInd individuals.
#'
#' @param pop an object of \code{\link{Pop-class}},
#'   \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#' @param nInd the number of individuals to select within a family
#' @param trait the trait for selection. Either a number indicating
#'   a single trait or a function returning a vector of length nInd.
#'   The function must work on a vector or matrix of \code{use} values as
#'   \code{trait(pop@use, ...)} - depending on what \code{use} is.
#'   See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values "gv", estimated breeding values "ebv", breeding values "bv",
#'   phenotypes "pheno", or randomly "rand") or
#'   a function returning a vector of length nInd.
#'   The function must work on \code{pop} as \code{use(pop, trait, ...)} or
#'   as \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param sex which sex to select. Use "B" for both, "F" for
#'   females and "M" for males. If the simulation is not using sexes,
#'   the argument is ignored.
#' @param famType which type of family to select. Use "B" for
#'   full-sib families, "F" for half-sib families on female side and "M"
#'   for half-sib families on the male side.
#' @param selectTop selects highest values if true.
#'   Selects lowest values if false.
#' @param returnPop should results be returned as a
#'   \code{\link{Pop-class}}. If FALSE, only the index of selected
#'   individuals is returned.
#' @param candidates an optional vector of eligible selection candidates.
#' @param onShortfall what to do when a family has fewer suitable
#'   candidates than \code{nInd}. "warning" (the default) warns and
#'   returns all of that family's candidates, "error" stops, and "none"
#'   returns them silently. Use "none" when \code{nInd} is meant as a cap
#'   on family size, so that small families are expected.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#'   \code{trait} and \code{use}
#'
#' @return Returns an object of \code{\link{Pop-class}},
#' \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#'
#' @family selection functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Create 3 biparental families with 10 progeny
#' pop2 = randCross(pop, nCrosses=3, nProgeny=10, simParam=SP)
#'
#' #Select best individual per family
#' pop3 = selectWithinFam(pop2, 1, simParam=SP)
#'
#' @export
selectWithinFam = function(pop,nInd,trait=1,use="pheno",sex="B",
                           famType="B",selectTop=TRUE,returnPop=TRUE,
                           candidates=NULL,simParam=NULL,nThreads=NULL,
                           onShortfall="warning",...){
  if(nInd<0) stop("nInd must be >= 0")
  onShortfall = checkShortfall(onShortfall)
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if(is(pop,"MultiPop")){
    if(!returnPop) stop("returnPop must be TRUE for a MultiPop")
    if(!is.null(candidates)) stop("candidates must be NULL for a MultiPop")
    pop@pops = lapply(pop@pops, function(x, ...){
      if(isEmptyMember(x)) return(x)
      selectWithinFam(x, nInd=nInd, trait=trait,
                      use=use, sex=sex, famType=famType,
                      selectTop=selectTop, returnPop=TRUE,
                      candidates=NULL, simParam=simParam,
                      nThreads=nThreads, onShortfall=onShortfall, ...)
    }, ...)
    return(pop)
  }
  checkNotEmpty(pop)
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
    candidates = getCandidates(pop=pop,candidates=candidates)
    eligible = eligible[eligible%in%candidates]
  }
  families = getFam(pop=pop,famType=famType)
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,nThreads=nThreads,...)
  if(is.matrix(response)){
    if(ncol(response)!=1) stop("response must have a single column")
  }
  nShort = 0L
  selInFam = function(selFam){
    index = which(families%in%selFam)
    y = response[index]
    index = index[order(y,decreasing=selectTop)]
    index = index[index%in%eligible]
    if(length(index)<nInd){
      nShort <<- nShort + 1L
      return(index)
    }else{
      return(index[0:nInd])
    }
  }
  take = unlist(sapply(unique(families),selInFam))
  if(nShort>0L){
    reportShortfall(paste0(nShort, " of ", length(unique(families)),
                           " families have fewer suitable candidates than ",
                           "nInd; all of their candidates are returned"),
                    onShortfall,
                    paste0(nShort, " of ", length(unique(families)),
                           " families have fewer suitable candidates than ",
                           "nInd = ", nInd))
  }
  if(returnPop){
    return(pop[take])
  }else{
    return(take)
  }
}

#' @title Select open pollinating plants
#'
#' @description
#' This function models selection in an open pollinating
#' plant population. It allows for varying the percentage of
#' selfing. The function also provides an option for modeling
#' selection as occurring before or after pollination.
#'
#' @param pop an object of \code{\link{Pop-class}}
#'   or \code{\link{MultiPop-class}}
#' @param nInd the number of plants to select
#' @param nSeeds number of seeds per plant
#' @param probSelf percentage of seeds expected from selfing.
#'   Value ranges from 0 to 1.
#' @param pollenControl are plants selected before pollination
#' @param trait the trait for selection. Either a number indicating
#'   a single trait or a function returning a vector of length nInd.
#'   The function must work on a vector or matrix of \code{use} values as
#'   \code{trait(pop@use, ...)} - depending on what \code{use} is.
#'   See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values "gv", estimated breeding values "ebv", breeding values "bv",
#'   phenotypes "pheno", or randomly "rand") or
#'   a function returning a vector of length nInd.
#'   The function must work on \code{pop} as \code{use(pop, trait, ...)} or
#'   as \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param selectTop selects highest values if true.
#'   Selects lowest values if false.
#' @param candidates an optional vector of eligible selection candidates.
#' @param onShortfall what to do when fewer suitable candidates are
#'   available than requested. "warning" (the default) warns and returns
#'   all of them, "error" stops, and "none" returns them silently.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#'   \code{trait} and \code{use}
#'
#' @return Returns an object of \code{\link{Pop-class}}
#' or \code{\link{MultiPop-class}}
#'
#' @family mating functions
#' @family selection functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Create new population by selecting the best 3 plant
#' #Assuming 50% selfing in plants and 10 seeds per plant
#' pop2 = selectOP(pop, nInd=3, nSeeds=10, probSelf=0.5, simParam=SP)
#'
#' @export
selectOP = function(pop,nInd,nSeeds,probSelf=0,
                    pollenControl=FALSE,trait=1,
                    use="pheno",selectTop=TRUE,
                    candidates=NULL,simParam=NULL,nThreads=NULL,
                    onShortfall="warning",...){
  if(nInd<0) stop("nInd must be >= 0")
  onShortfall = checkShortfall(onShortfall)
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if(is(pop,"MultiPop")){
    if(!is.null(candidates)) stop("candidates must be NULL for a MultiPop")
    pop@pops = lapply(pop@pops, function(x, ...){
      if(isEmptyMember(x)) return(x)
      selectOP(x, nInd=nInd, nSeeds=nSeeds,
               probSelf=probSelf,
               pollenControl=pollenControl, trait=trait, use=use,
               selectTop=selectTop, candidates=NULL,
               simParam=simParam, nThreads=nThreads,
               onShortfall=onShortfall, ...)
    }, ...)
    return(pop)
  }
  checkNotEmpty(pop)
  female = selectInd(pop=pop,nInd=nInd,trait=trait,
                     use=use,sex="B",selectTop=selectTop,
                     returnPop=FALSE,candidates=candidates,
                     simParam=simParam,nThreads=nThreads,
                     onShortfall=onShortfall,...)
  # selectInd returns at most as many individuals as there are eligible
  # candidates, which can be fewer than nInd. Everything below indexes
  # female, so the count has to come from female itself.
  nInd = length(female)
  nSelf = rbinom(n=nInd,prob=probSelf,size=nSeeds)
  if(pollenControl){
    male = female
  }else{
    male = seq_len(pop@nInd)
  }
  crossPlan = lapply(seq_len(nInd),function(x){
    male = male[!male==female[x]]
    if(length(male)==1){
      #Account for "convenience" feature of sample when length = 1
      cbind(rep(female[x],nSeeds),
            c(rep(female[x],nSelf[x]),
              rep(male,nSeeds-nSelf[x])))
    }else{
      cbind(rep(female[x],nSeeds),
            c(rep(female[x],nSelf[x]),
              sample(male,nSeeds-nSelf[x],replace=TRUE)))
    }
  })
  crossPlan = mergeMultIntMat(crossPlan,rep(nSeeds,nInd),2L)
  return(makeCross(pop=pop,crossPlan=crossPlan,simParam=simParam,
                   nThreads=nThreads))
}

#' @title Select populations
#'
#' @description
#' Selects a subset of \code{nPop} populations from a 
#' \code{\link{MultiPop-class}} object based on a summary value 
#' calculated from the selection criterion \code{use}.
#' 
#' For nested \code{MultiPop} objects (containing one or more 
#' \code{MultiPop} objects), selection is performed at the specified 
#' \code{level}, where \code{level=1} selects among direct children 
#' that are \code{\link{Pop-class}} objects, \code{level=2} descends 
#' one level into nested \code{MultiPop}s to select among their direct 
#' children, and so on.
#' 
#' Internally, \code{selectPop} executes the summary function \code{FUN} 
#' on the vector or matrix of \code{use} values for each population. The 
#' resulting vector of summary values is used to rank and select populations.
#'
#' @param x a \code{\link{MultiPop-class}} object
#' @param nPop the number of populations to select
#' @param level nesting depth at which selection is performed. 
#'   \code{level=1} selects among direct children, \code{level=2} 
#'   descends one level into nested \code{MultiPop}s, etc.
#' @param trait the trait for selection. Either a number indicating 
#'   a single trait or a character for a trait name, or a function 
#'   returning a vector of length \code{nInd(pop)} for each population 
#'   at \code{level}. The function must work on a vector or matrix of 
#'   \code{use} values as \code{trait(pop@use, ...)} - depending on what 
#'   \code{use} is. See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values \code{"gv"}, estimated breeding values \code{"ebv"}, 
#'   phenotypes \code{"pheno"}, or random values \code{"rand"}) or a function 
#'   returning a vector of length \code{nInd(pop)} for each population 
#'   at \code{level}. The function must work on each \code{pop} at 
#'   \code{level} in \code{x} as \code{use(pop, trait, ...)} or as 
#'   \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param selectTop selects highest values if \code{TRUE}. 
#'   Selects lowest values if \code{FALSE}.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#'   \code{trait} or \code{use}
#' @param FUN a summary function to be applied to the vector of 
#'   selection criterion values per population. Default is \code{\link{mean}}.
#' @param FUN.ARGS a list of additional arguments passed to \code{FUN}
#'
#' @details
#' The \code{level} argument controls at which nesting depth selection 
#' operates. At \code{level=1}, selection occurs among direct children 
#' of \code{x} (which must be \code{\link{Pop-class}} objects). At 
#' \code{level=2}, the function descends one level into nested 
#' \code{MultiPop}s and selects among their direct children, and so on. 
#' If \code{level} exceeds the nesting depth of \code{x}, the 
#' function stops with an error.
#' 
#' If the number of eligible populations is less than \code{nPop}, the
#' function returns all eligible populations and issues a warning.
#'
#' @return Returns a \code{\link{MultiPop-class}} object
#'   containing the selected populations. When \code{x} is a 
#'   \code{\link{Pop-class}}, the original object is returned unchanged.
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=20, nChr=1, segSites=10)
#' 
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' 
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' # Create multi-population with one level of nesting
#' mp1 = newMultiPop(pop[1:4], pop[5:8], pop[9:12], pop[13:16], pop[17:20])
#' 
#' # Select best 2 populations
#' selectPop(mp1, nPop=2, simParam=SP)
#' 
#' # Select 2 most deviating populations from an optima
#' squaredDeviation = function(x, optima=0) (x - optima)^2
#' selectPop(mp1, nPop=2, trait=squaredDeviation, simParam=SP)
#' 
#' # Create multi-population with two levels of nesting
#' pop@misc = list(smth=rnorm(20), smth2=rnorm(20))
#' mp2 = newMultiPop(
#'   do.call(newMultiPop, split(pop[1:10], rep(1:5, 2))),
#'   do.call(newMultiPop, split(pop[11:20], rep(1:5, 2)))
#' )
#' 
#' # Select best 2 populations at level 2
#' selectPop(mp2, nPop=2, level=2, simParam=SP)
#' 
#' # Select 2 populations at level 2 based on miscellaneous info with use function
#' useFunc = function(pop, trait=NULL) pop@misc$smth + pop@misc$smth2
#' selectPop(mp2, nPop=2, level=2, simParam=SP, use=useFunc)
#' 
#' # ... equivalent result with the use & trait function
#' useFunc2 = function(pop, trait=NULL) cbind(pop@misc$smth, pop@misc$smth2)
#' trtFunc = function(x) rowSums(x)
#' selectPop(mp2, nPop=2, level=2, use=useFunc2, trait=trtFunc, simParam=SP)
#'
#' @export
selectPop = function(
  x,
  nPop,
  level = 1,
  trait = 1,
  use = "pheno",
  selectTop = TRUE,
  simParam = NULL,
  ...,
  FUN = mean,
  FUN.ARGS = list()
) {

  stopifnot("nPop must be >= 0" = nPop >= 0)
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }

  if (isPop(x)) {
    return(x)
  }
  stopifnot("x must be a Pop or a MultiPop" = isMultiPop(x))

  if (
    length(level) != 1L ||
      !is.numeric(level) ||
      is.na(level) ||
      !is.finite(level) ||
      level < 1 ||
      level != trunc(level) ||
      abs(level) > .Machine$integer.max
  ) {
    stop("`level` must be a single non-NA integer value.")
  }
  if (is.character(use) && length(use) == 1L && tolower(use) == "bv") {
    stop("use='bv' is not currently supported for populations.")
  }
  if (!is.logical(selectTop) || length(selectTop) != 1L || is.na(selectTop)) {
    stop("`selectTop` must be TRUE or FALSE.")
  }
  if (!is.list(FUN.ARGS)) {
    stop("`FUN.ARGS` must be a list.")
  }

  is_multi = vapply(unname(x@pops), isMultiPop, logical(1L))
  is_pop = vapply(unname(x@pops), isPop, logical(1L))
  multi = which(is_multi)

  while (level > 1L) {
    level = level - 1L
    for (i in multi) {
      x@pops[[i]] = selectPop(
        x = x[[i]],
        nPop = nPop,
        level = level,
        trait = trait,
        use = use,
        selectTop = selectTop,
        simParam = simParam,
        FUN = FUN,
        FUN.ARGS = FUN.ARGS,
        ...
      )
    }
    return(x)
  }

  if (!identical(multi, integer(0))) {
    stop(paste(
      "This level contains",
      length(multi),
      "MultiPop-class objects.",
      "\nSelection can only be performed when all populations at this level are",
      "Pop-class objects.\nYou may want to increase the value of 'level'"
    ))
  }

  eligible = which(is_pop)

  if (length(eligible) < nPop) {
    nPop = length(eligible)
    warning(
      "Suitable candidate populations smaller than nPop, returning ",
      nPop,
      " populations"
    )
  }

  response = calcPopValue(
    x,
    FUN = function(pop) {
      getResponse(pop = pop, trait = trait, use = use, simParam = simParam, ...)
    },
    simplify = FALSE,
    level = 1L
  )

  popValues = vapply(response, function(res) {
    do.call(FUN, c(list(res), FUN.ARGS))
  }, numeric(1L))

  take = order(popValues, decreasing = selectTop)
  take = take[take %in% eligible]

  return(x[take[seq_len(nPop)]])
}

#' @title Select units from a MultiPop hierarchy
#'
#' @description
#' Selects \code{\link{Pop-class}} or \code{\link{MultiPop-class}}
#' objects at a specified level of a nested \code{MultiPop} hierarchy.
#' Selection is based on response summaries calculated \emph{within} \code{Pop}
#' objects and recursively aggregated \emph{across} \code{MultiPop} objects.
#'
#' @param x a \code{\link{MultiPop-class}} object
#' @param nMulti the number of population units to select
#' @param level an integer indicating the nesting level at which selection is
#'   performed. \code{level=1} selects among direct children, \code{level=2} 
#'   descends one level into nested \code{MultiPop}s, etc.
#' @param trait the trait for selection. Either a number indicating 
#'   a single trait or a character for a trait name, or a function 
#'   returning a vector of length \code{nInd(pop)} for each population 
#'   at \code{level}. The function must work on a vector or matrix of 
#'   \code{use} values as \code{trait(pop@use, ...)} - depending on what 
#'   \code{use} is. See \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values \code{"gv"}, estimated breeding values \code{"ebv"}, 
#'   phenotypes \code{"pheno"}, or random values \code{"rand"}) or a function 
#'   returning a vector of length \code{nInd(pop)} for each population 
#'   at \code{level}. The function must work on each \code{pop} at 
#'   \code{level} in \code{x} as \code{use(pop, trait, ...)} or as 
#'   \code{trait(pop@use, ...)} depending on what \code{trait} is.
#' @param selectTop selects highest values if \code{TRUE}. 
#'   Selects lowest values if \code{FALSE}.
#' @param withinPopFUN a function used to summarize individual responses
#'   within each \code{Pop} object. The function is applied to the vector of
#'   selection criterion values per population and returns one summary value.
#'   Default is \code{\link{mean}}.
#' @param withinPopArgs a list of additional arguments passed to
#'   \code{withinPopFUN}
#' @param acrossPopFUN a function used to aggregate population summaries
#'   across the immediate children of each \code{MultiPop} object. The
#'   function is applied to a vector of population summary values and returns
#'   one aggregated value. Default is \code{\link{mean}}.
#' @param acrossPopArgs a list of additional arguments passed to
#'   \code{acrossPopFUN}
#' @param simParam an object of \code{\link{SimParam}}. If \code{NULL},
#'   the function uses the object named \code{SP} from the global environment.
#' @param ... additional arguments passed to the functions supplied to
#'   \code{trait} or \code{use}
#'
#' @details
#' Individual responses are first obtained for every terminal \code{Pop}
#' object using the specified \code{trait} and \code{use}. These responses are
#' summarized within each \code{Pop} using \code{withinPopFUN}. The resulting
#' summaries are then recursively aggregated with \code{acrossPopFUN} until
#' the requested \code{level} is reached.
#'
#' Population units at the requested \code{level} are ranked using their
#' aggregated values. The highest \code{nMulti} units are retained when
#' \code{selectTop = TRUE}, and the lowest \code{nMulti} units are retained
#' when \code{selectTop = FALSE}.
#'
#' Aggregation with \code{acrossPopFUN} is unweighted with respect to the
#' number of individuals in each child population, unless weighting is
#' explicitly implemented by the user-supplied function.
#'
#' If a branch ends in a \code{Pop} object before the requested \code{level} is
#' reached, that terminal \code{Pop} is treated as a candidate at its deepest
#' available level. If candidates have equal selection values at the selection
#' boundary, ties are resolved according to their existing order in the
#' \code{MultiPop} hierarchy.
#'
#' @return
#' A \code{\link{MultiPop-class}} object containing the selected
#' \code{nMulti} population units at the requested \code{level}. Population
#' units not selected are removed, together with any intermediate
#' \code{MultiPop} objects left empty by their removal.
#'
#' If \code{nMulti} is greater than the number of available candidates, a
#' warning is issued and \code{x} is returned unchanged.
#'
#' @examples
#' founderPop = quickHaplo(nInd = 16, nChr = 1, segSites = 50)
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(nQtlPerChr = 20)
#' SP$setVarE(h2 = 0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' pop = newPop(founderPop, simParam = SP)
#'
#' # Create four groups, each divided into two replicates
#' multiPop = splitPop(
#'   pop,
#'   by = list(
#'     rep(LETTERS[1:4], each = 4),
#'     function(x) {
#'       rep(c("R1", "R2"), length.out = nInd(x))
#'     }
#'   )
#' )
#'
#' # Select the two groups with the highest mean phenotypic values
#' selected = selectMultiPop(
#'   multiPop,
#'   nMulti = 2,
#'   level = 1,
#'   trait = 1,
#'   use = "pheno",
#'   withinPopFUN = mean,
#'   acrossPopFUN = mean,
#'   simParam = SP
#' )
#' selected
#'
#' # Select the three replicates with the lowest sum of phenotypic values
#' selected = selectMultiPop(
#'   multiPop,
#'   nMulti = 3,
#'   level = 2,
#'   trait = 1,
#'   use = "pheno",
#'   selectTop = FALSE,
#'   withinPopFUN = sum,
#'   acrossPopFUN = mean,
#'   simParam = SP
#' )
#' selected
#'
#' @seealso \code{\link{calcPopSummary}}, \code{\link{selectPop}},
#'   \code{\link{selectInd}}, \code{\link{meanPPop}}
#'
#' @export
selectMultiPop = function(
  x,
  nMulti,
  level = 1,
  trait = 1,
  use = "pheno",
  selectTop = TRUE,
  withinPopFUN = mean,
  withinPopArgs = list(),
  acrossPopFUN = mean,
  acrossPopArgs = list(),
  simParam = NULL,
  ...
) {
  if (!isMultiPop(x)) {
    stop("`x` must be a MultiPop object.")
  }
  if (
    !is.numeric(nMulti) ||
      length(nMulti) != 1L ||
      is.na(nMulti) ||
      !is.finite(nMulti) ||
      nMulti != trunc(nMulti) ||
      nMulti < 1 ||
      nMulti > .Machine$integer.max
  ) {
    stop("`nMulti` must be a single positive integer.")
  }
  if (
    !is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      !is.finite(level) ||
      level != trunc(level) ||
      level < 1 ||
      level > .Machine$integer.max
  ) {
    stop("`level` must be a single positive integer.")
  }
  if (is.character(use) && length(use) == 1L && tolower(use) == "bv") {
    stop("use='bv' is not currently supported for populations.")
  }
  if (!is.logical(selectTop) || length(selectTop) != 1L || is.na(selectTop)) {
    stop("`selectTop` must be TRUE or FALSE.")
  }
  if (!is.function(withinPopFUN)) {
    stop("`withinPopFUN` must be a function.")
  }
  if (!is.list(withinPopArgs)) {
    stop("`withinPopArgs` must be a list.")
  }

  nMulti = as.integer(nMulti)
  level = as.integer(level)
  dots = list(...)

  # Calculate a response-based summary for each Pop object
  summarizeResponse = function(pop) {
    responseArgs = c(
      list(pop = pop, trait = trait, use = use, simParam = simParam),
      dots
    )
    response = do.call(getResponse, responseArgs)
    return(do.call(withinPopFUN, c(list(response), withinPopArgs)))
  }

  candidateScores = calcPopSummary(
    # Remove names so that source paths contain positional indices.
    # These paths are later used to prune the original MultiPop object.
    x = unnameMultiPop(x),
    withinPopFUN = summarizeResponse,
    acrossPopFUN = acrossPopFUN,
    acrossPopArgs = acrossPopArgs,
    level = level
  )

  src = attr(candidateScores, "source")

  if (!is.matrix(candidateScores) || ncol(candidateScores) != 1L) {
    stop(
      "`withinPopFUN` and `acrossPopFUN` must produce exactly one ",
      "selection value per candidate."
    )
  }
  if (anyNA(candidateScores) || any(!is.finite(candidateScores))) {
    stop("Candidate selection values must be finite and non-missing.")
  }

  nCandidates = nrow(candidateScores)
  if (nMulti > nCandidates) {
    warning(
      "`nMulti` exceeds the number of available candidates; ",
      "returning all candidates."
    )
    return(x)
  }
  if (nMulti == nCandidates) {
    return(x)
  }

  ranking = order(candidateScores, decreasing = selectTop)

  drop = ranking[-seq_len(nMulti)]
  dropPaths = src[drop, , drop = FALSE]

  x = .stripPops(x, dropPaths)

  return(x)
}

#' Helper function to collect leaf paths in a \code{MultiPop}
#'
#' @param x \code{\link{Pop-class}} or \code{\link{MultiPop-class}} object.
#' @param path Character vector of labels representing the current path.
#'
#' @keywords internal
.collectLeafPaths = function(x, path = list()) {
  if (isPop(x)) {
    return(list(path))
  }
  stopifnot(isMultiPop(x))

  out = list()
  nm = names(x@pops)
  if (!is.null(nm)) {
    valid_names = !is.na(nm) & nzchar(nm)
    duplicated_names = valid_names & duplicated(nm) |
      valid_names & duplicated(nm, fromLast = TRUE)
  } else {
    valid_names = duplicated_names = logical(length(x@pops))
  }

  for (i in seq_along(x@pops)) {
    child = x@pops[[i]]
    label = if (valid_names[i] && !duplicated_names[i]) nm[i] else i
    childPath = c(path, list(label))
    out = c(out, .collectLeafPaths(child, path = childPath))
  }

  return(out)
}

#' Helper function to format the source attribute for simplified \code{MultiPop}
#' output
#'
#' @param paths List of character vectors representing leaf paths.
#' @param nRows Vector of row counts for each path.
#' @param level_offset Integer scalar used to label source levels.
#'
#' @keywords internal
.formatPopSource = function(paths, nRows, level_offset = 0L) {
  if (length(paths) == 0L) {
    return(NULL)
  }

  maxDepth = max(lengths(paths))
  pathList = lapply(paths, function(p) {
    pad = rep(list(NA), maxDepth - length(p))
    c(p, pad)
  })

  pathDf = data.frame(sapply(pathList, `[[`, 1))
  for (i in seq_len(maxDepth)[-1]) {
    pathDf = cbind(pathDf, sapply(pathList, `[[`, i))
  }
  colnames(pathDf) = paste0("level", level_offset + seq_len(maxDepth))

  for (j in seq_len(ncol(pathDf))) {
    col = pathDf[[j]]
    first = col[!is.na(col)][1]
    if (!is.null(first) && is.numeric(first)) {
      pathDf[[j]] = as.integer(col)
    } else {
      pathDf[[j]] = as.character(col)
    }
  }

  idx = rep(seq_along(nRows), nRows)
  pathDf = pathDf[idx, , drop = FALSE]
  rownames(pathDf) = NULL

  return(pathDf)
}

#' @title Recursive worker for \code{calcPopSummary}
#'
#' @description
#' Recurses over a nested \code{\link{MultiPop-class}} object, applying a
#' summary function to each \code{\link{Pop-class}} object and aggregating 
#' the resulting summaries across \code{MultiPop} objects.
#'
#' @param node a \code{\link{Pop-class}} or
#'   \code{\link{MultiPop-class}} object
#' @param levelsRemaining an integer indicating the number of nesting levels
#'   remaining before aggregation stops
#' @param withinPopFUN a function applied to each \code{Pop} object. The
#'   function must return a vector or one-row matrix of summary values.
#' @param withinPopArgs a list of additional arguments passed to
#'   \code{withinPopFUN}
#' @param acrossPopFUN a function used to aggregate summaries across the
#'   immediate children of each \code{MultiPop} object. The function must
#'   accept a matrix of child summaries and return a vector or one-row matrix
#'   with the same columns.
#' @param acrossPopArgs a list of additional arguments passed to
#'   \code{acrossPopFUN}
#'
#' @return A vector or matrix of population summaries, aggregated until the
#'   requested nesting level is reached.
#'
#' @keywords internal
.calcPopSummary = function(
  node,
  levelsRemaining,
  withinPopFUN,
  withinPopArgs,
  acrossPopFUN,
  acrossPopArgs
) {
  # Normalize the output of withinPopFUN and acrossPopFUN.
  # Each Pop or aggregated MultiPop must contribute exactly one summary
  # row, represented internally as a numeric vector.
  normalizeSummary = function(value, functionName) {
    if (is.data.frame(value)) {
      value = as.matrix(value)
    }
    if (is.matrix(value) && nrow(value) != 1L) {
      stop("`", functionName, "` must return a vector or a one-row matrix.")
    }
    # Missing values are allowed, so that a summary such as colMeans() can
    # report NA for one trait without losing the others. Zero-length
    # summaries are allowed too, because a SimParam with no traits is legal
    # and colMeans() of its zero-column gv or pheno is numeric(0). Callers
    # that need one complete value, such as selectMultiPop(), check for it
    # themselves.
    if (!is.numeric(value)) {
      stop(
        "`",
        functionName,
        "` must return a numeric vector or one-row matrix."
      )
    }
    # A one-row matrix becomes a named vector, so that every summary is
    # compared and returned in the same form regardless of which shape the
    # user-supplied function chose.
    if (is.matrix(value)) {
      value = setNames(as.vector(value), colnames(value))
    }
    return(value)
  }

  if (isPop(node) || isHybridPop(node)) {
    if (node@nInd == 0L) {
      stop("`x` contains an empty Pop object.")
    }
    value = do.call(withinPopFUN, c(list(node), withinPopArgs))
    return(normalizeSummary(value, "withinPopFUN"))
  }
  if (length(node) == 0L) {
    stop("`x` contains an empty MultiPop object.")
  }

  childSummaries = lapply(
    node@pops,
    .calcPopSummary,
    levelsRemaining = levelsRemaining - 1L,
    withinPopFUN = withinPopFUN,
    withinPopArgs = withinPopArgs,
    acrossPopFUN = acrossPopFUN,
    acrossPopArgs = acrossPopArgs
  )

  same_columns = all(vapply(
    childSummaries,
    function(x) {
      reference = childSummaries[[1]]
      x_names = if (is.atomic(x) && is.null(dim(x))) names(x) else colnames(x)
      reference_names = if (is.atomic(reference) && is.null(dim(reference))) {
        names(reference)
      } else {
        colnames(reference)
      }
      identical(x_names, reference_names)
    },
    logical(1)
  ))
  if (!same_columns) {
    stop("Population summaries do not have consistent names.")
  }
  
  same_lengths = vapply(
    childSummaries,
    function(x) if (is.null(dim(x))) length(x) else ncol(x),
    integer(1)
  )
  if (length(unique(same_lengths)) != 1L) {
    stop("Population summaries do not have consistent numbers of columns.")
  }

  childSummaryMatrix = do.call(rbind, childSummaries)

  # Continue aggregating when the requested level has been reached.
  if (levelsRemaining < 1L) {
    aggregated = do.call(
      acrossPopFUN,
      c(list(childSummaryMatrix), acrossPopArgs)
    )
    
    aggregated = normalizeSummary(aggregated, "acrossPopFUN")

    # Aggregation must preserve the summary columns. This is important
    # for unbalanced hierarchies, where aggregated MultiPop summaries
    # may later be combined with terminal Pop summaries.
    if (length(aggregated) != ncol(childSummaryMatrix)) {
      stop(
        "`acrossPopFUN` must return one value for each column ",
        "produced by `withinPopFUN`."
      )
    }

    # Names returned by acrossPopFUN say which column each value belongs
    # to, so they are checked rather than overwritten. A permutation of the
    # columns is put back in order by name. Unnamed results, such as those
    # from mean(), are assumed to follow the column order and take the
    # column names.
    columnNames = colnames(childSummaryMatrix)
    aggregatedNames = names(aggregated)
    if (
      !is.null(columnNames) &&
        !is.null(aggregatedNames) &&
        !identical(aggregatedNames, columnNames)
    ) {
      if (
        anyDuplicated(columnNames) ||
          anyDuplicated(aggregatedNames) ||
          !setequal(aggregatedNames, columnNames)
      ) {
        stop(
          "`acrossPopFUN` returned names that do not match the columns ",
          "produced by `withinPopFUN`."
        )
      }
      aggregated = aggregated[match(columnNames, aggregatedNames)]
    }
    names(aggregated) = columnNames
    return(aggregated)
  }
  return(childSummaryMatrix)
}

#' Remove nested populations from a MultiPop
#'
#' Recursively removes one or more \code{\link{Pop-class}} or
#' \code{\link{MultiPop-class}} objects identified by positional paths through
#' a nested \code{MultiPop} object. Empty \code{MultiPop} objects created by
#' the removal are also removed.
#'
#' @param x a \code{\link{MultiPop-class}} object
#' @param paths an integer matrix containing one positional path per row.
#'   Columns correspond to successive nesting levels. Trailing \code{NA}
#'   values indicate that a path ends at an earlier level.
#'
#' @details
#' Population paths are interpreted relative to the original structure of
#' \code{x}. All requested objects are removed in one recursive traversal,
#' avoiding changes to positional indices while pruning is in progress.
#' The order and names of retained objects are preserved.
#'
#' @return A \code{\link{MultiPop-class}} object with the specified objects
#'   removed. Intermediate \code{MultiPop} objects left empty by the removal
#'   are also removed.
#'
#' @keywords internal
.stripPops = function(x, paths) {
  paths = as.matrix(paths)

  if (nrow(paths) == 0L) {
    return(x)
  }

  prune = function(node, paths) {
    children = node@pops
    keep = rep(TRUE, length(children))

    # Group the rows by child once, rather than scanning every row for each
    # child, which would be quadratic in the number of removed siblings.
    rowsByChild = split(seq_len(nrow(paths)), paths[, 1L])

    for (child in names(rowsByChild)) {
      idx = as.integer(child)
      childPaths = paths[rowsByChild[[child]], , drop = FALSE]

      # Remove the complete child if the path ends at this level.
      pathEnds = ncol(childPaths) == 1L || all(is.na(childPaths[, 2L]))

      if (pathEnds) {
        keep[idx] = FALSE
        next
      }

      remainingPaths = childPaths[, -1L, drop = FALSE]
      children[[idx]] = prune(children[[idx]], remainingPaths)

      # Remove parents left empty by recursive pruning.
      if (length(children[[idx]]) == 0L) {
        keep[idx] = FALSE
      }
    }

    # Subsetting a list preserves the order and names of retained children.
    children = children[keep]

    # Defer validation until the complete hierarchy has been pruned.
    slot(node, "pops", check = FALSE) = children

    return(node)
  }

  x = prune(x, paths)
  validObject(x)

  return(x)
}
