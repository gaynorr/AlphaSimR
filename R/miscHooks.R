#' @title Functions (hooks) for subsetting and merging individual and
#' population miscellaneous data
#'
#' @description Default functions (hooks) for subsetting and merging of
#' the individual-level \code{misc} and
#' the population-level \code{miscPop} slots of a \code{\link{Pop-class}}.
#' These functions can also be called inside custom functions.
#'
#' @param pop A \code{\link{Pop-class}} object.
#' @param popList A list of \code{\link{Pop-class}} object,
#'   with \code{NULL} entries removed.
#' @param i An index for subsetting. Character input is matched to
#' \code{pop@id}, while numeric and logical indices retain their usual R meaning,
#'   including repeated and negative indices.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#'
#' @return A list for the corresponding slot: \code{subsetMisc} and
#' \code{mergeMisc} return the new \code{misc}; \code{subsetMiscPop} and
#' \code{mergeMiscPop} return the new \code{miscPop}. Each element of
#' \code{misc} must have one item per resulting individual (one row for a matrix).
#'
#' @details
#' Set \code{SP$subsetMisc = function(pop, i, simParam) { code }} and
#' \code{SP$mergeMisc = function(popList, simParam) { code }} for individual data.
#' Set \code{SP$subsetMiscPop} and \code{SP$mergeMiscPop} with the same
#' respective signatures for population data.
#' These four slots and their behavior are independent.
#' They must be functions or \code{NULL}.
#' Each \code{simParam/SP} function is initialized to its corresponding default
#' function. Assigning \code{NULL} also selects the default.
#' The four functions are read at each operation.
#' If you are not naming \code{simParam}  as \code{SP} globally,
#' pass it explicitly with
#' \code{pop[i, simParam=SPX]}, \code{c(pop1, pop2, simParam=SPX)} or
#' \code{mergePops(popList, simParam=SPX)}.
#'
#' Both subset functions receive the original population and index.
#' Do not call \code{pop[i]} inside its own subset function, which would recurse.
#'
#' Both merge hooks are called when \code{\link{mergePops}} or \code{c}
#' combines \code{\link{Pop-class}} objects into a single population.
#' For nested \code{\link{MultiPop-class}} inputs, \code{mergePops}
#' merges the innermost populations first, then merges those results.
#' Each hook must therefore accept populations produced by an earlier merge.
#'
#' The defaults subset vectors and generic objects with \code{[},
#' subset matrices by rows without dropping dimensions, and
#' merge matching miscellaneous elements with \code{c} or \code{rbind}.
#' Mismatched element counts or names during merging warn and
#' clear \code{misc}.
#' 
#' The \code{miscPop} default functions return an empty list,
#' clearing \code{miscPop}, because the information may no longer
#' describe the resulting population correctly. You can change this
#' behavior by providing your own functions to populate \code{miscPop}.
#' 
#' Custom results are checked against the population's validity requirements.
#'
#' @examples
#' founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' # \dontshow{SP$nThreads = 1L}
#'
#' # Recompute a population total on subset and merge
#' SP$subsetMiscPop = function(pop, i, simParam=NULL){
#'   list(total=sum(pop@misc$measurement[i]))
#' }
#' SP$mergeMiscPop = function(popList, simParam=NULL){
#'   list(total=sum(vapply(popList, function(pop){
#'     sum(pop@misc$measurement)
#'   }, numeric(1))))
#' }
#' pop = newPop(founderPop)
#' pop@misc$measurement = c(2, 4, 6, 8)
#' pop@miscPop$total = sum(pop@misc$measurement)
#' pop@misc
#' pop@miscPop # 20
#' selected = pop[c(4, 2, 2)]
#' selected@misc
#' selected@miscPop # 16
#' combined = c(pop[1], pop[2])
#' combined@misc
#' combined@miscPop$total # 6
#' pop[0]@miscPop$total # 0
#'
#' @export
subsetMisc = function(pop, i, simParam=NULL){
  if(is.character(i)){
    i = match(i, pop@id)
  }
  misc = lapply(pop@misc, function(x){
    if(isPop(x)){
      x[i, simParam=simParam]
    }else if(is.matrix(x)){
      x[i, , drop=FALSE]
    }else{
      x[i]
    }
  })
  return(misc)
}

#' @describeIn subsetMisc Default function for merging individual-level miscellaneous data
#' @export
mergeMisc = function(popList, simParam=NULL){
  # vapply rather than sapply, which returns a list for an empty popList,
  # and the count is held separately because tmp is reused for the names
  # below. Unnamed misc entries have NULL names, so taking the count from
  # the names dropped every unnamed misc element.
  nMisc = vapply(popList, function(x) length(x@misc), integer(1))
  if(!all(nMisc==nMisc[1])){
    warning("number of misc elements differs - setting misc to an empty list!")
    misc = list()
  }else{
    if(nMisc[1]>0){
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
        misc = vector("list", length=nMisc[1])
        for(i in seq_len(nMisc[1])){
          miscTmp = lapply(popList, function(x) x@misc[[i]])
          if(isPop(miscTmp[[1]]) || isMultiPop(miscTmp[[1]])){
            misc[[i]] = do.call("c", c(miscTmp, list(simParam=simParam)))
          }else if(is.matrix(miscTmp[[1]])){
            misc[[i]] = do.call("rbind", miscTmp)
          }else{
            misc[[i]] = do.call("c", miscTmp)
          }
        }
        names(misc) = tmp[[1]]
      }else{
        warning(
          "misc element names do not match - setting misc to an empty list!"
        )
        misc = list()
      }
    }else{
      misc = list()
    }
  }
  return(misc)
}

#' @describeIn subsetMisc Default function for subsetting population-level miscellaneous data
#' @export
subsetMiscPop = function(pop, i, simParam=NULL){
  list()
}

#' @describeIn subsetMisc Default function for merging population-level miscellaneous data
#' @export
mergeMiscPop = function(popList, simParam=NULL){
  list()
}

#' Validate miscellaneous-data class
#'
#' Check that a miscellaneous-data hook returned a list.
#'
#' @param result The value returned by a miscellaneous-data hook.
#' @param name The hook name to include in an error message.
#' @return The supplied list, unchanged. An error is raised for a non-list result.
#' @keywords internal
isMiscAList = function(result, name){
  if(!is.list(result)){
    stop(paste0(name, " must return a list"))
  }
  return(result)
}
