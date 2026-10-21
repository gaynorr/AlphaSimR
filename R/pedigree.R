# Pedigree relationships ----

# Collects the recorded pedigree of a set of individuals and all of their
# ancestors, renumbered from 1 in iid order with 0 for an unknown parent.
# Parents always have smaller iids than their progeny, so the result has
# every parent before its progeny, which the C++ functions in pedigree.cpp
# require. Individuals that are not ancestors of the requested ones are
# left out, because they cannot change the relationships among them.
# Returns the pedigree and, in row, where each requested individual is in
# it.
.prunePed = function(iid, simParam){
  if(!simParam$isTrackPed){
    stop("Pedigree relationships can only be calculated with the trackPed option, see setTrackPed in ?SimParam",
         call.=FALSE)
  }
  ped = simParam$pedigree
  # The pedigree has one row per iid only when it has been tracked since
  # the first individual was made. Tracking turned on later leaves the
  # rows out of step with the iids.
  if(nrow(ped)!=simParam$lastId){
    stop("The pedigree does not cover every individual, which happens when setTrackPed is turned on after the simulation has started",
         call.=FALSE)
  }
  iid = as.integer(iid)
  if(any(is.na(iid) | iid<1L | iid>nrow(ped))){
    stop("Some individuals are not in the recorded pedigree",
         call.=FALSE)
  }

  need = logical(nrow(ped))
  cur = unique(iid)
  while(length(cur)>0L){
    cur = cur[cur>0L]
    cur = cur[!need[cur]]
    need[cur] = TRUE
    cur = unique(c(ped[cur,1:2]))
  }
  keep = which(need)

  newRow = integer(nrow(ped))
  newRow[keep] = seq_along(keep)
  mother = ped[keep,1]
  father = ped[keep,2]
  take = mother>0L
  mother[take] = newRow[mother[take]]
  take = father>0L
  father[take] = newRow[father[take]]

  return(list(iid=keep,
              mother=as.integer(mother),
              father=as.integer(father),
              isDH=as.integer(ped[keep,3]),
              row=newRow[iid]))
}

# Checks shared by the exported functions below
.checkPedPop = function(pop){
  if(!is(pop, "Pop")){
    stop("pop must be an object of class Pop", call.=FALSE)
  }
  if(pop@ploidy!=2L){
    stop("Pedigree relationships are only available for diploids",
         call.=FALSE)
  }
  invisible(NULL)
}

#' @title Pedigree relationship matrix
#'
#' @description
#' Calculates the numerator relationship matrix of a population from the
#' pedigree recorded in the \code{\link{SimParam}}. The pedigree must have
#' been tracked from the start of the simulation, see
#' \code{\link[=SimParam-setTrackPed]{setTrackPed}}.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#'
#' @details
#' The relationship matrix is calculated by the tabular method from the
#' pedigree of the population's individuals and all of their ancestors.
#' Individuals in the base population, those created by
#' \code{\link{newPop}}, are taken to be unrelated and not inbred, unless
#' the founders were inbred, as with \code{runMacs(inbred=TRUE)}, in which
#' case they are taken to be unrelated and fully inbred.
#'
#' Doubled haploids, made with \code{\link{makeDH}}, are fully inbred and
#' are related to every other individual exactly as their parent is,
#' because a doubled haploid's genome is two copies of a single gamete from
#' that parent. A selfed individual has the same relationships to others as
#' its parent, but is only partly inbred.
#'
#' The relationships assume a diploid pedigree. An error is given for a
#' population that is not diploid, but a diploid population descended from
#' polyploids, as through \code{\link{reduceGenome}}, cannot be detected,
#' and its relationships will be wrong.
#'
#' The matrix is dense and holds every ancestor of the population while it
#' is built, so its memory use grows with the square of the size of the
#' pedigree.
#'
#' @return Returns a symmetric matrix with a row and column for each
#' individual, named by their IDs. The diagonal is one plus the inbreeding
#' coefficient.
#'
#' @references
#' \insertRef{henderson1976Simple}{AlphaSimR}
#'
#' @seealso \code{\link{pedigreeInbreeding}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$setTrackPed(TRUE)
#'
#' #Create a pedigree with crosses and doubled haploids
#' pop = newPop(founderPop, simParam=SP)
#' pop2 = randCross(pop, nCrosses=2, nProgeny=2, simParam=SP)
#' pop3 = makeDH(pop2, nDH=2, simParam=SP)
#'
#' #Relationships among the doubled haploids
#' pedigreeA(pop3, simParam=SP)
#'
#' @export
pedigreeA = function(pop, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  .checkPedPop(pop)
  ped = .prunePed(pop@iid, simParam)
  A = calcPedA(ped$mother, ped$father, ped$isDH)
  A = A[ped$row, ped$row, drop=FALSE]
  dimnames(A) = list(pop@id, pop@id)
  return(A)
}

#' @title Pedigree inbreeding coefficients
#'
#' @description
#' Calculates the inbreeding coefficient of each individual in a population
#' from the pedigree recorded in the \code{\link{SimParam}}. The pedigree
#' must have been tracked from the start of the simulation, see
#' \code{\link[=SimParam-setTrackPed]{setTrackPed}}.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#'
#' @details
#' The inbreeding coefficients are calculated by the method of
#' \insertCite{meuwissen1992Computing;textual}{AlphaSimR}, which works from
#' the pedigree without building the relationship matrix, so it can be used
#' on large pedigrees. Doubled haploids and inbred founders have an
#' inbreeding coefficient of 1, and other founders of 0. The assumptions are
#' those described in \code{\link{pedigreeA}}, whose diagonal is one plus
#' these values.
#'
#' @return Returns a numeric vector named by the individuals' IDs.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{pedigreeA}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$setTrackPed(TRUE)
#'
#' #Two generations of selfing
#' pop = newPop(founderPop, simParam=SP)
#' pop2 = self(pop, simParam=SP)
#' pop3 = self(pop2, simParam=SP)
#'
#' pedigreeInbreeding(pop3, simParam=SP)
#'
#' @export
pedigreeInbreeding = function(pop, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  .checkPedPop(pop)
  ped = .prunePed(pop@iid, simParam)
  F = calcPedInbreeding(ped$mother, ped$father, ped$isDH)
  F = F[ped$row]
  names(F) = pop@id
  return(F)
}
