
# Pedigree ------------------------------------------------------------------

#' @title Pedigree
#' 
#' @description 
#' A pedigree
#' 
#' @slot nInd number of individuals
#' @slot ids ID of each individual
#' @slot mother position of mother in pedigree (0 for not in pedigree)
#' @slot father postiion of father in pedigree (0 for nor in pedigree) 
#' 
#' @export
setClass("Pedigree",
         slots=c(nInd="numeric",
                 ids="character",
                 mother="integer",
                 father="integer"))

setValidity("Pedigree",function(object){
  errors = character()
  if(object@nInd!=length(object@ids)){
    errors = c(errors,"nInd!=length(ids)")
  }
  if(object@nInd!=length(object@mother)){
    errors = c(errors,"nChr!=length(mother)")
  }
  if(object@nInd!=length(object@father)){
    errors = c(errors,"nChr!=length(father)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#' @title Sorts a pedigree
#' 
#' @description
#' Sorts a pedigree so parents appear before children
#'
#' @param x an object of code{\link{Pedigree-class}}
#' @param maxGen how many times should it loop when trying to sort
#'
#' @return Returns an object of \code{\link{Pedigree-class}}
#' 
#' @export
sortPed = function(x, maxGen=100){
  alldone = FALSE
  pedsize = x@nInd
  done = matrix(data=FALSE,nrow=pedsize)
  gen = matrix(nrow=pedsize)
  g = 0
  
  while (!alldone) {
    g = g + 1
    alldone = TRUE
    for (i in 1:pedsize) {
      if (!done[i]) {
        if (((x@mother[i] == 0) || done[x@mother[i]]) && ((x@father[i] == 0) || done[x@father[i]])) {
          gen[i] = g
          done[i] = TRUE
        }
        else {
          alldone = FALSE
        }
      }
    }
    if (g > maxGen)
    {
      stop("Can not sort pedigree - loops in pedigree?")
    }
  }
  
  ids=character(x@nInd)
  mother=integer(x@nInd)
  father=integer(x@nInd)
  maxg = g
  c = 0
  for (g in 1:maxg) {
    for (i in 1:pedsize) {
      if (gen[i] == g) {
        c = c + 1
        ids[c] = x@ids[i]
        father[c] = x@father[i]
        mother[c] = x@mother[i]
      }
    }
  }
  return(new("Pedigree",nInd=pedsize,ids=ids,mother=mother,father=father))
}