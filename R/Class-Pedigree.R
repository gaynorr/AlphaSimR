
# Pedigree ------------------------------------------------------------------

#' @title Pedigree
#' 
#' @description 
#' A pedigree
#' 
#' @slot nInd number of individuals

#' 
#' 
#' @export
setClass("Pedigree",
         slots=c(nInd="numeric",
                 ids="character",
                 mother="numeric",
                 father="numeric"))

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
  
  ids=character(ped@nInd)
  mother=numeric(ped@nInd)
  father=numeric(ped@nInd)
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