#' @title Loads a pedigree from file
#' 
#' @description
#' Loads a pedigree from file.  File should be space sperated with
#' three values per line: Individual ID, Father ID, Mother ID.
#' Father and Mother ID should be "0" if not in pedigree
#'
#' @param pedname file name
#'
#' @return Returns an object of \code{\link{Pedigree-class}}
#' 
#' @export
loadPedigreeFromFile = function(pedname) {
  fped = read.table(pedname,colClasses="character")
  
  pedsize = nrow(fped)
  
  ids = hashmap(keys=fped[,1], values = 1:pedsize)
  ids$insert("0",0)
  
  numped = matrix(nrow=pedsize,ncol=2)
  for (i in 1:pedsize) {
    numped[i,1] = ids[[fped[i,2]]]
    numped[i,2] = ids[[fped[i,3]]]
  }
  
  pedigree = new("Pedigree",nInd=nrow(fped),ids=fped[,1],father=numped[,1],mother=numped[,2])
  
  return(pedigree)
}