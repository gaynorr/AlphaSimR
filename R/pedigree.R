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
  pedigree = new("Pedigree",nInd=nrow(fped),ids=fped[,1],
                 mother=match(fped[,2],fped[,1],nomatch=0L),
                 father=match(fped[,3],fped[,1],nomatch=0L))
  return(pedigree)
}
