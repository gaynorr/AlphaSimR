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