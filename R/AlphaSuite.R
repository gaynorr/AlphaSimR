#' @export
writeAlphaGenotypes = function(pop,file,chips,simParam=SIMPARAM) {
  genotypes = AlphaSimR::pullMultipleSnpGeno(pop,chips,simParam=simParam)
  names = sprintf("%20i",1:pop@nInd)
  write.table(genotypes,file=file,quote=FALSE,row.names=names,col.names=FALSE)
}

#' @export
writeAlphaHaplotypes = function(pop,file,chips,simParam=SIMPARAM) {
  haplotypes = AlphaSimR::pullMultipleSnpHaplo(pop,chips,simParam=simParam)
  names = sprintf("%20i",rep(1:pop@nInd,each=2))
  write.table(haplotypes,file=file,quote=FALSE,row.names=names,col.names=FALSE)
}