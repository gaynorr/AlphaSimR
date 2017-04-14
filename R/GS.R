#' @title Write data records
#' 
#' @description
#' Saves a population's phenotypic and marker data to a directory.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param dir directory for saving output
#' @param snpChip which SNP chip genotype to save. If useQtl=TRUE, this 
#' value will indicate which trait's QTL genotype to save.
#' @param useQtl should QTL genotype be written instead of SNP chip 
#' genotypes.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
writeRecords = function(pop,dir,snpChip,useQtl=FALSE,simParam=SIMPARAM){
  stopifnot(dir.exists(dir))
  if(useQtl){
    nMarkers = simParam@traits[[snpChip]]@nLoci
    markerType = paste0("QTL_",snpChip)
  }else{
    nMarkers = simParam@snpChips[[snpChip]]@nLoci
    markerType = paste0("SNP_",snpChip)
  }
  #Check that the marker set isn't being changed
  markerInfo = data.frame(markerType=markerType,
                          nMarkers=nMarkers,
                          stringsAsFactors=FALSE)
  markerInfoPath = file.path(dir,"markerInfo.txt")
  if(file.exists(markerInfoPath)){
    tmp = read.table(markerInfoPath,stringsAsFactors=FALSE)
    stopifnot(identical(tmp,markerInfo))
  }else{
    write.table(markerInfo,markerInfoPath)
  }
  
  
}