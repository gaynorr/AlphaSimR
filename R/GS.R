#' @title Write data records
#' 
#' @description
#' Saves a population's phenotypic and marker data to a directory.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param dir path to a directory for saving output
#' @param snpChip which SNP chip genotype to save. If useQtl=TRUE, this 
#' value will indicate which trait's QTL genotype to save.
#' @param useQtl should QTL genotype be written instead of SNP chip 
#' genotypes.
#' @param reps number of reps for phenotypes. This values is used for modelling 
#' heterogenous error variance in genomic selection models. Leave value as 1 
#' unless using reps for phenotypes.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
writeRecords = function(pop,dir,snpChip,useQtl=FALSE,reps=1,simParam=SIMPARAM){
  stopifnot(dir.exists(dir))
  if(useQtl){
    nMarkers = simParam@traits[[snpChip]]@nLoci
    markerType = paste0("QTL_",snpChip)
  }else{
    nMarkers = simParam@snpChips[[snpChip]]@nLoci
    markerType = paste0("SNP_",snpChip)
  }
  #Check that the marker set isn't being changed
  nMarkerPath = file.path(dir,"nMarkers.txt")
  if(file.exists(nMarkerPath)){
    nMarkersDir = scan(nMarkerPath,integer(),quiet=TRUE)
    stopifnot(nMarkersDir==nMarkers)
  }else{
    writeLines(as.character(nMarkers),nMarkerPath)
  }
  markerTypePath = file.path(dir,"markerType.txt")
  if(file.exists(markerTypePath)){
    markerTypeDir = scan(markerTypePath,character(),quiet=TRUE)
    stopifnot(markerTypeDir==markerType)
  }else{
    writeLines(markerType,markerTypePath)
  }
  #Write info.txt
  info = data.frame(id=pop@id,mother=pop@mother,father=pop@father,
                    reps=rep(reps,pop@nInd),stringsAsFactors=FALSE)
  filePath = file.path(dir,"info.txt")
  if(file.exists(filePath)){
    write.table(info,filePath,append=TRUE,col.names=FALSE,
                row.names=FALSE)
  }else{
    write.table(info,filePath,row.names=FALSE)
  }
  #Write gv.txt
  write.table(pop@gv,file.path(dir,"gv.txt"),append=TRUE,
              col.names=FALSE,row.names=FALSE)
  #Write pheno.txt
  write.table(pop@pheno,file.path(dir,"pheno.txt"),append=TRUE,
              col.names=FALSE,row.names=FALSE)
  #Write markers.txt
  if(useQtl){
    write.table(pullQtlGeno(pop,snpChip,simParam=simParam),
                file.path(dir,"markers.txt"),append=TRUE,
                col.names=FALSE,row.names=FALSE)
  }else{
    write.table(pullSnpGeno(pop,snpChip,simParam=simParam),
                file.path(dir,"markers.txt"),append=TRUE,
                col.names=FALSE,row.names=FALSE)
  }
}