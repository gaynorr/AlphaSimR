#' @title Write data records
#'
#' @description
#' Saves a population's phenotypic and marker data to a directory.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param dir path to a directory for saving output
#' @param snpChip which SNP chip genotype to save. If useQtl=TRUE, this
#' value will indicate which trait's QTL genotype to save. A value of
#' 0 will skip writing a snpChip.
#' @param useQtl should QTL genotype be written instead of SNP chip
#' genotypes.
#' @param includeHaplo should markers be separated by female and male
#' haplotypes.
#' @param append if true, new records are added to any existing records.
#' If false, any existing records are deleted before writing new records.
#' Note that this will delete all files in the 'dir' directory.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @export
writeRecords = function(pop,dir,snpChip=1,useQtl=FALSE,
                        includeHaplo=FALSE,append=TRUE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  snpChip = as.integer(snpChip)
  dir = normalizePath(dir, mustWork=TRUE)
  if(!append){
    #Delete any existing files
    tmp = list.files(dir,full.names=TRUE)
    if(length(tmp)>0){
      unlink(tmp,recursive=TRUE)
    }
  }
  if(snpChip==0){
    nMarkers = 0
    markerType = "NULL"
  }else{
    if(useQtl){
      nMarkers = simParam$traits[[snpChip]]@nLoci
      markerType = paste("QTL",snpChip,sep="_")
    }else{
      nMarkers = simParam$snpChips[[snpChip]]@nLoci
      markerType = paste("SNP",snpChip,sep="_")
    }
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
                    fixEff=pop@fixEff,stringsAsFactors=FALSE)
  filePath = file.path(dir,"info.txt")
  if(file.exists(filePath)){
    write.table(info,filePath,append=TRUE,col.names=FALSE,
                row.names=FALSE,quote=FALSE)
  }else{
    write.table(info,filePath,row.names=FALSE,quote=FALSE)
  }
  #Write gv.txt
  write.table(pop@gv,file.path(dir,"gv.txt"),append=TRUE,
              col.names=FALSE,row.names=FALSE)
  #Write pheno.txt
  write.table(pop@pheno,file.path(dir,"pheno.txt"),append=TRUE,
              col.names=FALSE,row.names=FALSE)
  #Write genotype.txt, unless snpChip=0
  if(snpChip!=0){
    if(useQtl){
      writeGeno(pop@geno,simParam$traits[[snpChip]]@lociPerChr,
                simParam$traits[[snpChip]]@lociLoc,
                file.path(dir,"genotype.txt"),simParam$nThreads)
      if(includeHaplo){
        writeOneHaplo(pop@geno,simParam$traits[[snpChip]]@lociPerChr,
                      simParam$traits[[snpChip]]@lociLoc,1L,
                      file.path(dir,"haplotype1.txt"),simParam$nThreads)
        writeOneHaplo(pop@geno,simParam$traits[[snpChip]]@lociPerChr,
                      simParam$traits[[snpChip]]@lociLoc,2L,
                      file.path(dir,"haplotype2.txt"),simParam$nThreads)
      }
    }else{
      writeGeno(pop@geno,simParam$snpChips[[snpChip]]@lociPerChr,
                simParam$snpChips[[snpChip]]@lociLoc,
                file.path(dir,"genotype.txt"),simParam$nThreads)
      if(includeHaplo){
        writeOneHaplo(pop@geno,simParam$snpChips[[snpChip]]@lociPerChr,
                      simParam$snpChips[[snpChip]]@lociLoc,1L,
                      file.path(dir,"haplotype1.txt"),simParam$nThreads)
        writeOneHaplo(pop@geno,simParam$snpChips[[snpChip]]@lociPerChr,
                      simParam$snpChips[[snpChip]]@lociLoc,2L,
                      file.path(dir,"haplotype2.txt"),simParam$nThreads)
      }
    }
  }
}
