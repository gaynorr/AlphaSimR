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
#' @param reps number of reps for phenotypes. This values is used for modelling 
#' heterogenous error variance in genomic selection models. Leave value as 1 
#' unless using reps for phenotypes.
#' @param fixEff an integer indicating levels of fixed effect. Leave 
#' value as 1 if not using different levels of fixed effects.
#' @param includeHaplo should markers be seperated by female and male 
#' haplotypes.
#' @param append if true, new records are added to any existing records. 
#' If false, any existing records are deleted before writing new records. 
#' Note that this will delete all files in the 'dir' directory.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
writeRecords = function(pop,dir,snpChip,useQtl=FALSE,reps=1,fixEff=1,
                        includeHaplo=FALSE,append=TRUE,simParam=SIMPARAM){
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
      nMarkers = simParam@traits[[snpChip]]@nLoci
      markerType = paste("QTL",snpChip,sep="_")
    }else{
      nMarkers = simParam@snpChips[[snpChip]]@nLoci
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
                    reps=rep(reps,pop@nInd),fixEff=rep(fixEff,pop@nInd),
                    stringsAsFactors=FALSE)
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
      write.table(pullQtlGeno(pop,snpChip,simParam=simParam),
                  file.path(dir,"genotype.txt"),append=TRUE,
                  col.names=FALSE,row.names=FALSE)
      if(includeHaplo){
        write.table(pullQtlHaplo(pop,snpChip,haplo=1,simParam=simParam),
                    file.path(dir,"haplotype1.txt"),append=TRUE,
                    col.names=FALSE,row.names=FALSE)
        write.table(pullQtlHaplo(pop,snpChip,haplo=2,simParam=simParam),
                    file.path(dir,"haplotype2.txt"),append=TRUE,
                    col.names=FALSE,row.names=FALSE)
      }
    }else{
      write.table(pullSnpGeno(pop,snpChip,simParam=simParam),
                  file.path(dir,"genotype.txt"),append=TRUE,
                  col.names=FALSE,row.names=FALSE)
      if(includeHaplo){
        write.table(pullSnpHaplo(pop,snpChip,haplo=1,simParam=simParam),
                    file.path(dir,"haplotype1.txt"),append=TRUE,
                    col.names=FALSE,row.names=FALSE)
        write.table(pullSnpHaplo(pop,snpChip,haplo=2,simParam=simParam),
                    file.path(dir,"haplotype2.txt"),append=TRUE,
                    col.names=FALSE,row.names=FALSE)
      }
    }
  }
}

#' @title RR-BLUP Model
#' 
#' @description
#' Fits a typical RR-BLUP model for genomic predictions.
#'
#' @param dir path to a directory with output from /code{/link{writeRecords}}
#' @param traits an integer indicating the trait or traits to model, or a 
#' function of the traits returning a single value.
#' @param use train model using genetic value (\code{gv}) 
#' or phenotypes (\code{pheno}, default)
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
RRBLUP = function(dir, traits=1, use="pheno", simParam=SIMPARAM){
  dir = normalizePath(dir, mustWork=TRUE)
  #Read and calculate basic information
  markerInfo = read.table(file.path(dir,"info.txt"),header=TRUE,
                          comment.char="",stringsAsFactors=FALSE)
  nInd = nrow(markerInfo)
  nMarkers = scan(file.path(dir,"nMarkers.txt"),integer(),quiet=TRUE)
  markerType =scan(file.path(dir,"markerType.txt"),character(),quiet=TRUE)
  #Set trait/traits for genomic selection
  use = tolower(use)
  if(use == "gv"){
    y = scan(file.path(dir,"gv.txt"),numeric(),quiet=TRUE)
  }else if(use == "pheno"){
    y = scan(file.path(dir,"pheno.txt"),numeric(),quiet=TRUE)
  }else{
    stop(paste0("Use=",use," is not an option"))
  }
  y = matrix(y,nrow=nInd,ncol=length(y)/nInd,byrow=TRUE)
  if(is.function(traits)){
    y = apply(y,1,traits)
    y = as.matrix(y)
  }else{
    y = y[,traits,drop=FALSE]
  }
  #Fit model
  fixEff = as.integer(factor(markerInfo$fixEff))
  if(ncol(y)>1){
    ans = callRRBLUP_MV(y,fixEff,markerInfo$reps,
                           file.path(dir,"genotype.txt"),nMarkers)
  }else{
    ans = callRRBLUP(y,fixEff,markerInfo$reps,
                        file.path(dir,"genotype.txt"),nMarkers)
  }
  tmp = unlist(strsplit(markerType,"_"))
  if(tmp[1]=="SNP"){
    markers = simParam@snpChips[[as.integer(tmp[2])]]
  }else{
    markers = simParam@traits[[as.integer(tmp[2])]]
  }
  output = new("RRsol",
               nLoci=markers@nLoci,
               lociPerChr=markers@lociPerChr,
               lociLoc=markers@lociLoc,
               markerEff=ans$u,
               fixEff=ans$beta)
  return(output)
}

#' @title RR-BLUP GCA Model
#' 
#' @description
#' Fits an RR-BLUP model that estimates seperate marker effects for 
#' the female and male gametes. Used for predicting GCA of parents 
#' in single cross hybrids.
#'
#' @param dir path to a directory with output from /code{/link{writeRecords}}
#' @param traits an integer indicating the trait or traits to model, or a 
#' function of the traits returning a single value.
#' @param use train model using genetic value (\code{gv}) 
#' or phenotypes (\code{pheno}, default)
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
RRBLUP_GCA = function(dir, traits=1, use="pheno", simParam=SIMPARAM){
  dir = normalizePath(dir, mustWork=TRUE)
  #Read and calculate basic information
  markerInfo = read.table(file.path(dir,"info.txt"),header=TRUE,
                          comment.char="",stringsAsFactors=FALSE)
  nInd = nrow(markerInfo)
  nMarkers = scan(file.path(dir,"nMarkers.txt"),integer(),quiet=TRUE)
  markerType =scan(file.path(dir,"markerType.txt"),character(),quiet=TRUE)
  #Set trait/traits for genomic selection
  use = tolower(use)
  if(use == "gv"){
    y = scan(file.path(dir,"gv.txt"),numeric(),quiet=TRUE)
  }else if(use == "pheno"){
    y = scan(file.path(dir,"pheno.txt"),numeric(),quiet=TRUE)
  }else{
    stop(paste0("Use=",use," is not an option"))
  }
  y = matrix(y,nrow=nInd,ncol=length(y)/nInd,byrow=TRUE)
  if(is.function(traits)){
    y = apply(y,1,traits)
    y = as.matrix(y)
  }else{
    y = y[,traits,drop=FALSE]
  }
  stopifnot(ncol(y)==1)
  #Fit model
  fixEff = as.integer(factor(markerInfo$fixEff))
  ans = callRRBLUP_GCA(y,fixEff,markerInfo$reps,
                       file.path(dir,"haplotype1.txt"),
                       file.path(dir,"haplotype2.txt"),
                       nMarkers)
  tmp = unlist(strsplit(markerType,"_"))
  if(tmp[1]=="SNP"){
    markers = simParam@snpChips[[as.integer(tmp[2])]]
  }else{
    markers = simParam@traits[[as.integer(tmp[2])]]
  }
  output = new("GCAsol",
               nLoci=markers@nLoci,
               lociPerChr=markers@lociPerChr,
               lociLoc=markers@lociLoc,
               femaleEff=ans$u[[1]],
               maleEff=ans$u[[2]],
               fixEff=ans$beta)
  return(output)
}

#' @title RR-BLUP SCA Model
#' 
#' @description
#' Fits an RR-BLUP model that models seperate effects for both female 
#' and male gametes and dominance effects. Used for predicting single 
#' cross hybrid performance.
#'
#' @param dir path to a directory with output from /code{/link{writeRecords}}
#' @param traits an integer indicating the trait or traits to model, or a 
#' function of the traits returning a single value.
#' @param use train model using genetic value (\code{gv}) 
#' or phenotypes (\code{pheno}, default)
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
RRBLUP_SCA = function(dir, traits=1, use="pheno", simParam=SIMPARAM){
  dir = normalizePath(dir, mustWork=TRUE)
  #Read and calculate basic information
  markerInfo = read.table(file.path(dir,"info.txt"),header=TRUE,
                          comment.char="",stringsAsFactors=FALSE)
  nInd = nrow(markerInfo)
  nMarkers = scan(file.path(dir,"nMarkers.txt"),integer(),quiet=TRUE)
  markerType =scan(file.path(dir,"markerType.txt"),character(),quiet=TRUE)
  #Set trait/traits for genomic selection
  use = tolower(use)
  if(use == "gv"){
    y = scan(file.path(dir,"gv.txt"),numeric(),quiet=TRUE)
  }else if(use == "pheno"){
    y = scan(file.path(dir,"pheno.txt"),numeric(),quiet=TRUE)
  }else{
    stop(paste0("Use=",use," is not an option"))
  }
  y = matrix(y,nrow=nInd,ncol=length(y)/nInd,byrow=TRUE)
  if(is.function(traits)){
    y = apply(y,1,traits)
    y = as.matrix(y)
  }else{
    y = y[,traits,drop=FALSE]
  }
  stopifnot(ncol(y)==1)
  #Fit model
  fixEff = as.integer(factor(markerInfo$fixEff))
  ans = callRRBLUP_SCA(y,fixEff,markerInfo$reps,
                       file.path(dir,"haplotype1.txt"),
                       file.path(dir,"haplotype2.txt"),
                       nMarkers)
  tmp = unlist(strsplit(markerType,"_"))
  if(tmp[1]=="SNP"){
    markers = simParam@snpChips[[as.integer(tmp[2])]]
  }else{
    markers = simParam@traits[[as.integer(tmp[2])]]
  }
  output = new("SCAsol",
               nLoci=markers@nLoci,
               lociPerChr=markers@lociPerChr,
               lociLoc=markers@lociLoc,
               femaleEff=ans$u[[1]],
               maleEff=ans$u[[2]],
               scaEff=ans$u[[3]],
               fixEff=ans$beta)
  return(output)
}

#' @title Set EBV
#' 
#' @description
#' Sets a population's EBV with genomic estimated 
#' values from \code{\link{RRBLUP}}, \code{\link{RRBLUP_GCA}}, 
#' or \code{\link{RRBLUP_SCA}}.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param solution an object of \code{\link{RRsol-class}}, 
#' \code{\link{SCAsol-class}}, or \code{\link{GCAsol-class}}
#' @param gender either "male" or "female" if solution is 
#' \code{\link{GCAsol-class}}
#' @param append should EBVs be appended to existing EBVs
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
setEBV = function(pop, solution, gender=NULL, append=FALSE){
  if(class(solution)=="RRsol"){
    ebv = gebvRR(solution, pop)
  }else if(class(solution)=="GCAsol"){
    if(toupper(gender)=="FEMALE"){
      asFemale = TRUE
    }else if(toupper(gender)=="MALE"){
      asFemale = FALSE
    }else{
      stop("You must specify gender as 'male' or 'female' with class(solution)='GCAsol'")
    }
    ebv = gebvGCA(solution, pop, asFemale)
  }else if(class(solution)=="SCAsol"){
    ebv = gebvSCA(solution, pop)
  }else{
    stop("No method for class(solution)=",class(solution))
  }
  if(append){
    pop@ebv = cbind(pop@ebv,ebv)
  }else{
    pop@ebv = ebv
  }
  return(pop)
}

#' @title RRBLUP Memory Usage
#' 
#' @description
#' Estimates the amount of RAM needed to run the \code{\link{RRBLUP}} 
#' function for a given training population size. Note that this functions 
#' may underestimate total usage.
#' 
#' @param nInd the number of individuals in the training population
#' @param nMarker the number of markers per individual
#' 
#' @return Returns an estimate for the required gigabytes of RAM
#' 
#' @export
RRBLUPMemUse = function(nInd,nMarker){
  y = nInd
  X = nInd #times fixed effects, assuming 1 here
  Z = nInd*nMarker
  K = nMarker*nMarker
  S = nInd*nInd
  ZK = nInd*nMarker
  ZKZ = nInd*nInd
  eigval = nInd
  eigvec = nInd*nInd
  eta = nInd
  Hinv = nInd*nInd
  u = nMarker
  bytes = sapply(ls(),function(x) get(x))
  bytes = 8*sum(bytes)
  return(bytes*1e-9) #GB
}
