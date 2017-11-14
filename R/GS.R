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
                        includeHaplo=FALSE,append=TRUE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
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
      writeGeno(pop@geno,simParam@traits[[snpChip]]@lociPerChr,
                simParam@traits[[snpChip]]@lociLoc,
                file.path(dir,"genotype.txt"))
      if(includeHaplo){
        writeOneHaplo(pop@geno,simParam@traits[[snpChip]]@lociPerChr,
                      simParam@traits[[snpChip]]@lociLoc,1L,
                      file.path(dir,"haplotype1.txt"))
        writeOneHaplo(pop@geno,simParam@traits[[snpChip]]@lociPerChr,
                      simParam@traits[[snpChip]]@lociLoc,2L,
                      file.path(dir,"haplotype2.txt"))
      }
    }else{
      writeGeno(pop@geno,simParam@snpChips[[snpChip]]@lociPerChr,
                simParam@snpChips[[snpChip]]@lociLoc,
                file.path(dir,"genotype.txt"))
      if(includeHaplo){
        writeOneHaplo(pop@geno,simParam@snpChips[[snpChip]]@lociPerChr,
                      simParam@snpChips[[snpChip]]@lociLoc,1L,
                      file.path(dir,"haplotype1.txt"))
        writeOneHaplo(pop@geno,simParam@snpChips[[snpChip]]@lociPerChr,
                      simParam@snpChips[[snpChip]]@lociLoc,2L,
                      file.path(dir,"haplotype2.txt"))
      }
    }
  }
}

#' @title RR-BLUP Model
#'
#' @description
#' Fits a typical RR-BLUP model for genomic predictions.
#'
#' @param dir path to a directory with output from \code{\link{writeRecords}}
#' @param traits an integer indicating the trait or traits to model, or a
#' function of the traits returning a single value.
#' @param use train model using genetic value (\code{gv})
#' or phenotypes (\code{pheno}, default)
#' @param skip number of older records to skip
#' @param maxIter maximum number of iterations. Only used 
#' when number of traits is greater than 1.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
RRBLUP = function(dir, traits=1, use="pheno", 
                  skip=0, maxIter=1000, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  dir = normalizePath(dir, mustWork=TRUE)
  #Read and calculate basic information
  markerInfo = read.table(file.path(dir,"info.txt"),header=TRUE,
                          comment.char="",stringsAsFactors=FALSE)
  if(skip>0) markerInfo = markerInfo[-(1:skip),]
  nInd = nrow(markerInfo)
  nMarkers = scan(file.path(dir,"nMarkers.txt"),integer(),quiet=TRUE)
  markerType = scan(file.path(dir,"markerType.txt"),character(),quiet=TRUE)
  #Set trait/traits for genomic selection
  use = tolower(use)
  if(use == "gv"){
    y = scan(file.path(dir,"gv.txt"),numeric(),quiet=TRUE)
  }else if(use == "pheno"){
    y = scan(file.path(dir,"pheno.txt"),numeric(),quiet=TRUE)
  }else{
    stop(paste0("Use=",use," is not an option"))
  }
  y = matrix(y,nrow=nInd+skip,ncol=length(y)/(nInd+skip),byrow=TRUE)
  if(is.function(traits)){
    y = apply(y,1,traits)
    y = as.matrix(y)
  }else{
    y = y[,traits,drop=FALSE]
  }
  if(skip>0) y=y[-(1:skip),,drop=FALSE]
  #Fit model
  fixEff = as.integer(factor(markerInfo$fixEff))
  if(ncol(y)>1){
    ans = callRRBLUP_MV(y,fixEff,markerInfo$reps,
                        file.path(dir,"genotype.txt"),nMarkers,
                        skip,maxIter)
  }else{
    ans = callRRBLUP(y,fixEff,markerInfo$reps,
                     file.path(dir,"genotype.txt"),nMarkers,
                     skip)
  }
  tmp = unlist(strsplit(markerType,"_"))
  if(tmp[1]=="SNP"){
    markers = simParam@snpChips[[as.integer(tmp[2])]]
  }else{
    markers = simParam@traits[[as.integer(tmp[2])]]
  }
  markerEff=ans$u
  if(is.null(ans[["iter"]])){
    iter = 0
  }else{
    iter = ans$iter
  }
  output = new("RRsol",
               nLoci=markers@nLoci,
               lociPerChr=markers@lociPerChr,
               lociLoc=markers@lociLoc,
               markerEff=markerEff,
               fixEff=ans$beta,
               Vu=ans$Vu,
               Ve=ans$Ve,
               LL=ans$LL,
               iter=iter)
  return(output)
}

#' @title RR-BLUP Model with Dominance
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions that includes 
#' dominance effects.
#'
#' @param dir path to a directory with output from \code{\link{writeRecords}}
#' @param traits an integer indicating the trait to model, or a
#' function of the traits returning a single value.
#' @param use train model using genetic value (\code{gv})
#' or phenotypes (\code{pheno}, default)
#' @param skip number of older records to skip
#' @param maxIter maximum number of iterations.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
RRBLUP_D = function(dir, traits=1, use="pheno", 
                    skip=0, maxIter=1000, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  dir = normalizePath(dir, mustWork=TRUE)
  #Read and calculate basic information
  markerInfo = read.table(file.path(dir,"info.txt"),header=TRUE,
                          comment.char="",stringsAsFactors=FALSE)
  if(skip>0) markerInfo = markerInfo[-(1:skip),]
  nInd = nrow(markerInfo)
  nMarkers = scan(file.path(dir,"nMarkers.txt"),integer(),quiet=TRUE)
  markerType = scan(file.path(dir,"markerType.txt"),character(),quiet=TRUE)
  #Set trait/traits for genomic selection
  use = tolower(use)
  if(use == "gv"){
    y = scan(file.path(dir,"gv.txt"),numeric(),quiet=TRUE)
  }else if(use == "pheno"){
    y = scan(file.path(dir,"pheno.txt"),numeric(),quiet=TRUE)
  }else{
    stop(paste0("Use=",use," is not an option"))
  }
  y = matrix(y,nrow=nInd+skip,ncol=length(y)/(nInd+skip),byrow=TRUE)
  if(is.function(traits)){
    y = apply(y,1,traits)
    y = as.matrix(y)
  }else{
    y = y[,traits,drop=FALSE]
  }
  if(skip>0) y=y[-(1:skip),,drop=FALSE]
  stopifnot(ncol(y)==1)
  #Fit model
  fixEff = as.integer(factor(markerInfo$fixEff))
  ans = callRRBLUP_D(y,fixEff,markerInfo$reps,
                     file.path(dir,"genotype.txt"),nMarkers,
                     skip)
  tmp = unlist(strsplit(markerType,"_"))
  if(tmp[1]=="SNP"){
    markers = simParam@snpChips[[as.integer(tmp[2])]]
  }else{
    markers = simParam@traits[[as.integer(tmp[2])]]
  }
  output = new("RRDsol",
               nLoci=markers@nLoci,
               lociPerChr=markers@lociPerChr,
               lociLoc=markers@lociLoc,
               markerEff=ans$u[[1]],
               domEff=ans$u[[2]],
               fixEff=ans$beta,
               Vu=ans$Vu,
               Ve=ans$Ve,
               LL=ans$LL,
               iter=ans$iter)
  return(output)
}


#' @title RR-BLUP GCA Model
#'
#' @description
#' Fits an RR-BLUP model that estimates seperate marker effects for
#' the female and male gametes. Used for predicting GCA of parents
#' in single cross hybrids.
#'
#' @param dir path to a directory with output from \code{\link{writeRecords}}
#' @param traits an integer indicating the trait or traits to model, or a
#' function of the traits returning a single value.
#' @param use train model using genetic value (\code{gv})
#' or phenotypes (\code{pheno}, default)
#' @param skip number of older records to skip
#' @param maxIter maximum number of iterations for convergence.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
RRBLUP_GCA = function(dir, traits=1, use="pheno",
                      skip=0, maxIter=40, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  dir = normalizePath(dir, mustWork=TRUE)
  #Read and calculate basic information
  markerInfo = read.table(file.path(dir,"info.txt"),header=TRUE,
                          comment.char="",stringsAsFactors=FALSE)
  if(skip>0) markerInfo = markerInfo[-(1:skip),]
  nInd = nrow(markerInfo)
  nMarkers = scan(file.path(dir,"nMarkers.txt"),integer(),quiet=TRUE)
  markerType = scan(file.path(dir,"markerType.txt"),character(),quiet=TRUE)
  #Set trait/traits for genomic selection
  use = tolower(use)
  if(use == "gv"){
    y = scan(file.path(dir,"gv.txt"),numeric(),quiet=TRUE)
  }else if(use == "pheno"){
    y = scan(file.path(dir,"pheno.txt"),numeric(),quiet=TRUE)
  }else{
    stop(paste0("Use=",use," is not an option"))
  }
  y = matrix(y,nrow=nInd+skip,ncol=length(y)/(nInd+skip),byrow=TRUE)
  if(is.function(traits)){
    y = apply(y,1,traits)
    y = as.matrix(y)
  }else{
    y = y[,traits,drop=FALSE]
  }
  if(skip>0) y=y[-(1:skip),,drop=FALSE]
  stopifnot(ncol(y)==1)
  #Fit model
  fixEff = as.integer(factor(markerInfo$fixEff))
  ans = callRRBLUP_GCA(y,fixEff,markerInfo$reps,
                       file.path(dir,"haplotype1.txt"),
                       file.path(dir,"haplotype2.txt"),
                       nMarkers,skip,maxIter)
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
               fixEff=ans$beta,
               Vu=ans$Vu,
               Ve=ans$Ve,
               LL=ans$LL,
               iter=ans$iter)
  return(output)
}

#' @title RR-BLUP SCA Model
#'
#' @description
#' Fits an RR-BLUP model that models seperate effects for both female
#' and male gametes and dominance effects. Used for predicting single
#' cross hybrid performance.
#'
#' @param dir path to a directory with output from \code{\link{writeRecords}}
#' @param traits an integer indicating the trait or traits to model, or a
#' function of the traits returning a single value.
#' @param use train model using genetic value (\code{gv})
#' or phenotypes (\code{pheno}, default)
#' @param skip number of older records to skip
#' @param maxIter maximum number of iterations for convergence.
#' @param onFailGCA if true, \code{\link{RRBLUP_GCA}} is used if 
#' RRBLUP_SCA gives a variance component of zero
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @export
RRBLUP_SCA = function(dir, traits=1, use="pheno",
                      skip=0, maxIter=40, onFailGCA=TRUE, 
                      simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  dir = normalizePath(dir, mustWork=TRUE)
  #Read and calculate basic information
  markerInfo = read.table(file.path(dir,"info.txt"),header=TRUE,
                          comment.char="",stringsAsFactors=FALSE)
  if(skip>0) markerInfo = markerInfo[-(1:skip),]
  nInd = nrow(markerInfo)
  nMarkers = scan(file.path(dir,"nMarkers.txt"),integer(),quiet=TRUE)
  markerType = scan(file.path(dir,"markerType.txt"),character(),quiet=TRUE)
  #Set trait/traits for genomic selection
  use = tolower(use)
  if(use == "gv"){
    y = scan(file.path(dir,"gv.txt"),numeric(),quiet=TRUE)
  }else if(use == "pheno"){
    y = scan(file.path(dir,"pheno.txt"),numeric(),quiet=TRUE)
  }else{
    stop(paste0("Use=",use," is not an option"))
  }
  y = matrix(y,nrow=nInd+skip,ncol=length(y)/(nInd+skip),byrow=TRUE)
  if(is.function(traits)){
    y = apply(y,1,traits)
    y = as.matrix(y)
  }else{
    y = y[,traits,drop=FALSE]
  }
  if(skip>0) y=y[-(1:skip),,drop=FALSE]
  stopifnot(ncol(y)==1)
  #Fit model
  fixEff = as.integer(factor(markerInfo$fixEff))
  ans = callRRBLUP_SCA(y,fixEff,markerInfo$reps,
                       file.path(dir,"haplotype1.txt"),
                       file.path(dir,"haplotype2.txt"),
                       nMarkers,skip,maxIter)
  tmp = unlist(strsplit(markerType,"_"))
  if(tmp[1]=="SNP"){
    markers = simParam@snpChips[[as.integer(tmp[2])]]
  }else{
    markers = simParam@traits[[as.integer(tmp[2])]]
  }
  if(onFailGCA & any(ans$Vu<1e-10)){
    warning("using RRBLUP_GCA due to zero variance components")
    output = RRBLUP_GCA(dir=dir, traits=traits, use=use,
                        skip=skip, maxIter=maxIter, 
                        simParam=simParam)
    return(output)
  }
  output = new("SCAsol",
               nLoci=markers@nLoci,
               lociPerChr=markers@lociPerChr,
               lociLoc=markers@lociLoc,
               femaleEff=ans$u[[1]],
               maleEff=ans$u[[2]],
               scaEff=ans$u[[3]],
               fixEff=ans$beta,
               Vu=ans$Vu,
               Ve=ans$Ve,
               LL=ans$LL,
               iter=ans$iter)
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
#' @param gender either NULL, "male" or "female". If 
#' solution is \code{\link{GCAsol-class}} or 
#' \code{\link{SCAsol-class}} the EBV is the GCA if used in 
#' the corresponding pool
#' @param useD if model is \code{\link{RRDsol-class}}, should 
#' dominance be included in the EBV. If yes, the "EBV" is an 
#' estimate of genetic value and not an estimate of breeding value.
#' @param append should EBVs be appended to existing EBVs
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @export
setEBV = function(pop, solution, gender=NULL, useD=FALSE, 
                  append=FALSE){
  if(class(solution)=="RRsol"){
    ebv = gebvRR(solution, pop)
  }else if(class(solution)=="RRDsol"){
    if(useD){
      ebv = gebvRRD(solution, pop)
    }else{
      ebv = gebvRR(solution, pop)
    }
  }else if(class(solution)=="GCAsol"){
    if(is.null(gender)){
      ebv = gebvSCA(solution, pop, FALSE)
    }else if(toupper(gender)=="FEMALE"){
      ebv = gebvGCA(solution, pop, TRUE)
    }else if(toupper(gender)=="MALE"){
      ebv = gebvGCA(solution, pop, FALSE)
    }else{
      stop(paste0("gender=",gender," is not a valid option"))
    }
  }else if(class(solution)=="SCAsol"){
    if(is.null(gender)){
      ebv = gebvSCA(solution, pop)
    }else if(toupper(gender)=="FEMALE"){
      ebv = gebvGCA(solution, pop, TRUE, TRUE)
    }else if(toupper(gender)=="MALE"){
      ebv = gebvGCA(solution, pop, FALSE, TRUE)
    }else{
      stop(paste0("gender=",gender," is not a valid option"))
    }
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
#' and its related functions for a given training population size. 
#' Note that this functions may underestimate total usage.
#'
#' @param nInd the number of individuals in the training population
#' @param nMarker the number of markers per individual
#' @param model either "REG", "GCA", or "SCA" for \code{\link{RRBLUP}} 
#' \code{\link{RRBLUP_GCA}} and \code{\link{RRBLUP_SCA}} respectively.
#'
#' @return Returns an estimate for the required gigabytes of RAM
#'
#' @export
RRBLUPMemUse = function(nInd,nMarker,model="REG"){
  y = nInd
  X = nInd #times fixed effects, assuming 1 here
  M = nInd*nMarker
  u = nMarker
  if(toupper(model)=="REG"){
    S = nInd*nInd
    eigval = nInd
    eigvec = nInd*nInd
    eta = nInd
    Hinv = nInd*nInd
  }else if(toupper(model)=="GCA"){
    M = M*2
    V = nInd*nInd*3
    W = W0 = WQX= nInd*nInd
    WX = ee = nInd
    u = u*2
  }else if(toupper(model)=="SCA"){
    M = M*3
    V = nInd*nInd*4
    W = W0 = WQX= nInd*nInd
    WX = ee = nInd
    u = u*3
  }else{
    stop(paste0("model=",toupper(model)," not recognized"))
  }
  objects = ls()
  objects = objects[objects!="model"]
  bytes = sapply(objects,function(x) get(x))
  bytes = 8*sum(bytes)
  return(bytes*1e-9) #GB
}
