# A stop-gap solution until MaCS is integrated

#' @title Create founder genotypes using MaCS
#'
#' @description Uses an external programs MaCS and AlphaFormatter to produce initial founder genotypes.
#' 
#' @param macs path to MaCS
#' @param nInd number of individuals to simulate
#' @param nChr number of chromosomes to simulate
#' @param segSites number of segregating sites to keep per chromosome
#' @param inbred should founder individuals be inbred
#' @param species species history to simulate
#' @param split an optional historic population split in terms of generations ago.
#' @param manualCommand user provided MaCS options. For advanced users only.
#'
#' @export
runMacs = function(macs,nInd,nChr,segSites,inbred=TRUE,species="TEST",
                   split=NULL,manualCommand=NULL){
  ploidy = 2 #The only ploidy level currently supported
  if(inbred){
    popSize = nInd
  }else{
    popSize = ploidy*nInd
  }
  if(!is.null(manualCommand)){
    command = paste0(macs," ",popSize," ",manualCommand," -s ",sample.int(1e8,1)," 1>output.txt 2>/dev/null")
  }else{
    species = toupper(species)
    if(species=="WHEAT"){
      Ne = 50
      speciesParams = "800000000 -t 0.40E-06 -r 0.36E-06"
      speciesHist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 1.00 40 -eN 2.00 60 -eN 3.00 80 -eN 4.00 100 -eN 5.00 120 -eN 10.00 140 -eN 20.00 160 -eN 30.00 180 -eN 40.00 200 -eN 50.00 240 -eN 100.00 320 -eN 200.00 400 -eN 300.00 480 -eN 400.00 560 -eN 500.00 640"
    }else if(species=="MAIZE"){
      Ne = 100
      speciesParams = "200000000 -t 0.50E-05 -r 0.40E-05"
      speciesHist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 2.00 40 -eN 3.00 60 -eN 4.00 80 -eN 5.00 100"
    }else if(species=="TEST"){
      Ne = 100
      speciesParams = "100000000 -t 0.10E-04 -r 0.40E-05"
      speciesHist = ""
    }else{
      stop(paste("No rules for species",species))
    }
    if(is.null(split)){
      splitI = ""
      splitJ = ""
    }else{
      stopifnot(popSize%%2==0)
      splitI = paste(" -I 2",popSize%/%2,popSize%/%2)
      splitJ = paste(" -ej",split/(4*Ne)+0.000001,"2 1")
    }
    command = paste0(macs," ",popSize," ",speciesParams,splitI," ",speciesHist,splitJ," -s ",sample.int(1e8,1)," 1>output.txt 2>/dev/null")
  }
  currentDir = getwd()
  tmpDir = tempdir()
  setwd(tmpDir)
  errorHandler = function(e){
    setwd(currentDir)
    print(e)
    stop(paste("Output in directory",tmpDir))
  }
  genMaps = list()
  geno = list()
  for(chr in 1:nChr){
    cat("Making chomosome",chr,"of",nChr,"\n")
    cat("  Running MaCS...\n")
    system(command)
    cat("  Running AlphaFormatter...\n")
    errorInt = AlphaSimR:::AlphaFormatter()
    if(errorInt!=0){
      setwd(currentDir)
      cat("Error in AlphaFormatter","\n")
      stop(paste("Output in directory",tmpDir))
    }
    cat("  Reading output...\n")
    initSegSites = tryCatch(scan("SegSites.txt",what=integer(),quiet=TRUE),
                            error=errorHandler)
    if(length(initSegSites)!=1){
      setwd(currentDir)
      cat("length(intSegSites) =",length(initSegSites),"\n")
      stop(paste("Output in directory",tmpDir))
    }
    if(segSites>initSegSites){
      setwd(currentDir)
      stop(paste("Requested",segSites,"segSites but only created",initSegSites))
    }
    keep = sort(sample.int(initSegSites,segSites))
    map = tryCatch(scan("PhysicalMapInput.txt",what=numeric(),quiet=TRUE),
                   error=errorHandler)
    map = map[keep]
    map = map - map[1]
    tmpGeno = tryCatch(AlphaSimR:::readAF(nInd,initSegSites,ploidy,
                                          keep-1, #R to C++ index
                                          inbred),
                       error=errorHandler)
    genMaps[[chr]] = map
    geno[[chr]] = tmpGeno
  }
  setwd(currentDir)
  output = new("MapPop",nInd=as.integer(nInd),nChr=as.integer(nChr),
               ploidy=as.integer(ploidy),nLoci=as.integer(rep(segSites,nChr)),
               gender=rep("H",nInd),geno=geno,genMaps=genMaps)
  cat("Done\n")
  return(output)
}



