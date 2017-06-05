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
#' @param species species history to simulate. See details.
#' @param split an optional historic population split in terms of generations ago.
#' @param manualCommand user provided MaCS options. For advanced users only.
#' @param manualGenLen user provided genLen option for use with manual command.  For advanced users only.
#' 
#' @details
#' The current species histories are included: WHEAT, MAIZE, MAIZELANDRACE, CATTLE, 
#' PIG, CHICKEN, RABBIT and TEST. TEST uses MaCS's default history.
#'
#' @export
runMacs = function(macs,nInd,nChr,segSites,inbred=TRUE,species="TEST",
                   split=NULL,manualCommand=NULL,manualGenLen=NULL){
  ploidy = 2 #The only ploidy level currently supported
  if(inbred){
    popSize = nInd
  }else{
    popSize = ploidy*nInd
  }
  if(!is.null(manualCommand)){
    command = paste0(macs," ",popSize," ",manualCommand," -s ",sample.int(1e8,1))
    genLen = manualGenLen
  }else{
    species = toupper(species)
    if(species=="WHEAT"){ #WHEAT----
      genLen = 1.43
      Ne = 50
      speciesParams = "800000000 -t 0.40E-06 -r 0.36E-06"
      speciesHist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 1.00 40 -eN 2.00 60 -eN 3.00 80 -eN 4.00 100 -eN 5.00 120 -eN 10.00 140 -eN 20.00 160 -eN 30.00 180 -eN 40.00 200 -eN 50.00 240 -eN 100.00 320 -eN 200.00 400 -eN 300.00 480 -eN 400.00 560 -eN 500.00 640"
    }else if(species=="MAIZE"){ #MAIZE----
      genLen = 2.0
      Ne = 100
      speciesParams = "200000000 -t 0.50E-05 -r 0.40E-05"
      speciesHist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 2.00 40 -eN 3.00 60 -eN 4.00 80 -eN 5.00 100"
    }else if(species=="MAIZELANDRACE"){ #MAIZELANDRACE----
      genLen = 2.0
      Ne = 100
      speciesParams = "200000000 -t 0.50E-05 -r 0.40E-05"
      speciesHist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 2.00 40 -eN 3.00 60 -eN 4.00 80 -eN 5.00 100 -eN 6.00 120 -eN 7.00 140 -eN 8.00 160 -eN 9.00 180 -eN 10.00 200 -eN 12.50 400 -eN 15.00 600 -eN 17.50 800 -eN 20.00 1000 -eN 22.50 1200 -eN 25.00 1400 -eN 27.50 1600 -eN 30.00 2000"
    }else if(species=="CATTLE"){ #CATTLE----
      genLen = 1.0
      Ne = 100
      speciesParams = "100000000 -t 0.10E-04 -r 0.40E-05"
      speciesHist = "-eN 0.06 2.0 -eN 0.13 3.0 -eN 0.25 5.0 -eN 0.50 7.0 -eN 0.75 9.0 -eN 1.00 11.0 -eN 1.25 12.5 -eN 1.50 13.0 -eN 1.75 13.5 -eN 2.00 14.0 -eN 2.25 14.5 -eN 2.50 15.0 -eN 5.00 20.0 -eN 7.50 25.0 -eN 10.00 30.0 -eN 12.50 35.0 -eN 15.00 40.0 -eN 17.50 45.0 -eN 20.00 50.0 -eN 22.50 55.0 -eN 25.00 60.0 -eN 50.00 70.0 -eN 100.00 80.0 -eN 150.00 90.0 -eN 200.00 100.0 -eN 250.00 120.0 -eN 500.00 200.0 -eN 1000.00 400.0 -eN 1500.00 600.0 -eN 2000.00 800.0 -eN 2500.00 1000.0"
    }else if(species=="PIG"){ #PIG----
      genLen = 1.71
      Ne = 100
      speciesParams = "675000000 -t 0.95E-06 -r 0.10E-05"
      speciesHist = "-eN 25.00 100.0 -eN 50.00 200.0 -eN 75.00 300.0 -eN 100.00 400.0 -eN 125.00 500.0 -eN 150.00 600.0 -eN 175.00 700.0 -eN 200.00 800.0 -eN 225.00 900.0 -eN 250.00 1000.0 -eN 275.00 2000.0 -eN 300.00 3000.0 -eN 325.00 4000.0 -eN 350.00 5000.0 -eN 375.00 6000.0 -eN 400.00 7000.0 -eN 425.00 8000.0 -eN 450.00 9000.0 -eN 475.00 10000.0"
    }else if(species=="CHICKEN"){ #CHICKEN----
      genLen = 0.84
      Ne = 70
      speciesParams = "300000000 -t 0.23E-05 -r 0.78E-06"
      speciesHist = "-eN 0.18 0.71 -eN 0.36 1.43 -eN 0.54 2.14 -eN 0.71 2.86 -eN 0.89 3.57 -eN 1.07 4.29 -eN 1.25 5.00 -eN 1.43 5.71"
    }else if(species=="RABBIT"){ #RABBIT----
      genLen = 1.36
      Ne = 100
      speciesParams = "159000000 -t 0.44E-05 -r 0.34E-05"
      speciesHist = "-eN 0.05 1.25 -eN 0.08 1.50 -eN 0.10 1.75 -eN 0.13 2.00 -eN 0.15 2.25 -eN 0.18 2.50 -eN 0.20 2.75 -eN 0.23 3.00 -eN 0.25 3.25 -eN 0.50 4.00 -eN 1.00 5.00 -eN 1.50 6.00 -eN 2.00 7.00 -eN 2.50 8.00 -eN 3.00 90.00 -eN 3.50 10.00 -eN 4.00 11.00 -eN 4.50 12.00 -eN 5.00 1000.00"
    }else if(species=="TEST"){ #TEST----
      genLen = 1.0
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
    command = paste0(macs," ",popSize," ",speciesParams,splitI," ",speciesHist,splitJ," -s ",sample.int(1e8,1))
  }
  if(.Platform$OS.type=="windows"){
    command = paste0("powershell \"(",command,") 2>$null | out-file -filePath output.txt -encoding ASCII\"")
  }else{
    command = paste(command,"1>output.txt 2>/dev/null")
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
    errorInt = AlphaFormatter()
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
    map = map*genLen
    tmpGeno = tryCatch(readAF(nInd,initSegSites,ploidy,keep,inbred),
                       error=errorHandler)
    genMaps[[chr]] = map
    geno[[chr]] = tmpGeno
  }
  setwd(currentDir)
  output = new("MapPop",nInd=as.integer(nInd),nChr=as.integer(nChr),
               ploidy=as.integer(ploidy),nLoci=as.integer(rep(segSites,nChr)),
               gender=rep("H",nInd),geno=as.matrix(geno),genMaps=as.matrix(genMaps))
  cat("Done\n")
  return(output)
}



