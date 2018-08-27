#' @title New MapPop
#'
#' @description 
#' Creates a new \code{\link{MapPop-class}} from user supplied 
#' genetic maps and haplotypes.
#' 
#' @param genMap a list of genetic maps
#' @param haplotypes a list of matrices or data.frames that 
#' can be coerced to matrices. See details.
#' @param inbred are individuals fully inbred
#' 
#' @details
#' Each item of genMap must be a vector of ordered genetic lengths in 
#' Morgans. The first value must be zero. The length of the vector 
#' determines the number of segregating sites on the chromosome.
#' 
#' Each item of haplotypes must be coercible to a matrix. The columns 
#' of this matrix correspond to segregating sites and their number must 
#' match
#' 
#' @return an object of \code{\link{MapPop-class}}
#' 
#' @examples 
#' # Create genetic map for two chromosomes, each 1 Morgan long
#' # Each chromosome contains 11 equally spaced segregating sites
#' genMap = list(seq(0,1,length.out=11),
#'                seq(0,1,length.out=11))
#'                
#' # Create haplotypes for 10 outbred individuals
#' chr1 = sample(x=0:1,size=20*11,replace=TRUE)
#' chr1 = matrix(chr1,nrow=20,ncol=11)
#' chr2 = sample(x=0:1,size=20*11,replace=TRUE)
#' chr2 = matrix(chr2,nrow=20,ncol=11)
#' haplotypes = list(chr1,chr2)
#' 
#' founderPop = newMapPop(genMap=genMap,haplotypes=haplotypes)
#' 
#' @export
newMapPop = function(genMap,haplotypes,inbred=FALSE){
  stopifnot(length(genMap)==length(haplotypes))
  ploidy = 2 #The only ploidy level currently supported
  nRow = lapply(haplotypes,nrow)
  nRow = unlist(nRow)
  if(length(nRow)>1L){
    if(any(nRow[1]!=nRow)){
      stop("Number of rows must be equal in haplotypes")
    }
    nRow = nRow[1]
  }
  if(inbred){
    nInd = nRow
  }else{
    if(ploidy==2L){
      if(nRow%%2 == 1L){
        stop("Number of haplotypes must be divisible by 2")
      }
      nInd = nRow/2
    }
  }
  nCol = lapply(haplotypes,ncol)
  nCol = unlist(nCol)
  segSites = lapply(genMap,length)
  segSites = unlist(segSites)
  if(!all.equal(nCol,segSites)){
    stop("Number of segregating sites in haplotypes and genMap don't match")
  }
  output = vector("list",length(genMap))
  for(chr in 1:length(genMap)){
    geno = packHaplo(as.matrix(haplotypes[[chr]]),
                     ploidy=ploidy,inbred=inbred)
    output[[chr]] = new("MapPop",
                        nInd=as.integer(nInd),
                        nChr=1L,
                        ploidy=as.integer(ploidy),
                        nLoci=as.integer(segSites[chr]),
                        geno=as.matrix(list(geno)),
                        genMap=as.matrix(genMap[chr]))
  }
  output = do.call("c",output)
  return(output)
}

#' @title Haplotype tracking population
#' 
#' @description
#' Creates a population contain haplotypes numbered for 
#' identity be descent tracking.
#'
#' @param genMap a list of genetic maps
#' @param nInd number of individuals
#' @param inbred should individuals be fully inbred
#' 
#' @details
#' Each item of genMap must be a vector of ordered genetic lengths in 
#' Morgans. The first value must be zero. The length of the vector 
#' determines the number of segregating sites on the chromosome.
#' 
#' If inbred=FALSE, the value of nInd must be less than or equal to 
#' 128. Otherwise, it must be less than or equal to 256.
#' 
#' @examples
#' # Create genetic map for a single chromosome with 1 Morgan
#' # Chromosome contains 11 equally spaced segregating sites
#' genMap = list(seq(0,1,length.out=11))
#' founderPop = trackHaploPop(genMap=genMap,nInd=10)
#' 
#' @export
trackHaploPop = function(genMap,nInd,inbred=FALSE){
  stopifnot(is.list(genMap))
  if(inbred){
    stopifnot(nInd<=128)
  }else{
    stopifnot(nInd<=256)
  }
  nInd = as.integer(nInd)
  nChr = length(genMap)
  nLoci = unlist(lapply(genMap,length))
  geno = vector("list",nChr)
  for(i in 1:nChr){
    tmpGeno = as.raw(0:(2*nInd-1))
    tmpGeno = array(raw(),dim=c(nLoci[i],2L,nInd))
    tmp=-1
    for(j in 1:nInd){
      if(inbred){
        tmp=tmp+1
        tmpGeno[,1:2,j] = as.raw(tmp)
      }else{
        for(k in 1:2){
          tmp=tmp+1
          tmpGeno[,k,j] = as.raw(tmp)
        } 
      }
    }
    geno[[i]] = tmpGeno
  }
  output = new("MapPop",nInd=nInd,nChr=nChr,ploidy=2L,
               nLoci=nLoci,geno=as.matrix(geno),
               genMap=as.matrix(genMap))
  return(output)
}

#' @title Create founder haplotypes using MaCS
#'
#' @description Uses the MaCS software to produce founder haplotypes.
#' 
#' @param nInd number of individuals to simulate
#' @param nChr number of chromosomes to simulate
#' @param segSites number of segregating sites to keep per chromosome. A 
#' value of NULL results in all sites being retained.
#' @param inbred should founder individuals be inbred
#' @param species species history to simulate. See details.
#' @param split an optional historic population split in terms of generations ago.
#' @param manualCommand user provided MaCS options. For advanced users only.
#' @param manualGenLen user provided genetic length. This must be supplied if using 
#' manualCommand. If not using manualCommand, this value will replace the predefined 
#' genetic length for the species. However, this the genetic length is only used by 
#' AlphaSimR and is not passed to MaCS, so MaCS still uses the predefined genetic length. 
#' For advanced users only.
#' @param suppressMessages should messages on status be suppressed
#' 
#' @details
#' The current species histories are included: GENERIC, CATTLE, WHEAT, MAIZE,  
#' and EUROPEAN. 
#'
#' @return an object of \code{\link{MapPop-class}}
#' 
#' @examples 
#' # Creates a populations of 10 outbred individuals
#' # Their genome consists of 1 chromosome and 100 segregating sites
#' founderPop = runMacs(nInd=10,nChr=1,segSites=100)
#' 
#' @export
runMacs = function(nInd,nChr=1,segSites=NULL,inbred=FALSE,species="GENERIC",
                   split=NULL,manualCommand=NULL,manualGenLen=NULL,
                   suppressMessages=FALSE){
  nInd = as.integer(nInd)
  ploidy = 2L #The only ploidy level currently supported
  if(is.null(segSites)){
    segSites = rep(0,nChr)
  }else if(length(segSites)==1){
    segSites = rep(segSites,nChr)
  }
  popSize = ifelse(inbred,nInd,ploidy*nInd)
  if(!is.null(manualCommand)){
    if(is.null(manualGenLen)) stop("You must define manualGenLen")
    command = paste(popSize,manualCommand,"-s",sample.int(1e8,1))
    genLen = manualGenLen
  }else{
    species = toupper(species)
    if(species=="GENERIC"){ #GENERIC----
      genLen = 1.0
      Ne = 100
      speciesParams = "1E8 -t 1E-5 -r 4E-6"
      speciesHist = "-eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0"
    }else if(species=="CATTLE"){ #CATTLE----
      genLen = 1.0
      Ne = 90
      speciesParams = "1E8 -t 9E-6 -r 3.6E-6"
      speciesHist = "-eN 0.011 1.33 -eN 0.019 2.78 -eN 0.036 3.89 -eN 0.053 11.11 -eN 0.069 16.67 -eN 0.431 22.22 -eN 1.264 27.78 -eN 1.819 38.89 -eN 4.875 77.78 -eN 6.542 111.11 -eN 9.319 188.89 -eN 92.097 688.89 -eN 2592.097 688.89"
    }else if(species=="WHEAT"){ #WHEAT----
      genLen = 1.43
      Ne = 50
      speciesParams = "8E8 -t 4E-7 -r 3.6E-7"
      speciesHist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 1.00 40 -eN 2.00 60 -eN 3.00 80 -eN 4.00 100 -eN 5.00 120 -eN 10.00 140 -eN 20.00 160 -eN 30.00 180 -eN 40.00 200 -eN 50.00 240 -eN 100.00 320 -eN 200.00 400 -eN 300.00 480 -eN 400.00 560 -eN 500.00 640"
    }else if(species=="MAIZE"){ #MAIZE----
      genLen = 2.0
      Ne = 100
      speciesParams = "2E8 -t 5E-6 -r 4E-6"
      speciesHist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 2.00 40 -eN 3.00 60 -eN 4.00 80 -eN 5.00 100" 
    }else if(species=="EUROPEAN"){ #EUROPEAN----
      genLen = 1.3
      Ne = 512000
      speciesParams = "1.3E8 -t 0.0483328 -r 0.02054849"
      speciesHist = "-G 1.0195 -eG 0.0001000977 1.0031 -eN 0.0004492188 0.002015625 -eN 0.000449707 0.003634766"
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
    command = paste0(popSize," ",speciesParams,splitI," ",speciesHist,splitJ," -s ",sample.int(1e8,1))
  }
  if(!is.null(manualGenLen)){
    genLen = manualGenLen
  }
  output = vector("list",nChr)
  for(chr in 1:nChr){
    if(!suppressMessages){
      cat("Making chomosome",chr,"of",nChr,"\n")
    }
    macsOut = MaCS(command,segSites[chr])
    genMap = c(macsOut$genMap)
    genMap = genLen*(genMap-min(genMap))
    geno = packHaplo(macsOut$haplo,ploidy=ploidy,
                     inbred=inbred)
    output[[chr]] = new("MapPop",
                        nInd=nInd,
                        nChr=1L,
                        ploidy=ploidy,
                        nLoci=dim(geno)[1],
                        geno=as.matrix(list(geno)),
                        genMap=as.matrix(list(genMap)))
  }
  output = do.call("c",output)
  if(!suppressMessages){
    cat("Done\n")
  }
  return(output)
}

#' @title Alternative wrapper for MaCS
#'
#' @description 
#' A wrapper function for \code{\link{runMacs}}. This wrapper is designed 
#' to be easier to use than supply custom comands to manualCommand in 
#' \code{\link{runMacs}}. It effectively automates the creation of an 
#' appropriate manualCommand using user supplied variables, but only deals  
#' with a subset of the possibilities. The defaults were chosen to match 
#' species="GENERIC" in \code{\link{runMacs}}.
#' 
#' @param nInd number of individuals to simulate
#' @param nChr number of chromosomes to simulate
#' @param segSites number of segregating sites to keep per chromosome
#' @param Ne effective population size
#' @param bp base pair length of chromosome
#' @param genLen genetic length of chromosome in Morgans
#' @param mutRate per base pair mutation rate
#' @param histNe effective population size in previous 
#' generations
#' @param histGen number of generations ago for effective 
#' population sizes given in histNe
#' @param inbred should founder individuals be inbred
#' @param split an optional historic population split in terms of generations ago
#' @param returnCommand should the command passed to manualCommand in 
#' \code{\link{runMacs}} be returned. If TRUE, MaCS will not be called and 
#' the command is returned instead.
#' @param suppressMessages should messages on status be suppressed
#'
#' @return an object of \code{\link{MapPop-class}} or if 
#' returnCommand is true a string giving the MaCS command passed 
#' the manualCommand argument of \code{\link{runMacs}}.
#' 
#' @examples 
#' # Creates a populations of 10 outbred individuals
#' # Their genome consists of 1 chromosome and 100 segregating sites
#' # The command is equivalent to using species="GENERIC" in runMacs
#' founderPop = runMacs2(nInd=10,nChr=1,segSites=100)
#' 
#' @export
runMacs2 = function(nInd,nChr=1,segSites=NULL,Ne=100,
                    bp=1e8,genLen=1,mutRate=2.5e-8,
                    histNe=c(500,1500,6000,12000,100000),
                    histGen=c(100,1000,10000,100000,1000000),
                    inbred=FALSE,split=NULL,returnCommand=FALSE,
                    suppressMessages=FALSE){
  stopifnot(length(histNe)==length(histGen))
  speciesParams = paste(bp,"-t",4*Ne*mutRate,
                        "-r",4*Ne*genLen/bp)
  speciesHist = ""
  if(length(histNe)>0){
    histNe = histNe/Ne
    histGen = histGen/(4*Ne)
    for(i in 1:length(histNe)){
      speciesHist = paste(speciesHist,"-eN",
                      histGen[i],histNe[i])
    }
  }
  if(is.null(split)){
    command = paste(speciesParams,speciesHist)
  }else{
    popSize = ifelse(inbred,nInd,2*nInd)
    command = paste(speciesParams,
                    paste("-I 2",popSize%/%2,popSize%/%2),
                    speciesHist,
                    paste("-ej",split/(4*Ne)+0.000001,"2 1"))
  }
  if(returnCommand){
    return(command)
  }
  return(runMacs(nInd=nInd,nChr=nChr,segSites=segSites,
                 inbred=inbred,species="TEST",split=NULL,
                 manualCommand=command,manualGenLen=genLen,
                 suppressMessages=suppressMessages))
}

#' @title Sample haplotypes from a MapPop
#'
#' @description 
#' Creates a new \code{\link{MapPop-class}} from an existing 
#' \code{\link{MapPop-class}} by randomly sampling haplotypes.
#' 
#' @param mapPop the \code{\link{MapPop-class}} used to 
#' sample haplotypes
#' @param nInd the number of individuals to create
#' @param inbred should new individuals be fully inbred
#' @param replace should haplotypes be sampled with replacement
#' 
#' @return an object of \code{\link{MapPop-class}}
#' 
#' @examples 
#' # Create genetic map for a single chromosome with 1 Morgan
#' # Chromosome contains 11 equally spaced segregating sites
#' genMap = list(seq(0,1,length.out=11))
#' founderPop = trackHaploPop(genMap=genMap,nInd=2,inbred=TRUE)
#' founderPop = sampleHaplo(nInd=20,mapPop=founderPop)
#' 
#' @export
sampleHaplo = function(mapPop,nInd,inbred=FALSE,replace=TRUE){
  nHaplo = mapPop@nInd*mapPop@ploidy
  if(inbred){
    nSamp = nInd
  }else{
    nSamp = nInd*mapPop@ploidy
  }
  output = vector("list",mapPop@nChr)
  for(chr in 1:mapPop@nChr){
    haplo = sample.int(nHaplo,nSamp,replace=replace)
    geno = array(data=as.raw(0),
                 dim=c(mapPop@nLoci[chr],
                       mapPop@ploidy,nInd))
    outHap = 1L 
    outInd = 1L
    for(i in 1:length(haplo)){
      inHap = (haplo[i]-1L)%%mapPop@ploidy + 1L
      inInd = (haplo[i]-1L)%/%mapPop@ploidy + 1L
      if(inbred){
        for(outHap in 1:mapPop@ploidy){
          geno[,outHap,outInd] = 
            mapPop@geno[[chr]][,inHap,inInd]
        }
        outInd = outInd+1L
      }else{
        geno[,outHap,outInd] = 
          mapPop@geno[[chr]][,inHap,inInd]
        outHap = outHap%%mapPop@ploidy+1L
        if(outHap==1L){
          outInd = outInd+1L
        }
      }
    }
    output[[chr]] = new("MapPop",
                        nInd=as.integer(nInd),
                        nChr=1L,
                        ploidy=mapPop@ploidy,
                        nLoci=mapPop@nLoci[chr],
                        geno=as.matrix(list(geno)),
                        genMap=as.matrix(mapPop@genMap[chr]))
  }
  output = do.call("c",output)
  return(output)
}
