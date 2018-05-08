#' @title New MapPop
#'
#' @description 
#' Creates a new \code{\link{MapPop-class}} from user supplied 
#' genetic maps and haplotypes.
#' 
#' @param genMaps a list of genetic maps
#' @param haplotypes a list of matrices or data.frames that 
#' can be coerced to matrices. See details.
#' @param inbred are individuals fully inbred
#' 
#' @details
#' Each item of genMaps must be a vector of ordered genetic lengths in 
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
#' genMaps = list(seq(0,1,length.out=11),
#'                seq(0,1,length.out=11))
#'                
#' # Create haplotypes for 10 outbred individuals
#' chr1 = sample(x=0:1,size=20*11,replace=TRUE)
#' chr1 = matrix(chr1,nrow=20,ncol=11)
#' chr2 = sample(x=0:1,size=20*11,replace=TRUE)
#' chr2 = matrix(chr2,nrow=20,ncol=11)
#' haplotypes = list(chr1,chr2)
#' 
#' founderPop = newMapPop(genMaps=genMaps,haplotypes=haplotypes)
#' 
#' @export
newMapPop = function(genMaps,haplotypes,inbred=FALSE){
  stopifnot(length(genMaps)==length(haplotypes))
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
  segSites = lapply(genMaps,length)
  segSites = unlist(segSites)
  if(!all.equal(nCol,segSites)){
    stop("Number of segregating sites in haplotypes and genMaps don't match")
  }
  output = vector("list",length(genMaps))
  for(chr in 1:length(genMaps)){
    geno = packHaplo(as.matrix(haplotypes[[chr]]),
                     ploidy=ploidy,inbred=inbred)
    output[[chr]] = new("MapPop",
                        nInd=as.integer(nInd),
                        nChr=1L,
                        ploidy=as.integer(ploidy),
                        nLoci=as.integer(segSites[chr]),
                        geno=as.matrix(list(geno)),
                        genMaps=as.matrix(genMaps[chr]))
  }
  output = do.call("c",output)
  return(output)
}

#' @title Haplotype tracking population
#' 
#' @description
#' Creates a population for tracking haplotypes.
#'
#' @param genMaps a list of genetic maps
#' @param nInd number of individuals
#' @param inbred should individuals be fully inbred
#' 
#' @details
#' Each item of genMaps must be a vector of ordered genetic lengths in 
#' Morgans. The first value must be zero. The length of the vector 
#' determines the number of segregating sites on the chromosome.
#' 
#' If inbred=FALSE, the value of nInd must be less than or equal to 
#' 128. Otherwise, it must be less than or equal to 256.
#' 
#' @examples
#' # Create genetic map for a single chromosome with 1 Morgan
#' # Chromosome contains 11 equally spaced segregating sites
#' genMaps = list(seq(0,1,length.out=11))
#' founderPop = trackHaploPop(genMaps=genMaps,nInd=10)
#' 
#' @export
trackHaploPop = function(genMaps,nInd,inbred=FALSE){
  stopifnot(is.list(genMaps))
  if(inbred){
    stopifnot(nInd<=128)
  }else{
    stopifnot(nInd<=256)
  }
  nInd = as.integer(nInd)
  nChr = length(genMaps)
  nLoci = unlist(lapply(genMaps,length))
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
               genMaps=as.matrix(genMaps))
  return(output)
}

#' @title Create founder genotypes using MaCS
#'
#' @description Uses an external programs MaCS and AlphaFormatter to produce initial founder genotypes.
#' 
#' @param nInd number of individuals to simulate
#' @param nChr number of chromosomes to simulate
#' @param segSites number of segregating sites to keep per chromosome. A 
#' value of NULL results in all sites being retained.
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
#' @return an object of \code{\link{MapPop-class}}
#' 
#' @examples 
#' # Creates a populations of 10 outbred individuals
#' # Their genome consists of 1 chromosome and 100 segregating sites
#' founderPop = runMacs(nInd=10,nChr=1,segSites=100)
#' 
#' @export
runMacs = function(nInd,nChr=1,segSites=NULL,inbred=FALSE,species="TEST",
                   split=NULL,manualCommand=NULL,manualGenLen=NULL){
  nInd = as.integer(nInd)
  ploidy = 2L #The only ploidy level currently supported
  if(is.null(segSites)){
    segSites = rep(0,nChr)
  }else if(length(segSites)==1){
    segSites = rep(segSites,nChr)
  }
  if(inbred){
    popSize = nInd
  }else{
    popSize = ploidy*nInd
  }
  if(!is.null(manualCommand)){
    if(is.null(manualGenLen)) stop("You must define manualGenLen")
    command = paste(popSize,manualCommand,"-s",sample.int(1e8,1))
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
    command = paste0(popSize," ",speciesParams,splitI," ",speciesHist,splitJ," -s ",sample.int(1e8,1))
  }
  output = vector("list",nChr)
  for(chr in 1:nChr){
    cat("Making chomosome",chr,"of",nChr,"\n")
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
                        genMaps=as.matrix(list(genMap)))
  }
  output = do.call("c",output)
  cat("Done\n")
  return(output)
}

#' @title Alternative wrapper for MaCS
#'
#' @description 
#' A wrapper function for \code{\link{runMacs}}. This wrapper 
#' is an alternative to directly using manualCommand in 
#' \code{\link{runMacs}}. It automatically creates an appropriate 
#' manualCommand based on user supplied variables.
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
#'
#' @return an object of \code{\link{MapPop-class}}
#' 
#' @examples 
#' # Creates a populations of 10 outbred individuals
#' # Their genome consists of 1 chromosome and 100 segregating sites
#' # The command is equivalent to using species="TEST" in runMacs
#' founderPop = runMacs2(nInd=10,nChr=1,segSites=100)
#' 
#' @export
runMacs2 = function(nInd,nChr,segSites,Ne=100,
                    bp=1e8,genLen=1,mutRate=2.5e-8,
                    histNe=NULL,histGen=NULL,
                    inbred=FALSE){
  stopifnot(length(histNe)==length(histGen))
  command = paste(bp,"-t",4*Ne*mutRate,
                  "-r",4*Ne*genLen/bp)
  if(length(histNe)>0){
    histNe = histNe/Ne
    histGen = histGen/(4*Ne)
    for(i in 1:length(histNe)){
      command = paste(command,"-eN",
                      histGen[i],histNe[i])
    }
  }
  return(runMacs(nInd=nInd,nChr=nChr,segSites=segSites,
                 inbred=inbred,species="TEST",split=NULL,
                 manualCommand=command,manualGenLen=genLen))
}

#' @title Sample haplotypes from a MapPop
#'
#' @description 
#' Creates a new \code{\link{MapPop-class}} from an existing 
#' \code{\link{MapPop-class}} by randomly sampling haplotypes.
#' 
#' @param nInd the number of individuals to create
#' @param mapPop the \code{\link{MapPop-class}} used to 
#' sample haplotypes
#' @param inbred should new individuals be fully inbred
#' @param replace should haplotypes be sampled with replacement
#' 
#' @return an object of \code{\link{MapPop-class}}
#' 
#' @examples 
#' # Create genetic map for a single chromosome with 1 Morgan
#' # Chromosome contains 11 equally spaced segregating sites
#' genMaps = list(seq(0,1,length.out=11))
#' founderPop = trackHaploPop(genMaps=genMaps,nInd=2,inbred=TRUE)
#' founderPop = sampleHaplo(nInd=20,mapPop=founderPop)
#' 
#' @export
sampleHaplo = function(nInd,mapPop,inbred=FALSE,replace=TRUE){
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
                        genMaps=as.matrix(mapPop@genMaps[chr]))
  }
  output = do.call("c",output)
  return(output)
}
