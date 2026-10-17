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
#' @param ploidy ploidy level of the organism
#' @param genLen an optional vector of chromosome genetic lengths in
#' Morgans, either one value or one per chromosome. Each must be at
#' least the position of the chromosome's last locus. If NULL, each
#' chromosome ends at its last locus.
#'
#' @details
#' Each item of genMap must be a vector of ordered genetic positions in
#' Morgans. Positions are absolute, measured from the start of the
#' chromosome, so they must be non-negative but the first need not be
#' zero. The length of the vector determines the number of segregating
#' sites on the chromosome.
#'
#' A chromosome runs from position 0 to its genetic length, which can
#' extend past the last locus. Crossovers are simulated over the whole
#' chromosome. The centromere is placed at half of the genetic length.
#'
#' Each item of haplotypes must be coercible to a matrix. The columns
#' of this matrix correspond to segregating sites. The number of rows
#' must match the number of individuals times the ploidy if using
#' inbred=FALSE. If using inbred=TRUE, the number of rows must equal
#' the number of individuals. The haplotypes can be stored as numeric,
#' integer or raw. The underlying C++ function will use raw.
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
#' founderPop = newMapPop(genMap=genMap, haplotypes=haplotypes)
#'
#' @export
newMapPop = function(genMap, haplotypes, inbred=FALSE,
                     ploidy=2L, genLen=NULL){
  stopifnot("genMap and haplotypes must have one element per chromosome each"=
              length(genMap)==length(haplotypes))
  .checkGenMap(genMap)
  genLen = .checkGenLen(genLen, genMap)
  nRow = lapply(haplotypes,nrow)
  nRow = unlist(nRow)
  ploidy = as.integer(ploidy)
  if(length(nRow)>1L){
    if(any(nRow[1]!=nRow)){
      stop("Number of rows must be equal in haplotypes")
    }
    nRow = nRow[1]
  }
  if(inbred){
    nInd = nRow
  }else{
    if(nRow%%ploidy != 0L){
      stop("Number of haplotypes must be divisible by ploidy")
    }
    nInd = nRow/ploidy
  }
  nCol = lapply(haplotypes,ncol)
  nCol = unlist(nCol)
  segSites = lapply(genMap,length)
  segSites = unlist(segSites)
  if(!all(nCol == segSites)){
    stop("Number of segregating sites in haplotypes and genMap don't match")
  }
  output = vector("list",length(genMap))
  for(chr in seq_len(length(genMap))){
    geno = packHaplo(as.matrix(haplotypes[[chr]]),
                     ploidy=ploidy,inbred=inbred)
    output[[chr]] = new("MapPop",
                        nInd=as.integer(nInd),
                        nChr=1L,
                        ploidy=ploidy,
                        nLoci=as.integer(segSites[chr]),
                        geno=list(geno),
                        genMap=genMap[chr],
                        centromere=genLen[chr]/2,
                        inbred=inbred,
                        genLen=genLen[chr])
    if(is.null(names(output[[chr]]@genMap))){
      names(output[[chr]]@genMap) = as.character(chr)
    }
    if(is.null(names(output[[chr]]@genMap[[1]]))){
      names(output[[chr]]@genMap[[1]]) = paste(chr,1:segSites[chr],sep="_")
    }
  }
  output = do.call("cChr",output)
  return(output)
}

# Checks a list of chromosome genetic maps. Positions are absolute, measured
# from the start of the chromosome, so they are validated rather than shifted
# to start at zero: on every chromosome they must be finite, non-negative and
# in non-decreasing order.
.checkGenMap = function(genMap){
  for(i in seq_along(genMap)){
    x = genMap[[i]]
    if(!is.numeric(x) || !all(is.finite(x))){
      stop("Genetic map positions must be finite numbers")
    }
    if(any(x<0)){
      stop("Genetic map positions must be non-negative")
    }
    if(is.unsorted(x)){
      stop("Genetic map positions must be in non-decreasing order")
    }
  }
  invisible(NULL)
}

# Checks the chromosome genetic lengths given for a genetic map and returns
# one per chromosome. A chromosome runs from position 0 to its genetic
# length, which can lie past its last locus but not before it. NULL means
# each chromosome ends at its last locus, and a single value is used for
# every chromosome.
.checkGenLen = function(genLen, genMap){
  lastLocus = .lastLocus(genMap)
  if(is.null(genLen)){
    return(lastLocus)
  }
  if(length(genLen)==1L){
    genLen = rep(genLen, length(genMap))
  }
  if(!is.numeric(genLen) || length(genLen)!=length(genMap) ||
     !all(is.finite(genLen))){
    stop("genLen must be finite numbers, one value or one per chromosome")
  }
  if(any(genLen<lastLocus)){
    stop("genLen must be at least the position of the last locus")
  }
  return(unname(as.numeric(genLen)))
}

#' @title Create founder haplotypes using MaCS
#'
#' @description Uses the `MaCS` software to produce founder haplotypes
#' \insertCite{chen2009Fast}{AlphaSimR}.
#'
#' @param nInd number of individuals to simulate.
#' @param nChr number of chromosomes to simulate (using independent simulations)
#' @param segSites number of segregating sites to keep per chromosome. A
#' value of \code{NULL} results in retaining all simulated segregating sites.
#' Requesting more sites than were generated raises an error.
#' @param inbred should founder individuals be inbred by replicating their
#' whole-chromosome haplotypes (see also \code{ploidy}).
#' This creates fully homozygous founders without simulating a history of selfing.
#' @param species demographic model to simulate (see details).
#' @param split an optional historic population split in terms of generations ago.
#' @param ploidy level of species to simulate \code{nInd * ploidy}
#' whole-chromosome haplotypes per chromosome unless \code{inbred=TRUE},
#' where \code{nInd} haplotypes are simulated and each is copied across
#' all homologues of an individual. Unlike \code{\link{runMacs2}}, this
#' does not rescale the demographic model's effective population sizes
#' when ploidy changes.
#' @param manualCommand user provided `MaCS` options. For advanced users only.
#' See the \code{vignette("DemographicModels", package = "AlphaSimR")} for details.
#' @param manualGenLen user provided genetic length. This must be supplied if using
#' \code{manualCommand}. If not using \code{manualCommand},
#' this value will replace the predefined genetic length for the species.
#' However, this genetic length is only used by AlphaSimR and is not passed to `MaCS`,
#' so `MaCS` still uses the predefined genetic length. For advanced users only.
#' @param nThreads if OpenMP is available, this will allow for simulating chromosomes in parallel.
#' If the value is NULL, the number of threads is automatically detected.
#' @param mutRate optional mutation rate per base pair per generation
#' to reduce or increase the number of segregating sites;
#' NULL uses the preset mutation rate in species demographic models
#' (GENERIC: 1e-8; CATTLE: 9.4e-9; WHEAT: 2e-9; MAIZE: 1.25e-8).
#' Cannot be combined with manualCommand.
#'
#' @details
#' There are currently four demographic models available:
#' GENERIC, CATTLE, WHEAT, and MAIZE.
#'
#' The GENERIC demographic model is meant to be a reasonable all-purpose choice
#' for breeding simulations. It runs quickly and models a population
#' with an effective population size that has gone through several historic bottlenecks.
#' This demographic model has the same parameters as the default arguments
#' in the \code{\link{runMacs2}} function, so examine that function for the details.
#'
#' The CATTLE demographic model is based on estimated demographic parameters
#' from real genome sequence data \insertCite{macleod2013Inferring}{AlphaSimR}.
#'
#' The WHEAT \insertCite{gaynor2017Two}{AlphaSimR} and
#' MAIZE \insertCite{hickey2014Evaluation}{AlphaSimR} demographic models
#' have been included due to their use in previous simulations.
#' Note that neither faithfully represent demographic histories of these species.
#' This is apparent by the lower number of segregating sites
#' simulated by each demographic models relative to their real-world analogs.
#' Adjusting these demographic models to better represent their real-world analogs
#' would result in a drastic increase to runtime with `MaCS`.
#'
#' For a fixed sample size, chromosome length, and demographic history,
#' expected site counts scale with \code{mutRate}. Reducing the rate reduces
#' sequence diversity and memory usage and may save simulation time,
#' but ancestry simulation limits the runtime saving.
#' The number of sites remains random, so allow a margin above \code{segSites}
#' when choosing a rate. Consider whether reduced diversity suits the
#' simulation's goals and report the mutation rate used.
#'
#' See the \code{vignette("DemographicModels", package = "AlphaSimR")} for
#' further information on \code{runMacs} and demographic models.
#'
#' Segregating sites are placed at absolute genetic map positions, scaled
#' from their positions on the simulated chromosome by its genetic length.
#' Each chromosome keeps that genetic length, which extends past its last
#' site, and its centromere is placed at half of it, whichever sites are
#' returned.
#'
#' @return an object of \code{\link{MapPop-class}}
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Creates a population of 10 outbred individuals
#' # Their genome consists of 1 chromosome and retains 100 segregating sites
#' \dontrun{
#' founderPop = runMacs(nInd=10,nChr=1,segSites=100)
#' }
#' @export
runMacs = function(nInd, nChr=1, segSites=NULL, inbred=FALSE, species="GENERIC",
                   split=NULL, ploidy=2L, manualCommand=NULL, manualGenLen=NULL,
                   nThreads=NULL, mutRate=NULL){
  nInd = as.integer(nInd)
  nChr = as.integer(nChr)
  ploidy = as.integer(ploidy)
  if(is.null(nThreads)){
    nThreads = getNumThreads()
  }
  if(nChr<nThreads){
    nThreads = nChr
  }

  if(is.null(segSites)){
    segSites = rep(0L,nChr)
  }else if(length(segSites)==1L){
    segSites = rep(as.integer(segSites),nChr)
  }else if(length(segSites)!=nChr){
    # One seed is drawn per chromosome below, and the C++ code reads a seed
    # for every element of segSites
    stop("segSites must have length 1 or nChr")
  }

  popSize = ifelse(inbred,nInd,ploidy*nInd)

  if(!is.null(mutRate)){
    if(!is.numeric(mutRate) || length(mutRate)!=1L ||
       !is.finite(mutRate) || mutRate<=0){
      stop("mutRate must be a single positive finite number")
    }
    if(!is.null(manualCommand)){
      stop("mutRate cannot be combined with manualCommand; set -t 4*Ne*mutRate in manualCommand instead")
    }
  }

  if(!is.null(manualCommand)){
    if(is.null(manualGenLen)) stop("You must define manualGenLen")
    command = paste0(popSize," ",manualCommand," -s ")
    genLen = manualGenLen
  }else{
    species = toupper(species)
    if(species=="GENERIC"){
      #GENERIC----
      phyLen = 1e8
      genLen = 1.0
      Ne = 100
      presetMutRate = 1e-8
      recRate = 1e-8
      histNe  = c(500, 1500, 6000, 12000, 100000)
      histGen = c(100, 1000, 10000, 100000, 1000000)
    }else if(species=="CATTLE"){
      #CATTLE----
      # Macleod et al. (2013) https://doi.org/10.1093/molbev/mst125
      phyLenGenome = 2.8e9 # https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/
      phyLen = round(phyLenGenome / 30)
      presetMutRate = 9.4e-9
      recRate = 9.26e-9
      genLen = recRate * phyLen
      Ne = 90
      histNe  = c(120, 250, 350, 1000, 1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000)
      histGen = c(  3,   6,  12,   18,   24,  154,  454,  654, 1754,  2354,  3354, 33154)
    }else if(species=="WHEAT"){
      #WHEAT----
      phyLen = 8e8
      genLen = 1.43
      Ne = 50
      presetMutRate = 2e-9
      recRate = 1.8e-9
      histNe  = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 12000, 16000, 20000, 24000, 28000, 32000)
      histGen = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 400, 600, 800, 1000, 2000, 4000, 6000, 8000, 10000, 20000, 40000, 60000, 80000, 100000)
    }else if(species=="MAIZE"){
      #MAIZE----
      phyLen = 2e8
      genLen = 2.0
      Ne = 100
      presetMutRate = 1.25e-08
      recRate = 1e-8
      histNe  = c(200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 4000, 6000, 8000, 10000)
      histGen = c(20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 800, 1200, 1600, 2000)
    }else{
      stop(paste("No rules for species",species))
    }

    effectiveMutRate = if(is.null(mutRate)) presetMutRate else mutRate
    scaledMutRate = 4 * Ne * effectiveMutRate
    scaledRecRate = 4 * Ne * recRate
    speciesParams = paste(phyLen, "-t", scaledMutRate, "-r", scaledRecRate)
    # GENERIC: 1E8 -t 4e-06 -r 4e-06
    # CATTLE: 93333333 -t 3.384e-06 -r 3.3336e-06
    # WHEAT: 8E8 -t 4E-7 -r 3.6E-7
    # MAIZE: 2E8 -t 5E-6 -r 4E-6
    # MaCS times are in coalescent units (4*Ne generations) and
    # all past Ne are relative to the present-time Ne
    speciesHist = NULL
    for(i in seq_along(histNe)){
      speciesHist = paste(speciesHist,"-eN", histGen[i] / (4*Ne), histNe[i] / Ne)
    }

    if(is.null(split)){
      splitI = ""
      splitJ = ""
    }else{
      stopifnot("popSize must be even to split it into two populations"=
                popSize%%2==0)
      splitI = paste(" -I 2",popSize%/%2,popSize%/%2)
      splitJ = paste(" -ej",split/(4*Ne)+0.000001,"2 1")
    }
    command = paste0(popSize," ",speciesParams,splitI," ",speciesHist,splitJ," -s ")
  }

  if(!is.null(manualGenLen)){
    genLen = manualGenLen
  }

  if(length(genLen)==1){
    genLen = rep(genLen, nChr)
  }

  # MaCS uses its own RNG, so runMacs() creates explicit chromosome seeds from
  # R's RNG and passes them through the R/C++ bridge. The same chromosome seed
  # is also reused for possible post-MaCS site subsampling in C++.
  seed = sample.int(n = 1e8, size = nChr)

  # Run MaCS
  macsOut = MaCS(command, segSites, inbred, ploidy,
                 nThreads, seed)
  dim(macsOut$geno) = NULL # Account for matrix bug in RcppArmadillo

  # Check if desired number of loci were obtained
  nLoci = sapply(macsOut$genMap,length)
  isLimited = segSites>0
  if(any(nLoci[isLimited] != segSites[isLimited])){
    stop("MaCS did not return enough segSites, use segSites=NULL to return all sites generated by MaCS")
  }

  # MaCS reports each site as a fraction of the chromosome's length, so
  # scaling by genLen gives an absolute map position
  genMap = vector("list",nChr)
  for(i in seq_len(nChr)){
    genMap[[i]] = genLen[i]*c(macsOut$genMap[[i]])
    names(genMap[[i]]) = paste(i,1:length(genMap[[i]]),sep="_")
  }
  names(genMap) = as.character(1:nChr)

  # The chromosome is as long as MaCS was asked to simulate, and the
  # centromere sits in its middle, so neither depends on which sites
  # happened to be returned
  output = new("MapPop",
               nInd=nInd,
               nChr=nChr,
               ploidy=ploidy,
               nLoci=nLoci,
               geno=macsOut$geno,
               genMap=genMap,
               centromere=as.numeric(genLen)/2,
               inbred=inbred,
               genLen=as.numeric(genLen))
  return(output)
}

#' @title Alternative wrapper for MaCS
#'
#' @description
#' A wrapper function for \code{\link{runMacs}}. This wrapper is designed
#' to provide a more intuitive interface for writing custom commands
#' in `MaCS` \insertCite{chen2009Fast}{AlphaSimR}. It effectively automates the creation
#' of an appropriate line for the manualCommand argument in \code{\link{runMacs}}
#' using user supplied variables, but only allows for a subset of the functionality
#' offered by this argument. The default demographic model matches
#' \code{species="GENERIC"} in \code{\link{runMacs}}.
#'
#' @param nInd number of individuals to simulate
#' @param nChr number of chromosomes to simulate (using independent simulations)
#' @param segSites number of segregating sites to keep per chromosome. A
#' value of \code{NULL} results in retaining all simulated segregating sites.
#' Requesting more sites than were generated raises an error.
#' @param Ne present effective population size
#' @param bp base pair length of chromosome
#' @param genLen genetic length of chromosome in Morgans
#' @param mutRate mutation rate per base pair per generation
#' @param histNe historical effective population sizes
#' Each size applies backward in time from the corresponding boundary in
#' \code{histGen} until the next boundary, or indefinitely for the oldest size.
#' @param histGen generation boundaries measured in generations ago,
#' in increasing order, corresponding to \code{histNe}. The present size
#' \code{Ne} applies from the present to the first boundary.
#' @param inbred should founder individuals be inbred by replicating their
#' whole-chromosome haplotypes (see also \code{ploidy}).
#' This creates fully homozygous founders without simulating a history of selfing.
#' @param split an optional historic population split in terms of generations ago
#' @param ploidy level of species to simulate \code{nInd * ploidy}
#' whole-chromosome haplotypes per chromosome unless \code{inbred=TRUE},
#' where \code{nInd} haplotypes are simulated and each is copied across
#' all homologues of an individual.
#' Also multiplies \code{Ne} and \code{histNe} by \code{ploidy / 2},
#' regardless of \code{inbred}, consistent with the autotetraploid
#' coalescent scaling of \insertCite{arnold2012Extending}{AlphaSimR}.
#' @param returnCommand should the command passed to manualCommand in
#' \code{\link{runMacs}} be returned. If TRUE, `MaCS` will not be called and
#' the command is returned instead.
#' @param nThreads if OpenMP is available, this will allow for simulating chromosomes in parallel.
#' If the value is NULL, the number of threads is automatically detected.
#'
#' @details
#' See \code{vignette("DemographicModels", package = "AlphaSimR")} for
#' further information on \code{runMacs2} and demographic models.
#'
#' @return an object of \code{\link{MapPop-class}} or if
#' returnCommand is true a string giving the MaCS command passed to
#' the manualCommand argument of \code{\link{runMacs}}.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Creates a population of 10 outbred individuals
#' # Their genome consists of 1 chromosome and retains 100 segregating sites
#' # The command is equivalent to species="GENERIC" in runMacs
#' \dontrun{
#' founderPop = runMacs2(nInd=10,nChr=1,segSites=100)
#'
#' # runMacs() Implementation of the CATTLE demographic model following
#' #  Macleod et al. (2013) https://doi.org/10.1093/molbev/mst125
#' cattleGenome = 2.8e9 # https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/
#' (cattleChr = round(cattleGenome / 30))
#' recRate = 9.26e-09
#' (cattleGenLen = recRate * cattleChr)
#' cattleMutRate = 9.4e-9
#' runMacs2(nInd = 10, nChr = 1, Ne = 90, bp = cattleChr,
#'          genLen = cattleGenLen, mutRate = cattleMutRate,
#'          histNe  = c(120, 250, 350, 1000, 1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000),
#'          histGen = c(  3,   6,  12,   18,   24,  154,  454,  654, 1754,  2354,  3354, 33154),
#'          returnCommand = TRUE)
#' }
#' @export
runMacs2 = function(nInd,nChr=1,segSites=NULL,Ne=100,
                    bp=1e8,genLen=1,mutRate=1e-8,
                    histNe=c(500,1500,6000,12000,100000),
                    histGen=c(100,1000,10000,100000,1000000),
                    inbred=FALSE,split=NULL,ploidy=2L,returnCommand=FALSE,
                    nThreads=NULL){
  stopifnot("histNe and histGen must have the same length"=
              length(histNe)==length(histGen))
  # Adjust Ne according to ploidy level
  Ne = Ne*(ploidy/2L)
  if(!is.null(histNe)){
    histNe = histNe*(ploidy/2L)
  }
  speciesParams = paste(bp,"-t",4*Ne*mutRate,
                        "-r",4*Ne*genLen/bp)
  speciesHist = ""
  if(length(histNe)>0){
    histNe = histNe/Ne
    histGen = histGen/(4*Ne)
    for(i in seq_len(length(histNe))){
      speciesHist = paste(speciesHist,"-eN",
                      histGen[i],histNe[i])
    }
  }
  if(is.null(split)){
    command = paste(speciesParams,speciesHist)
  }else{
    popSize = ifelse(inbred,nInd,ploidy*nInd)
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
                 ploidy=ploidy,manualCommand=command,
                 manualGenLen=genLen,nThreads=nThreads))
}

#' Predict moments of segregating-site counts
#'
#' Calculates tree-length and segregating-site moments
#' under piecewise constant demography by numerically tracking
#' the expected total branch length and its second moment
#' through successive population-size epochs.
#'
#' @param n Number of sampled haplotypes (an integer of at least two).
#' For outbred diploid individuals, use twice the number of individuals.
#' @param Ne Current diploid effective population size.
#' @param bp Length of the region in base pairs. Use 1 for a single base pair.
#' @param mutRate Mutation rate per base pair per generation.
#' @param histNe Historical diploid effective population sizes, ordered
#' backward in time, as in \code{\link{runMacs2}}.
#' @param histGen Increasing numbers of generations ago at which the sizes
#' in \code{histNe} take effect. Must have the same length as \code{histNe}.
#' The oldest size applies indefinitely. Set both historical vectors to
#' \code{numeric(0)} or \code{NULL} for a constant population size.
#' @param tolerance Positive numerical integration tolerance, smaller than 1.
#' This controls local numerical error in the scaled first and second
#' tree-length moments; it is not a guaranteed relative error bound
#' on the returned variance. The default, \code{1e-5}, favors speed.
#' Lower this value for higher accuracy.
#'
#' @details
#' Theory for segregating-site counts under piecewise constant demography
#' is described in \insertCite{zivkovic2008Second}{AlphaSimR}.
#' This function implements and extends the theory by
#' approximately solving differential equations of
#' how the first and second moments of total tree length
#' change over time under a piecewise constant population size model.
#' First and second tree-length moments are calculated numerically from
#' backward coalescent reward equations, including covariance between
#' lineage waiting times.
#' 
#' Total tree length \eqn{T_c} is the sum of branch lengths
#' over all sampled lineages until their most recent common ancestor.
#' With the moments of total tree length,
#' the expected number of segregating sites are calculated using
#' \eqn{E(S) = (\theta/2) E(T_c)} and
#' \eqn{Var(S) = E(S) + (\theta^2/4)Var(T_c)},
#' where \eqn{\theta = 4 N_0 \mu L}, with
#' \eqn{\mu} mutation rate per base pair per generation,
#' \eqn{N_0} present-time effective population size, and
#' \eqn{L} the length of the region in base pairs.
#' In these identities, \eqn{T_c} is in coalescent units
#' (time in units of \eqn{2N_0} generations).
#' 
#' These results follow from Poisson mutation counts conditional on tree length
#' under an infinite-sites model with a single genealogy of the genome region.
#' Recombination changes the variance by averaging correlated genealogies
#' along the region, but not the expected site count under the same demography.
#' The returned variance assumes no recombination and should not be used as
#' the chromosome variance for simulations with recombination.
#'
#' @return A named list with four numeric components:
#' \describe{
#'   \item{expectedTc}{Mean total tree length over the region in generations
#'     (divide by \eqn{2N_0} to obtain coalescent units)}
#'   \item{varianceTc}{Variance of total tree length over the region in generations squared
#'     (divide by \eqn{(2N_0)^2} to obtain squared coalescent units)}
#'   \item{expectedS}{Mean number of segregating sites in the region}
#'   \item{varianceS}{Variance of the number of segregating sites in the region}
#' }
#'
#' Similar approaches can be used to calculate moments
#' of other summary statistics under the same demographic model,
#' including pairwise nucleotide diversity and the allele frequency spectrum.
#' These extensions require additional development.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # GENERIC demographic model, with the same defaults as runMacs2()
#' predictSegSitesMoments(n=20, Ne=100, bp=1e8, mutRate=1e-8,
#'   histNe=c(500,1500,6000,12000,100000),
#'   histGen=c(100,1000,10000,100000,1000000))
#'
#' # Constant population size
#' predictSegSitesMoments(n=20, Ne=100, bp=1e8, mutRate=1e-8,
#'   histNe=numeric(0), histGen=numeric(0))
#' @export
predictSegSitesMoments = function(n,Ne,bp,mutRate,histNe,histGen,
                                  tolerance=1e-5){
  if(!is.numeric(n) || length(n)!=1L || !is.finite(n) ||
     n<2 || n!=floor(n) || n>=.Machine$integer.max){
    stop("n must be an integer of at least two")
  }
  if(!is.numeric(Ne) || length(Ne)!=1L || !is.finite(Ne) || Ne<=0){
    stop("Ne must be a finite positive number")
  }
  if(!is.numeric(bp) || length(bp)!=1L || !is.finite(bp) || bp<=0){
    stop("bp must be a finite positive number")
  }
  if(!is.numeric(mutRate) || length(mutRate)!=1L ||
     !is.finite(mutRate) || mutRate<0){
    stop("mutRate must be a finite non-negative number")
  }
  if(is.null(histNe)) histNe = numeric(0)
  if(is.null(histGen)) histGen = numeric(0)
  if(!is.numeric(histNe) || any(!is.finite(histNe)) || any(histNe<=0)){
    stop("histNe must contain finite positive population sizes")
  }
  if(!is.numeric(histGen) || length(histGen)!=length(histNe) ||
     any(!is.finite(histGen)) || any(histGen<=0) || any(diff(histGen)<=0)){
    stop("histGen must contain increasing positive epoch boundaries, one per histNe value")
  }
  if(!is.numeric(tolerance) || length(tolerance)!=1L ||
     !is.finite(tolerance) || tolerance<=0 || tolerance>=1){
    stop("tolerance must be a finite number between zero and one")
  }
  referenceNe = Ne
  theta = 4*Ne*mutRate*bp
  if(!is.finite(theta)) stop("The mutation parameter 4*Ne*mutRate*bp must be finite")
  relativeSizes = c(Ne,histNe)/referenceNe
  scaledTimes = histGen/(2*referenceNe)
  if(any(!is.finite(relativeSizes)) || any(relativeSizes<=0) ||
     any(!is.finite(scaledTimes)) || any(scaledTimes<=0) ||
     any(diff(scaledTimes)<=0)){
    stop("Population sizes and times cannot be represented on the reference scale")
  }
  # Integrate first and second reward moments backward through all epochs.
  # The oldest epoch has the exact constant-population terminal moments.
  moments = .segSitesTreeMoments(as.integer(n),relativeSizes,scaledTimes,tolerance)
  expectedS = (theta/2)*moments[1]
  list(expectedTc=moments[1]*(2*referenceNe),
       varianceTc=moments[2]*(2*referenceNe)^2,
       expectedS=expectedS,
       varianceS=expectedS+(theta/2)^2*moments[2])
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
#' @param ploidy new ploidy level for organism. If NULL,
#' the ploidy level of the mapPop is used.
#' @param replace should haplotypes be sampled with replacement
#'
#' @return an object of \code{\link{MapPop-class}}
#'
#' @examples
#' founderPop = quickHaplo(nInd=2,nChr=1,segSites=11,inbred=TRUE)
#' founderPop = sampleHaplo(mapPop=founderPop,nInd=20)
#'
#' @export
sampleHaplo = function(mapPop,nInd,inbred=FALSE,ploidy=NULL,replace=TRUE){
  if(is.null(ploidy)) ploidy = mapPop@ploidy
  nHaplo = mapPop@nInd*mapPop@ploidy
  if(inbred){
    nSamp = nInd
  }else{
    nSamp = nInd*ploidy
  }
  nBin = mapPop@nLoci%/%8L + (mapPop@nLoci%%8L > 0L)
  if(!replace) stopifnot("mapPop has too few haplotypes to sample without replacement; use replace=TRUE"=
                            nHaplo>=nSamp)
  output = vector("list",mapPop@nChr)
  for(chr in seq_len(mapPop@nChr)){
    haplo = sample.int(nHaplo,nSamp,replace=replace)
    geno = array(data=as.raw(0),
                 dim=c(nBin[chr],
                       ploidy,nInd))
    outHap = 1L
    outInd = 1L
    for(i in seq_len(length(haplo))){
      inHap = (haplo[i]-1L)%%mapPop@ploidy + 1L
      inInd = (haplo[i]-1L)%/%mapPop@ploidy + 1L
      if(inbred){
        for(outHap in seq_len(ploidy)){
          geno[,outHap,outInd] =
            mapPop@geno[[chr]][,inHap,inInd]
        }
        outInd = outInd+1L
      }else{
        geno[,outHap,outInd] =
          mapPop@geno[[chr]][,inHap,inInd]
        outHap = outHap%%ploidy+1L
        if(outHap==1L){
          outInd = outInd+1L
        }
      }
    }
    output[[chr]] = new("MapPop",
                        nInd=as.integer(nInd),
                        nChr=1L,
                        ploidy=as.integer(ploidy),
                        nLoci=mapPop@nLoci[chr],
                        geno=list(geno),
                        genMap=mapPop@genMap[chr],
                        centromere=mapPop@centromere[chr],
                        inbred=inbred,
                        genLen=.getGenLen(mapPop)[chr])
  }
  output = do.call("cChr",output)
  return(output)
}

#' @title Quick founder haplotype simulation
#'
#' @description Rapidly simulates founder haplotypes by randomly
#' sampling 0s and 1s. This is equivalent to having all loci with
#' allele frequency 0.5 and being in linkage equilibrium.
#'
#' @param nInd number of individuals to simulate
#' @param nChr number of chromosomes to simulate
#' @param segSites number of segregating sites per chromosome
#' @param genLen genetic length of chromosomes
#' @param ploidy ploidy level of organism
#' @param inbred should founder individuals be inbred
#'
#' @details
#' See \code{\link{runMacs}} and \code{\link{runMacs2}} for
#' more realistic founder haplotype simulation.
#'
#' @return an object of \code{\link{MapPop-class}}
#'
#' @examples
#' # Creates a populations of 10 outbred individuals
#' # Their genome consists of 1 chromosome and 100 segregating sites
#' founderPop = quickHaplo(nInd=10,nChr=1,segSites=100)
#'
#' @export
quickHaplo = function(nInd,nChr,segSites,genLen=1,ploidy=2L,inbred=FALSE){
  ploidy = as.integer(ploidy)
  nInd = as.integer(nInd)
  nChr = as.integer(nChr)
  segSites = as.integer(segSites)
  if(length(segSites)==1) segSites = rep(segSites,nChr)
  if(length(genLen)==1) genLen = rep(genLen,nChr)
  nBins = segSites%/%8L + (segSites%%8L > 0L)
  centromere = genLen/2

  genMap = vector("list",nChr)
  geno = vector("list",nChr)
  for(i in seq_len(nChr)){
    genMap[[i]] = seq(0,genLen[i],length.out=segSites[i])
    names(genMap[[i]]) = paste(i, 1:segSites[i], sep="_")
    geno[[i]] = array(sample(as.raw(0:255),
                             nInd*ploidy*nBins[i],
                             replace=TRUE),
                      dim = c(nBins[i],ploidy,nInd))

    if(inbred){
      if(ploidy>1){
        for(j in 2:ploidy){
          geno[[i]][,j,] = geno[[i]][,1,]
        }
      }
    }
  }
  names(genMap) = as.character(1:nChr)
  return(new("MapPop",
             nInd=nInd,
             nChr=nChr,
             ploidy=ploidy,
             nLoci=segSites,
             geno=geno,
             genMap=genMap,
             centromere=centromere,
             inbred=inbred,
             genLen=as.numeric(genLen)))
}

#' @title Add segregating site to MapPop
#'
#' @description This function allows for adding a new
#' segregating site with user supplied genotypes to a MapPop.
#' The position of the site is set using a genetic map position.
#'
#' @param mapPop an object of \code{\link{MapPop-class}}
#' @param siteName name to give the segregating site
#' @param chr which chromosome to add the site
#' @param mapPos absolute genetic map position of site in Morgans
#' @param haplo haplotypes for the site
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{getNumThreads()}.
#' @param extendGenLen should the chromosome's genetic length be
#' extended to \code{mapPos} when the site lies past it. If FALSE, a
#' site past the end of the chromosome is an error. The centromere is
#' not moved.
#'
#' @return an object of \code{\link{MapPop-class}}
#'
#' @examples
#' # Creates a populations of 10 outbred individuals
#' # Their genome consists of 1 chromosome and 2 segregating sites
#' founderPop = quickHaplo(nInd=10,nChr=1,segSites=2)
#'
#' # Add a locus a the 0.5 Morgan map position
#' haplo = matrix(sample(x=0:1, size=20, replace=TRUE), ncol=1)
#'
#' founderPop2 = addSegSite(founderPop, siteName="x", chr=1, mapPos=0.5,
#'                          haplo=haplo, nThreads=1L)
#'
#' pullSegSiteHaplo(founderPop2)
#'
#' @export
addSegSite = function(mapPop, siteName, chr, mapPos, haplo, nThreads=NULL,
                      extendGenLen=FALSE){
  # Check validity of input data
  stopifnot("mapPop must be a MapPop"=is(mapPop, "MapPop"))
  stopifnot("haplo must have one value per haplotype, nInd*ploidy"=
              length(haplo)==(mapPop@nInd*mapPop@ploidy))
  if(is.null(nThreads)){
    nThreads = getNumThreads()
  }else{
    nThreads = as.integer(nThreads)
  }

  # Map positions are absolute, so the new site cannot sit before position 0
  if(length(mapPos)!=1L || !is.numeric(mapPos) || !is.finite(mapPos) ||
     mapPos<0){
    stop("mapPos must be a single non-negative number")
  }

  # Check that name isn't already present
  allSiteNames = unlist(unname(lapply(mapPop@genMap, names)))
  if(siteName%in%allSiteNames){
    stop("The siteName '",siteName, "' is already in use.")
  }

  # Coerce haplo to a raw matrix
  haplo = matrix(as.raw(haplo), ncol=1)

  # Convert chr to a number, if needed
  if(is.character(chr)){
    chr = match(chr, names(mapPop@genMap))
  }

  # A site past the end of the chromosome changes how long the chromosome
  # is, so it is only allowed when asked for
  genLen = .getGenLen(mapPop)
  if(mapPos>genLen[chr]){
    if(isTRUE(extendGenLen)){
      genLen[chr] = mapPos
    }else{
      stop("mapPos is past the end of the chromosome (genLen=", genLen[chr],
           "), use extendGenLen=TRUE to extend it")
    }
  }
  mapPop@genLen = genLen

  # Extract temporary map and geno for target chromosome
  chrMap = mapPop@genMap[[chr]]
  markerNames = names(chrMap)
  chrGeno = mapPop@geno[chr] # Leaving in list format for getHaplo

  # Extract haplotype matrix from the bit array
  haploMat = getHaplo(chrGeno,
                      mapPop@nLoci[chr],
                      1:mapPop@nLoci[chr],
                      nThreads)

  # Find position of insertion
  pos = findInterval(x=mapPos, vec=chrMap)

  # Insert site
  if(pos==length(chrMap)){
    # Inserting site at the end
    chrMap = c(chrMap, mapPos)
    names(chrMap) = c(markerNames, siteName)
    haploMat = cbind(haploMat, haplo)
  }else if(pos==0){
    # Inserting site at the beginning
    chrMap = c(mapPos, chrMap)
    names(chrMap) = c(siteName, markerNames)
    haploMat = cbind(haplo, haploMat)
  }else{
    # Inserting site in the middle
    chrMap = c(chrMap[1:pos], mapPos,
               chrMap[(pos+1L):length(chrMap)])
    names(chrMap) = c(markerNames[1:pos], siteName,
                      markerNames[(pos+1L):length(markerNames)])
    haploMat = cbind(haploMat[,1:pos,drop=FALSE], haplo,
                     haploMat[,(pos+1L):ncol(haploMat),drop=FALSE])
  }

  # Convert haplotype matrix back to a bit array
  chrGeno = packHaplo(haploMat, ploidy=mapPop@ploidy, inbred=FALSE)

  # Return temporary map and geno
  mapPop@genMap[[chr]] = chrMap
  mapPop@geno[[chr]] = chrGeno

  # Add one to locus count
  mapPop@nLoci[chr] = mapPop@nLoci[chr] + 1L

  # Check if an inbred population has become outbred
  if(mapPop@inbred){
    M = pullMarkerGeno(mapPop, markers=siteName, nThreads=nThreads)
    isHomo0 = M==0L
    isHomo1 = M==mapPop@ploidy
    isHet = !(isHomo0 | isHomo1)
    if(any(isHet)){
      mapPop@inbred = FALSE
    }
  }

  return(mapPop)
}
