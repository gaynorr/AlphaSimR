convToImat = function(X){
  return(matrix(as.integer(X),nrow=nrow(X),ncol=ncol(X)))
}

selectLoci = function(chr,inLociPerChr,inLociLoc){
  if(is.null(chr)){
    return(list(lociPerChr=inLociPerChr,
                lociLoc=inLociLoc))
  }
  nChr = length(inLociPerChr)
  stopifnot(any(chr%in%(1:nChr)),
            max(chr)<=nChr)
  outLociPerChr = numeric(nChr)
  outLociPerChr[chr] = inLociPerChr[chr]
  outLociLoc = numeric(sum(outLociPerChr))
  inStart = outStart = inEnd = outEnd = 0L
  for(i in 1:nChr){
    inStart = inStart + 1L
    inEnd = inEnd + inLociPerChr[i]
    if(outLociPerChr[i]>0){
      outStart = outStart + 1L
      outEnd = outEnd + inLociPerChr[i]
      outLociLoc[outStart:outEnd] = inLociLoc[inStart:inEnd]
      outStart = outEnd
    }
    inStart = inEnd
  }
  return(list(lociPerChr=outLociPerChr,
              lociLoc=outLociLoc))
}

#' @title Get SNP genetic map
#' 
#' @description Retrieves the genetic map for a 
#' given SNP chip.
#' 
#' @param snpChip an integer. Indicates which SNP
#' chip's map to retrieve.
#' @param sex determines which sex specific map 
#' is returned. Options are "A" for average map, "F" 
#' for female map, and "M" for male map. All options are 
#' equivalent if not using sex specific maps.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a data.frame with:
#' \describe{
#'   \item{id}{Unique identifier for the SNP}
#'   \item{chr}{Chromosome containing the SNP}
#'   \item{site}{Segregating site on the chromosome}
#'   \item{pos}{Genetic map position}
#' }
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addSnpChip(5)
#' 
#' #Pull SNP map
#' getSnpMap(snpChip=1, simParam=SP)
#' 
#' @export
getSnpMap = function(snpChip=1, sex="A", simParam=NULL){
  
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  #Extract genetic map and SNP positions
  if(sex=="A"){
    genMap = simParam$genMap
  }else if(sex=="F"){
    genMap = simParam$femaleMap
  }else if(sex=="M"){
    genMap = simParam$maleMap
  }else{
    stop(paste("sex =",sex,"is not a valid option"))
  }
  snp = simParam$snpChips[[snpChip]] #SNP positions
  
  #Create a list of SNP positions on the genetic map
  #Each list element corresponds to a chromosome
  snpMap = lapply(1:simParam$nChr, function(x){
    if(snp@lociPerChr[x]==0){
      #No SNPs on chromosome
      return(NULL)
    }else{
      if(x==1){
        #First chromosome, start at position 1
        take = 1:snp@lociPerChr[x]
      }else{
        #All other chromosomes
        take = (sum(snp@lociPerChr[1:(x-1)])+1):sum(snp@lociPerChr[1:x])
      }
      return(genMap[[x]][snp@lociLoc[take]])
    }
  })
  
  #Create a data.frame with SNP postions on genetic map
  output = data.frame(id=1:snp@nLoci,
                      chr=rep(1:simParam$nChr,snp@lociPerChr),
                      site=snp@lociLoc,
                      pos=do.call("c",snpMap))
  return(output)
}

#' @title Get QTL genetic map
#' 
#' @description Retrieves the genetic map for the 
#' QTL of a given trait.
#' 
#' @param trait an integer for the 
#' @param sex determines which sex specific map 
#' is returned. Options are "A" for average map, "F" 
#' for female map, and "M" for male map. All options are 
#' equivalent if not using sex specific maps.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a data.frame with:
#' \describe{
#'   \item{id}{Unique identifier for the QTL}
#'   \item{chr}{Chromosome containing the QTL}
#'   \item{site}{Segregating site on the chromosome}
#'   \item{pos}{Genetic map position}
#' }
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(5)
#' 
#' #Pull SNP map
#' getQtlMap(trait=1, simParam=SP)
#' 
#' @export
getQtlMap = function(trait=1, sex="A", simParam=NULL){
  
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  #Extract genetic map and SNP positions
  if(sex=="A"){
    genMap = simParam$genMap
  }else if(sex=="F"){
    genMap = simParam$femaleMap
  }else if(sex=="M"){
    genMap = simParam$maleMap
  }else{
    stop(paste("sex =",sex,"is not a valid option"))
  }
  qtl = simParam$traits[[trait]] #QTL positions
  
  #Create a list of QTL positions on the genetic map
  #Each list element corresponds to a chromosome
  qtlMap = lapply(1:simParam$nChr, function(x){
    if(qtl@lociPerChr[x]==0){
      #No QTL on chromosome
      return(NULL)
    }else{
      if(x==1){
        #First chromosome, start at position 1
        take = 1:qtl@lociPerChr[x]
      }else{
        #All other chromosomes
        take = (sum(qtl@lociPerChr[1:(x-1)])+1):sum(qtl@lociPerChr[1:x])
      }
      return(genMap[[x]][qtl@lociLoc[take]])
    }
  })
  
  #Create a data.frame with QTL positions on genetic map
  output = data.frame(id=1:qtl@nLoci,
                      chr=rep(1:simParam$nChr,qtl@lociPerChr),
                      site=qtl@lociLoc,
                      pos=do.call("c",qtlMap))
  return(output)
}

#' @title Pull SNP genotype
#'
#' @description Retrieves SNP genotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param snpChip an integer. Indicates which SNP
#' chip's genotypes to retrieve.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of SNP genotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullSnpGeno(pop, simParam=SP)
#' 
#' @export
pullSnpGeno = function(pop, snpChip=1, chr=NULL, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  tmp = selectLoci(chr,
                   simParam$snpChips[[snpChip]]@lociPerChr,
                   simParam$snpChips[[snpChip]]@lociLoc)
  output = getGeno(pop@geno,tmp$lociPerChr,tmp$lociLoc,simParam$nThreads)
  output = convToImat(output)
  if(class(pop)=="Pop"){
    rownames(output) = pop@id
  }else{
    rownames(output) = as.character(1:pop@nInd)
  }
  colnames(output) = paste("SNP",1:ncol(output),sep="_")
  return(output)
}

#' @title Pull QTL genotype
#'
#' @description Retrieves QTL genotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param trait an integer. Indicates which trait's
#' QTL genotypes to retrieve.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of QTL genotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullQtlGeno(pop, simParam=SP)
#' 
#' @export
pullQtlGeno = function(pop, trait=1, chr=NULL, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  tmp = selectLoci(chr,
                   simParam$traits[[trait]]@lociPerChr,
                   simParam$traits[[trait]]@lociLoc)
  output = getGeno(pop@geno,tmp$lociPerChr,tmp$lociLoc,simParam$nThreads)
  output = convToImat(output)
  if(class(pop)=="Pop"){
    rownames(output) = pop@id
  }else{
    rownames(output) = as.character(1:pop@nInd)
  }
  colnames(output) = paste("QTL",1:ncol(output),sep="_")
  return(output)
}

#' @title Pull seg site genotypes
#'
#' @description
#' Retrieves genotype data for all segregating sites
#'
#' @param pop an object of \code{\link{Pop-class}} or
#' \code{\link{RawPop-class}}
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of genotypes
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullSegSiteGeno(pop, simParam=SP)
#' 
#' @export
pullSegSiteGeno = function(pop, chr=NULL, simParam=NULL){
  if(class(pop)=="MapPop" | class(pop)=="NamedMapPop"){
    allLoci = unlist(c(sapply(pop@nLoci, function(x) 1:x)))
    lociTot = pop@nLoci
    nThreads = getNumThreads()
    map = pop@genMap
  }else{
    if(is.null(simParam)){
      simParam = get("SP",envir=.GlobalEnv)
    }
    allLoci = unlist(c(sapply(simParam$segSites, function(x) 1:x)))
    lociTot = simParam$segSites
    nThreads = simParam$nThreads
    map = simParam$genMap
  }
  tmp = selectLoci(chr,lociTot,allLoci)
  output = getGeno(pop@geno,tmp$lociPerChr,tmp$lociLoc,nThreads)
  output = convToImat(output)
  if(class(pop)=="Pop" | class(pop)=="NamedMapPop"){
    rownames(output) = pop@id
  }else{
    rownames(output) = as.character(1:pop@nInd)
  }
  colnames(output) = unlist(lapply(map[chr], names))
  return(output)
}

#' @title Pull SNP haplotypes
#'
#' @description Retrieves SNP haplotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param snpChip an integer. Indicates which SNP
#' chip's haplotypes to retrieve.
#' @param haplo either "all" for all haplotypes or an integer
#' for a single set of haplotypes. Use a value of 1 for female
#' haplotypes and a value of 2 for male haplotypes.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of SNP haplotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullSnpHaplo(pop, simParam=SP)
#' 
#' @export
pullSnpHaplo = function(pop, snpChip=1, haplo="all",
                        chr=NULL, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  tmp = selectLoci(chr,
                   simParam$snpChips[[snpChip]]@lociPerChr,
                   simParam$snpChips[[snpChip]]@lociLoc)
  lociPerChr = tmp$lociPerChr
  lociLoc = tmp$lociLoc
  if(haplo=="all"){
    output = getHaplo(pop@geno,lociPerChr,lociLoc,simParam$nThreads)
    output = convToImat(output)
    if(class(pop)=="Pop"){
      rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(rep(1:pop@nInd,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }
  }else{
    stopifnot(haplo%in%c(1,2))
    output = getOneHaplo(pop@geno,lociPerChr,lociLoc,
                         as.integer(haplo),simParam$nThreads)
    output = convToImat(output)
    if(class(pop)=="Pop"){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  colnames(output) = paste("SNP",1:ncol(output),sep="_")
  return(output)
}

#' @title Pull QTL haplotypes
#'
#' @description Retrieves QTL haplotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param trait an integer. Indicates which trait's
#' QTL haplotypes to retrieve.
#' @param haplo either "all" for all haplotypes or an integer
#' for a single set of haplotypes. Use a value of 1 for female
#' haplotypes and a value of 2 for male haplotypes.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of QTL haplotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullQtlHaplo(pop, simParam=SP)
#' 
#' @export
pullQtlHaplo = function(pop, trait=1, haplo="all",
                        chr=NULL, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  tmp = selectLoci(chr,
                   simParam$traits[[trait]]@lociPerChr,
                   simParam$traits[[trait]]@lociLoc)
  lociPerChr = tmp$lociPerChr
  lociLoc = tmp$lociLoc
  if(haplo=="all"){
    output = getHaplo(pop@geno,lociPerChr,lociLoc,simParam$nThreads)
    output = convToImat(output)
    if(class(pop)=="Pop"){
      rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(rep(1:pop@nInd,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }
  }else{
    stopifnot(haplo%in%c(1,2))
    output = getOneHaplo(pop@geno,lociPerChr,lociLoc,
                         as.integer(haplo),simParam$nThreads)
    output = convToImat(output)
    if(class(pop)=="Pop"){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  colnames(output) = paste("QTL",1:ncol(output),sep="_")
  return(output)
}

#' @title Pull seg site haplotypes
#'
#' @description
#' Retrieves haplotype data for all segregating sites
#'
#' @param pop an object of \code{\link{Pop-class}} or
#' \code{\link{RawPop-class}}
#' @param haplo either "all" for all haplotypes or an integer
#' for a single set of haplotypes. Use a value of 1 for female
#' haplotypes and a value of 2 for male haplotypes.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of haplotypes
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullSegSiteHaplo(pop, simParam=SP)
#' 
#' @export
pullSegSiteHaplo = function(pop, haplo="all",
                            chr=NULL, simParam=NULL){
  if(class(pop)=="MapPop" | class(pop)=="NamedMapPop"){
    allLoci = unlist(c(sapply(pop@nLoci, function(x) 1:x)))
    lociTot = pop@nLoci
    nThreads = getNumThreads()
    map = pop@genMap
  }else{
    if(is.null(simParam)){
      simParam = get("SP",envir=.GlobalEnv)
    }
    allLoci = unlist(c(sapply(simParam$segSites, function(x) 1:x)))
    lociTot = simParam$segSites
    nThreads = simParam$nThreads
    map = simParam$genMap
  }
  if(!is.null(chr)){
    tmp = selectLoci(chr,lociTot,allLoci)
    lociTot = tmp$lociPerChr
    allLoci = tmp$lociLoc
  }else{
    chr = 1:pop@nChr
  }
  if(haplo=="all"){
    output = getHaplo(pop@geno,
                      lociTot,
                      allLoci,
                      nThreads)
    output = convToImat(output)
    if(class(pop)=="Pop" | class(pop)=="NamedMapPop"){
      rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(rep(1:pop@nInd,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }
  }else{
    stopifnot(haplo%in%c(1,2))
    output = getOneHaplo(pop@geno,
                         lociTot,
                         allLoci,
                         as.integer(haplo),
                         nThreads)
    output = convToImat(output)
    if(class(pop)=="Pop" | class(pop)=="NamedMapPop"){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  colnames(output) = unlist(lapply(map[chr], names))
  return(output)
}


#' @title Pull IBD haplotypes
#'
#' @description Retrieves IBD haplotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of SNP haplotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' SP$setTrackRec(TRUE)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullIbdHaplo(pop, simParam=SP)
#' 
#' @export
pullIbdHaplo = function(pop, chr=NULL, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(!simParam$isTrackRec){
    stop(
      "pullIbdHaplo can only be used with the trackRec option, see trackRec in ?SimParam"
      )
  }
  
  if(is.null(chr)){
    chr = 1:pop@nChr
  }
  
  # Retrieve IBD data
  ibd = simParam$ibdHaplo(pop@iid)
  
  # Fill in output matrix
  output = createIbdMat(ibd=ibd, chr=chr, 
                        nLoci=pop@nLoci, 
                        ploidy=pop@ploidy,
                        nThreads=simParam$nThreads)
  
  rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                           rep(1:pop@ploidy,pop@nInd),sep="_")
  
  colnames(output) = unlist(lapply(simParam$genMap[chr], names))
  
  return(output)
}
