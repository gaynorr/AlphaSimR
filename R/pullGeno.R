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
#' @param gender determines which gender specific map 
#' is returned. Options are "A" for average map, "F" 
#' for female map, and "M" for male map. All options are 
#' equivalent if not using gender specific maps.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a data.frame for the SNP map.
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
getSnpMap = function(snpChip=1, gender="A", simParam=NULL){
  
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  #Extract genetic map and SNP positions
  if(gender=="A"){
    genMap = simParam$genMap
  }else if(gender=="F"){
    genMap = simParam$femaleMap
  }else if(gender=="M"){
    genMap = simParam$maleMap
  }else{
    stop(paste("gender =",gender,"is not a valid option"))
  }
  snp = simParam$snpChips[[snpChip]] #SNP positions
  
  #Create a list of SNP postions on the genetic map
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
                      pos=do.call("c",snpMap))
  return(output)
}

#' @title Get QTL genetic map
#' 
#' @description Retrieves the genetic map for the 
#' QTL of a given trait.
#' 
#' @param trait an integer for the 
#' @param gender determines which gender specific map 
#' is returned. Options are "A" for average map, "F" 
#' for female map, and "M" for male map. All options are 
#' equivalent if not using gender specific maps.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a data.frame for the SNP map.
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
getQtlMap = function(trait=1, gender="A", simParam=NULL){
  
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  #Extract genetic map and SNP positions
  if(gender=="A"){
    genMap = simParam$genMap
  }else if(gender=="F"){
    genMap = simParam$femaleMap
  }else if(gender=="M"){
    genMap = simParam$maleMap
  }else{
    stop(paste("gender =",gender,"is not a valid option"))
  }
  qtl = simParam$traits[[trait]] #QTL positions
  
  #Create a list of QTL postions on the genetic map
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
  
  #Create a data.frame with QTL postions on genetic map
  output = data.frame(id=1:qtl@nLoci,
                      chr=rep(1:simParam$nChr,qtl@lociPerChr),
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
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
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

#' @title Pull SNP genotype for multiple snp chips
#'
#' @description Retrieves SNP genotype data for multiple snp chips
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param chips a vector. For each animal indicates what snp
#' chip to use
#' @param missing What value to use for missing
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @return Returns a matrix of SNP genotypes.
#' 
#' @export
pullMultipleSnpGeno = function(pop, chips,
                               missing=9, simParam=NULL) {
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  stopifnot(length(chips) == pop@nInd)
  # I feel like the next line shouldn't be needed but I don't know
  # enough R! (dmoney)
  missing = as.integer(missing)
  allSnps = numeric(0)
  uniqueChips = unique(chips)
  for (c in uniqueChips){
    allSnps = sort(union(allSnps,simParam$snpChips[[c]]@lociLoc))
  }

  output = matrix(pop@nInd,length(allSnps),data=missing)
  if(class(pop)=="Pop"){
    rownames(output) = pop@id
  }else{
    rownames(output) = as.character(1:pop@nInd)
  }

  for (snpChip in uniqueChips){
    mask = allSnps %in% simParam$snpChips[[snpChip]]@lociLoc
    one = getGeno(pop@geno,
                  simParam$snpChips[[snpChip]]@lociPerChr,
                  simParam$snpChips[[snpChip]]@lociLoc,
                  simParam$nThreads)
    one = convToImat(one)
    for (i in 1:pop@nInd){
      if (chips[i] == snpChip) {
        output[i,mask] = one[i,]
        output[i,mask] = one[i,]
      }
    }
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
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
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
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
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
  if(class(pop)=="MapPop"){
    allLoci = unlist(c(sapply(pop@nLoci, function(x) 1:x)))
    lociTot = pop@nLoci
    nThreads = getNumThreads()
  }else{
    if(is.null(simParam)){
      simParam = get("SP",envir=.GlobalEnv)
    }
    allLoci = unlist(c(sapply(simParam$segSites, function(x) 1:x)))
    lociTot = simParam$segSites
    nThreads = simParam$nThreads
  }
  tmp = selectLoci(chr,lociTot,allLoci)
  output = getGeno(pop@geno,tmp$lociPerChr,tmp$lociLoc,nThreads)
  output = convToImat(output)
  if(class(pop)=="Pop"){
    rownames(output) = pop@id
  }else{
    rownames(output) = as.character(1:pop@nInd)
  }
  colnames(output) = paste("SITE",1:ncol(output),sep="_")
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
#' haplotyes and a value of 2 for male haplotypes.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of SNP haplotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
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

#' @title Pull SNP haplotypes for multiple chips
#'
#' @description Retrieves SNP haplotype data for multiple snp
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param chips a vector. For each animal indicates what snp
#' chip to use
#' @param haplo either "all" for all haplotypes or an integer
#' for a single set of haplotypes. Use a value of 1 for female
#' haplotyes and a value of 2 for male haplotypes.
#' @param missing What value to use for missing
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of SNP haplotypes.
#' 
#' @export
pullMultipleSnpHaplo = function(pop, chips, haplo="all",
                                missing=9, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  stopifnot(length(chips) == pop@nInd)
  # I feel like the next line shouldn't be needed but I don't know
  # enough R! (dmoney)
  missing = as.integer(missing)
  allSnps = numeric(0)
  uniqueChips = unique(chips)
  for (c in uniqueChips){
    allSnps = sort(union(allSnps,simParam$snpChips[[c]]@lociLoc))
  }

  if (haplo == "all") {
    output = matrix(pop@nInd*2,length(allSnps),data=missing)
    if(class(pop)=="Pop"){
      rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(rep(1:pop@nInd,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }
  }else{
    output = matrix(pop@nInd,length(allSnps),data=missing)
    if(class(pop)=="Pop"){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  for (snpChip in uniqueChips){
    mask = allSnps %in% simParam$snpChips[[snpChip]]@lociLoc
    if (haplo == "all") {
      one = getHaplo(pop@geno,
                     simParam$snpChips[[snpChip]]@lociPerChr,
                     simParam$snpChips[[snpChip]]@lociLoc,
                     simParam$nThreads)
      one = convToImat(one)
      for (i in 1:pop@nInd){
        if (chips[i] == snpChip) {
          output[i*2-1,mask] = one[i*2-1,]
          output[i*2,mask] = one[i*2,]
        }
      }
    }
    else {
      one = getOneHaplo(pop@geno,
                        simParam$snpChips[[snpChip]]@lociPerChr,
                        simParam$snpChips[[snpChip]]@lociLoc,
                        as.integer(haplo),
                        simParam$nThreads)
      one = convToImat(one)
      for (i in 1:pop@nInd){
        if (chips[i] == snpChip) {
          output[i,mask] = one[i,]
          output[i,mask] = one[i,]
        }
      }
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
#' haplotyes and a value of 2 for male haplotypes.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of QTL haplotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
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
#' haplotyes and a value of 2 for male haplotypes.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosome are retrieved.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of haplotypes
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
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
  if(class(pop)=="MapPop"){
    allLoci = unlist(c(sapply(pop@nLoci, function(x) 1:x)))
    lociTot = pop@nLoci
    nThreads = getNumThreads()
  }else{
    if(is.null(simParam)){
      simParam = get("SP",envir=.GlobalEnv)
    }
    allLoci = unlist(c(sapply(simParam$segSites, function(x) 1:x)))
    lociTot = simParam$segSites
    nThreads = simParam$nThreads
  }
  if(!is.null(chr)){
    tmp = selectLoci(chr,lociTot,allLoci)
    lociTot = tmp$lociPerChr
    allLoci = tmp$lociLoc
  }
  if(haplo=="all"){
    output = getHaplo(pop@geno,
                      lociTot,
                      allLoci,
                      nThreads)
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
    output = getOneHaplo(pop@geno,
                         lociTot,
                         allLoci,
                         as.integer(haplo),
                         nThreads)
    output = convToImat(output)
    if(class(pop)=="Pop"){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  colnames(output) = paste("SITE",1:ncol(output),sep="_")
  return(output)
}

#' @title Pull Identity By Descent (IBD) haplotypes
#' 
#' @description 
#' Retrieves Identity By Descent (IBD) haplotype data
#'
#' @param pop an object of \code{\link{Pop-class}} or 
#' \code{\link{RawPop-class}}. If NULL, haplotypes for the whole
#' ancestral pedigree are retreived. Otherwise, haplotypes just for
#' the \code{pop} individuals are retreived. In both cases the base
#' population is controlled by \code{pedigree}.
#' @param chr a vector of chromosomes to retrieve. If NULL, 
#' all chromosomes are retrieved.
#' @param snpChip an integer. Indicates which SNP array loci
#' are retrieved. If NULL, all sites are retrieved.
#' @param pedigree a matrix with ancestral pedigree to set a base
#' population. It should be of the same form as \code{simParam$pedigree} 
#' (see \code{\link{SimParam_setTrackPed}}), i.e., two columns (mother
#' and father) and the same number of rows as \code{simParam$pedigree}.
#' Base population can be set by setting parents as 0. If NULL, pedigree
#' from \code{\link{SimParam}} is taken.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns a matrix of haplotypes with Identity By Descent
#' (IBD) coding of locus alleles. The matrix colnames reflect whether
#' all segregagting loci (sites) are retreived or only SNP array loci.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
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
pullIbdHaplo = function(pop = NULL, chr = NULL, snpChip = NULL, pedigree = NULL, simParam = NULL) {
  
  # ---- Setup -----
  
  if (is.null(simParam)) {
    simParam = get(x = "SP", envir = .GlobalEnv)
  }
  if (pop@ploidy != 2L) {
    stop("pullIbdHaplo() works (currently) only with diploids!")
  }
  if (!simParam$isTrackRec) {
    stop("To use pullIbdHaplo(), simParam must hold ancestral recombination data! See ?SimParam_setTrackRec")
  }
  if (is.null(pedigree) & !simParam$isTrackPed) {
    stop("To use pullIbdHaplo() with pedigree = NULL, simParam must hold ancestral pedigree data! See ?SimParam_setTrackPed")
  }
  lociLoc = c(sapply(X = simParam$segSites, FUN = function(x) 1L:x))
  lociPerChr = simParam$segSites
  if (is.null(chr)) {
    chr = 1L:simParam$nChr
  } else {
    tmp = selectLoci(chr = chr, inLociPerChr = lociPerChr, inLociLoc = lociLoc)
    lociLoc = tmp$lociLoc
    lociPerChr = tmp$lociPerChr
  }
  if (is.null(pedigree)) {
    pedigree = simParam$pedigree
  } else {
    if (nrow(pedigree) != length(simParam$recHist)) {
      stop("pedigree must have the same number of rows as simParam$recHist!")
    }
  }
  
  # ---- Get IBD recombinations -----
  
  ibdRecHist = getIbdRecHist(recHist     = simParam$recHist,
                             pedigree    = pedigree,
                             nLociPerChr = lociPerChr)$ibdRecHist

  # ---- Get IBD haplotypes -----
  
  if (!is.null(pop)) {
    nInd = pop@nInd
    individuals = as.integer(pop@id)
  } else {
    nInd = nrow(pedigree)
    individuals = 1L:nInd
  }
  output = getIbdHaplo(ibdRecHist  = ibdRecHist,
                       individuals = individuals,
                       nLociPerChr = lociPerChr)
  rownames(output) = paste(rep(x = individuals,        each  = pop@ploidy),
                           rep(x = 1L:pop@ploidy, times = nInd), sep = "_")
  
  # ---- Subset loci -----
  
  if (!is.null(snpChip)) {
    Sel = integer(length = sum(simParam$snpChips[[snpChip]]@lociPerChr[chr]))
    ArrayEnd = 0L
    ChrStart = 0L
    for (Chr in chr) {
      # Chr = 1L
      # Chr = 2L
      Tmp = selectLoci(chr = Chr,
                       inLociPerChr = simParam$snpChips[[snpChip]]@lociPerChr,
                       inLociLoc    = simParam$snpChips[[snpChip]]@lociLoc)
      ArrayStart = ArrayEnd + 1L
      ArrayEnd   = ArrayStart + Tmp$lociPerChr[Chr] - 1L
      # cat(ArrayStart, ArrayEnd, "\n")
      Sel[ArrayStart:ArrayEnd] = ChrStart + Tmp$lociLoc
      ChrStart = ChrStart + lociPerChr[Chr]
    }
    output = output[, Sel]
    colnames(output) = paste("SNP",  1L:ncol(output), sep = "_")
  } else {
    colnames(output) = paste("SITE", 1L:ncol(output), sep = "_")
  }
  
  return(output)
}