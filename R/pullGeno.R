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
#' @export
pullSnpGeno = function(pop, snpChip=1, chr=NULL, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  tmp = selectLoci(chr,
                   simParam$snpChips[[snpChip]]@lociPerChr,
                   simParam$snpChips[[snpChip]]@lociLoc)
  output = getGeno(pop@geno,tmp$lociPerChr,tmp$lociLoc)
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
#' @return Returns a matrix of SNP genotypes.
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
                  simParam$snpChips[[snpChip]]@lociLoc)
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
#' @export
pullQtlGeno = function(pop, trait=1, chr=NULL, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  tmp = selectLoci(chr,
                   simParam$traits[[trait]]@lociPerChr,
                   simParam$traits[[trait]]@lociLoc)
  output = getGeno(pop@geno,tmp$lociPerChr,tmp$lociLoc)
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
#' @export
pullSegSiteGeno = function(pop, chr=NULL, simParam=NULL){
  if(class(pop)=="MapPop"){
    allLoci = unlist(sapply(pop@nLoci,
                            function(x) 1:x))
    lociTot = pop@nLoci
  }else{
    if(is.null(simParam)){
      simParam = get("SP",envir=.GlobalEnv)
    }
    allLoci = unlist(sapply(simParam$segSites,
                            function(x) 1:x))
    lociTot = simParam$segSites
  }
  tmp = selectLoci(chr,lociTot,allLoci)
  output = getGeno(pop@geno,tmp$lociPerChr,tmp$lociLoc)
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
    output = getHaplo(pop@geno,lociPerChr,lociLoc)
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
                         as.integer(haplo))
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
                     simParam$snpChips[[snpChip]]@lociLoc)
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
                        as.integer(haplo))
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
    output = getHaplo(pop@geno,lociPerChr,lociLoc)
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
                         as.integer(haplo))
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
#' @export
pullSegSiteHaplo = function(pop, haplo="all", 
                            chr=NULL, simParam=NULL){
  if(class(pop)=="MapPop"){
    allLoci = unlist(sapply(pop@nLoci,
                            function(x) 1:x))
    lociTot = pop@nLoci
  }else{
    if(is.null(simParam)){
      simParam = get("SP",envir=.GlobalEnv)
    }
    allLoci = unlist(sapply(simParam$segSites,
                            function(x) 1:x))
    lociTot = simParam$segSites
  }
  if(!is.null(chr)){
    tmp = selectLoci(chr,lociTot,allLoci)
    lociTot = tmp$lociPerChr
    allLoci = tmp$lociLoc
  }
  if(haplo=="all"){
    output = getHaplo(pop@geno,
                      lociTot,
                      allLoci)
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
                         as.integer(haplo))
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
