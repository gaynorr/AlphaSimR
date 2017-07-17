#' @title Pull SNP genotype
#' 
#' @description Retrieves SNP genotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param snpChip an integer. Indicates which SNP 
#' chip's genotypes to retrieve.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of SNP genotypes.
#' @export
pullSnpGeno = function(pop, snpChip=1, simParam=SIMPARAM){
  output = getGeno(pop@geno,
                   simParam@snpChips[[snpChip]]@lociPerChr,
                   simParam@snpChips[[snpChip]]@lociLoc)
  output = convToImat(output)
  rownames(output) = pop@id
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
#' @param simParam an object of \code{\link{SimParam-class}}#'
#' @return Returns a matrix of SNP genotypes.
#' @export
pullMultipleSnpGeno = function(pop, chips,
                              missing = 9, simParam=SIMPARAM) {
  stopifnot(length(chips) == pop@nInd)
  # I feel like the next line shouldn't be needed but I don't know
  # enough R! (dmoney)
  missing = as.integer(missing)
  allSnps = numeric(0)
  uniqueChips = unique(chips)
  for (c in uniqueChips){
    allSnps = sort(union(allSnps,simParam@snpChips[[c]]@lociLoc))
  }
  
  output = matrix(pop@nInd,length(allSnps),data=missing)
  rownames(output) = pop@id
  
  for (snpChip in uniqueChips){
    mask = allSnps %in% simParam@snpChips[[snpChip]]@lociLoc
    one = getGeno(pop@geno,
                      simParam@snpChips[[snpChip]]@lociPerChr,
                      simParam@snpChips[[snpChip]]@lociLoc)
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
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of QTL genotypes.
#' @export
pullQtlGeno = function(pop, trait=1, simParam=SIMPARAM){
  output = getGeno(pop@geno,
                   simParam@traits[[trait]]@lociPerChr,
                   simParam@traits[[trait]]@lociLoc)
  output = convToImat(output)
  rownames(output) = pop@id
  colnames(output) = paste("QTL",1:ncol(output),sep="_")
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
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of SNP haplotypes.
#' @export
pullSnpHaplo = function(pop, snpChip=1, haplo="all", 
                        simParam=SIMPARAM){
  if(haplo=="all"){
    output = getHaplo(pop@geno,
                      simParam@snpChips[[snpChip]]@lociPerChr,
                      simParam@snpChips[[snpChip]]@lociLoc)
    output = convToImat(output)
    rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                             rep(1:pop@ploidy,pop@nInd),sep="_")
  }else{
    stopifnot(haplo%in%c(1,2))
    output = getOneHaplo(pop@geno,
                         simParam@snpChips[[snpChip]]@lociPerChr,
                         simParam@snpChips[[snpChip]]@lociLoc,
                         as.integer(haplo))
    output = convToImat(output)
    rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    
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
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of SNP haplotypes.
#' @export
pullMultipleSnpHaplo = function(pop, chips, haplo="all", 
                                missing = 9, simParam=SIMPARAM) {
  stopifnot(length(chips) == pop@nInd)
  # I feel like the next line shouldn't be needed but I don't know
  # enough R! (dmoney)
  missing = as.integer(missing)
  allSnps = numeric(0)
  uniqueChips = unique(chips)
  for (c in uniqueChips){
    allSnps = sort(union(allSnps,simParam@snpChips[[c]]@lociLoc))
  }
  
  if (haplo == "all") {
    output = matrix(pop@nInd*2,length(allSnps),data=missing)
    rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                             rep(1:pop@ploidy,pop@nInd),sep="_")
  }
  else {
    output = matrix(pop@nInd,length(allSnps),data=missing)
    rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
  }
  for (snpChip in uniqueChips){
    mask = allSnps %in% simParam@snpChips[[snpChip]]@lociLoc
    if (haplo == "all") {
      one = getHaplo(pop@geno,
                     simParam@snpChips[[snpChip]]@lociPerChr,
                     simParam@snpChips[[snpChip]]@lociLoc)
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
                     simParam@snpChips[[snpChip]]@lociPerChr,
                     simParam@snpChips[[snpChip]]@lociLoc,
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
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of QTL haplotypes.
#' @export
pullQtlHaplo = function(pop, trait=1, haplo="all", 
                        simParam=SIMPARAM){
  if(haplo=="all"){
    output = getHaplo(pop@geno,
                      simParam@traits[[trait]]@lociPerChr,
                      simParam@traits[[trait]]@lociLoc)
    output = convToImat(output)
    rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                             rep(1:pop@ploidy,pop@nInd),sep="_")
  }else{
    stopifnot(haplo%in%c(1,2))
    output = getOneHaplo(pop@geno,
                         simParam@traits[[trait]]@lociPerChr,
                         simParam@traits[[trait]]@lociLoc,
                         as.integer(haplo))
    output = convToImat(output)
    rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    
  }
  colnames(output) = paste("QTL",1:ncol(output),sep="_")
  return(output)
}
