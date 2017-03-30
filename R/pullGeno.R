#' @title Pull SNP genotype
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param snpChip an integer. Indicates which SNP 
#' chip's genotypes to retrieve.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of SNP genotypes.
#' @export
pullSnpGeno = function(pop, snpChip=1, simParam=SIMPARAM){
  output = AlphaSimR:::getGeno(pop@geno,
                               simParam@snpChips[[snpChip]]@lociPerChr,
                               simParam@snpChips[[snpChip]]@lociLoc)
  output = AlphaSimR:::convToImat(output)
  rownames(output) = pop@id
  colnames(output) = paste("SNP",1:ncol(output),sep="_")
  return(output)
}

#' @title Pull QTL genotype
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param trait an integer. Indicates which trait's
#' QTL genotypes to retrieve.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of QTL genotypes.
#' @export
pullQtlGeno = function(pop, trait=1, simParam=SIMPARAM){
  output = AlphaSimR:::getGeno(pop@geno,
                               simParam@traits[[trait]]@lociPerChr,
                               simParam@traits[[trait]]@lociLoc)
  output = AlphaSimR:::convToImat(output)
  rownames(output) = pop@id
  colnames(output) = paste("QTL",1:ncol(output),sep="_")
  return(output)
}

#' @title Pull SNP haplotypes
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param snpChip an integer. Indicates which SNP 
#' chip's haplotypes to retrieve.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of SNP haplotypes.
#' @export
pullSnpHaplo = function(pop, snpChip=1, simParam=SIMPARAM){
  output = AlphaSimR:::getHaplo(pop@geno,
                                simParam@snpChips[[snpChip]]@lociPerChr,
                                simParam@snpChips[[snpChip]]@lociLoc)
  output = AlphaSimR:::convToImat(output)
  rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                           rep(1:pop@ploidy,pop@nInd),sep="_")
  colnames(output) = paste("SNP",1:ncol(output),sep="_")
  return(output)
}

#' @title Pull QTL haplotypes
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param trait an integer. Indicates which trait's
#' QTL haplotypes to retrieve.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of QTL haplotypes.
#' @export
pullQtlHaplo = function(pop, trait=1, simParam=SIMPARAM){
  output = AlphaSimR:::getHaplo(pop@geno,
                                simParam@traits[[trait]]@lociPerChr,
                                simParam@traits[[trait]]@lociLoc)
  output = AlphaSimR:::convToImat(output)
  rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                           rep(1:pop@ploidy,pop@nInd),sep="_")
  colnames(output) = paste("QTL",1:ncol(output),sep="_")
  return(output)
}
