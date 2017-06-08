#' @title Calculate G matrix
#' 
#' @description Calculates the genomic relationship matrix. Note 
#' this function uses single precision and \code{\link{calcG}} 
#' uses double precision.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param useSnp use a SNP chip to construct G matrix. If 
#' false, QTL genotypes are used to construct the G matrix.
#' @param snpChip an integer. Indicates which SNP 
#' chip's genotypes are used to construct G, if useSnp=TRUE.
#' @param trait an integer. Indicates which trait's 
#' genotypes are used to contruct G, if useSnp=FALSE.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of genomic relationships
#' 
#' @export
calcPopG = function(pop, useSnp=TRUE, snpChip=1, trait=1, 
                    simParam=SIMPARAM){
  if(useSnp){
    return(.calcPopG(pop@geno, 
                     simParam@snpChips[[snpChip]]@lociPerChr,
                     simParam@snpChips[[snpChip]]@lociLoc))
  }else{
    return(.calcPopG(pop@geno, 
                     simParam@traits[[trait]]@lociPerChr,
                     simParam@traits[[trait]]@lociLoc))
  }
}

#' @title Calculate IBS G matrix
#' 
#' @description Calculates an identity-by-state genomic 
#' relationship matrix based on simple matching. Note that 
#' this function uses single precision and \code{\link{calcGIbs}} 
#' uses double precision.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param useSnp use a SNP chip to construct G matrix. If 
#' false, QTL genotypes are used to construct the G matrix.
#' @param snpChip an integer. Indicates which SNP 
#' chip's genotypes are used to construct G, if useSnp=TRUE.
#' @param trait an integer. Indicates which trait's 
#' genotypes are used to contruct G, if useSnp=FALSE.
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns a matrix of genomic relationships
#' 
#' @export
calcPopGIbs = function(pop, useSnp=TRUE, snpChip=1, trait=1, 
                       simParam=SIMPARAM){
  if(useSnp){
    return(.calcPopGIbs(pop@geno, 
                        simParam@snpChips[[snpChip]]@lociPerChr,
                        simParam@snpChips[[snpChip]]@lociLoc))
  }else{
    return(.calcPopGIbs(pop@geno, 
                        simParam@traits[[trait]]@lociPerChr,
                        simParam@traits[[trait]]@lociLoc))
  }
}