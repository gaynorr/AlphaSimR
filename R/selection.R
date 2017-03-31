#' @title Select individuals
#' 
#' @description Selects a subset of nInd individuals from a 'Pop' superclass using various types of traits.
#' 
#' @param pop and object of superclass 'Pop' or 'HybridPop'
#' @param nInd the number of individuals to select
#' @param trait the trait for selection. Either a number for one of the traits or an object of superclass 'SelIndex'
#' @param useGv should genetic value be used instead of phenotypes
#' @param selectTop selects highest values if true. Selects lowest if false.
#' 
#' @export
selectPop = function(pop,nInd,trait=1,useGv=FALSE,
                     selectTop=TRUE){
  if(class(pop)=="Pop") stop("Must call addGv first")
  stopifnot(nInd<=pop@nInd)
  if(class(trait)=="SelIndex"){
    stop("Not currently implemented")
  }else{
    stopifnot(length(trait)==1,trait<=pop@nTraits)
    if(useGv){
      response = pop@gv[,trait]
    }else{
      response = pop@pheno[,trait]
    }
  }
  take = order(response,decreasing=selectTop)
  return(pop[take[1:nInd]])
}
