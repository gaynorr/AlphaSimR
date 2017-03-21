#Pop----
#' @title Population
#' 
#' @description Population class
#' 
#' @slot nInd number of individuals
#' @slot nChr number of chromosomes
#' @slot ploidy level of ploidy
#' @slot gender gender of individuals
#' @slot geno list containing chromosome genotypes
#' 
#' @export
setClass("Pop",
         slots=c(nInd="integer",
                 nChr="integer",
                 ploidy="integer",
                 gender="character",
                 geno="list"))

setValidity("Pop",function(object){
  errors = character()
  if(object@nInd!=length(object@gender)){
    errors = c(errors,"nInd!=length(gender)")
  }
  if(object@nChr!=length(object@geno)){
    errors = c(errors,"nChr!=length(geno)")
  }
  if(object@ploidy!=length(object@geno[[1]])){
    errors = c(errors,"ploidy!=length(geno[[1]])")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

setMethod("[",
          signature(x = "Pop"),
          function(x, i, j=NULL, ..., drop = TRUE){
            x@gender = x@gender[i]
            x@nInd = length(x@gender)
            for(chr in 1:x@nChr){
              for(chrI in 1:x@ploidy){
                x@geno[[chr]][[chrI]] = 
                  matrix(x@geno[[chr]][[chrI]][i,],nrow=x@nInd)
              }
            }
            validObject(x)
            return(x)
          }
)

#TraitPop----
#' @title Population with traits
#' 
#' @description Extends \code{\link{Pop-class}}
#' 
#' @slot gv matrix of genetic values
#' @slot pheno matrix of phenotypic values
#' @slot nTraits number of traits
#'
#' @export
setClass("TraitPop",
         slots=c(gv="matrix",
                 pheno="matrix",
                 nTraits="integer"),
         contains="Pop")

setValidity("TraitPop",function(object){
  errors = character()
  if(object@nInd!=nrow(object@gv)){
    errors = c(errors,"nInd!=nrow(gv)")
  }
  if(object@nInd!=nrow(object@pheno)){
    errors = c(errors,"nInd!=nrow(pheno)")
  }
  if(ncol(object@gv)!=object@nTraits){
    errors = c(errors,"ncol(gv)!=nTraits")
  }
  if(ncol(object@pheno)!=object@nTraits){
    errors = c(errors,"ncol(pheno)!=nTraits")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

setMethod("[",
          signature(x = "TraitPop"),
          function(x, i, j=NULL, ..., drop = TRUE){
            x@gv = matrix(x@gv[i,],ncol=x@nTraits)
            x@pheno = matrix(x@pheno[i,],ncol=x@nTraits)
            x@gender = x@gender[i]
            x@nInd = length(x@gender)
            for(chr in 1:x@nChr){
              for(chrI in 1:x@ploidy){
                x@geno[[chr]][[chrI]] = 
                  matrix(x@geno[[chr]][[chrI]][i,],nrow=x@nInd)
              }
            }
            validObject(x)
            return(x)
          }
)

#PedPop----
#' @title Population with pedigree
#' 
#' @description Extends \code{\link{TraitPop-class}}
#' 
#' @slot id individual's identifier
#' @slot par1 individual's female parent
#' @slot par2 individual's male parent
#'
#' @export
setClass("PedPop",
         slots=c(id="character",
                 par1="character",
                 par2="character"),
         contains="TraitPop")

setValidity("PedPop",function(object){
  errors = character()
  if(object@nInd!=length(object@id)){
    errors = c(errors,"nInd!=length(id)")
  }
  if(object@nInd!=length(object@par1)){
    errors = c(errors,"nInd!=length(par1)")
  }
  if(object@nInd!=length(object@par2)){
    errors = c(errors,"nInd!=length(par2)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

setMethod("[",
          signature(x = "PedPop"),
          function(x, i, j=NULL, ..., drop = TRUE){
            x@id = x@id[i]
            x@par1 = x@par1[i]
            x@par2 = x@par2[i]
            x@gv = matrix(x@gv[i,],ncol=x@nTraits)
            x@pheno = matrix(x@pheno[i,],ncol=x@nTraits)
            x@gender = x@gender[i]
            x@nInd = length(x@gender)
            for(chr in 1:x@nChr){
              for(chrI in 1:x@ploidy){
                x@geno[[chr]][[chrI]] = 
                  matrix(x@geno[[chr]][[chrI]][i,],nrow=x@nInd)
              }
            }
            validObject(x)
            return(x)
          }
)

#' @title Add genetic values
#' 
#' @description Promotes class 'Pop' to class 'TraitPop'
#' 
#' @param pop an object of class 'Pop'
#' @param simParam an object of class 'SimParam'
#' 
#' @export
addGv = function(pop, simParam=SIMPARAM){
  stopifnot(class(pop)=="Pop")
  gv = lapply(simParam@traits,getGv,pop=pop,w=0)
  gv = do.call("cbind",gv)
  pheno = matrix(NA_real_,nrow=nrow(gv),ncol=ncol(gv))
  pop = new("TraitPop",pop,gv=gv,pheno=pheno,
            nTraits=simParam@nTraits)
  return(pop)
}

#' @title Add pedigree
#' 
#' @description Promotes class 'Pop' or 'TraitPop' to 'PedPop'
#' 
#' @param pop an object of class 'Pop' or 'TraitPop'
#' @param id a unique name for each individual
#' @param par1 the first/female parent for each individual
#' @param par2 the second/male parent for each individual
#' @param simParam an object of class 'SimParam'
#' 
#' @export
addPed = function(pop, id, par1, par2, simParam=SIMPARAM){
  if(class(pop)=="Pop"){
    pop = addGv(pop,simParam=simParam)
  }
  stopifnot(class(pop)=="TraitPop")
  pop = new("PedPop",pop,id=as.character(id),
            par1=as.character(par1),par2=as.character(par2))
  return(pop)
}