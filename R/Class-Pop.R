
# RawPop ------------------------------------------------------------------

#' @title Raw Population
#' 
#' @description 
#' The raw population class contains only genotype data. 
#' 
#' @param x a 'RawPop'
#' @param i index of individuals
#' @param ... additional 'RawPop' objects
#' 
#' @slot nInd number of individuals
#' @slot nChr number of chromosomes
#' @slot ploidy level of ploidy
#' @slot nLoci number of loci per chromosome
#' @slot geno "matrix" containing chromosome genotypes. The "matrix" 
#' has dimensions nChr by 1 and each element is a three dimensional
#' array of raw values. The array dimensions are nLoci by ploidy by nInd.
#' 
#' 
#' @export
setClass("RawPop",
         slots=c(nInd="integer",
                 nChr="integer",
                 ploidy="integer",
                 nLoci="integer",
                 geno="matrix"))

setValidity("RawPop",function(object){
  errors = character()
  if(object@nChr!=length(object@geno)){
    errors = c(errors,"nChr!=length(geno)")
  }
  if(object@nChr!=length(object@nLoci)){
    errors = c(errors,"nChr!=length(nLoci)")
  }
  for(i in 1:object@nChr){
    if(object@nLoci[i]!=dim(object@geno[[i]])[1]){
      errors = c(errors,
                 paste0("nLoci[",i,
                        "]!=dim(geno[[",i,
                        "]][1]"))
    }
    if(object@ploidy!=dim(object@geno[[i]])[2]){
      errors = c(errors,
                 paste0("ploidy!=dim(geno[[",i,
                        "]][2]"))
    }
    if(object@nInd!=dim(object@geno[[i]])[3]){
      errors = c(errors,
                 paste0("nInd!=dim(geno[[",i,
                        "]][3]"))
    }
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})


#' @describeIn RawPop Extract RawPop by index
setMethod("[",
          signature(x = "RawPop"),
          function(x, i){
            for(chr in 1:x@nChr){
              x@geno[[chr]] = x@geno[[chr]][,,i,drop=FALSE]
            }
            x@nInd = dim(x@geno[[1]])[3]
            validObject(x)
            return(x)
          }
)


#' @describeIn RawPop Combine multiple RawPops
setMethod("c",
          signature(x = "RawPop"),
          function (x, ...){
            for(y in list(...)){
              stopifnot(class(y)=="RawPop",
                        x@nChr==y@nChr,
                        x@ploidy==y@ploidy,
                        x@nLoci==y@nLoci)
              x@nInd = x@nInd+y@nInd
              x@geno = mergeGeno(x@geno,y@geno)
            }
            validObject(x)
            return(x)
          }
)

# MapPop ------------------------------------------------------------------

#' @title Raw population with genetic map
#' 
#' @description 
#' Extends \code{\link{RawPop-class}} to add a genetic map. 
#' This is the first object created in a simulation. It is used
#' for creating initial populations and setting traits in the 
#' \code{\link{SimParam-class}}.
#' 
#' @param x a 'MapPop'
#' @param i index of chromosomes
#' @param ... aditional 'MapPop' objects
#' 
#' @slot genMaps "matrix" of chromsome genetic maps
#' 
#' @export
setClass("MapPop",
         slots=c(genMaps="matrix"),
         contains="RawPop")

setValidity("MapPop",function(object){
  errors = character()
  if(object@nChr!=length(object@genMaps)){
    errors = c(errors,"nInd!=length(id)")
  }
  for(i in 1:object@nChr){
    if(object@nLoci[i]!=length(object@genMaps[[i]])){
      errors = c(errors,
                 paste0("nLoci[",i,"]!=length(genMaps[[",i,"]]"))
    }
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#' @describeIn MapPop Extract \code{\link{RawPop-class}} by index
setMethod("[",
          signature(x = "MapPop"),
          function(x, i){
            for(chr in 1:x@nChr){
              x@geno[[chr]] = x@geno[[chr]][,,i,drop=FALSE]
            }
            x@nInd = dim(x@geno[[1]])[3]
            class(x) = "RawPop"
            validObject(x)
            return(x)
          }
)

#' @describeIn MapPop Combine MapPop chromosomes
setMethod("c",
          signature(x = "MapPop"),
          function (x, ...){
            for(y in list(...)){
              stopifnot(class(y)=="MapPop",
                        x@nInd==y@nInd,
                        x@ploidy==y@ploidy)
              x@nChr = x@nChr+y@nChr
              x@geno = rbind(x@geno,y@geno)
              x@genMaps = rbind(x@genMaps,y@genMaps)
              x@nLoci = c(x@nLoci,y@nLoci)
            }
            validObject(x)
            return(x)
          }
)

# Pop ---------------------------------------------------------------------

#' @title Population
#' 
#' @description 
#' Extends \code{\link{RawPop-class}} to add gender, genetic values, 
#' phenotypes, and pedigrees.
#' 
#' @param x a 'Pop'
#' @param i index of individuals
#' @param ... additional 'Pop' objects
#' 
#' @slot id an individual's identifier
#' @slot mother the identifier of the individual's mother
#' @slot father the identifier of the individual's father
#' @slot gender gender of individuals
#' @slot nTraits number of traits
#' @slot gv matrix of genetic values. When using GxE traits,
#' gv reflects gv when w=0. Dimensions are nInd by nTraits.
#' @slot pheno matrix of phenotypic values. Dimensions are
#' nInd by nTraits.
#' @slot ebv matrix of estimated breeding values. Dimensions 
#' are nInd rows and a variable number of columns.
#' @slot gxe list containing GxE slopes for GxE traits
#' 
#' @export
setClass("Pop",
         slots=c(id="character",
                 mother="character",
                 father="character",
                 gender="character",
                 nTraits="integer",
                 gv="matrix",
                 pheno="matrix",
                 ebv="matrix",
                 gxe="list"),
         contains="RawPop")

setValidity("Pop",function(object){
  errors = character()
  if(any(grepl(" ",object@id,fixed=TRUE))){
    errors = c(errors,"id can not contain spaces")
  }
  if(any(grepl(" ",object@mother,fixed=TRUE))){
    errors = c(errors,"mother can not contain spaces")
  }
  if(any(grepl(" ",object@father,fixed=TRUE))){
    errors = c(errors,"father can not contain spaces")
  }
  if(object@nInd!=length(object@gender)){
    errors = c(errors,"nInd!=length(gender)")
  }
  if(object@nInd!=length(object@id)){
    errors = c(errors,"nInd!=length(id)")
  }
  if(object@nInd!=length(object@mother)){
    errors = c(errors,"nInd!=length(mother)")
  }
  if(object@nInd!=length(object@father)){
    errors = c(errors,"nInd!=length(father)")
  }
  if(object@nInd!=nrow(object@gv)){
    errors = c(errors,"nInd!=nrow(gv)")
  }
  if(object@nInd!=nrow(object@pheno)){
    errors = c(errors,"nInd!=nrow(pheno)")
  }
  if(object@nInd!=nrow(object@ebv)){
    errors = c(errors,"nInd!=nrow(ebv)")
  }
  if(!is.numeric(object@gv)){
    errors = c(errors,"!is.numeric(gv)")
  }
  if(!is.numeric(object@pheno)){
    errors = c(errors,"!is.numeric(pheno)")
  }
  if(!is.numeric(object@ebv)){
    errors = c(errors,"!is.numeric(ebv)")
  }
  if(object@nTraits!=ncol(object@gv)){
    errors = c(errors,"nTraits!=ncol(gv)")
  }
  if(object@nTraits!=ncol(object@pheno)){
    errors = c(errors,"nTraits!=ncol(pheno)")
  }
  if(object@nTraits!=length(object@gxe)){
    errors = c(errors,"nTraits!=length(gxe)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#' @describeIn Pop Extract Pop by index or id
setMethod("[",
          signature(x = "Pop"),
          function(x, i){
            n = x@nInd
            if(is.character(i)){
              i = x@id%in%i
            }
            x@id = x@id[i]
            x@mother = x@mother[i]
            x@father = x@father[i]
            x@gv = x@gv[i,,drop=FALSE]
            x@pheno = x@pheno[i,,drop=FALSE]
            x@ebv = x@ebv[i,,drop=FALSE]
            x@gender = x@gender[i]
            x@nInd = length(x@gender)
            if(x@nTraits>=1){
              for(trait in 1:x@nTraits){
                if(!is.null(x@gxe[[trait]])){
                  x@gxe[[trait]] = x@gxe[[trait]][i]
                }
              }
            }
            for(chr in 1:x@nChr){
              x@geno[[chr]] = x@geno[[chr]][,,i,drop=FALSE]
            }
            validObject(x)
            return(x)
          }
)

#' @describeIn Pop Combine multiple Pops
setMethod("c",
          signature(x = "Pop"),
          function (x, ...){
            for(y in list(...)){
              stopifnot(class(y)=="Pop",
                        x@nChr==y@nChr,
                        x@ploidy==y@ploidy,
                        x@nLoci==y@nLoci)
              x@nInd = x@nInd+y@nInd
              x@id = c(x@id,y@id)
              x@mother = c(x@mother,y@mother)
              x@father = c(x@father,y@father)
              x@gv = rbind(x@gv,y@gv)
              x@pheno = rbind(x@pheno,y@pheno)
              x@gender = c(x@gender,y@gender)
              x@geno = mergeGeno(x@geno,y@geno)
              if(x@nTraits>=1){
                for(trait in 1:x@nTraits){
                  if(!is.null(x@gxe[[trait]])){
                    x@gxe[[trait]] = c(x@gxe[[trait]],y@gxe[[trait]])
                  }
                }
              }
              #Account for variable number of ebv columns
              tmp = try({
                x@ebv = rbind(x@ebv,y@ebv)
                },silent=TRUE)
              if(class(tmp)=="try-error"){
                x@ebv = matrix(NA_real_,nrow=x@nInd,ncol=0)
              }
            }
            validObject(x)
            return(x)
          }
)

#' @title Reset population
#' 
#' @description
#' Recalculates a population's genetic values and 
#' resets phenotypes and EBVs.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return an object of \code{\link{Pop-class}}
#' 
#' @export
resetPop = function(pop,simParam){
  pop@nTraits = simParam@nTraits
  pop@pheno = matrix(NA_real_,
                     nrow=pop@nInd,
                     ncol=simParam@nTraits)
  pop@ebv = matrix(NA_real_,
                   nrow=pop@nInd,
                   ncol=0)
  pop@gxe = vector("list",simParam@nTraits)
  pop@gv = matrix(NA_real_,nrow=pop@nInd,
                  ncol=simParam@nTraits)
  if(simParam@nTraits>=1){
    for(i in 1:simParam@nTraits){
      tmp = getGv(simParam@traits[[i]],pop)
      pop@gv[,i] = tmp[[1]]
      if(length(tmp)>1){
        pop@gxe[[i]] = tmp[[2]]
      }
    }
  }
  validObject(pop)
  return(pop)
}
