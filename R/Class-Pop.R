# RawPop ------------------------------------------------------------------

#' @title Raw Population
#' 
#' @description 
#' The raw population class contains only genotype data. 
#' 
#' @param object a 'RawPop' object
#' @param x a 'RawPop' object
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
    DIM1 = object@nLoci[i]%/%8L + (object@nLoci[i]%%8L > 0L)
    if(DIM1!=dim(object@geno[[i]])[1]){
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
            if(any(abs(i)>x@nInd)){
              stop("Trying to select invalid individuals")
            }
            for(chr in 1:x@nChr){
              x@geno[[chr]] = x@geno[[chr]][,,i,drop=FALSE]
            }
            x@nInd = dim(x@geno[[1]])[3]
            return(x)
          }
)


#' @describeIn RawPop Combine multiple RawPops
setMethod("c",
          signature(x = "RawPop"),
          function (x, ...){
            for(y in list(...)){
              if(class(y)=="NULL"){
                # Do nothing
              }else{
                stopifnot(class(y)=="RawPop",
                          x@nChr==y@nChr,
                          x@ploidy==y@ploidy,
                          x@nLoci==y@nLoci)
                x@nInd = x@nInd+y@nInd
                x@geno = mergeGeno(x@geno,y@geno)
              }
            }
            return(x)
          }
)

#' @describeIn RawPop Show population summary
setMethod("show",
          signature(object = "RawPop"),
          function (object){
            cat("An object of class", 
                classLabel(class(object)), "\n")
            cat("Ploidy:", object@ploidy,"\n")
            cat("Individuals:", object@nInd,"\n")
            cat("Chromosomes:", object@nChr,"\n")
            cat("Loci:", sum(object@nLoci),"\n")
            invisible()
          }
)

# MapPop ------------------------------------------------------------------

#' @title Raw population with genetic map
#' 
#' @description 
#' Extends \code{\link{RawPop-class}} to add a genetic map. 
#' This is the first object created in a simulation. It is used
#' for creating initial populations and setting traits in the 
#' \code{\link{SimParam}}.
#' 
#' @param x a 'MapPop' object
#' @param i index of individuals
#' @param ... additional 'MapPop' objects
#' 
#' @slot genMap "matrix" of chromosome genetic maps
#' @slot centromere vector of centromere positions
#' @slot inbred indicates whether the individuals are fully inbred
#' 
#' @export
setClass("MapPop",
         slots=c(genMap="list",
                 centromere="numeric",
                 inbred="logical"),
         contains="RawPop")

setValidity("MapPop",function(object){
  errors = character()
  if(object@nChr!=length(object@genMap)){
    errors = c(errors,"nInd!=length(id)")
  }
  for(i in 1:object@nChr){
    if(object@nLoci[i]!=length(object@genMap[[i]])){
      errors = c(errors,
                 paste0("nLoci[",i,"]!=length(genMap[[",i,"]]"))
    }
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#' @describeIn MapPop Extract MapPop by index
setMethod("[",
          signature(x = "MapPop"),
          function(x, i){
            if(any(abs(i)>x@nInd)){
              stop("Trying to select invalid individuals")
            }
            for(chr in 1:x@nChr){
              x@geno[[chr]] = x@geno[[chr]][,,i,drop=FALSE]
            }
            x@nInd = dim(x@geno[[1]])[3]
            return(x)
          }
)

#' @describeIn MapPop Combine multiple MapPops 
setMethod("c",
          signature(x = "MapPop"),
          function (x, ...){
            for(y in list(...)){
              if(class(y)=="NULL"){
                # Do nothing
              }else{
                stopifnot(class(y)=="MapPop",
                          x@nChr==y@nChr,
                          x@ploidy==y@ploidy,
                          x@nLoci==y@nLoci,
                          all.equal(x@genMap, y@genMap))
                x@nInd = x@nInd+y@nInd
                x@geno = mergeGeno(x@geno,y@geno)
                x@inbred = x@inbred & y@inbred
              }
            }
            return(x)
          }
)

# NamedMapPop ------------------------------------------------------------------

#' @title Raw population with genetic map and id
#' 
#' @description 
#' Extends \code{\link{MapPop-class}} to add id, mother and father. 
#' 
#' @param x a 'NamedMapPop' object
#' @param i index of individuals
#' @param ... additional 'NamedMapPop' objects
#' 
#' @slot id an individual's identifier
#' @slot mother the identifier of the individual's mother
#' @slot father the identifier of the individual's father
#' 
#' @export
setClass("NamedMapPop",
         slots=c(id="character",
                 mother="character",
                 father="character"),
         contains="MapPop")

setValidity("NamedMapPop",function(object){
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
  if(object@nInd!=length(object@id)){
    errors = c(errors,"nInd!=length(id)")
  }
  if(object@nInd!=length(object@mother)){
    errors = c(errors,"nInd!=length(mother)")
  }
  if(object@nInd!=length(object@father)){
    errors = c(errors,"nInd!=length(father)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#' @describeIn NamedMapPop Extract NamedMapPop by index
setMethod("[",
          signature(x = "NamedMapPop"),
          function(x, i){
            if(any(abs(i)>x@nInd)){
              stop("Trying to select invalid individuals")
            }
            for(chr in 1:x@nChr){
              x@geno[[chr]] = x@geno[[chr]][,,i,drop=FALSE]
            }
            x@nInd = dim(x@geno[[1]])[3]
            x@id = x@id[i]
            x@mother = x@mother[i]
            x@father = x@father[i]
            return(x)
          }
)

#' @describeIn NamedMapPop Combine multiple NamedMapPops 
setMethod("c",
          signature(x = "NamedMapPop"),
          function (x, ...){
            for(y in list(...)){
              if(class(y)=="NULL"){
                # Do nothing
              }else{
                stopifnot(class(y)=="NamedMapPop",
                          x@nChr==y@nChr,
                          x@ploidy==y@ploidy,
                          x@nLoci==y@nLoci,
                          all.equal(x@genMap, y@genMap))
                x@nInd = x@nInd+y@nInd
                x@id = c(x@id, y@id)
                x@mother = c(x@mother, y@mother)
                x@father = c(x@father, y@father)
                x@geno = mergeGeno(x@geno,y@geno)
                x@inbred = x@inbred & y@inbred
              }
            }
            return(x)
          }
)

#' @title Combine MapPop chromosomes
#' 
#' @description
#' Merges the chromosomes of multiple \code{\link{MapPop-class}} or 
#' \code{\link{NamedMapPop-class}} objects. 
#' Each MapPop must have the same number of chromosomes
#'
#' @param ... \code{\link{MapPop-class}} or \code{\link{NamedMapPop-class}} 
#' objects to be combined
#' 
#' @return Returns an object of \code{\link{MapPop-class}}
#' 
#' @examples 
#' pop1 = quickHaplo(nInd=10, nChr=1, segSites=10)
#' pop2 = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' combinedPop = cChr(pop1, pop2)
#'
#' @export
cChr = function(...){
  for(y in list(...)){
    if(class(y)=="NULL"){
      #Do nothing
    }else{
      stopifnot(class(y)=="MapPop" | class(y)=="NamedMapPop")
      if(!exists("x",inherits=FALSE)){
        x = y
      }else{
        stopifnot(x@nInd==y@nInd,
                  x@ploidy==y@ploidy)
        x@nChr = x@nChr+y@nChr
        x@geno = rbind(x@geno,y@geno)
        x@genMap = c(x@genMap,y@genMap)
        x@centromere = c(x@centromere,y@centromere)
        x@nLoci = c(x@nLoci,y@nLoci)
        x@inbred = x@inbred & y@inbred
      }
    }
  }
  return(x)
}

# Pop ---------------------------------------------------------------------

#' @title Population
#' 
#' @description 
#' Extends \code{\link{RawPop-class}} to add sex, genetic values, 
#' phenotypes, and pedigrees.
#' 
#' @param object a 'Pop' object
#' @param x a 'Pop' object
#' @param i index of individuals
#' @param ... additional 'Pop' objects
#' 
#' @slot id an individual's identifier
#' @slot iid an individual's internal identifier
#' @slot mother the identifier of the individual's mother
#' @slot father the identifier of the individual's father
#' @slot sex sex of individuals: "M" for males, "F" for females,
#' and "H" for hermaphrodites
#' @slot nTraits number of traits
#' @slot gv matrix of genetic values. When using GxE traits,
#' gv reflects gv when p=0.5. Dimensions are nInd by nTraits.
#' @slot pheno matrix of phenotypic values. Dimensions are
#' nInd by nTraits.
#' @slot ebv matrix of estimated breeding values. Dimensions 
#' are nInd rows and a variable number of columns.
#' @slot gxe list containing GxE slopes for GxE traits
#' @slot fixEff a fixed effect relating to the phenotype. 
#' Used by genomic selection models but otherwise ignored.
#' @slot reps the number of replications used to measure the 
#' phenotype. Used by genomic selection models, but otherwise ignored.
#' @slot misc a list whose elements correspond to individuals in the 
#' population. This list is normally empty and exists solely as an 
#' open slot available for uses to store extra information about 
#' individuals.
#' 
#' @export
setClass("Pop",
         slots=c(id="character",
                 iid="integer",
                 mother="character",
                 father="character",
                 sex="character",
                 nTraits="integer",
                 gv="matrix",
                 pheno="matrix",
                 ebv="matrix",
                 gxe="list",
                 fixEff="integer",
                 reps="numeric",
                 misc="list"),
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
  if(object@nInd!=length(object@sex)){
    errors = c(errors,"nInd!=length(sex)")
  }
  if(object@nInd!=length(object@id)){
    errors = c(errors,"nInd!=length(id)")
  }
  if(object@nInd!=length(object@iid)){
    errors = c(errors,"nInd!=length(iid)")
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
  if(object@nInd!=length(object@fixEff)){
    errors = c(errors,"nInd!=length(fixEff)")
  }
  if(object@nInd!=length(object@reps)){
    errors = c(errors,"nInd!=length(reps)")
  }
  if(object@nInd!=length(object@misc)){
    errors = c(errors,"nInd!=length(misc)")
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
            if(is.character(i)){
              i = match(i, x@id)
              if(any(is.na(i))){
                stop("Trying to select invalid individuals")
              }
              if(any(is.null(i))){
                stop("Not valid ids")
              }
            }else{
              if(any(abs(i)>x@nInd)){
                stop("Trying to select invalid individuals")
              }
            }
            x@id = x@id[i]
            x@iid = x@iid[i]
            x@mother = x@mother[i]
            x@father = x@father[i]
            x@fixEff = x@fixEff[i]
            x@reps = x@reps[i]
            x@misc = x@misc[i]
            x@gv = x@gv[i,,drop=FALSE]
            x@pheno = x@pheno[i,,drop=FALSE]
            x@ebv = x@ebv[i,,drop=FALSE]
            x@sex = x@sex[i]
            x@nInd = length(x@sex)
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
            return(x)
          }
)

#' @describeIn Pop Combine multiple Pops
setMethod("c",
          signature(x = "Pop"),
          function (x, ...){
            # Uses mergePops for increased speed
            x = mergePops(c(list(x),list(...)))
            return(x)
          }
)

#' @describeIn Pop Show population summary
setMethod("show",
          signature(object = "Pop"),
          function (object){
            cat("An object of class", 
                classLabel(class(object)), "\n")
            cat("Ploidy:", object@ploidy,"\n")
            cat("Individuals:", object@nInd,"\n")
            cat("Chromosomes:", object@nChr,"\n")
            cat("Loci:", sum(object@nLoci),"\n")
            cat("Traits:", object@nTraits,"\n")
            invisible()
          }
)

#' @title Create new Population
#' 
#' @description
#' Creates a new \code{\link{Pop-class}} from an object of 
#' \code{\link{MapPop-class}} or \code{\link{RawPop-class}}. 
#' The function is intended for creating initial populations from 
#' 'FOUNDERPOP' created by \code{\link{runMacs}}.
#'
#' @param rawPop an object of \code{\link{MapPop-class}} or 
#' \code{\link{RawPop-class}}
#' @param id optional id for new individuals.
#' @param mother optional id for mothers.
#' @param father optional id for fathers.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments used internally
#'
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' @export
newPop = function(rawPop,id=NULL,mother=NULL,father=NULL,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  args = list(...)
  stopifnot(sapply(simParam$genMap,length)==rawPop@nLoci)
  lastId = simParam$lastId
  iid = (1:rawPop@nInd) + lastId
  lastId = max(iid)
  if(is.null(id)){
    if(class(rawPop)=="NamedMapPop"){
      id = rawPop@id
    }else{
      id = as.character(iid)
    }
  }else{
    id = as.character(id)
    stopifnot(length(id)==rawPop@nInd)
  }
  if(any(names(args)=="iMother")){
    iMother = args$iMother
  }else{
    iMother = rep(0L, rawPop@nInd)
  }
  if(any(names(args)=="iFather")){
    iFather = args$iFather
  }else{
    iFather = rep(0L, rawPop@nInd)
  }
  if(any(names(args)=="isDH")){
    isDH = args$isDH
  }else{
    if(class(rawPop)=="MapPop" | class(rawPop)=="NamedMapPop"){
      isDH = rawPop@inbred
    }else{
      isDH = FALSE
    }
  }
  if(is.null(mother)){
    if(class(rawPop)=="NamedMapPop"){
      mother = rawPop@mother
    }else{
      mother = rep("0", rawPop@nInd)
    }
  }else{
    mother = as.character(mother)
  }
  if(is.null(father)){
    if(class(rawPop)=="NamedMapPop"){
      father = rawPop@father
    }else{
      father = rep("0", rawPop@nInd)
    }
  }else{
    father = as.character(father)
  }
  stopifnot(length(id)==length(mother),
            length(id)==length(father))
  if(simParam$sexes=="no"){
    sex = rep("H", rawPop@nInd)
  }else if(simParam$sexes=="yes_rand"){
    sex = sample(c("M","F"), rawPop@nInd, replace=TRUE)
  }else if(simParam$sexes=="yes_sys"){
    sex = rep_len(c("M","F"), rawPop@nInd)
  }else{
    stop(paste("no rules for sex type", simParam$sexes))
  }
  gxe = vector("list",simParam$nTraits)
  gv = matrix(NA_real_,nrow=rawPop@nInd,
              ncol=simParam$nTraits)
  if(simParam$nTraits>=1){
    for(i in 1:simParam$nTraits){
      tmp = getGv(simParam$traits[[i]], rawPop, simParam$nThreads)
      gv[,i] = tmp[[1]]
      if(length(tmp)>1){
        gxe[[i]] = tmp[[2]]
      }
    }
  }
  output = new("Pop",
               nInd=rawPop@nInd,
               nChr=rawPop@nChr,
               ploidy=rawPop@ploidy,
               nLoci=rawPop@nLoci,
               sex=sex,
               geno=rawPop@geno,
               id=id,
               iid=iid,
               mother=mother,
               father=father,
               fixEff=rep(1L,rawPop@nInd),
               reps=rep(1,rawPop@nInd),
               nTraits=simParam$nTraits,
               gv=gv,
               gxe=gxe,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam$nTraits),
               ebv=matrix(NA_real_,
                          nrow=rawPop@nInd,
                          ncol=0),
               misc=vector("list",rawPop@nInd))
  if(simParam$nTraits>=1){
    output = setPheno(output, varE=NULL, reps=1, 
                      fixEff=1L, p=NULL, onlyPheno=FALSE, 
                      simParam=simParam)
  }
  output = simParam$finalizePop(output,...)
  if(simParam$isTrackPed){
    if(simParam$isTrackRec){
      if(any(names(args)=="hist")){
        hist = args$hist
      }else{
        hist = NULL
      }
      simParam$addToRec(lastId,id,iMother,iFather,isDH,hist,output@ploidy)
    }else{
      simParam$addToPed(lastId,id,iMother,iFather,isDH)
    }
  }else{
    simParam$updateLastId(lastId)
  }
  return(output)
}

#' @title Reset population
#' 
#' @description
#' Recalculates a population's genetic values and 
#' resets phenotypes and EBVs.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Rescale to set mean to 1
#' SP$rescaleTraits(mean=1)
#' pop = resetPop(pop, simParam=SP)
#' 
#' @export
resetPop = function(pop,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  pop@nTraits = simParam$nTraits
  pop@pheno = matrix(NA_real_,
                     nrow=pop@nInd,
                     ncol=simParam$nTraits)
  pop@ebv = matrix(NA_real_,
                   nrow=pop@nInd,
                   ncol=0)
  pop@gxe = vector("list",simParam$nTraits)
  pop@gv = matrix(NA_real_,nrow=pop@nInd,
                  ncol=simParam$nTraits)
  pop@fixEff = rep(1L,pop@nInd)
  pop@reps = rep(1,pop@nInd) 
  if(simParam$nTraits>=1){
    for(i in 1:simParam$nTraits){
      tmp = getGv(simParam$traits[[i]],pop,simParam$nThreads)
      pop@gv[,i] = tmp[[1]]
      if(length(tmp)>1){
        pop@gxe[[i]] = tmp[[2]]
      }
    }
  }
  return(pop)
}

# MegaPop ------------------------------------------------------------------

#' @title Mega-Population
#' 
#' @description 
#' The mega-population represents a population of populations. 
#' It is designed to behave like a list of populations. 
#' 
#' @param x a 'MegaPop' object
#' @param i index of populations or mega-populations
#' @param ... additional 'MegaPop' or 'Pop' objects
#' 
#' @slot pops list of \code{\link{Pop-class}} and/or 
#' \code{MegaPop-class}
#' 
#' 
#' @export
setClass("MegaPop",
         slots=c(pops="list"))

setValidity("MegaPop",function(object){
  errors = character()
    # Check that all populations are valid
    for(i in 1:length(object@pops)){
      if(!validObject(object@pops[[i]]) & 
         (class(object@pops[[i]])=="Pop" | 
                class(object@pops[[i]])=="MegaPop")){
        errors = c(errors,paste("object",i,"is not a valid pop"))
      }
    }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#' @describeIn MegaPop Extract MegaPop by index
setMethod("[",
          signature(x = "MegaPop"),
          function(x, i){
            x@pops = x@pops[i]
            return(x)
          }
)

#' @describeIn MegaPop Extract Pop by index
setMethod("[[",
          signature(x = "MegaPop"),
          function (x, i){
            return(x@pops[[i]])
          }
)

#' @describeIn MegaPop Combine multiple MegaPops
setMethod("c",
          signature(x = "MegaPop"),
          function (x, ...){
            for(y in list(...)){
              if(class(y)=="NULL"){
                # Do nothing
              }else{
                if(class(y)=="Pop"){
                  x@pops = c(x@pops, y)
                }else{
                  stopifnot(class(y)=="MegaPop")
                  x@pops = c(x@pops, y@pops)
                }
              }
            }
            return(x)
          }
)

#' @title Create new Mega Population
#' 
#' @description
#' Creates a new \code{\link{MegaPop-class}} from one or more
#' \code{\link{Pop-class}} and/or \code{\link{MegaPop-class}} 
#' objects.  
#'
#' @param ... one or more \code{\link{Pop-class}} and/or 
#' \code{\link{MegaPop-class}} objects.
#'
#' @return Returns an object of \code{\link{MegaPop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' megaPop = newMegaPop(pop=pop)
#' 
#' @export
newMegaPop = function(...){
  input = list(...)
  class = sapply(input, "class")
  stopifnot(all(class=="Pop" | class=="MegaPop"))
  output = new("MegaPop", pops=input)
  return(output)
}

