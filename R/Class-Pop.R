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
#' @slot geno list of nChr length containing chromosome genotypes.
#' Each element is a three dimensional array of raw values.
#' The array dimensions are nLoci by ploidy by nInd.
#'
#' @export
setClass("RawPop",
         slots=c(nInd="integer",
                 nChr="integer",
                 ploidy="integer",
                 nLoci="integer",
                 geno="list"))

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
              if(is(y,"NULL")){
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

#' @describeIn RawPop Test if object is of a RawPop class
#' @export
isRawPop = function(x) {
  ret = is(x, class2 = "RawPop")
  return(ret)
}

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
#' @slot genMap list of chromosome genetic maps
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
              if(is(y,"NULL")){
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

#' @describeIn MapPop Test if object is of a MapPop class
#' @export
isMapPop = function(x) {
  ret = is(x, class2 = "MapPop")
  return(ret)
}

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
              if(is(y,"NULL")){
                # Do nothing
              }else{
                stopifnot(is(y,"NamedMapPop"),
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
    if(is(y,"NULL")){
      #Do nothing
    }else{
      stopifnot(is(y,"MapPop"))
      if(!exists("x",inherits=FALSE)){
        x = y
      }else{
        stopifnot(x@nInd==y@nInd,
                  x@ploidy==y@ploidy)
        x@nChr = x@nChr+y@nChr
        x@geno = c(x@geno,y@geno)
        x@genMap = c(x@genMap,y@genMap)
        x@centromere = c(x@centromere,y@centromere)
        x@nLoci = c(x@nLoci,y@nLoci)
        x@inbred = x@inbred & y@inbred
      }
    }
  }
  return(x)
}

#' @describeIn NamedMapPop Test if object is a NamedMapPop class
#' @export
isNamedMapPop = function(x) {
  ret = is(x, class2 = "NamedMapPop")
  return(ret)
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
#' @slot misc a list whose elements correspond to additional miscellaneous
#' nodes with the items for individuals in the population (see example in
#' \code{\link{newPop}}).
#' This list is normally empty and exists solely as an
#' open slot available for uses to store extra information about
#' individuals.
#' @slot miscPop a list of any length containing optional meta data for the
#' population (see example in \code{\link{newPop}}).
#' This list is empty unless information is supplied by the user.
#' Note that the list is emptied every time the population is subsetted or
#' combined because the meta data for old population might not be valid anymore.
#'
#' @seealso \code{\link{newPop}}, \code{\link{newEmptyPop}}, \code{\link{resetPop}}
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
                 misc="list",
                 miscPop="list"),
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
  if(any(object@nInd!=sapply(object@misc, length))){
    errors = c(errors,"any(nInd!=sapply(misc, length))")
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
            x@misc = lapply(x@misc, FUN = function(z) z[i])
            x@miscPop = list()
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

#' @title Create new population
#'
#' @description
#' Creates an initial \code{\link{Pop-class}} from an object of
#' \code{\link{MapPop-class}} or \code{\link{NamedMapPop-class}}.
#' The function is intended for us with output from functions such
#' as \code{\link{runMacs}}, \code{\link{newMapPop}}, or
#' \code{\link{quickHaplo}}.
#'
#' @param rawPop an object of \code{\link{MapPop-class}} or
#' \code{\link{NamedMapPop-class}}
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
#' isPop(pop)
#'
#' #Misc
#' pop@misc$tmp1 = rnorm(n=2)
#' pop@misc$tmp2 = rnorm(n=2)
#'
#' #MiscPop
#' pop@miscPop$tmp1 = sum(pop@misc$tmp1)
#' pop@miscPop$tmp2 = sum(pop@misc$tmp2)
#' @export
newPop = function(rawPop,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  return(.newPop(rawPop=rawPop,simParam=simParam,...))
}

#' @title Create new population (internal)
#'
#' @description
#' Creates a new \code{\link{Pop-class}} from an object of
#' of the Pop superclass.
#'
#' @param rawPop an object of the pop superclass
#' @param id optional id for new individuals
#' @param mother optional id for mothers
#' @param father optional id for fathers
#' @param iMother optional internal id for mothers
#' @param iFather optional internal id for fathers
#' @param isDH optional indicator for DH/inbred individuals
#' @param femaleParentPop optional population of female parents
#' @param maleParentPop optional population of male parents
#' @param hist optional recombination history
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments passed to the finalizePop
#' function in simParam
#'
#' @return Returns an object of \code{\link{Pop-class}}
.newPop = function(rawPop, id=NULL, mother=NULL, father=NULL,
                   iMother=NULL, iFather=NULL, isDH=NULL,
                   femaleParentPop=NULL, maleParentPop=NULL,
                   hist=NULL, simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }

  stopifnot(sapply(simParam$genMap,length)==rawPop@nLoci)

  lastId = simParam$lastId
  iid = (1:rawPop@nInd) + lastId
  lastId = max(iid)

  if(is.null(id)){
    if(is(rawPop, "NamedMapPop")){
      id = rawPop@id
    }else{
      id = as.character(iid)
    }
  }

  if(is.null(iMother)){
    iMother = rep(0L, rawPop@nInd)
  }

  if(is.null(iFather)){
    iFather = rep(0L, rawPop@nInd)
  }

  if(is.null(isDH)){
    if(is(rawPop, "MapPop")){
      isDH = rawPop@inbred
    }else{
      isDH = FALSE
    }
  }

  if(is.null(mother)){
    if(is(rawPop, "NamedMapPop")){
      mother = rawPop@mother
    }else{
      mother = rep("0", rawPop@nInd)
    }
  }

  if(is.null(father)){
    if(is(rawPop, "NamedMapPop")){
      father = rawPop@father
    }else{
      father = rep("0", rawPop@nInd)
    }
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

  gxe = vector("list", simParam$nTraits)

  gv = matrix(NA_real_,nrow=rawPop@nInd,
              ncol=simParam$nTraits)
  colnames(gv) = rep(NA_character_, simParam$nTraits)
  pheno = gv

  if(simParam$nTraits>=1){
    for(i in 1:simParam$nTraits){
      tmp = getGv(simParam$traits[[i]], rawPop, simParam$nThreads)
      gv[,i] = tmp[[1]]

      colnames(gv)[i] = simParam$traits[[i]]@name

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
               misc=list(),
               miscPop=list(),
               nTraits=simParam$nTraits,
               gv=gv,
               gxe=gxe,
               pheno=pheno,
               ebv=matrix(NA_real_,
                          nrow=rawPop@nInd,
                          ncol=0))
  if(simParam$nTraits>=1){
    output = setPheno(output, varE=NULL, reps=1,
                      fixEff=1L, p=NULL, onlyPheno=FALSE,
                      simParam=simParam)
  }

  output = simParam$finalizePop(output,...)

  if(simParam$isTrackPed){
    if(simParam$isTrackRec){
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

  # Extract names to add back at the end
  traitNames = colnames(pop@gv)

  # Create empty slots for traits
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

  # Calculate genetic values
  if(simParam$nTraits>=1){
    for(i in 1:simParam$nTraits){
      tmp = getGv(simParam$traits[[i]],pop,simParam$nThreads)
      pop@gv[,i] = tmp[[1]]
      if(length(tmp)>1){
        pop@gxe[[i]] = tmp[[2]]
      }
    }
  }

  # Add back trait names
  colnames(pop@pheno) = colnames(pop@gv) = traitNames

  return(pop)
}


#' @title Test if object is of a Population class
#'
#' @description Utilify function to test if object is of a Population class
#'
#' @param x \code{\link{Pop-class}}
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
#' isPop(pop)
#' isPop(SP)
#'
#' @export
isPop = function(x) {
  ret = is(x, class2 = "Pop")
  return(ret)
}

#' @title Creates an empty population
#'
#' @description
#' Creates an empty \code{\link{Pop-class}} object with user
#' defined ploidy and other parameters taken from simParam.
#'
#' @param ploidy the ploidy of the population
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns an object of \code{\link{Pop-class}} with
#' zero individuals
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#'
#' #Create empty population
#' pop = newEmptyPop(simParam=SP)
#' isPop(pop)
#'
#' @export
newEmptyPop = function(ploidy=2L, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP", envir=.GlobalEnv)
  }

  # Create 0 x nTrait matrix with trait names
  # For pheno and gv slots
  traitMat = matrix(NA_real_,
                    nrow = 0L,
                    ncol = simParam$nTraits)

  traitNames = character(simParam$nTraits)

  if(simParam$nTraits > 0L){
    # Get trait names
    for(i in 1:simParam$nTraits){
      traitNames[i] = simParam$traits[[i]]@name
    }
  }

  colnames(traitMat) = traitNames

  # Create empty geno list
  nLoci = unname(sapply(simParam$genMap, length))
  geno = vector("list", simParam$nChr)
  for(i in 1:simParam$nChr){
    DIM1 = nLoci[i]%/%8L + (nLoci[i]%%8L > 0L)
    geno[[i]] = array(as.raw(0), dim=c(DIM1, ploidy, 0))
  }

  output = new("Pop",
               nInd = 0L,
               nChr = simParam$nChr,
               ploidy = as.integer(ploidy),
               nLoci = nLoci,
               sex = character(),
               geno = geno,
               id = character(),
               iid = integer(),
               mother = character(),
               father = character(),
               fixEff = integer(),
               misc = list(),
               miscPop = list(),
               nTraits = simParam$nTraits,
               gv = traitMat,
               gxe = vector("list", simParam$nTraits),
               pheno = traitMat,
               ebv = matrix(NA_real_,
                            nrow=0L,
                            ncol=0L))
  return(output)
}

# MultiPop ------------------------------------------------------------------

#' @title Multi-Population
#'
#' @description
#' The mega-population represents a population of populations.
#' It is designed to behave like a list of populations.
#'
#' @param x a 'MultiPop' object
#' @param i index of populations or mega-populations
#' @param ... additional 'MultiPop' or 'Pop' objects
#'
#' @slot pops list of \code{\link{Pop-class}} and/or
#' \code{MultiPop-class}
#'
#'
#' @export
setClass("MultiPop",
         slots=c(pops="list"))

setValidity("MultiPop",function(object){
  errors = character()
    # Check that all populations are valid
    for(i in 1:length(object@pops)){
      if(!validObject(object@pops[[i]]) &
         (is(object@pops[[i]], "Pop") |
                is(object@pops[[i]],"MultiPop"))){
        errors = c(errors,paste("object",i,"is not a valid pop"))
      }
    }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#' @describeIn MultiPop Extract MultiPop by index
setMethod("[",
          signature(x = "MultiPop"),
          function(x, i){
            x@pops = x@pops[i]
            return(x)
          }
)

#' @describeIn MultiPop Extract Pop by index
setMethod("[[",
          signature(x = "MultiPop"),
          function (x, i){
            return(x@pops[[i]])
          }
)

#' @describeIn MultiPop Combine multiple MultiPops
setMethod("c",
          signature(x = "MultiPop"),
          function (x, ...){
            for(y in list(...)){
              if(is(y,"NULL")){
                # Do nothing
              }else{
                if(is(y,"Pop")){
                  x@pops = c(x@pops, y)
                }else{
                  stopifnot(is(y,"MultiPop"))
                  x@pops = c(x@pops, y@pops)
                }
              }
            }
            return(x)
          }
)

#' @title Create new Multi Population
#'
#' @description
#' Creates a new \code{\link{MultiPop-class}} from one or more
#' \code{\link{Pop-class}} and/or \code{\link{MultiPop-class}}
#' objects.
#'
#' @param ... one or more \code{\link{Pop-class}} and/or
#' \code{\link{MultiPop-class}} objects.
#'
#' @return Returns an object of \code{\link{MultiPop-class}}
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
#' megaPop = newMultiPop(pop=pop)
#' isMultiPop(megaPop)
#'
#' @export
newMultiPop = function(...){
  input = list(...)
  class = sapply(input, "class")
  stopifnot(all(class=="Pop" | class=="MultiPop"))
  output = new("MultiPop", pops=input)
  return(output)
}

#' @describeIn MultiPop Test if object is of a MultiPop class
#' @export
isMultiPop = function(x) {
  ret = is(x, class2 = "MultiPop")
  return(ret)
}
