#' Returns a vector response from a population
#'
#' @param pop an object of class Pop or HybirdPop
#' @param trait a vector or custom function
#' @param use a character ("rand", "gv", "ebv", "pheno", or "bv"; 
#' note that "bv" doesn't work on class HybridPop)
#' @param simParam simulation parameters are only used when use="bv"
#' @param ... are additional arguments passed to trait when trait is a function
#'
#' @keywords internal
getResponse = function(pop,trait,use,simParam=NULL,...){
  if(is(trait,"function")){
    if(is.character(use)){
      use = tolower(use)
      if(use=="rand"){
        return(rnorm(pop@nInd))
      }else if(use=="gv"){
        response = trait(pop@gv,...)
      }else if(use=="ebv"){
        response = trait(pop@ebv,...)
      }else if(use=="pheno"){
        response = trait(pop@pheno,...)
      }else if(use=="bv"){
        if(is(pop,"HybridPop")){
          stop("Use='bv' is not a valid option for HybridPop")
        }
        response = genParam(pop,simParam=simParam)$bv
        response = trait(response,...)
      }else{
        stop(paste0("Use=",use," is not an option"))
      }
    }else if(is(use,"function")){
      response = trait(use(pop, ...), ...)
    }else{
      stop("use must be a character or a function")
    }
  }else{ # trait is not a function, so must be numeric or character
    if(is.character(trait)){ # Suspect trait is a name
      take = match(trait, simParam$traitNames)
      if(any(is.na(take))){
        stop("'",trait[is.na(take)],"' did not match any trait names")
      }
      trait = take
    }
    if(is.character(use)){
      use = tolower(use)
      if(use=="rand"){
        return(rnorm(pop@nInd))
      }else if(use == "gv"){
        response = pop@gv[,trait,drop=FALSE]
      }else if(use=="ebv"){
        response = pop@ebv[,trait,drop=FALSE]
      }else if(use=="pheno"){
        response = pop@pheno[,trait,drop=FALSE]
      }else if(use=="bv"){
        if(is(pop,"HybridPop")){
          stop("Use='bv' is not a valid option for HybridPop")
        }
        response = genParam(pop,simParam=simParam)$bv[,trait,drop=FALSE]
      }else{
        stop(paste0("Use=",use," is not an option"))
      }
    }else if(is(use,"function")){
      response = use(pop, trait=trait, ...)
    }else{
      stop("use must be a character or a function")
    }
  }
  if(any(is.na(response))){
    stop("selection trait has missing values, phenotype may need to be set")
  }
  return(response)
}

#' Identify candidate individuals
#' 
#' This function handles indexing by id and negative value indexing
#'
#' @param pop a population object
#' @param candidates a vector of candidates (numeric or character)
#'
#' @keywords internal
getCandidates = function(pop, candidates){
  if(is.character(candidates)){
    candidates = match(candidates, pop@id)
    if(any(is.na(candidates))){
      stop("Trying to select invalid individuals")
    }
    if(any(is.null(candidates))){
      stop("Not valid ids")
    }
  }else{
    if(any(abs(candidates)>pop@nInd)){
      stop("Trying to select invalid individuals")
    }
    candidates = (1:pop@nInd)[candidates]
  }
  return(candidates)
}

#' Find individuals of desired sex
#' 
#' This function ignores sex if sexes is "no" in SimParam.
#'
#' @param pop a population
#' @param sex the desired sex (M or F)
#' @param simParam simulation parameters object
#' @param ...  captures use of old gender argument
#'
#' @keywords internal
checkSexes = function(pop,sex,simParam,...){
  sex = toupper(sex)
  eligible = 1:pop@nInd
  if(simParam$sexes=="no"){
    return(eligible)
  }else{
    # Check in gender is incorrectly being used
    args = list(...)
    if(any(names(args)=="gender")){
      stop("The discontinued 'gender' argument appears to be in use. This argument was renamed as 'sex' in AlphaSimR version 0.13.0.")
    }
    if(sex=="B"){
      return(eligible)
    }else if(sex=="F" | sex=="M"){
      return(eligible[pop@sex%in%sex])
    }else{
      stop("Invalid option supplied for 'sex'")
    }
  }
}

#' Determine families
#'
#' @param pop a population object
#' @param famType type of family (fullsib "B", maternal "M", or paternal "F")
#'
#' @keywords internal
getFam = function(pop,famType){
  famType = toupper(famType)
  if(famType=="B"){
    return(paste(pop@mother,pop@father,sep="_"))
  }else if(famType=="F"){
    return(pop@mother)
  }else if(famType=="M"){
    return(pop@father)
  }else{
    stop(paste0("famType=",famType," is not a valid option"))
  }
}

#' @title Select individuals
#'
#' @description Selects a subset of nInd individuals from a
#' population.
#'
#' @param pop and object of \code{\link{Pop-class}},
#'   \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#' @param nInd the number of individuals to select
#' @param trait the trait for selection. Either a number indicating
#'   a single trait or a function returning a vector of length nInd.
#'   The function must work on a vector or matrix of \code{use} values as
#'   \code{trait(pop@use, ...)} - depending on what \code{use} is.
#'   See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values "gv", estimated breeding values "ebv", breeding values "bv",
#'   phenotypes "pheno", or randomly "rand") or
#'   a function returning a vector of length nInd.
#'   The function must work on \code{pop} as \code{use(pop, trait, ...)} or
#'   as \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param sex which sex to select. Use "B" for both, "F" for
#'   females and "M" for males. If the simulation is not using sexes,
#'   the argument is ignored.
#' @param selectTop selects highest values if true.
#'   Selects lowest values if false.
#' @param returnPop should results be returned as a
#'   \code{\link{Pop-class}}. If FALSE, only the index of selected
#'   individuals is returned.
#' @param candidates an optional vector of eligible selection candidates.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for
#'   \code{trait} or \code{use}
#'
#' @return Returns an object of \code{\link{Pop-class}},
#' \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Select top 5 (directional selection)
#' pop2 = selectInd(pop, 5, simParam=SP)
#' hist(pop@pheno); abline(v=pop@pheno, lwd=2)
#' abline(v=pop2@pheno, col="red", lwd=2)
#'
#' #Select 5 most deviating from an optima (disruptive selection)
#' squaredDeviation = function(x, optima=0) (x - optima)^2
#' pop3 = selectInd(pop, 5, trait=squaredDeviation, selectTop=TRUE, simParam=SP)
#' hist(pop@pheno); abline(v=pop@pheno, lwd=2)
#' abline(v=pop3@pheno, col="red", lwd=2)
#'
#' #Select 5 least deviating from an optima (stabilising selection)
#' pop4 = selectInd(pop, 5, trait=squaredDeviation, selectTop=FALSE, simParam=SP)
#' hist(pop@pheno); abline(v=pop@pheno, lwd=2)
#' abline(v=pop4@pheno, col="red", lwd=2)
#'
#' #Select 5 individuals based on miscelaneous information with use function
#' pop@misc = list(smth=rnorm(10), smth2=rnorm(10))
#' useFunc = function(pop, trait=NULL) pop@misc$smth + pop@misc$smth2
#' pop5 = selectInd(pop, 5, use=useFunc, simParam=SP)
#' pop5@id
#'
#' #... equivalent result with the use & trait function
#' useFunc2 = function(pop, trait=NULL) cbind(pop@misc$smth, pop@misc$smth2)
#' trtFunc = function(x) rowSums(x)
#' pop6 = selectInd(pop, 5, trait=trtFunc, use=useFunc2, simParam=SP)
#' pop6@id
#'
#' @export
selectInd = function(pop,nInd,trait=1,use="pheno",sex="B",
                     selectTop=TRUE,returnPop=TRUE,
                     candidates=NULL,simParam=NULL,...){
  stopifnot(nInd>=0)
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is(pop,"MultiPop")){
    stopifnot(returnPop, is.null(candidates))
    pop@pops = lapply(pop@pops, selectInd, nInd=nInd, trait=trait,
                      use=use, sex=sex, selectTop=selectTop,
                      returnPop=TRUE, candidates=NULL,
                      simParam=simParam, ...)
    return(pop)
  }
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
    candidates = getCandidates(pop=pop,candidates=candidates)
    eligible = eligible[eligible%in%candidates]
  }
  if(length(eligible)<nInd){
    nInd = length(eligible)
    warning("Suitable candidates smaller than nInd, returning ",nInd," individuals")
  }
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,...)
  if(is.matrix(response)){
    stopifnot(ncol(response)==1)
  }
  take = order(response,decreasing=selectTop)
  take = take[take%in%eligible]
  if(returnPop){
    return(pop[take[0:nInd]])
  }else{
    return(take[0:nInd])
  }
}

#' @title Select families
#'
#' @description Selects a subset of full-sib families from a
#' population.
#'
#' @param pop and object of \code{\link{Pop-class}},
#'   \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#' @param nFam the number of families to select
#' @param trait the trait for selection. Either a number indicating
#'   a single trait or a function returning a vector of length nInd.
#'   The function must work on a vector or matrix of \code{use} values as
#'   \code{trait(pop@use, ...)} - depending on what \code{use} is.
#'   See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values "gv", estimated breeding values "ebv", breeding values "bv",
#'   phenotypes "pheno", or randomly "rand") or
#'   a function returning a vector of length nInd.
#'   The function must work on \code{pop} as \code{use(pop, trait, ...)} or
#'   as \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param sex which sex to select. Use "B" for both, "F" for
#'   females and "M" for males. If the simulation is not using sexes,
#'   the argument is ignored.
#' @param famType which type of family to select. Use "B" for
#'   full-sib families, "F" for half-sib families on female side and "M"
#'   for half-sib families on the male side.
#' @param selectTop selects highest values if true.
#'   Selects lowest values if false.
#' @param returnPop should results be returned as a
#'   \code{\link{Pop-class}}. If FALSE, only the index of selected
#'   individuals is returned.
#' @param candidates an optional vector of eligible selection candidates.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for
#'   \code{trait} and \code{use}
#'
#' @return Returns an object of \code{\link{Pop-class}},
#' \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Create 3 biparental families with 10 progeny
#' pop2 = randCross(pop, nCrosses=3, nProgeny=10, simParam=SP)
#'
#' #Select best 2 families
#' pop3 = selectFam(pop2, 2, simParam=SP)
#'
#' @export
selectFam = function(pop,nFam,trait=1,use="pheno",sex="B",
                     famType="B",selectTop=TRUE,returnPop=TRUE,
                     candidates=NULL,simParam=NULL,...){
  stopifnot(nFam>=0)
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is(pop,"MultiPop")){
    stopifnot(returnPop, is.null(candidates))
    pop@pops = lapply(pop@pops, selectFam, nFam=nFam, trait=trait,
                      use=use, sex=sex, famType=famType,
                      selectTop=selectTop, returnPop=TRUE,
                      candidates=NULL, simParam=simParam, ...)
    return(pop)
  }
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
    candidates = getCandidates(pop=pop,candidates=candidates)
    eligible = eligible[eligible%in%candidates]
  }
  allFam = getFam(pop=pop,famType=famType)
  availFam = allFam[eligible]
  if(nFam>length(unique(availFam))){
    nFam = length(unique(availFam))
    warning("Suitable families smaller than nFam, returning ", nFam, " families")
  }
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,...)
  if(is.matrix(response)){
    stopifnot(ncol(response)==1)
  }
  response = response[eligible]
  #Calculate family means
  famMeans = aggregate(response,list(families=availFam),mean)
  response = famMeans$x
  #Select families
  bestFam = order(response,decreasing=selectTop)[0:nFam]
  bestFam = famMeans$families[bestFam]
  take = which(allFam%in%bestFam)
  take = take[take%in%eligible]
  if(returnPop){
    return(pop[take])
  }else{
    return(take)
  }
}

#' @title Select individuals within families
#'
#' @description Selects a subset of nInd individuals from each
#' full-sib family within a population. Will return all individuals
#' from a full-sib family if it has less than or equal to nInd individuals.
#'
#' @param pop and object of \code{\link{Pop-class}},
#'   \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#' @param nInd the number of individuals to select within a family
#' @param trait the trait for selection. Either a number indicating
#'   a single trait or a function returning a vector of length nInd.
#'   The function must work on a vector or matrix of \code{use} values as
#'   \code{trait(pop@use, ...)} - depending on what \code{use} is.
#'   See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values "gv", estimated breeding values "ebv", breeding values "bv",
#'   phenotypes "pheno", or randomly "rand") or
#'   a function returning a vector of length nInd.
#'   The function must work on \code{pop} as \code{use(pop, trait, ...)} or
#'   as \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param sex which sex to select. Use "B" for both, "F" for
#'   females and "M" for males. If the simulation is not using sexes,
#'   the argument is ignored.
#' @param famType which type of family to select. Use "B" for
#'   full-sib families, "F" for half-sib families on female side and "M"
#'   for half-sib families on the male side.
#' @param selectTop selects highest values if true.
#'   Selects lowest values if false.
#' @param returnPop should results be returned as a
#'   \code{\link{Pop-class}}. If FALSE, only the index of selected
#'   individuals is returned.
#' @param candidates an optional vector of eligible selection candidates.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for
#'   \code{trait} and \code{use}
#'
#' @return Returns an object of \code{\link{Pop-class}},
#' \code{\link{HybridPop-class}} or \code{\link{MultiPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Create 3 biparental families with 10 progeny
#' pop2 = randCross(pop, nCrosses=3, nProgeny=10, simParam=SP)
#'
#' #Select best individual per family
#' pop3 = selectWithinFam(pop2, 1, simParam=SP)
#'
#' @export
selectWithinFam = function(pop,nInd,trait=1,use="pheno",sex="B",
                           famType="B",selectTop=TRUE,returnPop=TRUE,
                           candidates=NULL,simParam=NULL,...){
  stopifnot(nInd>=0)
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is(pop,"MultiPop")){
    stopifnot(returnPop, is.null(candidates))
    pop@pops = lapply(pop@pops, selectWithinFam, nInd=nInd, trait=trait,
                      use=use, sex=sex, selectTop=selectTop,
                      returnPop=TRUE, candidates=NULL,
                      simParam=simParam, ...)
    return(pop)
  }
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
    candidates = getCandidates(pop=pop,candidates=candidates)
    eligible = eligible[eligible%in%candidates]
  }
  families = getFam(pop=pop,famType=famType)
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,...)
  if(is.matrix(response)){
    stopifnot(ncol(response)==1)
  }
  warn = FALSE
  selInFam = function(selFam){
    index = which(families%in%selFam)
    y = response[index]
    index = index[order(y,decreasing=selectTop)]
    index = index[index%in%eligible]
    if(length(index)<nInd){
      warn <<- TRUE
      return(index)
    }else{
      return(index[0:nInd])
    }
  }
  take = unlist(sapply(unique(families),selInFam))
  if(warn){
    warning("One or more families are smaller than nInd")
  }
  if(returnPop){
    return(pop[take])
  }else{
    return(take)
  }
}

#' @title Select open pollinating plants
#'
#' @description
#' This function models selection in an open pollinating
#' plant population. It allows for varying the percentage of
#' selfing. The function also provides an option for modeling
#' selection as occuring before or after pollination.
#'
#' @param pop and object of \code{\link{Pop-class}}
#'   or \code{\link{MultiPop-class}}
#' @param nInd the number of plants to select
#' @param nSeeds number of seeds per plant
#' @param probSelf percentage of seeds expected from selfing.
#'   Value ranges from 0 to 1.
#' @param pollenControl are plants selected before pollination
#' @param trait the trait for selection. Either a number indicating
#'   a single trait or a function returning a vector of length nInd.
#'   The function must work on a vector or matrix of \code{use} values as
#'   \code{trait(pop@use, ...)} - depending on what \code{use} is.
#'   See the examples and \code{\link{selIndex}}.
#' @param use the selection criterion. Either a character
#'   (genetic values "gv", estimated breeding values "ebv", breeding values "bv",
#'   phenotypes "pheno", or randomly "rand") or
#'   a function returning a vector of length nInd.
#'   The function must work on \code{pop} as \code{use(pop, trait, ...)} or
#'   as \code{trait(pop@use, ...)} depending on what \code{trait} is.
#'   See the examples.
#' @param selectTop selects highest values if true.
#'   Selects lowest values if false.
#' @param candidates an optional vector of eligible selection candidates.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for
#'   \code{trait} and \code{use}
#'
#' @return Returns an object of \code{\link{Pop-class}}
#' or \code{\link{MultiPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Create new population by selecting the best 3 plant
#' #Assuming 50% selfing in plants and 10 seeds per plant
#' pop2 = selectOP(pop, nInd=3, nSeeds=10, probSelf=0.5, simParam=SP)
#'
#' @export
selectOP = function(pop,nInd,nSeeds,probSelf=0,
                    pollenControl=FALSE,trait=1,
                    use="pheno",selectTop=TRUE,
                    candidates=NULL,simParam=NULL,...){
  stopifnot(nInd>=0)
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is(pop,"MultiPop")){
    stopifnot(is.null(candidates))
    pop@pops = lapply(pop@pops, selectOP, nInd=nInd, nSeeds=nSeeds,
                      pollenControl=pollenControl, trait=trait, use=use,
                      selectTop=selectTop, candidates=NULL,
                      simParam=simParam, ...)
    return(pop)
  }
  female = selectInd(pop=pop,nInd=nInd,trait=trait,
                     use=use,sex="B",selectTop=selectTop,
                     returnPop=FALSE,candidates=candidates,
                     simParam=simParam,...)
  nSelf = rbinom(n=nInd,prob=probSelf,size=nSeeds)
  if(pollenControl){
    male = female
  }else{
    male = 1:pop@nInd
  }
  crossPlan = lapply(1:nInd,function(x){
    male = male[!male==female[x]]
    if(length(male)==1){
      #Account for "convenience" feature of sample when length = 1
      cbind(rep(female[x],nSeeds),
            c(rep(female[x],nSelf[x]),
              rep(male,nSeeds-nSelf[x])))
    }else{
      cbind(rep(female[x],nSeeds),
            c(rep(female[x],nSelf[x]),
              sample(male,nSeeds-nSelf[x],replace=TRUE)))
    }
  })
  crossPlan = mergeMultIntMat(crossPlan,rep(nSeeds,nInd),2L)
  return(makeCross(pop=pop,crossPlan=crossPlan,simParam=simParam))
}
