# Returns a vector response from a population
# pop is an object of class Pop or HybridPop
# trait is a vector of traits or a function
# use is "rand", "gv", "ebv", "pheno", or "bv"
#  "bv" doesn't work on class HybridPop
# simParam is an object of class SimParam, it is only called when use="bv"
# ... are additional arguments passed to trait when trait is a function
getResponse = function(pop,trait,use,simParam=NULL,...){
  use = tolower(use)
  if(use=="rand"){
    return(rnorm(pop@nInd))
  }
  if(class(trait)=="function"){
    if(use=="gv"){
      response = trait(pop@gv,...)
    }else if(use=="ebv"){
      response = trait(pop@ebv,...)
    }else if(use=="pheno"){
      response = trait(pop@pheno,...)
    }else if(use=="bv"){
      if(class(pop)=="HybridPop"){
        stop("Use='bv' is not a valid option for HybridPop")
      }
      response = genParam(pop,simParam=simParam)$bv
      response = trait(response,...)
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }else{
    if(use == "gv"){
      response = pop@gv[,trait,drop=FALSE]
    }else if(use=="ebv"){
      response = pop@ebv[,trait,drop=FALSE]
    }else if(use=="pheno"){
      response = pop@pheno[,trait,drop=FALSE]
    }else if(use=="bv"){
      if(class(pop)=="HybridPop"){
        stop("Use='bv' is not a valid option for HybridPop")
      }
      response = genParam(pop,simParam=simParam)$bv[,trait,drop=FALSE]
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }
  if(any(is.na(response))){
    stop("selection trait has missing values, phenotype may need to be set")
  }
  return(response)
}

# Returns a vector of individuals in a population with the required sex
checkSexes = function(pop,sex,simParam,...){
  sex = toupper(sex)
  eligible = 1:pop@nInd
  if(simParam$sexes=="no"){
    return(eligible)
  }else{
    if(sex=="B"){
      # Check in gender is incorrectly being used
      args = list(...)
      if(any(names(args)=="gender")){
        stop("The discontinued 'gender' argument appears to be in use. This argument was renamed as 'sex' in AlphaSimR version 0.13.0.")
      }
      return(eligible)
    }else{
      return(eligible[pop@sex%in%sex])
    }
  }
}

# Returns a vector of families
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
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param nInd the number of individuals to select
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values "gv", estimated
#' breeding values "ebv", breeding values "bv", phenotypes "pheno", 
#' or randomly "rand"
#' @param sex which sex to select. Use "B" for both, "F" for 
#' females and "M" for males. If the simulation is not using sexes, 
#' the argument is ignored.
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param candidates an optional vector of eligible selection candidates. 
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Select best 5
#' pop2 = selectInd(pop, 5, simParam=SP)
#' 
#' @export
selectInd = function(pop,nInd,trait=1,use="pheno",sex="B",
                     selectTop=TRUE,returnPop=TRUE,
                     candidates=NULL,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
    eligible = eligible[eligible%in%candidates]
  }
  if(length(eligible)<nInd){
    stop("Not enough suitable candidates, check request value and sex")
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
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param nFam the number of families to select
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values "gv", estimated
#' breeding values "ebv", breeding values "bv", phenotypes "pheno", 
#' or randomly "rand"
#' @param sex which sex to select. Use "B" for both, "F" for 
#' females and "M" for males. If the simulation is not using sexes, 
#' the argument is ignored.
#' @param famType which type of family to select. Use "B" for 
#' full-sib families, "F" for half-sib families on female side and "M" 
#' for half-sib families on the male side.
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param candidates an optional vector of eligible selection candidates.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
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
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
    eligible = eligible[eligible%in%candidates]
  }
  allFam = getFam(pop=pop,famType=famType)
  availFam = allFam[eligible]
  if(nFam>length(unique(availFam))){
    stop(paste(nFam,"families requested but only",length(unique(availFam)),
               "families are available"))
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
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param nInd the number of individuals to select within a family
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values "gv", estimated
#' breeding values "ebv", breeding values "bv", phenotypes "pheno", 
#' or randomly "rand"
#' @param sex which sex to select. Use "B" for both, "F" for 
#' females and "M" for males. If the simulation is not using sexes, 
#' the argument is ignored.
#' @param famType which type of family to select. Use "B" for 
#' full-sib families, "F" for half-sib families on female side and "M" 
#' for half-sib families on the male side.
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param candidates an optional vector of eligible selection candidates.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
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
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  eligible = checkSexes(pop=pop,sex=sex,simParam=simParam,...)
  if(!is.null(candidates)){
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
#' @param pop an object of \code{\link{Pop-class}}
#' @param nInd the number of plants to select
#' @param nSeeds number of seeds per plant
#' @param probSelf percentage of seeds expected from selfing. 
#' Value ranges from 0 to 1.
#' @param pollenControl are plants selected before pollination
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values "gv", estimated
#' breeding values "ebv", breeding values "bv", phenotypes "pheno", 
#' or randomly "rand"
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param candidates an optional vector of eligible selection candidates.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
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
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
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
