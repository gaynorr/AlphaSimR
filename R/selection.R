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
      response = genParam(pop,TRUE,simParam=simParam)$bv
      response = trait(response,...)
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }else{
    stopifnot(length(trait)==1)
    if(use == "gv"){
      response = pop@gv[,trait]
    }else if(use=="ebv"){
      response = pop@ebv[,trait]
    }else if(use=="pheno"){
      response = pop@pheno[,trait]
    }else if(use=="bv"){
      if(class(pop)=="HybridPop"){
        stop("Use='bv' is not a valid option for HybridPop")
      }
      response = genParam(pop,TRUE,simParam=simParam)$bv[,trait]
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }
  if(any(is.na(response))){
    stop("selection trait has missing values, phenotype may need to be set")
  }
  return(response)
}

# Returns a vector of individuals in a population with the required gender
checkGender = function(pop,gender,simParam){
  gender = toupper(gender)
  eligible = 1:pop@nInd
  if(simParam$gender=="no"){
    return(eligible)
  }else{
    if(gender=="B"){
      return(eligible)
    }else{
      return(eligible[pop@gender%in%gender])
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
#' @param gender which gender to select. Use "B" for both, "F" for 
#' females and "M" for males. If the simulation is not using gender, 
#' the argument is ignored.
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectInd = function(pop,nInd,trait=1,use="pheno",gender="B",
                     selectTop=TRUE,returnPop=TRUE,
                     simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  eligible = checkGender(pop=pop,gender=gender,simParam=simParam)
  if(length(eligible)<nInd){
    stop("Not enough suitable candidates, check request value and gender")
  }
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,...)
  take = order(response,decreasing=selectTop)
  take = take[take%in%eligible]
  if(returnPop){
    return(pop[take[1:nInd]])
  }else{
    return(take[1:nInd])
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
#' @param gender which gender to select. Use "B" for both, "F" for 
#' females and "M" for males. If the simulation is not using gender, 
#' the argument is ignored.
#' @param famType which type of family to select. Use "B" for 
#' full-sib families, "F" for half-sib families on female side and "M" 
#' for half-sib families on the male side.
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectFam = function(pop,nFam,trait=1,use="pheno",gender="B",
                     famType="B",selectTop=TRUE,returnPop=TRUE,
                     simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  eligible = checkGender(pop=pop,gender=gender,simParam=simParam)
  allFam = getFam(pop=pop,famType=famType)
  availFam = allFam[eligible]
  if(nFam>length(unique(availFam))){
    stop(paste(nFam,"families requested but only",length(unique(availFam)),
               "families are available"))
  }
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,...)[eligible]
  #Calculate family means
  famMeans = aggregate(response,list(families=availFam),mean)
  response = famMeans$x
  #Select families
  bestFam = order(response,decreasing=selectTop)[1:nFam]
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
#' @param gender which gender to select. Use "B" for both, "F" for 
#' females and "M" for males. If the simulation is not using gender, 
#' the argument is ignored.
#' @param famType which type of family to select. Use "B" for 
#' full-sib families, "F" for half-sib families on female side and "M" 
#' for half-sib families on the male side.
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectWithinFam = function(pop,nInd,trait=1,use="pheno",gender="B",
                           famType="B",selectTop=TRUE,returnPop=TRUE,
                           simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  eligible = checkGender(pop=pop,gender=gender,simParam=simParam)
  families = getFam(pop=pop,famType=famType)
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,...)
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
      return(index[1:nInd])
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
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
selectOP = function(pop,nInd,nSeeds,probSelf=0,
                    pollenControl=FALSE,trait=1,
                    use="pheno",selectTop=TRUE,
                    simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  female = selectInd(pop=pop,nInd=nInd,trait=trait,
                     use=use,gender="B",selectTop=selectTop,
                     returnPop=FALSE,simParam=simParam,...)
  nSelf = rbinom(n=nInd,prob=probSelf,size=nSeeds)
  if(pollenControl){
    male = female
  }else{
    male = 1:nInd
  }
  crossPlan = lapply(1:nInd,function(x){
    cbind(rep(female[x],nSeeds),
          c(rep(female[x],nSelf[x]),
            sample(male[!male==female[x]],nSeeds-nSelf[x],replace=TRUE)))
  })
  crossPlan = mergeMultIntMat(crossPlan,rep(nSeeds,nInd),2L)
  return(makeCross(pop=pop,crossPlan=crossPlan,
                   rawPop=FALSE,simParam=simParam))
}
