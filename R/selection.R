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
#' @param use select on genetic values (\code{gv}), estimated
#' breeding values (\code{ebv}), breeding values (\code{bv}), 
#' or phenotypes (\code{pheno}, default)
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectInd = function(pop,nInd,trait=1,use="pheno",selectTop=TRUE,
                     returnPop=TRUE,simParam=NULL,...){
  if(is.null(simParam) & use=="bv"){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  stopifnot(nInd<=pop@nInd)
  use = tolower(use)
  if(class(trait)=="function"){
    if(use == "gv"){
      response = trait(pop@gv,...)
    }else if(use == "ebv"){
      response = trait(pop@ebv,...)
    }else if(use == "pheno"){
      response = trait(pop@pheno,...)
    }else if(use == "bv"){
      response = varAD(pop,retGenParam=TRUE,simParam=simParam)$bv
      response = trait(response,...)
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }else{
    stopifnot(length(trait)==1,trait<=pop@nTraits)
    if(use == "gv"){
      response = pop@gv[,trait]
    }else if(use == "ebv"){
      response = pop@ebv[,trait]
    }else if(use == "pheno"){
      response = pop@pheno[,trait]
    }else if(use == "bv"){
      response = varAD(pop,retGenParam=TRUE,simParam=simParam)$bv[,trait]
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }
  if(any(is.na(response))){
    stop("selection trait has missing values, phenotype may need to be set")
  }
  take = order(response,decreasing=selectTop)
  if(returnPop){
    return(pop[take[1:nInd]])
  }else{
    return(take[1:nInd])
  }
}

#' @title Select males
#' 
#' @description Selects a subset of nInd males from a 
#' population.
#' 
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param nInd the number of individuals to select
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values (\code{gv}), estimated
#' breeding values (\code{ebv}), breeding values (\code{bv}), 
#' or phenotypes (\code{pheno}, default)
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectMale = function(pop,nInd,trait=1,use="pheno",selectTop=TRUE,
                      returnPop=TRUE,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  pop = pop[which(pop@gender=="M")]
  if(nInd>pop@nInd){
    stop(paste("the population only contains",pop@nInd,"males"))
  }
  output = selectInd(pop=pop,nInd=nInd,trait=trait,use=use,
                     selectTop=selectTop,returnPop=returnPop,
                     simParam=simParam,...)
  return(output)
}

#' @title Select females
#' 
#' @description Selects a subset of nInd females from a 
#' population.
#' 
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param nInd the number of individuals to select
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values (\code{gv}), estimated
#' breeding values (\code{ebv}), breeding values (\code{bv}), 
#' or phenotypes (\code{pheno}, default)
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectFemale = function(pop,nInd,trait=1,use="pheno",selectTop=TRUE,
                        returnPop=TRUE,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  pop = pop[which(pop@gender=="F")]
  if(nInd>pop@nInd){
    stop(paste("the population only contains",pop@nInd,"females"))
  }
  output = selectInd(pop=pop,nInd=nInd,trait=trait,use=use,
                     selectTop=selectTop,returnPop=returnPop,
                     simParam=simParam,...)
  return(output)
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
#' @param use select on genetic values (\code{gv}), estimated
#' breeding values (\code{ebv}), breeding values (\code{bv}), 
#' or phenotypes (\code{pheno}, default)
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectFam = function(pop,nFam,trait=1,use="pheno",selectTop=TRUE,
                     returnPop=TRUE,simParam=NULL,...){
  if(is.null(simParam) & use=="bv"){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  families = paste(pop@mother,pop@father,sep="_")
  availFam = length(unique(families))
  if(nFam>availFam){
    stop(paste(nFam,"families requested but only",availFam,
               "families are available"))
  }
  use = tolower(use)
  if(class(trait)=="function"){
    if(use == "gv"){
      response = trait(pop@gv,...)
    }else if(use == "ebv"){
      response = trait(pop@ebv,...)
    }else if(use == "pheno"){
      response = trait(pop@pheno,...)
    }else if(use == "bv"){
      response = varAD(pop,retGenParam=TRUE,simParam=simParam)$bv
      response = trait(response,...)
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }else{
    stopifnot(length(trait)==1,trait<=pop@nTraits)
    if(use == "gv"){
      response = pop@gv[,trait]
    }else if(use == "ebv"){
      response = pop@ebv[,trait]
    }else if(use == "pheno"){
      response = pop@pheno[,trait]
    }else if(use == "bv"){
      response = varAD(pop,retGenParam=TRUE,simParam=simParam)$bv[,trait]
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }
  if(any(is.na(response))){
    stop("selection trait has missing values, phenotype may need to be set")
  }
  #Calculate family means
  famMeans=aggregate(response,list(families=families),mean)
  response = famMeans$x
  #Select families
  take = order(response,decreasing=selectTop)[1:nFam]
  take = families%in%(famMeans$families[take])
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
#' @param use select on genetic values (\code{gv}), estimated
#' breeding values (\code{ebv}), breeding values (\code{bv}), 
#' or phenotypes (\code{pheno}, default)
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param returnPop should results be returned as a 
#' \code{\link{Pop-class}}. If FALSE, only the index of selected 
#' individuals is returned.
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectWithinFam = function(pop,nInd,trait=1,use="pheno",
                           selectTop=TRUE,returnPop=TRUE,
                           simParam=NULL,...){
  if(is.null(simParam) & use=="bv"){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  families = paste(pop@mother,pop@father,sep="_")
  use = tolower(use)
  if(class(trait)=="function"){
    if(use == "gv"){
      response = trait(pop@gv,...)
    }else if(use == "ebv"){
      response = trait(pop@ebv,...)
    }else if(use == "pheno"){
      response = trait(pop@pheno,...)
    }else if(use == "bv"){
      response = varAD(pop,retGenParam=TRUE,simParam=simParam)$bv
      response = trait(response,...)
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }else{
    stopifnot(length(trait)==1,trait<=pop@nTraits)
    if(use == "gv"){
      response = pop@gv[,trait]
    }else if(use == "ebv"){
      response = pop@ebv[,trait]
    }else if(use == "pheno"){
      response = pop@pheno[,trait]
    }else if(use == "bv"){
      response = varAD(pop,retGenParam=TRUE,simParam=simParam)$bv[,trait]
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }
  if(any(is.na(response))){
    stop("selection trait has missing values, phenotype may need to be set")
  }
  selInFam = function(selFam){
    index = which(families%in%selFam)
    y = response[index]
    index = index[order(y,decreasing=selectTop)]
    index = index[1:min(nInd,length(index))]
    return(index)
  }
  take = unlist(sapply(unique(families),selInFam))
  if(returnPop){
    return(pop[take])
  }else{
    return(take)
  }
}
