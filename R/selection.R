#' @title Select individuals
#' 
#' @description Selects a subset of nInd individuals from a 
#' population.
#' 
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param nInd the number of individuals to select
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a single value.
#' @param use select on genetic value (\code{gv}), estimated
#' genetic values (\code{ebv}) or phenotypes (\code{pheno}, default)
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectInd = function(pop,nInd,trait=1,use="pheno",selectTop=TRUE){
  stopifnot(nInd<=pop@nInd)
  use = tolower(use)
  if(class(trait)=="function"){
    if(use == "gv"){
      response = apply(pop@gv,1,trait)
    }else if(use == "ebv"){
      response = apply(pop@ebv,1,trait)
    }else if(use == "pheno"){
      response = apply(pop@pheno,1,trait)
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
    }else{
      stop(paste0("Use=",use," is not an option"))
    }
  }
  if(any(is.na(response))){
    stop("selection trait has missing values, phenotype may need to be set")
  }
  take = order(response,decreasing=selectTop)
  return(pop[take[1:nInd]])
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
#' a single trait or a function returning a single value.
#' @param use select on genetic value (\code{gv}), estimated
#' genetic values (\code{ebv}) or phenotypes (\code{pheno})
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectMale = function(pop,nInd,trait=1,use="pheno",selectTop=TRUE){
  pop = pop[which(pop@gender=="M")]
  if(nInd>pop@nInd){
    stop(paste("the population only contains",pop@nInd,"males"))
  }
  pop = selectInd(pop=pop,nInd=nInd,trait=trait,use=use,
                  selectTop=selectTop)
  return(pop)
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
#' a single trait or a function returning a single value.
#' @param use select on genetic value (\code{gv}), estimated
#' genetic values (\code{ebv}) or phenotypes (\code{pheno})
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectFemale = function(pop,nInd,trait=1,use="pheno",selectTop=TRUE){
  pop = pop[which(pop@gender=="F")]
  if(nInd>pop@nInd){
    stop(paste("the population only contains",pop@nInd,"females"))
  }
  pop = selectInd(pop=pop,nInd=nInd,trait=trait,use=use,
                  selectTop=selectTop)
  return(pop)
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
#' a single trait or a function returning a single value.
#' @param use select on genetic value (\code{gv}), estimated
#' genetic values (\code{ebv}) or phenotypes (\code{pheno}, default)
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectFam = function(pop,nFam,trait=1,use="pheno",selectTop=TRUE){
  families = paste(pop@mother,pop@father,sep="_")
  availFam = length(unique(families))
  if(nFam>availFam){
    stop(paste(nFam,"families requested but only",availFam,
               "families are available"))
  }
  use = tolower(use)
  if(class(trait)=="function"){
    if(use == "gv"){
      response = apply(pop@gv,1,trait)
    }else if(use == "ebv"){
      response = apply(pop@ebv,1,trait)
    }else if(use == "pheno"){
      response = apply(pop@pheno,1,trait)
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
  return(pop[take])
}

#' @title Select individuals within families
#' 
#' @description Selects a subset of nInd individuals from each  
#' full-sib family within a population.
#' 
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param nInd the number of individuals to select within a family
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a single value.
#' @param use select on genetic value (\code{gv}), estimated
#' genetic values (\code{ebv}) or phenotypes (\code{pheno}, default)
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
selectWithinFam = function(pop,nInd,trait=1,use="pheno",
                           selectTop=TRUE){
  families = paste(pop@mother,pop@father,sep="_")
  if(any(table(families)<nInd)){
    stop("some families have less than nInd individuals")
  }
  use = tolower(use)
  if(class(trait)=="function"){
    if(use == "gv"){
      response = apply(pop@gv,1,trait)
    }else if(use == "ebv"){
      response = apply(pop@ebv,1,trait)
    }else if(use == "pheno"){
      response = apply(pop@pheno,1,trait)
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
    index = index[1:nInd]
    return(index)
  }
  take = c(sapply(unique(families),selInFam))
  return(pop[take])
}

