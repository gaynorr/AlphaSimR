#' @title Set phenotype
#' 
#' @description Calculates phenotypes for all traits by adding random error from a multivariate normal distribution.
#' 
#' @param pop an object of superclass 'pop'
#' @param eVar error variances for phenotype. A vector of length nTraits for independent error or a square matrix of dimensions nTraits for correlated errors.
#' @param simParam an object of class 'SimParam'
#' 
#' @export
setPheno = function(pop,eVar,simParam=SIMPARAM){
  if(is.matrix(eVar)){
    stopifnot(nrow(eVar)==simParam@nTraits,
              ncol(eVar)==simParam@nTraits)
  }else{
    stopifnot(length(eVar)==simParam@nTraits)
    if(length(eVar)==1){
      eVar = matrix(eVar)
    }else{
      eVar = diag(eVar)
    }
  }
  if(class(pop)=="Pop"){
    pop = addGv(pop,simParam=simParam)
  }
  pop@pheno = pop@gv + MASS::mvrnorm(pop@nInd,
                                     mu=rep(0,simParam@nTraits),
                                     Sigma=eVar)
  return(pop)
}

#' @title Select individuals
#' 
#' @description Selects a subset of nInd individuals from a 'Pop' superclass using various types of traits.
#' 
#' @param pop and object of superclass 'Pop'
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

#' @title Population summary
#' 
#' @description Generates summary statistics for an object of 'Pop' superclass
#' 
#' @param pop an object of 'Pop' superclass
#' @param simParam an object of class 'SimParam'
#' @param w environmental covariate
#' 
#' @export
popSummary = function(pop,simParam=SIMPARAM,w=0){
  if(class(pop)=="Pop"){
    pop = addGv(pop,simParam=simParam)
  }
  output = list(bv=NULL,dd=NULL)
  #Loop through bv and dd calculations
  for(i in 1:simParam@nTraits){
    trait = simParam@traits[[i]]
    if(class(trait)=="TraitA"){
      tmp = list()
      tmp$bv = scale(pop@gv[,i],scale=F)
      tmp$dd = rep(0,pop@nInd)
      tmp$alpha = trait@addEff
    }else if(class(trait)=="TraitAD"){
      tmp = calcGenParam(trait,pop,trait@addEff,trait@domEff)
    }else if(class(trait)=="TraitAG"){
      
    }else if(class(trait)=="TraitADG"){
      
    }else{
      stop("No method for trait class",class(simParam@traits[[i]]))
    }
    output$bv = cbind(output$bv,tmp$bv)
    output$dd = cbind(output$dd,tmp$dd)
  }
  output$meanG = colMeans(pop@gv)
  output$meanP = colMeans(pop@pheno)
  output$varA = var(output$bv)
  output$varD = var(output$dd)
  output$varG = var(pop@gv)
  output$varP = var(pop@pheno)
  return(output)
}


