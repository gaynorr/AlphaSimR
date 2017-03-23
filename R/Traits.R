#' @title Set phenotype
#' 
#' @description Calculates phenotypes for all traits by adding random error from a multivariate normal distribution.
#' 
#' @param pop an object of superclass 'pop'
#' @param varE error variances for phenotype. A vector of length nTraits for independent error or a square matrix of dimensions nTraits for correlated errors.
#' @param simParam an object of class 'SimParam'
#' 
#' @export
setPheno = function(pop,varE,w=0,simParam=SIMPARAM){
  if(is.matrix(varE)){
    stopifnot(nrow(varE)==simParam@nTraits,
              ncol(varE)==simParam@nTraits)
  }else{
    stopifnot(length(varE)==simParam@nTraits)
    if(length(varE)==1){
      eVar = matrix(varE)
    }else{
      eVar = diag(varE)
    }
  }
  if(class(pop)=="Pop"){
    pop = addGv(pop,simParam=simParam)
  }
  gv = pop@gv
  for(i in 1:simParam@nTraits){
    traitClass = class(simParam@traits[[i]])
    if(traitClass=="TraitAG" | traitClass=="TraitADG"){
      gv[,i] = getGv(simParam@traits[[i]],pop=pop,w=w)
    }
  }
  pop@pheno = gv + MASS::mvrnorm(pop@nInd,
                                 mu=rep(0,simParam@nTraits),
                                 Sigma=varE)
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
  output$varA = popVar(output$bv)
  output$varD = popVar(output$dd)
  output$varG = popVar(pop@gv)
  output$varP = popVar(pop@pheno)
  return(output)
}

#' @title Mean genetic values
#' 
#' @description Returns the mean genetic values for all traits
#' 
#' @param pop an object of class 'TraitPop' or 'PedPop'
#' 
#' @export
meanG = function(pop){
  stopifnot(class(pop)=="TraitPop" | class(pop)=="PedPop")
  colMeans(pop@gv)
}

#' @title Mean phenotypic values
#' 
#' @description Returns the mean phenotypic values for all traits
#' 
#' @param pop an object of class 'TraitPop' or 'PedPop'
#' 
#' @export
meanP = function(pop){
  stopifnot(class(pop)=="TraitPop" | class(pop)=="PedPop")
  colMeans(pop@pheno)
}

#' @title Total genetic variance
#' 
#' @description Returns total genetic variance for all traits
#' 
#' @param pop an object of class 'TraitPop' or 'PedPop'
#' 
#' @export
varG = function(pop){
  stopifnot(class(pop)=="TraitPop" | class(pop)=="PedPop")
  popVar(pop@gv)
}

#' @title Phenotypic variance
#' 
#' @description Returns phenotypic variance for all traits
#' 
#' @param pop an object of class 'TraitPop' or 'PedPop'
#' 
#' @export
varP = function(pop){
  stopifnot(class(pop)=="TraitPop" | class(pop)=="PedPop")
  popVar(pop@pheno)
}
