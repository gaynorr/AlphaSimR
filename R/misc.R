#' @title Selection Intensity
#' 
#' @description 
#' Calculates the standardized selection intensity
#' 
#' @param p the proportion of individuals selected
#' 
#' @export
selInt = function(p){
  return(dnorm(qnorm(1-p))/p)
}

#' @title Calculate Smith-Hazel Weights
#' 
#' @description
#' Calculates weights for Smith-Hazel index given economice weights 
#' and phenotypic and genotypic variance-covariance matrices.
#' 
#' @param econWt vector of economic weights
#' @param varG the genetic variance-covariance matrix
#' @param varP the phenotypic variance-covariance matrix
#' 
#' @return a vector of weight for calculating index values
#' 
#' @export
smithHazel = function(econWt,varG,varP){
  return(solve(varP)%*%varG%*%econWt)
}

#' @title Selection Index
#' 
#' @description
#' Calculates values of a selection index given trait values and 
#' weights. This function is intended to be used in combination with 
#' selection functions working on populations such as 
#' \code{\link{selectInd}}.
#' 
#' @param Y a matrix of trait values
#' @param b a vector of weights
#' @param scale should Y be scaled and centered
#' 
#' @export
selIndex = function(Y,b,scale=FALSE){
  if(scale){
    return(scale(Y)%*%b)
  }
  return(Y%*%b)
}

#' @title Edit Genome
#' 
#' @description
#' Edits selected loci of selected individuals to a homozygous 
#' state for either the 1 or 0 allele. The gv slot is recalculated to 
#' reflect the any changes due to editing, but other slots remain the same.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param ind a vector of individuals to edit
#' @param chr a vector of chromosomes to edit. Length must match 
#' length of segSites.
#' @param segSites a vector of segregating sites to edit. Length must 
#' match length of chr.
#' @param allele either 0 or 1 for desired allele
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
editGenome = function(pop,ind,chr,segSites,allele,
                      simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  ind = unique(as.integer(ind))
  stopifnot(all(ind%in%(1:pop@nInd)))
  chr = as.integer(chr)
  segSites = as.integer(segSites)
  stopifnot(length(chr)==length(segSites))
  allele = as.integer(allele)
  stopifnot(allele==0L | allele==1L)
  allele = as.raw(allele)
  for(selChr in unique(chr)){
    selSegSites = segSites[chr==selChr]
    pop@geno[[selChr]][selSegSites,,ind] = allele
  }
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

#' @title Correlated variable
#' 
#' @description
#' Creates a correlated vector by adding random error.
#'
#' @param x a numeric vector
#' @param rho desired correlation. Must be greater than 
#' 0 and less than or equal to 1.
#'
#' @return a numeric vector
#'
#' @export
corVar = function(x,rho){
  stopifnot(rho>0, rho<=1)
  varX = var(x)
  varY = varX/(rho^2)-varX
  return(x+rnorm(length(x),sd=sqrt(varY)))
}

#' @title Usefulness Criterion
#' 
#' @description Calculates the usefulness criterion
#' 
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values (\code{gv}, default), estimated
#' breeding values (\code{ebv}), breeding values (\code{bv}), 
#' or phenotypes (\code{pheno})
#' @param p the proportion of individuals selected
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns a numeric value
#' 
#' @export
usefulness = function(pop,trait=1,use="gv",p=0.1,
                      selectTop=TRUE,simParam=NULL,...){
  if(is.null(simParam) & use=="bv"){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
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
  response = sort(response,decreasing=selectTop)
  response = response[1:ceiling(p*length(response))]
  return(mean(response))
}
