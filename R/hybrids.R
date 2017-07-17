#A wrapper for calling getHybridGv
#This function uses chunking to reduce RAM usage
getHybridGvByChunk = function(trait,fPop,fPar,
                              mPop,mPar,w,chunkSize){
  nOut = length(fPar)
  if(nOut<=chunkSize){
    output = getHybridGv(trait,fPop,fPar,mPop,mPar,w)
  }else{
    Chunks = split(1:nOut,ceiling(seq_along(1:nOut)/chunkSize))
    output = matrix(NA_real_,nrow=nOut,ncol=1)
    for(chunk in Chunks){
      output[chunk,] = getHybridGv(trait,fPop,fPar[chunk],
                                   mPop,mPar[chunk],w)
    }
  }
  return(output)
}


#' @title Hybrid crossing
#' 
#' @description 
#' A convience function for hybrid plant breeding simulations. Allows for 
#' easy specification of a test cross scheme and/or creation of an object 
#' of \code{\link{HybridPop-class}}. Note that the \code{\link{HybridPop-class}} 
#' should only be used if the parents were created using the \code{\link{makeDH}} 
#' function or \code{\link{newPop}} using inbred founders. The id for 
#' new individuals is [mother_id]_[father_id]
#' 
#' @param fPop female population, an object of \code{\link{Pop-class}}
#' @param mPop male population, an object of \code{\link{Pop-class}}
#' @param crossPlan either "testcross" for all possible combinantions 
#' or a matrix with two columns for designed crosses
#' @param varE error variance for phenotypes. If NULL, phenotypes 
#' aren't calculated.
#' @param w environmental covariate for phenotypes
#' @param reps number of replications for phenotype. See 
#' \code{\link{setPheno}} for details.
#' @param returnHybridPop should results be returned as 
#' \code{\link{HybridPop-class}}. If false returns results as 
#' \code{\link{Pop-class}}. Population must be fully inbred if TRUE.
#' @param chunkSize when using returnHybridPop=TRUE, this 
#' parameter determines the maximum number of hybrids created 
#' at one time. Smaller values reduce RAM usage, but may take 
#' more time.
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @export
hybridCross = function(fPop,mPop,crossPlan="testcross",varE=NULL,
                       w=0.5,reps=1,returnHybridPop=FALSE,
                       chunkSize=1000,
                       simParam=SIMPARAM){
  stopifnot(simParam@gender=="no")
  if(simParam@ploidy!=2){
    stop("Only works with diploids")
  }
  #crossPlan for test cross
  if(crossPlan=="testcross"){
    crossPlan = cbind(rep(1:fPop@nInd,each=mPop@nInd),
                      rep(1:mPop@nInd,fPop@nInd))
  }
  #Set id
  fPar = fPop@id[crossPlan[,1]]
  mPar = mPop@id[crossPlan[,2]]
  id = paste(fPar,mPar,sep="_")
  
  #Return Pop-class
  if(!returnHybridPop){
    output = makeCross2(fPop,mPop,crossPlan,id,simParam)
    if(is.null(varE)){
      return(output)
    }else{ #Calculate phenotype
      output@pheno = calcPheno(pop=output,varE=varE,reps=reps,
                               w=w,simParam=simParam)
      return(output)
    }
  }
  
  #Return HybridPop-class
  #Calculate gv and pheno
  gv = NULL
  pheno = NULL
  for(trait in simParam@traits){
    tmp = getHybridGvByChunk(trait=trait,fPop=fPop,fPar=crossPlan[,1],
                             mPop=mPop,mPar=crossPlan[,2],w=0.5,
                             chunkSize=chunkSize)
    gv = cbind(gv, tmp)
    #Will a phenotype be calculated
    if(is.null(varE)){
      pheno = cbind(pheno,matrix(rep(NA_real_,length(fPar)),ncol=1))
    }else{
      #Does GxE matter
      if(class(trait)=="TraitAG" | class(trait)=="TraitADG"){
        pheno = cbind(pheno, 
                      getHybridGvByChunk(trait=trait,fPop=fPop,fPar=crossPlan[,1],
                                         mPop=mPop,mPar=crossPlan[,2],w=w,
                                         chunkSize=chunkSize))
      }else{
        pheno = cbind(pheno,tmp)
      }
    }
  }
  if(!is.null(varE)) pheno = addError(gv=pheno,varE=varE,reps=reps)
  
  output = new("HybridPop",
               nInd=length(id),
               id=id,
               mother=fPar,
               father=mPar,
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=pheno)
  return(output)
}

#' @title Calculate GCA
#' 
#' @description 
#' Calculate general combining ability of test crosses. Intended for 
#' output from hybridCross using the "testcross" option, but will work 
#' for any population.
#' 
#' @param pop an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param use true genetic value (\code{gv}) or phenotypes (\code{pheno}, default)
#' 
#' @export
calcGCA = function(pop,use="pheno"){
  use = tolower(use)
  if(use == "gv"){
    y=pop@gv
  }else if(use == "pheno"){
    y=pop@pheno
    if(any(is.na(y))){
      stop("Missing values in pop@pheno")
    }
  }else{
    stop(paste0("Use=",use," is not an option"))
  }
  colnames(y) = paste0("Trait",1:pop@nTraits)
  female = factor(pop@mother,
                  levels=unique(pop@mother))
  male = factor(pop@father,
                levels=unique(pop@father))
  females = matrix(NA_real_,nrow=length(unique(female)),ncol=pop@nTraits)
  males = matrix(NA_real_,nrow=length(unique(male)),ncol=pop@nTraits)
  colnames(females) = colnames(males) = colnames(y)
  for(i in 1:pop@nTraits){
    #Calculate female GCA
    if(length(unique(male))==1){
      females[,i] = y[,i]
    }else if(length(unique(female))==1){
      females[,i] = mean(y[,i])
    }else{
      ans = lm(y[,i]~female+male-1,contrasts=list(male="contr.sum"))
      females[,i] = coef(ans)[1:length(unique(female))]
    }
    #Calculate male GCA
    if(length(unique(female))==1){
      males[,i] = y[,i]
    }else if(length(unique(male))==1){
      males[,i] = mean(y[,i])
    }else{
      ans = lm(y[,i]~male+female-1,contrasts=list(female="contr.sum"))
      males[,i] = coef(ans)[1:length(unique(male))]
    }
  }
  #Create output
  output = list()
  output$females = data.frame(female=unique(female),females)
  output$males = data.frame(male=unique(male),males)
  output$females$female = as.character(output$females$female)
  output$males$male = as.character(output$males$male)
  #SCA
  output$SCA=aggregate(y,list(female=female,male=male),mean)
  output$SCA$female = as.character(output$SCA$female)
  output$SCA$male = as.character(output$SCA$male)
  return(output)
}

#' @title Set GCA as phenotype
#' 
#' @description 
#' Calculates general combining ability from a set of testers and 
#' returns these values as phenotypes for a population.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param testers an object of \code{\link{Pop-class}}
#' @param use true genetic value (\code{gv}) or phenotypes (\code{pheno}, default)
#' @param varE error variances for phenotype if \code{use="pheno"}. A vector
#' of length nTraits for independent error or a square matrix of 
#' dimensions nTraits for correlated errors.
#' @param reps number of replications for phenotype. See details.
#' @param w the environmental covariate used by GxE traits.
#' @param inbred are both pop and testers fully inbred. They are only 
#' fully inbred if created by \code{\link{newPop}} using inbred founders 
#' or by the \code{\link{makeDH}} function
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @details
#' The reps parameter is for convient representation of replicated data. 
#' It was intended for representation of replicated yield trials in plant 
#' breeding programs. In this case, varE is set to the plot error and 
#' reps is set to the number plots per entry. The resulting phenotype 
#' would reflect the mean of all replications.
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
setPhenoGCA = function(pop,testers,use="pheno",varE=NULL,reps=1,
                       w=0.5,inbred=FALSE,simParam=SIMPARAM){
  stopifnot(class(pop)=="Pop",class(testers)=="Pop")
  use = tolower(use)
  if(use == "pheno"){
    if(is.null(varE)){
      stop("varE must be specified if use=\"pheno\"")
    }
  }
  tmp = hybridCross(fPop=pop,mPop=testers,crossPlan="testcross",
                    varE=varE,w=w,reps=reps,returnHybridPop=inbred,
                    simParam=simParam)
  tmp = calcGCA(pop=tmp,use=use)
  pop@pheno = as.matrix(tmp$females[,-1])
  return(pop)
}

