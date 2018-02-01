#A wrapper for calling getHybridGv
#This function uses chunking to reduce RAM usage
getHybridGvByChunk = function(trait,females,femaleParents,
                              males,maleParents,chunkSize){
  nOut = length(femaleParents)
  if(nOut<=chunkSize){
    output = getHybridGv(trait,females,femaleParents,males,maleParents)
  }else{
    Chunks = split(1:nOut,ceiling(seq_along(1:nOut)/chunkSize))
    output = list()
    output[[1]] = matrix(NA_real_,nrow=nOut,ncol=1)
    for(chunk in Chunks){
      tmp = getHybridGv(trait,females,femaleParents[chunk],
                        males,maleParents[chunk])
      output[[1]][chunk,] = tmp[[1]]
      if(length(tmp)==2){
        if(length(output)==2){
          output[[2]] = c(output[[2]],tmp[[2]])
        }else{
          output[[2]] = tmp[[2]]
        }
      }
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
#' @param females female population, an object of \code{\link{Pop-class}}
#' @param males male population, an object of \code{\link{Pop-class}}
#' @param crossPlan either "testcross" for all possible combinantions 
#' or a matrix with two columns for designed crosses
#' @param returnHybridPop should results be returned as 
#' \code{\link{HybridPop-class}}. If false returns results as 
#' \code{\link{Pop-class}}. Population must be fully inbred if TRUE.
#' @param chunkSize when using returnHybridPop=TRUE, this 
#' parameter determines the maximum number of hybrids created 
#' at one time. Smaller values reduce RAM usage, but may take 
#' more time.
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
hybridCross = function(females,males,crossPlan="testcross",
                       returnHybridPop=FALSE,chunkSize=10000,
                       simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(simParam$ploidy!=2){
    stop("Only works with diploids")
  }
  #crossPlan for test cross
  if(crossPlan=="testcross"){
    crossPlan = cbind(rep(1:females@nInd,each=males@nInd),
                      rep(1:males@nInd,females@nInd))
  }
  #Set id
  femaleParents = females@id[crossPlan[,1]]
  maleParents = males@id[crossPlan[,2]]
  id = paste(femaleParents,maleParents,sep="_")
  
  #Return Pop-class
  if(!returnHybridPop){
    return(makeCross2(females=females,males=males,
                      crossPlan=crossPlan,
                      simParam=simParam))
  }
  
  #Return HybridPop-class
  gv = pheno = matrix(NA_real_,nrow=length(id),
                      ncol=simParam$nTraits)
  gxe = vector("list",simParam$nTraits)
  i = 0L
  for(trait in simParam$traits){
    i = i+1L
    tmp = getHybridGvByChunk(trait=trait,females=females,femaleParents=crossPlan[,1],
                             males=males,maleParents=crossPlan[,2],chunkSize=chunkSize)
    gv[,i] = tmp[[1]]
    if(length(tmp)==2){
      gxe[[i]] = tmp[[2]]
    }
  }
  
  output = new("HybridPop",
               nInd=length(id),
               id=id,
               mother=femaleParents,
               father=maleParents,
               nTraits=simParam$nTraits,
               gv=gv,
               pheno=pheno,
               gxe=gxe)
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
#' @param use true genetic value "gv" or phenotypes "pheno" (default)
#' 
#' @export
calcGCA = function(pop,use="pheno"){
  if(use=="pheno"){
    y = pop@pheno
  }else if(use=="gv"){
    y = pop@gv
  }else{
    stop(paste0("use=",use," is not a valid option"))
  }
  female = factor(pop@mother,
                  levels=unique(pop@mother))
  male = factor(pop@father,
                levels=unique(pop@father))
  sca = paste(as.character(female),as.character(male),sep="_")
  sca = factor(sca,levels=unique(sca))
  # Female GCA
  if(length(unique(female))==1){
    GCAf = matrix(colMeans(y),nrow=1)
  }else{
    if(length(unique(male))==1){
      GCAf = y
    }else{
      X = model.matrix(~female+male-1,contrasts=list(male="contr.sum"))
      GCAf = calcCoef(X,y)[1:length(unique(female)),,drop=FALSE]
    }
  }
  GCAf = data.frame(as.character(unique(female)),
                    GCAf,stringsAsFactors=FALSE)
  names(GCAf) = c("id",paste0("Trait",1:pop@nTraits))
  # Male GCA
  if(length(unique(male))==1){
    GCAm = matrix(colMeans(y),nrow=1)
  }else{
    if(length(unique(female))==1){
      GCAm = y
    }else{
      X = model.matrix(~male+female-1,contrasts=list(female="contr.sum"))
      GCAm = calcCoef(X,y)[1:length(unique(male)),,drop=FALSE]
    }
  }
  GCAm = data.frame(as.character(unique(male)),
                    GCAm,stringsAsFactors=FALSE)
  names(GCAm) = c("id",paste0("Trait",1:pop@nTraits))
  # SCA
  if(length(unique(sca))==1){
    SCA = y
  }else{
    X = model.matrix(~sca-1)
    SCA = calcCoef(X,y)
  }
  SCA = data.frame(as.character(unique(sca)),
                   SCA,stringsAsFactors=FALSE)
  names(SCA) = c("id",paste0("Trait",1:pop@nTraits))
  return(list(GCAf=GCAf,
              GCAm=GCAm,
              SCA=SCA))
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
#' @param p the p-value for the environmental covariate 
#' @param inbred are both pop and testers fully inbred. They are only 
#' fully inbred if created by \code{\link{newPop}} using inbred founders 
#' or by the \code{\link{makeDH}} function
#' @param chunkSize when using inbred=TRUE, this 
#' parameter determines the maximum number of hybrids created 
#' at one time. Smaller values reduce RAM usage, but may take 
#' more time.
#' @param onlyPheno should only the phenotype be returned
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @details
#' The reps parameter is for convient representation of replicated data. 
#' It was intended for representation of replicated yield trials in plant 
#' breeding programs. In this case, varE is set to the plot error and 
#' reps is set to the number plots per entry. The resulting phenotype 
#' would reflect the mean of all replications.
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' a matrix if onlyPheno=TRUE
#' 
#' @export
setPhenoGCA = function(pop,testers,use="pheno",varE=NULL,reps=1,
                       p=0.5,inbred=FALSE,chunkSize=10000,
                       onlyPheno=FALSE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  stopifnot(class(pop)=="Pop",class(testers)=="Pop")
  use = tolower(use)
  if(use == "pheno"){
    if(is.null(varE)){
      stop("varE must be specified if use=\"pheno\"")
    }
  }
  tmp = hybridCross(females=pop,males=testers,crossPlan="testcross",
                    returnHybridPop=inbred,chunkSize=chunkSize,simParam=simParam)
  if(use=="pheno"){
    tmp = setPheno(tmp,varE=varE,p=p,reps=reps,simParam=simParam)
  }
  tmp = calcGCA(pop=tmp,use=use)
  if(onlyPheno){
    return(as.matrix(tmp$GCAf[,-1]))
  }
  pop@pheno = as.matrix(tmp$GCAf[,-1])
  return(pop)
}
