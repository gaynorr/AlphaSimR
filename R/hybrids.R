#' @title Test cross
#' 
#' @description Creates test crosses of all possible combinations
#' 
#' @param fPop female population, an object of 'Pop' superclass
#' @param mPop male population, an object of 'Pop' superclass
#' @param crossPlan either "test" for a testcross or matrix with two columns for designed crosses
#' @param varE error variance for phenotypes, if NULL phenotypes aren't calculated
#' @param w environmental covariate for phenotype
#' @param returnHybridPop should result be returned as a class 'HybridPop'
#' @param simParam an object of 'SimParam' class
#' 
#' @export
hybridCross = function(fPop,mPop,crossPlan="test",varE=NULL,w=0,
                       returnHybridPop=TRUE,simParam=SIMPARAM){
  if(simParam@ploidy!=2){
    stop("Only works with diploids")
  }
  if(crossPlan=="test"){
    fPar = rep(1:fPop@nInd,each=mPop@nInd)
    mPar = rep(1:mPop@nInd,fPop@nInd)
  }else{
    fPar = crossPlan[,1]
    mPar = crossPlan[,2]
  }
  if(returnHybridPop){
    gv = NULL
    pheno = NULL
    for(trait in simParam@traits){
      tmp = getHybridGv(trait,fPop,fPar,mPop,mPar,w=0)
      gv = cbind(gv, tmp)
      #Will a phenotype be calculated
      if(is.null(varE)){
        pheno = cbind(pheno,rep(NA_real_,length(fPar)))
      }else{
        #Does GxE matter
        if(class(trait)=="TraitAG" | class(trait)=="TraitADG"){
          pheno = cbind(pheno, getHybridGv(trait,fPop,fPar,mPop,mPar,w=w))
        }else{
          pheno = cbind(pheno,tmp)
        }
      }
    }
    if(!is.null(varE)) pheno = addError(pheno,varE)
  }else{
    geno = cross2(fPop@geno,fPar,
                  mPop@geno,mPar,
                  simParam@genMaps)
    output = new("Pop")
    output@nInd=as.integer(length(fPar))
    output@nChr=fPop@nChr
    output@ploidy=fPop@ploidy
    output@gender=rep("H",length(fPar))
    output@geno=geno
    validObject(output)
  }
  if(class(fPop)=="PedPop"){
    fPar = fPop@id[fPar]
  }
  if(class(mPop)=="PedPop"){
    mPar = mPop@id[mPar] 
  }
  id = paste(fPar,mPar,sep="_")
  if(returnHybridPop){
    output = new("HybridPop",
                 nInd=length(id),
                 nTraits=simParam@nTraits,
                 gv=gv,
                 pheno=pheno,
                 id=id,
                 par1=fPar,
                 par2=mPar)
  }else{
    output = addPed(pop=output, id=id, par1=fPar, 
                    par2=mPar, simParam=simParam)
    if(!is.null(varE)){
      output = setPheno(pop=output, varE=varE,
                        w=w, simParam=simParam)
    }
  }
  return(output)
}

#' @title Calculate GCA
#' 
#' @description Calculate general combining ability of test crosses
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param useGv should genetic values be used insead of phenotypes
#' 
#' @export
calcGCA = function(pop,useGv=FALSE){
  if(useGv){
    y=pop@gv
  }else{
    y=pop@pheno
  }
  colnames(y) = paste0("Trait",1:pop@nTraits)
  output = list()
  output$females=aggregate(y,list(female=factor(pop@par1,
                                                levels=unique(pop@par1))),
                           mean)
  output$females$female = as.character(output$females$female)
  output$males=aggregate(y,list(male=factor(pop@par2,
                                            levels=unique(pop@par2))),
                         mean)
  output$males$male = as.character(output$males$male)
  output$SCA=aggregate(y,list(female=factor(pop@par1,
                                            levels=unique(pop@par1)),
                              male=factor(pop@par2,
                                          levels=unique(pop@par2))),
                       mean)
  output$SCA$female = as.character(output$SCA$female)
  output$SCA$male = as.character(output$SCA$male)
  return(output)
}

