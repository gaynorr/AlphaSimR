
#' @title Make designed crosses
#'
#' @param pop an object of superclass 'Pop'
#' @param crossPlan a matrix with two column representing female and male parents
#' @param id optional ids to give to progeny if pop is class 'PedPop'
#' @param simParam an object of class 'SimParam'
#'
#' @export
makeCross = function(pop,crossPlan,id=NULL,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  geno = cross2(pop@geno,crossPlan[,1],
                pop@geno,crossPlan[,2],
                simParam@genMaps)
  output = new("Pop")
  output@nInd=as.integer(nrow(crossPlan))
  output@nChr=pop@nChr
  output@ploidy=pop@ploidy
  output@gender=rep("H",nrow(crossPlan))
  output@geno=geno
  validObject(output)
  if(class(pop)=="PedPop"){
    if(is.null(id)) id = 1:nrow(crossPlan)
    output = addPed(pop=output,id=id,par1=pop@id[crossPlan[,1]],
                    par2=pop@id[crossPlan[,2]],simParam=simParam)
  }
  return(output)
}

#' @title Make random crosses
#' 
#' @description Randomly selects parental combinations for crossing. Does not consider gender.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross
#' @param id optional id to assign to F1s
#' @param simParam an object of 'SimParam' class
#' 
#' @export
randCross = function(pop,nCrosses,nProgeny=1,
                     id=NULL,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  nInd = pop@nInd
  stopifnot(nCrosses<=(nInd*(nInd-1)/2))
  parComb = t(combn(1:nInd,2))
  parComb = parComb[sample.int(nInd*(nInd-1)/2,nCrosses),]
  femalePar = rep(parComb[,1],nProgeny)
  malePar = rep(parComb[,2],nProgeny)
  geno = cross2(pop@geno,femalePar,
                  pop@geno,malePar,
                  simParam@genMaps)
  output = new("Pop")
  output@nInd=as.integer(nCrosses*nProgeny)
  output@nChr=pop@nChr
  output@ploidy=pop@ploidy
  output@gender=rep("H",nCrosses*nProgeny)
  output@geno=geno
  validObject(output)
  if(class(pop)=="PedPop"){
    femalePar = pop@id[as.integer(femalePar)]
    malePar = pop@id[as.integer(malePar)]
  }
  if(is.null(id)){
    id=1:output@nInd
  }
  output = addPed(output,id=id,
                  par1=femalePar,par2=malePar,
                  simParam=simParam)
  return(output)
}

#' @title Generates DH lines
#' 
#' @description Creates DH lines from each individual in a population.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nDH total number of DH lines per individual
#' @param id optional id to assign to DH lines
#' @param simParam an object of 'SimParam' class
#' 
#' @export
makeDH = function(pop,nDH,id=NULL,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  #Should be replaced with something more efficient
  geno = cross2(pop@geno,rep(1:pop@nInd,nDH),
                  pop@geno,rep(1:pop@nInd,nDH),
                  simParam@genMaps)
  for(i in 1:pop@nChr){
    geno[[i]][[2]] = geno[[i]][[1]]
  }
  output = new("Pop")
  output@nInd=as.integer(pop@nInd*nDH)
  output@nChr=pop@nChr
  output@ploidy=pop@ploidy
  output@gender=rep("H",pop@nInd*nDH)
  output@geno=geno
  validObject(output)
  if(class(pop)=="PedPop"){
    if(is.null(id)) id=1:pop@nInd
    output = addPed(output,id=id,
                    par1=rep(pop@par1,nDH),
                    par2=rep(pop@par2,nDH),
                    simParam=simParam)
  }
  return(output)
}

#' @title Self individuals
#' 
#' @description Creates selfed progeny from each individual in a population.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nProgeny total number of selfed progeny per individual
#' @param id optional id to give lines is pop is class 'PedPop'
#' @param simParam an object of 'SimParam' class
#' 
#' @export
self = function(pop,nProgeny,id=NULL,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  geno = cross2(pop@geno,rep(1:pop@nInd,nProgeny),
                pop@geno,rep(1:pop@nInd,nProgeny),
                simParam@genMaps)
  output = new("Pop")
  output@nInd=as.integer(pop@nInd*nProgeny)
  output@nChr=pop@nChr
  output@ploidy=pop@ploidy
  output@gender=rep("H",pop@nInd*nProgeny)
  output@geno=geno
  validObject(output)
  if(class(pop)=="PedPop"){
    if(is.null(id)) id = 1:output@nInd
    output = addPed(pop=output,id=id,
                    par1=rep(pop@par1,each=nProgeny),
                    par2=rep(pop@par2,each=nProgeny))
  }
  return(output)
}

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
#' @param pop an object of class 'HybridPop' or 'PedPop'
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

