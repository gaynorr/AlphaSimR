# Functions for starting a new simulation
# These functions cover running MaCS to setting founder populations

#' @title Create new simulation
#'
#' @description Starts the process of building a new simulation
#' 
#' @param maxQtl the maximum number of segSites for QTLs. 
#' Can be a single value or values for each chromosome.
#' @param maxSnp the maximum number of segSites for SNPs. 
#' Can be a single value or values for each chromosome.
#' @param snpQtlOverlap should SNPs and QTLs be allowed to overlap
#' @param minSnpFreq minimum allowable frequency for SNPs
#' @param gender should gender be considered in the simulation. 
#' Options are: no, yes_rand, or yes_sys. See details.
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @details
#' There are three options for gender. Option "no" means gender is not 
#' used. All individuals will have "H" for gender, which represents 
#' hermaphrodite. Option "yes_rand" means that all new individuals are 
#' randomly assigned either an "M" or "F" for gender, which represents 
#' male and female, respectively. Option "yes_sys" means gender assigned 
#' systematically when new individuals are created. Odd individuals recieve 
#' an "M" and even individuals recieve an "F". If gender is used, it will 
#' prevent crossing between same gender individuals using the \code{\link{randCross}} 
#' or \code{\link{randCross2}} functions.
#' 
#' @return Returns an object of \code{\link{SimParam-class}}
#' 
#' @export
createSimulation = function(maxQtl=0,maxSnp=0,snpQtlOverlap=FALSE,
                            minSnpFreq=NULL,gender="no",
                            founderPop=FOUNDERPOP){
  stopifnot(class(founderPop)=="MapPop")
  if(length(maxSnp)==1){
    maxSnp = rep(maxSnp,founderPop@nChr)
  }
  if(length(maxQtl)==1){
    maxQtl = rep(maxQtl,founderPop@nChr)
  }
  stopifnot(length(maxSnp)==founderPop@nChr,
            length(maxQtl)==founderPop@nChr)
  potSnp = list()
  potQtl = list()
  for(chr in 1:founderPop@nChr){
    if(snpQtlOverlap){
      stopifnot(founderPop@nLoci[chr]>=maxSnp[chr],
                founderPop@nLoci[chr]>=maxQtl[chr])
      if(is.null(minSnpFreq)){
        potSnp[[chr]] = sort(sample.int(founderPop@nLoci[chr],
                                        maxSnp[chr]))
      }else{
        q = calcChrMinorFreq(founderPop@geno[[chr]],
                             founderPop@ploidy)
        potSnp[[chr]] = sort(sample(which(q>=minSnpFreq),maxSnp[chr]))
      }
      potQtl[[chr]] = sort(sample.int(founderPop@nLoci[chr],
                                      maxQtl[chr]))
    }else{
      stopifnot(founderPop@nLoci[chr]>=sum(maxSnp[chr],maxQtl[chr]))
      if(is.null(minSnpFreq)){
        tmp = sample.int(founderPop@nLoci[chr],sum(maxSnp[chr],maxQtl[chr]))
        potSnp[[chr]] = sort(tmp[1:maxSnp[chr]])
        potQtl[[chr]] = sort(tmp[(maxSnp[chr]+1):length(tmp)])
      }else{
        q = calcChrMinorFreq(founderPop@geno[[chr]],
                             founderPop@ploidy)
        potSnp[[chr]] = sort(sample(which(q>=minSnpFreq),maxSnp[chr]))
        potQtl[[chr]] = sort(sample(which(!((1:founderPop@nLoci[chr])%in%potSnp[[chr]])),
                                    maxQtl[chr]))
      }
    }
  }
  output = new("SimParam",
               ploidy=founderPop@ploidy,
               nChr=founderPop@nChr,
               nTraits=0L,
               nSnpChips=0L,
               segSites=founderPop@nLoci,
               gender=gender,
               genMaps=founderPop@genMaps,
               traits=list(),
               snpChips=list(),
               potQtl=potQtl,
               potSnp=potSnp)
  assign("LASTID",0L,envir=.GlobalEnv)
  return(output)
}

#' @title Add SNP chip
#' 
#' @description 
#' Randomly assigns eligble SNPs to a SNP chip
#' 
#' @param nSnpPerChr number of SNPs per chromosome. 
#' Can be a single value or nChr values.
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object \code{\link{SimParam-class}}
#' 
#' @export
addSnpChip = function(nSnpPerChr,simParam=SIMPARAM){
  if(length(nSnpPerChr)==1){
    nSnpPerChr = rep(nSnpPerChr,simParam@nChr)
  }
  stopifnot(length(nSnpPerChr)==simParam@nChr)
  stopifnot(sapply(simParam@potSnp,length)>=nSnpPerChr)
  lociLoc = lapply(1:simParam@nChr,function(x){
    sort(sample(simParam@potSnp[[x]],nSnpPerChr[x]))
  })
  lociLoc = do.call("c",lociLoc)
  snpChip = new("LociMap",
                nLoci=as.integer(sum(nSnpPerChr)),
                lociPerChr=as.integer(nSnpPerChr),
                lociLoc=as.integer(lociLoc))
  simParam@nSnpChips = simParam@nSnpChips + 1L
  simParam@snpChips[[simParam@nSnpChips]] = snpChip
  validObject(simParam)
  return(simParam)
}

#' @title Add SNP chips
#' 
#' @description 
#' Randomly selects the number of snps in structure and then
#' assigns them to chips based on structure
#' 
#' @param nSnpPerChr number of SNPs per chromosome. 
#' Can be a single value or nChr values.
#' @param structure a matrix.  Rows are snp chips, columns are chips.
#' If value is true then that snp is on that chip.
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object \code{\link{SimParam-class}}
#' 
#' @export
addStructuredSnpChips = function(nSnpPerChr,structure,simParam=SIMPARAM){
  if(length(nSnpPerChr)==1){
    nSnpPerChr = rep(nSnpPerChr,simParam@nChr)
  }
  stopifnot(length(nSnpPerChr)==simParam@nChr)
  stopifnot(sapply(simParam@potSnp,length)>=nSnpPerChr)
  stopifnot(dim(structure)[2]==sum(nSnpPerChr))
  lociLoc = lapply(1:simParam@nChr,function(x){
    sort(sample(simParam@potSnp[[x]],nSnpPerChr[x]))
  })
  lociLoc = do.call("c",lociLoc)
  
  for (i in 1:nrow(structure)){
    snps = lociLoc[structure[i,]]
    start = 1
    numChr = numeric(length(nSnpPerChr))
    for (j in 1:length(nSnpPerChr)){
      end = start + nSnpPerChr[j] - 1
      numChr[j] = sum(structure[i,start:end])
      start = end + 1
    }
    snpChip = new("LociMap",
                  nLoci = length(snps),
                  lociPerChr = as.integer(numChr),
                  lociLoc = as.integer(snps))
    simParam@nSnpChips = simParam@nSnpChips + 1L
    simParam@snpChips[[simParam@nSnpChips]] = snpChip
  }
  
  validObject(simParam)
  return(simParam)
}

#Function for selecting QTL loci called by all addTrait functions
pickQtlLoci = function(nQtlPerChr, simParam){
  if(length(nQtlPerChr)==1){
    nQtlPerChr = rep(nQtlPerChr,simParam@nChr)
  }
  stopifnot(length(nQtlPerChr)==simParam@nChr)
  stopifnot(sapply(simParam@potQtl,length)>=nQtlPerChr)
  lociLoc = lapply(1:simParam@nChr,function(x){
    sort(sample(simParam@potQtl[[x]],nQtlPerChr[x]))
  })
  lociLoc = do.call("c",lociLoc)
  qtlLoci = new("LociMap",
                nLoci=as.integer(sum(nQtlPerChr)),
                lociPerChr=as.integer(nQtlPerChr),
                lociLoc=as.integer(lociLoc))
  return(qtlLoci)
}

#' @title Add additive traits
#' 
#' @description 
#' Randomly assigns eligble QTLs for one ore more additive traits. 
#' If simulating more than one trait, all traits will be pleiotrophic 
#' with correlated additive effects.
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param meanG a vector of mean genetic values for one or more traits
#' @param varG a vector of total genetic variances for one or more traits
#' @param corr a matrix of correlations between traits
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape value of the shape parameter if using a gamma distribution
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @return Returns an object of \code{\link{SimParam-class}}
#' 
#' @export
addTraitA = function(nQtlPerChr,meanG,varG,corr=matrix(1),
                     gamma=FALSE,shape=1,simParam=SIMPARAM,
                     founderPop=FOUNDERPOP){
  stopifnot(length(meanG)==length(varG),
            nrow(corr)==ncol(corr),
            length(meanG)==length(diag(corr)))
  nTraits = length(meanG)
  qtlLoci = pickQtlLoci(nQtlPerChr,simParam=simParam)
  addEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%chol(corr)
  if(gamma){
    addEff = qgamma(pnorm(addEff),shape=shape)
  }
  geno = getGeno(founderPop@geno,
                 qtlLoci@lociPerChr,
                 qtlLoci@lociLoc)
  for(i in 1:nTraits){
    tmp = tuneTraitA(geno,addEff[,i],varG[i])
    intercept = tmp$output$intercept
    addEff[,i] = addEff[,i]*tmp$parameter
    trait = new("TraitA",
                qtlLoci,
                addEff=addEff[,i],
                intercept=meanG[i]-intercept)
    simParam@nTraits = simParam@nTraits + 1L
    simParam@traits[[simParam@nTraits]] = trait
  }
  validObject(simParam)
  return(simParam)
}

#' @title Add additive and dominance traits
#' 
#' @description 
#' Randomly assigns eligble QTLs for one or more traits with dominance. 
#' If simulating more than one trait, all traits will be pleiotrophic 
#' with correlated additive effects.
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param meanG a vector of mean genetic values for one or more traits
#' @param varG a vector of total genetic variances for one or more traits
#' @param domDegree the dominance degree of individual loci. 
#' Can be a single value or nLoci values.
#' @param corr a matrix of correlations between traits
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape value of the shape parameter if using a gamma distribution
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @return Returns an object of \code{\link{SimParam-class}}
#'  
#' @export
addTraitAD = function(nQtlPerChr,meanG,varG,domDegree,corr=matrix(1),
                      gamma=FALSE,shape=1,simParam=SIMPARAM,
                      founderPop=FOUNDERPOP){
  stopifnot(length(meanG)==length(varG),
            nrow(corr)==ncol(corr),
            length(meanG)==length(diag(corr)))
  nTraits = length(meanG)
  qtlLoci = pickQtlLoci(nQtlPerChr,simParam=simParam)
  addEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%chol(corr)
  if(gamma){
    addEff = qgamma(pnorm(addEff),shape=shape)
  }
  if(length(domDegree)==1){
    domDegree = rep(domDegree,qtlLoci@nLoci)
  }else{
    stopifnot(length(domDegree)==qtlLoci@nLoci) 
  }
  domEff = sweep(abs(addEff),1,domDegree,"*")
  geno = getGeno(founderPop@geno,
                 qtlLoci@lociPerChr,
                 qtlLoci@lociLoc)
  for(i in 1:nTraits){
    tmp = tuneTraitAD(geno,addEff[,i],domEff[,i],varG[i])
    intercept = tmp$output$intercept
    addEff[,i] = addEff[,i]*tmp$parameter
    domEff[,i] = domEff[,i]*tmp$parameter
    trait = new("TraitAD",
                qtlLoci,
                addEff=addEff[,i],
                domEff=domEff[,i],
                intercept=meanG[i]-intercept)
    simParam@nTraits = simParam@nTraits + 1L
    simParam@traits[[simParam@nTraits]] = trait
  }
  validObject(simParam)
  return(simParam)
}

#' @title Add additive GxE traits
#' 
#' @description 
#' Randomly assigns eligble QTLs for one ore more additive GxE traits. 
#' If simulating more than one trait, all traits will be pleiotrophic 
#' with correlated effects.
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param meanG a vector of mean genetic values for one or more traits
#' @param varG a vector of total genetic variances for one or more traits
#' @param varGE a vector of total genotype-by-environment variances for the traits
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape value of the shape parameter if using a gamma distribution
#' @param corr a matrix of correlations between traits
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @return Returns an object of \code{\link{SimParam-class}}
#' 
#' @export
addTraitAG = function(nQtlPerChr,meanG,varG,varGE,corr=matrix(1),
                      gamma=FALSE,shape=1,simParam=SIMPARAM,
                      founderPop=FOUNDERPOP){
  stopifnot(length(meanG)==length(varG),
            nrow(corr)==ncol(corr),
            length(meanG)==length(diag(corr)),
            length(meanG)==length(varGE))
  nTraits = length(meanG)
  qtlLoci = pickQtlLoci(nQtlPerChr,simParam=simParam)
  addEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%chol(corr)
  if(gamma){
    addEff = qgamma(pnorm(addEff),shape=shape)
  }
  geno = getGeno(founderPop@geno,
                 qtlLoci@lociPerChr,
                 qtlLoci@lociLoc)
  for(i in 1:nTraits){
    tmp = tuneTraitA(geno,addEff[,i],varG[i])
    intercept = tmp$output$intercept
    addEff[,i] = addEff[,i]*tmp$parameter
    varGxeLoc = sqrt(popVar(addEff[,i,drop=FALSE])*varGE[i]/varG[i])
    trait = new("TraitAG",
                qtlLoci,
                addEff=addEff[,i],
                intercept=meanG[i]-intercept,
                gxeEff = rnorm(qtlLoci@nLoci,sd=sqrt(varGxeLoc)),
                varGxeLoci = c(varGxeLoc))
    simParam@nTraits = simParam@nTraits + 1L
    simParam@traits[[simParam@nTraits]] = trait
  }
  validObject(simParam)
  return(simParam)
}

#' @title Add an additive and dominance GxE trait
#' 
#' @description 
#' Randomly assigns eligble QTLs for a trait with dominance and GxE. 
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single 
#' value or nChr values.
#' @param meanG a vector of mean genetic values for one or more traits
#' @param varG a vector of total genetic variances for one or more traits
#' @param domDegree the dominance degree of individual loci. 
#' Can be a single value or nLoci values.
#' @param varGE a vector of total genotype-by-environment variances for the traits
#' @param corr a matrix of correlations between traits
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape value of the shape parameter if using a gamma distribution
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @return Returns an object of \code{\link{SimParam-class}}
#'  
#' @export
addTraitADG = function(nQtlPerChr,meanG,varG,domDegree,varGE,corr=matrix(1),
                       gamma=FALSE,shape=1,simParam=SIMPARAM,
                       founderPop=FOUNDERPOP){
  stopifnot(length(meanG)==length(varG),
            nrow(corr)==ncol(corr),
            length(meanG)==length(diag(corr)),
            length(meanG)==length(varGE))
  nTraits = length(meanG)
  qtlLoci = pickQtlLoci(nQtlPerChr,simParam=simParam)
  addEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%chol(corr)
  if(gamma){
    addEff = qgamma(pnorm(addEff),shape=shape)
  }
  if(length(domDegree)==1){
    domDegree = rep(domDegree,qtlLoci@nLoci)
  }else{
    stopifnot(length(domDegree)==qtlLoci@nLoci) 
  }
  domEff = sweep(abs(addEff),1,domDegree,"*")
  geno = getGeno(founderPop@geno,
                 qtlLoci@lociPerChr,
                 qtlLoci@lociLoc)
  for(i in 1:nTraits){
    tmp = tuneTraitAD(geno,addEff[,i],domEff[,i],varG[i])
    intercept = tmp$output$intercept
    addEff[,i] = addEff[,i]*tmp$parameter
    domEff[,i] = domEff[,i]*tmp$parameter
    varGxeLoc = sqrt(popVar(addEff[,i,drop=FALSE])*varGE[i]/varG[i])
    trait = new("TraitADG",
                qtlLoci,
                addEff=addEff[,i],
                domEff=domEff[,i],
                intercept=meanG[i]-intercept,
                gxeEff = rnorm(qtlLoci@nLoci,sd=sqrt(varGxeLoc)),
                varGxeLoci = c(varGxeLoc))
    simParam@nTraits = simParam@nTraits + 1L
    simParam@traits[[simParam@nTraits]] = trait
  }
  validObject(simParam)
  return(simParam)
}

#' @title Create new Population
#' 
#' @description
#' Creates a new \code{\link{Pop-class}} from an object of 
#' \code{\link{MapPop-class}} or \code{\link{RawPop-class}}. 
#' The function is intended for creating initial populations from 
#' 'FOUNDERPOP' created by \code{\link{runMacs}}.
#'
#' @param rawPop an object of \code{\link{MapPop-class}} or 
#' \code{\link{RawPop-class}}
#' @param id optional ids to assign individuals
#' @param simParam an object of \code{\link{SimParam-class}}
#'
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
newPop = function(rawPop, id=NULL, simParam=SIMPARAM){
  stopifnot(class(rawPop)=="RawPop" | class(rawPop)=="MapPop")
  if(is.null(id)){
    lastId = get("LASTID",envir=.GlobalEnv)
    id = (1:rawPop@nInd) + lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  if(simParam@gender=="no"){
    gender = rep("H",rawPop@nInd)
  }else if(simParam@gender=="yes_rand"){
    gender = sample(c("M","F"),rawPop@nInd,replace=TRUE)
  }else if(simParam@gender=="yes_sys"){
    gender = rep_len(c("M","F"),rawPop@nInd)
  }else{
    stop(paste("no rules for gender type",simParam@gender))
  }
  if(simParam@nTraits==0){
    gv = matrix(NA_real_,
                nrow=rawPop@nInd,
                ncol=simParam@nTraits)
  }else{
    gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
    gv = do.call("cbind",gv)
  }
  output = new("Pop",
               nInd=rawPop@nInd,
               nChr=rawPop@nChr,
               ploidy=rawPop@ploidy,
               nLoci=rawPop@nLoci,
               gender=gender,
               geno=rawPop@geno,
               id=as.character(id),
               mother=rep("0",rawPop@nInd),
               father=rep("0",rawPop@nInd),
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits),
               ebv=matrix(NA_real_,
                          nrow=rawPop@nInd,
                          ncol=0))
  if(updateId){
    assign("LASTID",lastId,envir=.GlobalEnv)
  }
  return(output)
}


