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
#' @param useGender should gender be considered in the simulation
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @return Returns an object of \code{\link{SimParam-class}}
#' 
#' @export
createSimulation = function(maxQtl,maxSnp,snpQtlOverlap=FALSE,
                            minSnpFreq=NULL,useGender=FALSE,
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
        q = AlphaSimR:::calcChrMinorFreq(founderPop@geno[[chr]],
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
        q = AlphaSimR:::calcChrMinorFreq(founderPop@geno[[chr]],
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
               useGender=useGender,
               genMaps=founderPop@genMaps,
               traits=list(),
               snpChips=list(),
               potQtl=potQtl,
               potSnp=potSnp,
               lastId=0L)
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

#' @title Add an additive trait
#' 
#' @description 
#' Randomly assigns eligble QTLs for an additive trait. 
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param meanG the mean genetic value for the trait
#' @param varG the total genetic variance for the trait
#' @param simParm an object of \code{\link{SimParam-class}}
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @return Returns an object of \code{\link{SimParam-class}}
#' 
#' @export
addTraitA = function(nQtlPerChr,meanG,varG,simParam=SIMPARAM,
                     founderPop=FOUNDERPOP){
  qtlLoci = pickQtlLoci(nQtlPerChr,simParam=simParam)
  addEff = rnorm(qtlLoci@nLoci)
  geno = AlphaSimR:::getGeno(founderPop@geno,
                             qtlLoci@lociPerChr,
                             qtlLoci@lociLoc)
  tmp = AlphaSimR:::tuneTraitA(geno,addEff,varG)
  intercept = tmp$output$intercept
  addEff = addEff*tmp$parameter
  trait = new("TraitA",
              qtlLoci,
              addEff=addEff,
              intercept=meanG-intercept)
  simParam@nTraits = simParam@nTraits + 1L
  simParam@traits[[simParam@nTraits]] = trait
  validObject(simParam)
  return(simParam)
}

#' @title Add an additive and dominance
#' 
#' @description 
#' Randomly assigns eligble QTLs for a trait with dominance. 
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param meanG the mean genetic value for the trait
#' @param varG the total genetic variance for the trait
#' @param domDegree the dominance degree of individual loci. Can be a single value or nLoci values.
#' @param simParm an object of \code{\link{SimParam-class}}
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @return Returns an object of \code{\link{SimParam-class}}
#'  
#' @export
addTraitAD = function(nQtlPerChr,meanG,varG,domDegree,simParam=SIMPARAM,
                      founderPop=FOUNDERPOP){
  qtlLoci = pickQtlLoci(nQtlPerChr,simParam=simParam)
  addEff = rnorm(qtlLoci@nLoci)
  if(length(domDegree)==1){
    domDegree = rep(domDegree,qtlLoci@nLoci)
  }else{
    stopifnot(length(domDegree)==qtlLoci@nLoci) 
  }
  addEff = rnorm(qtlLoci@nLoci)
  domEff = abs(addEff)*domDegree
  geno = AlphaSimR:::getGeno(founderPop@geno,
                             qtlLoci@lociPerChr,
                             qtlLoci@lociLoc)
  tmp = AlphaSimR:::tuneTraitAD(geno,addEff,domEff,varG)
  intercept = tmp$output$intercept
  addEff = addEff*tmp$parameter
  domEff = domEff*tmp$parameter
  trait = new("TraitAD",
              qtlLoci,
              addEff=addEff,
              domEff=domEff,
              intercept=meanG-intercept)
  simParam@nTraits = simParam@nTraits + 1L
  simParam@traits[[simParam@nTraits]] = trait
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
#' @export
newPop = function(rawPop, id=NULL, simParam=SIMPARAM){
  stopifnot(class(rawPop)=="RawPop" | class(rawPop)=="MapPop")
  if(is.null(id)){
    id = (1:rawPop@nInd) + simParam@lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  #Loop through all traits
  gv = lapply(simParam@traits,getGv,pop=rawPop,w=0)
  gv = do.call("cbind",gv)
  output = new("Pop",
               nInd=rawPop@nInd,
               nChr=rawPop@nChr,
               ploidy=rawPop@ploidy,
               nLoci=rawPop@nLoci,
               gender=rawPop@gender,
               geno=rawPop@geno,
               id=as.character(id),
               mother=rep("0",rawPop@nInd),
               father=rep("0",rawPop@nInd),
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits))
  validObject(output)
  if(updateId){
    #Update simParam@lastId
    AlphaSimR:::assignInt(simParam@lastId,lastId)
  }
  return(output)
}


#ToDo----
addTraitAG = function(nQtlPerChr,meanG,varG,varGE,simParam=SIMPARAM,
                      founderPop=FOUNDERPOP){
  qtlLoci = pickQtlLoci(nQtlPerChr,simParam=simParam)
  
}

addTraitADG = function(nQtlPerChr,meanG,varG,domDegree,varGE,
                       simParam=SIMPARAM,founderPop=FOUNDERPOP){
  qtlLoci = pickQtlLoci(nQtlPerChr,simParam=simParam)
  
}

