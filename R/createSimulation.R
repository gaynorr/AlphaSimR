# Functions for starting a new simulation
# Requires InitialHaplo
# Returns a list containing SimParam, FounderPop, snpLoc, and qtlLoc

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
#' @param founderPop an object of class 'MapPop"
#' 
#' @export
createSimulation = function(maxQtl,maxSnp,snpQtlOverlap=FALSE,
                            minSnpFreq=NULL,founderPop = FOUNDERPOP){
  stopifnot(class(founderPop)=="MapPop")
  if(length(maxSnp)==1){
    maxSnp = rep(maxSnp,founderPop@nChr)
  }
  if(length(maxQtl)==1){
    maxQtl = rep(maxQtl,founderPop@nChr)
  }
  stopifnot(length(maxSnp)==nChr,length(maxQtl)==nChr)
  
  if(inbred){
    nInd = sum(InitialHaplo@nHaplo)
  }else{
    nInd = sum(InitialHaplo@nHaplo%/%ploidy)
  }
  genMaps = list()
  segSites = NULL
  geno = list()
  snpLoc = list()
  qtlLoc = list()
  for(chr in 1:nChr){
    genMaps[[chr]] = InitialHaplo@chrData[[chr]]$map
    segSites = c(segSites,length(genMaps[[chr]]))
    tmpGeno = list()
    for(i in 1:ploidy){
      if(inbred){
        tmpGeno[[i]] = InitialHaplo@chrData[[chr]]$geno
      }else{
        take = seq(from=i,to=nInd*ploidy,by=ploidy)
        tmpGeno[[i]] = InitialHaplo@chrData[[chr]]$geno[take,]
      }
    }
    geno[[chr]] = tmpGeno
    if(snpQtlOverlap){
      stopifnot(segSites[chr]>=maxSnp[chr],segSites[chr]>=maxQtl[chr])
      if(is.null(minSnpFreq)){
        snpLoc[[chr]] = sort(sample.int(segSites[chr],maxSnp[chr]))
      }else{
        q = calcQ2(geno[[chr]])
        snpLoc[[chr]] = sort(sample(which(q>=minSnpFreq),maxSnp[chr]))
      }
      qtlLoc[[chr]] = sort(sample.int(segSites[chr],maxQtl[chr]))
    }else{
      stopifnot(segSites[chr]>=sum(maxSnp[chr],maxQtl[chr]))
      if(is.null(minSnpFreq)){
        tmp = sample.int(segSites[chr],sum(maxSnp[chr],maxQtl[chr]))
        snpLoc[[chr]] = sort(tmp[1:maxSnp[chr]])
        qtlLoc[[chr]] = sort(tmp[(maxSnp[chr]+1):length(tmp)])
      }else{
        q = calcQ2(geno[[chr]])
        snpLoc[[chr]] = sort(sample(which(q>=minSnpFreq),maxSnp[chr]))
        qtlLoc[[chr]] = sort(sample(which(!((1:segSites[chr])%in%snpLoc[[chr]])),
                                    maxQtl[chr]))
      }
    }
  }
  output$FounderPop = new("Pop",
                          nInd=as.integer(nInd),
                          nChr=as.integer(nChr),
                          ploidy=as.integer(ploidy),
                          gender=rep("H",nInd),
                          geno=geno)
  output$SimParam = new("SimParam",
                        ploidy=as.integer(ploidy),
                        nChr=as.integer(nChr),
                        nTraits=0L,
                        nSnpChips=0L,
                        segSites=as.integer(segSites),
                        useGender=F,
                        genMaps=genMaps,
                        traits=list(),
                        snpChips=list())
  output$snpLoc = snpLoc
  output$qtlLoc = qtlLoc
  class(output) = "SimTempList"
  return(output)
}

#' @title Assign gender
#' 
#' @description Sets gender of founder individuals
#' 
#' @param simInfo an object of class 'SimTempList'
#' @param gender the gender to assign to each individual. Either 'H', 'M', or 'F'.
#' 
#' @export
assignGender = function(simInfo,gender){
  nInd = simInfo$FounderPop@nInd
  stopifnot(length(gender)==nInd)
  simInfo$FounderPop@gender = gender
  return(simInfo)
}

#' @title Add SNP chip
#' 
#' @description Randomly assigns eligble SNPs to a SNP chip
#' 
#' @param simInfo an object of class 'SimTempList'
#' @param nSnpPerChr number of SNPs per chromosome. Can be a single value or nChr values.
#' 
#' @export
addSnpChip = function(simInfo,nSnpPerChr){
  if(length(nSnpPerChr)==1){
    nSnpPerChr = rep(nSnpPerChr,simInfo$SimParam@nChr)
  }
  stopifnot(length(nSnpPerChr)==simInfo$SimParam@nChr)
  stopifnot(sapply(simInfo$snpLoc,length)>=nSnpPerChr)
  lociLoc = lapply(1:simInfo$SimParam@nChr,function(x){
    sort(sample(simInfo$snpLoc[[x]],nSnpPerChr[x]))
  })
  lociLoc = do.call("c",lociLoc)
  snpChip = new("LociMap",
                nLoci=as.integer(sum(nSnpPerChr)),
                lociPerChr=as.integer(nSnpPerChr),
                lociLoc=as.integer(lociLoc))
  simInfo$SimParam@nSnpChips = simInfo$SimParam@nSnpChips + 1L
  simInfo$SimParam@snpChips[[simInfo$SimParam@nSnpChips]] = snpChip
  return(simInfo)
}

#Function for selecting QTL loci called by all addTrait functions
pickQtlLoci = function(simInfo, nQtlPerChr){
  if(length(nQtlPerChr)==1){
    nQtlPerChr = rep(nQtlPerChr,simInfo$SimParam@nChr)
  }
  stopifnot(length(nQtlPerChr)==simInfo$SimParam@nChr)
  stopifnot(sapply(simInfo$qtlLoc,length)>=nQtlPerChr)
  lociLoc = lapply(1:simInfo$SimParam@nChr,function(x){
    sort(sample(simInfo$qtlLoc[[x]],nQtlPerChr[x]))
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
#' @description Randomly assigns eligble QTLs for an additive trait. 
#' 
#' @param simInfo an object of class 'SimTempList'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param meanG the mean genetic value for the trait
#' @param varG the total genetic variance for the trait
#' 
#' @export
addTraitA = function(simInfo,nQtlPerChr,meanG,varG){
  qtlLoci = pickQtlLoci(simInfo,nQtlPerChr)
  addEff = rnorm(qtlLoci@nLoci)
  geno = getGeno(pop = simInfo$FounderPop,
                 lociMap = qtlLoci)
  tmp = tuneTraitA(geno,addEff,varG)
  intercept = tmp$output$intercept
  addEff = addEff*tmp$parameter
  trait = new("TraitA",
              nLoci=qtlLoci@nLoci,
              lociPerChr=qtlLoci@lociPerChr,
              lociLoc=qtlLoci@lociLoc,
              addEff=addEff,
              intercept=meanG-intercept)
  simInfo$SimParam@nTraits = simInfo$SimParam@nTraits + 1L
  simInfo$SimParam@traits[[simInfo$SimParam@nTraits]] = trait
  return(simInfo)
}

#' @title Add an additive and dominance
#' 
#' @description Randomly assigns eligble QTLs for a trait with dominance. 
#' 
#' @param simInfo an object of class 'SimTempList'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param meanG the mean genetic value for the trait
#' @param varG the total genetic variance for the trait
#' @param domDegree the dominance degree of individual loci. Can be a single value or nLoci values.
#' 
#' @export
addTraitAD = function(simInfo,nQtlPerChr,meanG,varG,domDegree){
  qtlLoci = pickQtlLoci(simInfo,nQtlPerChr)
  if(length(domDegree)==1){
    domDegree = rep(domDegree,qtlLoci@nLoci)
  }else{
   stopifnot(length(domDegree)==qtlLoci@nLoci) 
  }
  addEff = rnorm(qtlLoci@nLoci)
  domEff = abs(addEff)*domDegree
  geno = getGeno(pop = simInfo$FounderPop,
                 lociMap = qtlLoci)
  if(simInfo$inbred){
    tmp = tuneTraitA(geno,addEff,varG)
  }else{
    tmp = tuneTraitAD(geno,addEff,domEff,varG)
  }
  intercept = tmp$output$intercept
  addEff = addEff*tmp$parameter
  domEff = domEff*tmp$parameter
  trait = new("TraitAD",
              nLoci=qtlLoci@nLoci,
              lociPerChr=qtlLoci@lociPerChr,
              lociLoc=qtlLoci@lociLoc,
              addEff=addEff,
              domEff=domEff,
              intercept=meanG-intercept)
  simInfo$SimParam@nTraits = simInfo$SimParam@nTraits + 1L
  simInfo$SimParam@traits[[simInfo$SimParam@nTraits]] = trait
  return(simInfo)
}

addTraitAG = function(simInfo,nQtlPerChr,meanG,varG,varGE){
  qtlLoci = pickQtlLoci(simInfo,nQtlPerChr)
  
}

addTraitADG = function(simInfo,nQtlPerChr,meanG,varG,domDegree,varGE){
  qtlLoci = pickQtlLoci(simInfo,nQtlPerChr)
  
}

#' @title Get founder individuals
#' 
#' @description Returns an object of class 'Pop' or 'TraitPop' from a 'SimTempList' class.
#' 
#' @param simInfo an object of class 'SimTempList'
#' 
#' @export
getFounders = function(simInfo){
  pop = simInfo$FounderPop
  if(simInfo$SimParam@nTraits==0){
    return(pop)
  }
  pop = addGv(pop,simInfo$SimParam)
  return(pop)
}

#' @title Get SimParam object
#' 
#' @description Returns an object of 'SimParam' class from a 'SimTempList' class.
#' 
#' @param simInfo an object of class 'SimTempList'
#' 
#' @export
getSimParam = function(simInfo){
  output = simInfo$SimParam
  validObject(output)
  return(output)
}
