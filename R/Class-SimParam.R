# SimParam ----
#' @title Simulation parameters
#' 
#' @description 
#' Container for global simulation parameters. Saving this object 
#' as SP will allow it to be accessed by function defaults.
#' 
#' @field ploidy ploidy level of species
#' @field nChr number of chromosomes
#' @field nTraits number of traits
#' @field nSnpChips number of SNP chips
#' @field segSites segregating sites per chromosome
#' @field gender is gender used for mating
#' @field genMap "matrix" of chromsome genetic maps
#' @field femaleMap "matrix" of chromsome genetic maps for 
#' females
#' @field maleMap "matrix" of chromsome genetic maps for 
#' males
#' @field sepMap are there seperate genetic maps for 
#' males and females
#' @field recombRatio ratio of genetic recombination in 
#' females relative to male
#' @field traits list of trait
#' @field snpChips list of SNP chips
#' @field potQtl list of potential QTL segregating sites
#' @field potSnp list of potential SNP segregating sites
#' @field lastId last ID number assigned
#' @field isTrackPed is pedigree being tracked
#' @field pedigree pedigree matrix for all individuals
#' @field isTrackRec is recombination being tracked
#' @field recHist list of historic recombination events
#' @field varA additive genetic variance in founderPop
#' @field varG total genetic variance in founderPop
#' @field varE default error variance
#' @field founderPop the founder population used for scaling traits
#'
#' @export
SimParam = R6Class(
  "SimParam",
  public = list(),
  private = list(
    .ploidy="integer",
    .nChr="integer",
    .nTraits="integer",
    .nSnpChips="integer",
    .segSites="integer",
    .gender="character",
    .femaleMap="matrix",
    .maleMap="matrix",
    .sepMap="logical",
    .recombRatio="numeric",
    .traits="list",
    .snpChips="list",
    .potQtl="list",
    .potSnp="list",
    .lastId="integer",
    .isTrackPed="logical",
    .pedigree="matrix",
    .isTrackRec="logical",
    .recHist="list",
    .varA="numeric",
    .varG="numeric",
    .varE="numeric",
    .founderPop="MapPop"
  ),
  active = list(
    ploidy=function(value){
      if(missing(value)){
        private$.ploidy
      }else{
        stop("`$ploidy` is read only",call.=FALSE)
      }
    },
    nChr=function(value){
      if(missing(value)){
        private$.nChr
      }else{
        stop("`$nChr` is read only",call.=FALSE)
      }
    },
    nTraits=function(value){
      if(missing(value)){
        private$.nTraits
      }else{
        stop("`$nTraits` is read only",call.=FALSE)
      }
    },
    nSnpChips=function(value){
      if(missing(value)){
        private$.nSnpChips
      }else{
        stop("`$nSnpChips` is read only",call.=FALSE)
      }
    },
    segSites=function(value){
      if(missing(value)){
        private$.segSites
      }else{
        stop("`$segSites` is read only",call.=FALSE)
      }
    },
    gender=function(value){
      if(missing(value)){
        private$.gender
      }else{
        stop("`$gender` is read only",call.=FALSE)
      }
    },
    sepMap=function(value){
      if(missing(value)){
        private$.sepMap
      }else{
        stop("`$sepMap` is read only",call.=FALSE)
      }
    },
    genMap=function(value){
      if(missing(value)){
        if(private$.sepMap){
          genMap = vector("list",private$.nChr)
          for(i in 1:private$.nChr){
            genMap[[i]] = (private$.femaleMap[[i]]+
              private$.maleMap[[i]])/2
          }
          as.matrix(genMap)
        }else{
          private$.femaleMap
        }
      }else{
        stop("`$genMap` is read only",call.=FALSE)
      }
    },
    femaleMap=function(value){
      if(missing(value)){
        private$.femaleMap
      }else{
        stop("`$femaleMap` is read only",call.=FALSE)
      }
    },
    maleMap=function(value){
      if(missing(value)){
        if(private$.sepMap){
          private$.maleMap
        }else{
          private$.femaleMap
        }
      }else{
        stop("`$maleMap` is read only",call.=FALSE)
      }
    },
    traits=function(value){
      if(missing(value)){
        private$.traits
      }else{
        stop("`$traits` is read only",call.=FALSE)
      }
    },
    snpChips=function(value){
      if(missing(value)){
        private$.snpChips
      }else{
        stop("`$snpChips` is read only",call.=FALSE)
      }
    },
    potQtl=function(value){
      if(missing(value)){
        private$.potQtl
      }else{
        stop("`$potQtl` is read only",call.=FALSE)
      }
    },
    potSnp=function(value){
      if(missing(value)){
        private$.potSnp
      }else{
        stop("`$potSnp` is read only",call.=FALSE)
      }
    },
    lastId=function(value){
      if(missing(value)){
        private$.lastId
      }else{
        stop("`$lastId` is read only",call.=FALSE)
      }
    },
    isTrackPed=function(value){
      if(missing(value)){
        private$.isTrackPed
      }else{
        stop("`$isTrackPed` is read only",call.=FALSE)
      }
    },
    pedigree=function(value){
      if(missing(value)){
        private$.pedigree
      }else{
        stop("`$pedigree` is read only",call.=FALSE)
      }
    },
    isTrackRec=function(value){
      if(missing(value)){
        private$.isTrackRec
      }else{
        stop("`$isTrackRec` is read only",call.=FALSE)
      }
    },
    recHist=function(value){
      if(missing(value)){
        private$.recHist
      }else{
        stop("`$recHist` is read only",call.=FALSE)
      }
    },
    varA=function(value){
      if(missing(value)){
        private$.varA
      }else{
        stop("`$varA` is read only",call.=FALSE)
      }
    },
    varG=function(value){
      if(missing(value)){
        private$.varG
      }else{
        stop("`$varG` is read only",call.=FALSE)
      }
    },
    varE=function(value){
      if(missing(value)){
        private$.varE
      }else{
        stop("`$varE` is read only",call.=FALSE)
      }
    },
    founderPop=function(value){
      if(missing(value)){
        private$.founderPop
      }else{
        stop("`$founderPop` is read only",call.=FALSE)
      }
    }
  ),
  cloneable=FALSE
)

#' @title Create new simulation
#'
#' @description Starts the process of building a new simulation 
#' by creating a new SimParam object and assigning a founder 
#' population to the class. It is recommended that you save the 
#' object with the name "SP", because subsequent functions will 
#' check your global enviroment for an object of this name if 
#' their simParam arguments are NULL. This allows you to call 
#' these functions without explicitly supplying a simParam 
#' argument with every call.
#' 
#' @section Usage: SimParam$new(founderPop)
#' 
#' @param founderPop an object of \code{\link{MapPop-class}}
#' 
#' @name SimParam_new
NULL
# new ----
SimParam$set(
  "public",
  "initialize",
  function(founderPop){
    stopifnot(class(founderPop)=="MapPop")
    private$.ploidy = founderPop@ploidy
    private$.nChr = founderPop@nChr
    private$.nTraits = 0L
    private$.nSnpChips = 0L
    private$.segSites = founderPop@nLoci
    private$.gender = "no"
    private$.femaleMap = founderPop@genMap
    private$.maleMap = NULL
    private$.sepMap = FALSE
    private$.traits = list()
    private$.snpChips = list()
    private$.potQtl = lapply(
      founderPop@nLoci,
      function(x) 1:x
    )
    private$.potSnp = private$.potQtl
    private$.lastId = 0L
    private$.isTrackPed = FALSE
    private$.pedigree = matrix(NA_integer_,nrow=0,ncol=3)
    private$.isTrackRec = FALSE
    private$.recHist = list()
    private$.varA = numeric()
    private$.varG = numeric()
    private$.varE = numeric()
    private$.founderPop = founderPop
    invisible(self)
  }
)

# .isRunning ----
SimParam$set(
  "private",
  ".isRunning",
  function(){
    if(private$.lastId==0L){
      invisible(self)
    }else{
      stop("lastId doesn't equal 0, you must run resetPed to proceed")
    }
  }
)

# .addTrait ----
SimParam$set(
  "private",
  ".addTrait",
  function(lociMap,varA=NA_real_,varG=NA_real_){
    private$.nTraits = private$.nTraits+1L
    private$.traits[[private$.nTraits]] = lociMap
    private$.varA[private$.nTraits] = varA
    private$.varG[private$.nTraits] = varG
    private$.varE[private$.nTraits] = NA_real_
    invisible(self)
  }
)

# updateLastId ----
SimParam$set(
  "public",
  "updateLastId",
  function(lastId){
    lastId = as.integer(lastId)
    stopifnot(lastId>=private$.lastId)
    private$.lastId = lastId
    invisible(self)
  }
)

#' @title Set pedigree tracking
#'
#' @description Sets pedigree tracking for the simulation. 
#' By default pedigree tracking is turned off. When turned on, 
#' the pedigree of all individuals created will be tracked, 
#' except those created by \code{\link{hybridCross}}. Turning 
#' off pedigree tracking will turn off recombination tracking 
#' if it is turned on.
#' 
#' @section Usage: SP$setTrackPed(isTrackPed, force = FALSE)
#' 
#' @param isTrackPed should pedigree tracking be on.
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_setTrackPed
NULL
# setTrackPed ----
SimParam$set(
  "public",
  "setTrackPed",
  function(isTrackPed, force=FALSE){
    stopifnot(is.logical(isTrackPed))
    if(!force){
      private$.isRunning()
    }
    private$.isTrackPed = isTrackPed
    if(!isTrackPed){
      private$.isTrackRec = FALSE
    }
    invisible(self)
  }
)

#' @title Set recombination tracking
#'
#' @description Sets recombination tracking for the simulation. 
#' By default recombination tracking is turned off. When turned 
#' on recombination tracking will also turn on pedigree tracking. 
#' Recombination tracking keeps records of all individuals created, 
#' except those created by \code{\link{hybridCross}}, because their 
#' pedigree is not tracked.
#' 
#' @section Usage: SimParam$setTrackRec(isTrackRec, force = FALSE)
#' 
#' @param isTrackRec should recombination tracking be on.
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_setTrackRec
NULL
# setTrackRec ----
SimParam$set(
  "public",
  "setTrackRec",
  function(isTrackRec, force=FALSE){
    stopifnot(is.logical(isTrackRec))
    if(!force){
      private$.isRunning()
    }
    private$.isTrackRec = isTrackRec
    if(isTrackRec){
      private$.isTrackPed = TRUE
    }
    invisible(self)
  }
)

#' @title Reset pedigree
#'
#' @description Resets the internal lastId, the pedigree 
#' and recombination tracking, if it is being used, to the 
#' supplied lastId. Be careful using this function because 
#' it may introduce bug if you supsequently use individuals 
#' that come from a portion the pedigree that is being reset.
#' 
#' @param lastId last ID to include in pedigree
#' 
#' @section Usage: SP$resetPed(lastId = 0L)
#' 
#' @name SimParam_resetPed
NULL
# resetPed ----
SimParam$set(
  "public",
  "resetPed",
  function(lastId=0L){
    private$.lastId = lastId
    private$.pedigree = private$.pedigree[0:lastId,,drop=FALSE]
    if(private$.isTrackRec){
      private$.recHist = private$.recHist[0:lastId]
    }
    invisible(self)
  }
)

#' @title Restrict segregating sites
#'
#' @description Sets restrictions on which segregating sites 
#' can serve and SNP and/or QTL loci.
#' 
#' @section Usage: SP$restrSegSites(maxQtl = 0, maxSnp = 0, 
#' snpQtlOverlap = FALSE, minSnpFreq = NULL, force = FALSE)
#' 
#' @param maxQtl the maximum number of segSites for QTLs. 
#' Can be a single value or a vector values for each 
#' chromosome.
#' @param maxSnp the maximum number of segSites for SNPs. 
#' Can be a single value or a vector values for each 
#' chromosome.
#' @param snpQtlOverlap should SNP and QTL loci be allowed 
#' to overlap.
#' @param minSnpFreq minimum allowable frequency for SNP loci. 
#' No minimum SNP frequency is used if value is NULL.
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_restrSegSites
NULL
# restrSegSites ----
SimParam$set(
  "public",
  "restrSegSites",
  function(maxQtl=0,maxSnp=0,snpQtlOverlap=FALSE,
           minSnpFreq=NULL, force=FALSE){
    if(!force){
      private$.isRunning()
    }
    if(length(maxSnp)==1){
      maxSnp = rep(maxSnp,private$.nChr)
    }
    if(length(maxQtl)==1){
      maxQtl = rep(maxQtl,private$.nChr)
    }
    stopifnot(length(maxSnp)==private$.nChr,
              length(maxQtl)==private$.nChr)
    potSnp = list()
    potQtl = list()
    for(chr in 1:private$.nChr){
      if(snpQtlOverlap){
        stopifnot(private$.segSites[chr]>=maxSnp[chr],
                  private$.segSites[chr]>=maxQtl[chr])
        if(is.null(minSnpFreq)){
          potSnp[[chr]] = sort(sample.int(private$.segSites[chr],
                                          maxSnp[chr]))
        }else{
          q = calcChrFreq(private$.founderPop@geno[[chr]])
          q = 0.5-abs(q-0.5) #Convert to minor allele frequency
          potSnp[[chr]] = sort(sample(which(q>=minSnpFreq),maxSnp[chr]))
        }
        potQtl[[chr]] = sort(sample.int(private$.segSites[chr],
                                        maxQtl[chr]))
      }else{
        stopifnot(private$.segSites[chr]>=sum(maxSnp[chr],maxQtl[chr]))
        if(is.null(minSnpFreq)){
          tmp = sample.int(private$.segSites[chr],sum(maxSnp[chr],maxQtl[chr]))
          if(maxSnp[chr]>0){
            potSnp[[chr]] = sort(tmp[1:maxSnp[chr]])
          }else{
            potSnp[[chr]] = numeric()
          }
          if(maxQtl[chr]>0){
            potQtl[[chr]] = sort(tmp[(maxSnp[chr]+1):length(tmp)])
          }else{
            maxQtl[[chr]] = numeric()
          }
        }else{
          q = calcChrFreq(private$.founderPop@geno[[chr]])
          q = 0.5-abs(q-0.5)
          potSnp[[chr]] = sort(sample(which(q>=minSnpFreq),maxSnp[chr]))
          potQtl[[chr]] = sort(sample(which(!((1:private$.segSites[chr])%in%potSnp[[chr]])),
                                      maxQtl[chr]))
        }
      }
    }
    private$.potSnp = potSnp
    private$.potQtl = potQtl
    invisible(self)
  }
)

#' @title Set gender in simulation
#'
#' @description Changes how gender is used in the simulation. 
#' The default gender of a simulation is "no". To add gender 
#' to the simulation, run this function with "yes_sys" or 
#' "yes_rand". The value "yes_sys" will systematically assign 
#' gender to newly created individuals as first male, then female. 
#' Thus, odd numbers of individuals will have one more male than 
#' female. The value "yes_rand" will randomly assign gender to 
#' individuals.
#' 
#' @section Usage: SP$setGender(gender, force = FALSE)
#' 
#' @param gender acceptable value are "no", "yes_sys", or 
#' "yes_rand"
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_setGender
NULL
# setGender ----
SimParam$set(
  "public",
  "setGender",
  function(gender, force=FALSE){
    if(!force){
      private$.isRunning()
    }
    gender = tolower(gender)
    if(gender=="no"){
      private$.gender="no"
    }else if(gender=="yes_sys"){
      private$.gender="yes_sys"
    }else if(gender=="yes_rand"){
      private$.gender="yes_rand"
    }else{
      stop(paste0("gender=",gender," is not a valid option"))
    }
    invisible(self)
  }
)

#' @title Set gender specific recombination
#'
#' @description Defines a gender specific recombination ratio. 
#' The ratio is defined as the amount of recombination in 
#' females relative to male. Thus, a value of 1 (default) 
#' specifies equal recombination rates in both males and females. 
#' A value of 2 specifies twice as much recombination in females 
#' and a value of 0.5 specifies half as much recombination in 
#' females.
#' 
#' @section Usage: SP$setRecRatio(ratio, force = FALSE)
#' 
#' @param ratio any value greater than 0
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_setRecRatio
NULL
# setRecRatio ----
SimParam$set(
  "public",
  "setRecRatio",
  function(ratio, force=FALSE){
    if(!force){
      private$.isRunning()
    }
    stopifnot(ratio>0)
    genMap = self$genMap
    private$.sepMap = TRUE
    feSc = 2/(1/ratio+1)
    maSc = 2/(ratio+1)
    private$.femaleMap = as.matrix(
      lapply(genMap,
             function(x){
               feSc*x
             })
    )
    private$.maleMap = as.matrix(
      lapply(genMap,
             function(x){
               maSc*x
             })
    )
    invisible(self)
  }
)

#' @title Set simulation error variance
#'
#' @description Defines a default value for error 
#' variances in the simulation.
#' 
#' @section Usage: SP$setVarE(h2 = NULL, H2 = NULL, varE = NULL)
#' 
#' @param h2 a vector of desired narrow-sense heritabilities
#' @param H2 a vector of desired broad-sense heritabilities
#' @param varE a vector of error variances
#' 
#' @name SimParam_setVarE
NULL
# setVarE ----
SimParam$set(
  "public",
  "setVarE",
  function(h2=NULL,H2=NULL,varE=NULL){
    if(!is.null(h2)){
      stopifnot(length(h2)==private$.nTraits,
                all(private$.varG>0),
                all(private$.varA>0))
      varE = numeric(private$.nTraits)
      for(i in 1:length(h2)){
        tmp = private$.varA[i]/h2[i]-private$.varG[i]
        if(tmp<0){
          stop(paste0("h2=",h2[i]," is not possible for trait ",i))
        }
        varE[i] = tmp
      }
      private$.varE = varE
    }else if(!is.null(H2)){
      stopifnot(length(H2)==private$.nTraits)
      varE = numeric(private$.nTraits)
      for(i in 1:length(h2)){
        tmp = private$.varG[i]/H2[i]-private$.varG[i]
        varE[i] = tmp
      }
      private$.varE = varE
    }else if(!is.null(varE)){
      stopifnot(length(varE)==private$.nTraits)
      private$.varE = varE
    }else{
      private$.varE = rep(NA_real_,private$.nTraits)
    }
    invisible(self)
  }
)

#' @title Set correlated error variance
#'
#' @description Defines a correlation structure for default 
#' error variances. You must call \code{\link{SimParam_setVarE}} 
#' first to define the default error variances.
#' 
#' @section Usage: SP$setCorrVarE(corr)
#' 
#' @param corr a correlation matrix for the error variances
#' 
#' @name SimParam_setCorrVarE
NULL
# setCorrVarE ----
SimParam$set(
  "public",
  "setCorrVarE",
  function(corr){
    stopifnot(isSymmetric(corr),
              nrow(corr)==private$.nTraits)
    varE = diag(sqrt(private$.varE),
                nrow=private$.nTraits,
                ncol=private$.nTraits)
    varE = varE%*%corr%*%varE
    private$.varE = varE
    invisible(self)
  }
)


#' @title Add SNP chip
#' 
#' @description 
#' Randomly assigns eligble SNPs to a SNP chip
#' 
#' @section Usage: SP$addSnpChip(nSnpPerChr, force = FALSE)
#' 
#' @param nSnpPerChr number of SNPs per chromosome. 
#' Can be a single value or nChr values.
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_addSnpChip
NULL
# addSnpChip ----
SimParam$set(
  "public",
  "addSnpChip",
  function(nSnpPerChr, force=FALSE){
    if(!force){
      private$.isRunning()
    }
    if(length(nSnpPerChr)==1){
      nSnpPerChr = rep(nSnpPerChr,private$.nChr)
    }
    stopifnot(length(nSnpPerChr)==private$.nChr)
    stopifnot(sapply(private$.potSnp,length)>=nSnpPerChr)
    lociLoc = lapply(1:private$.nChr,function(x){
      sort(sample(private$.potSnp[[x]],nSnpPerChr[x]))
    })
    lociLoc = do.call("c",lociLoc)
    snpChip = new("LociMap",
                  nLoci=as.integer(sum(nSnpPerChr)),
                  lociPerChr=as.integer(nSnpPerChr),
                  lociLoc=as.integer(lociLoc))
    private$.nSnpChips = private$.nSnpChips + 1L
    private$.snpChips[[private$.nSnpChips]] = snpChip
    invisible(self)
  }
)

#' @title Add Structured SNP chips
#' 
#' @description 
#' Randomly selects the number of snps in structure and then
#' assigns them to chips based on structure
#' 
#' @section Usage: SP$addStructuredSnpChip(nSnpPerChr, structure, force = FALSE)
#' 
#' @param nSnpPerChr number of SNPs per chromosome. 
#' Can be a single value or nChr values.
#' @param structure a matrix.  Rows are snp chips, columns are chips.
#' If value is true then that snp is on that chip.
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_addStructuredSnpChips
NULL
# addStructuredSnpChip ----
SimParam$set(
  "public",
  "addStructuredSnpChip",
  function(nSnpPerChr,structure,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    if(length(nSnpPerChr)==1){
      nSnpPerChr = rep(nSnpPerChr,private$.nChr)
    }
    stopifnot(length(nSnpPerChr)==private$.nChr)
    stopifnot(sapply(simParam$potSnp,length)>=nSnpPerChr)
    stopifnot(dim(structure)[2]==sum(nSnpPerChr))
    lociLoc = lapply(1:private$.nChr,function(x){
      sort(sample(private$.potSnp[[x]],nSnpPerChr[x]))
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
      private$.nSnpChips = private$.nSnpChips + 1L
      private$.snpChips[[private$.nSnpChips]] = snpChip
    }
    invisible(self)
  }
)

#' @title Remove SNP chip
#' 
#' @description 
#' Removes designated SNP chip(s).
#' 
#' @section Usage: SP$removeSnpChip(chips, force = FALSE)
#' 
#' @param chips a vector of SNP chips to remove
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_removeSnpChip
NULL
# removeSnpChip ----
SimParam$set(
  "public",
  "removeSnpChip",
  function(chips, force=FALSE){
    if(!force){
      private$.isRunning()
    }
    chips = as.integer(chips)
    stopifnot(max(chips)<=private$.nSnpChips)
    private$.snpChips = private$.snpChips[-chips]
    private$.nSnpChips = length(private$.snpChips)
    invisible(self)
  }
)

#' @title Switch SNP chip
#' 
#' @description 
#' Replaces the \code{\link{LociMap-class}} for a SNP chip.
#' 
#' @section Usage: SP$switchSnpChips(lociMap, chip, force = FALSE)
#' 
#' @param lociMap a new \code{\link{LociMap-class}}
#' @param chip an integer indicating which chip to replace
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_switchSnpChip
NULL
# switchSnpChip ----
SimParam$set(
  "public",
  "switchSnpChip",
  function(lociMap,chip,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    stopifnot(length(chip)==1,chip<=private$.nSnpChips)
    private$.snpChips[[chip]] = lociMap
    invisible(self)
  }
)

#' @title Manually add SNP chip
#' 
#' @description 
#' Adds a new \code{\link{LociMap-class}} for a SNP chip.
#' 
#' @section Usage: SP$manAddSnpChips(lociMap, force = FALSE)
#' 
#' @param lociMap a new \code{\link{LociMap-class}}
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_manAddSnpChip
NULL
# manAddSnpChip ----
SimParam$set(
  "public",
  "manAddSnpChip",
  function(lociMap,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    private$.nSnpChips = private$.nSnpChips+1L
    private$.snpChips[[private$.nSnpChips]] = lociMap
    invisible(self)
  }
)

# pickQtlLoci ----
SimParam$set(
  "private",
  ".pickQtlLoci",
  function(nQtlPerChr){
    if(length(nQtlPerChr)==1){
      nQtlPerChr = rep(nQtlPerChr,private$.nChr)
    }
    stopifnot(length(nQtlPerChr)==private$.nChr)
    stopifnot(sapply(private$.potQtl,length)>=nQtlPerChr)
    lociLoc = lapply(1:private$.nChr,function(x){
      if(nQtlPerChr[x]==0){
        return(NULL)
      }else{
        return(sort(sample(private$.potQtl[[x]],nQtlPerChr[x])))
      }
    })
    lociLoc = do.call("c",lociLoc)
    qtlLoci = new("LociMap",
                  nLoci=as.integer(sum(nQtlPerChr)),
                  lociPerChr=as.integer(nQtlPerChr),
                  lociLoc=as.integer(lociLoc))
    return(qtlLoci)
  }
)

sampAddEff = function(qtlLoci,nTraits,corr,gamma,shape){
  addEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%chol(corr)
  if(any(gamma)){
    for(i in which(gamma)){
      x = (pnorm(addEff[,i])-0.5)*2
      addEff[,i] = sign(x)*qgamma(abs(x),shape=shape)
    }
  }
  return(addEff)
}

sampDomEff = function(qtlLoci,nTraits,addEff,corDD,
                      meanDD,varDD){
  domEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%chol(corDD)
  domEff = sweep(domEff,2,sqrt(varDD),"*")
  domEff = sweep(domEff,2,meanDD,"+")
  domEff = abs(addEff)*domEff
  return(domEff)
}

#' @title Add additive traits
#' 
#' @description 
#' Randomly assigns eligble QTLs for one ore more additive traits. 
#' If simulating more than one trait, all traits will be pleiotrophic 
#' with correlated additive effects.
#' 
#' @section Usage: SP$addTraitA(nQtlPerChr, mean = 0, var = 1, corr = NULL, 
#' gamma = FALSE, shape = 1, force = FALSE)
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of desired mean genetic values for one or more traits
#' @param var a vector of desired genetic variances for one or more traits
#' @param corr a matrix of correlations between additive effects
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_addTraitA
NULL
# addTraitA ----
SimParam$set(
  "public",
  "addTraitA",
  function(nQtlPerChr,mean=0,var=1,corr=NULL,
           gamma=FALSE,shape=1,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    nTraits = length(mean)
    if(length(gamma)==1) gamma = rep(gamma,nTraits)
    if(length(shape)==1) shape = rep(shape,nTraits)
    if(is.null(corr)) corr=diag(nTraits)
    stopifnot(length(mean)==length(var),
              isSymmetric(corr),
              length(mean)==nrow(corr))
    qtlLoci = private$.pickQtlLoci(nQtlPerChr)
    addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                        corr=corr,gamma=gamma,shape=shape)
    geno = getGeno(private$.founderPop@geno,
                   qtlLoci@lociPerChr,
                   qtlLoci@lociLoc)
    for(i in 1:nTraits){
      tmp = tuneTraitA(geno,addEff[,i],var[i])
      intercept = tmp$output$intercept
      addEff[,i] = addEff[,i]*tmp$parameter
      trait = new("TraitA",
                  qtlLoci,
                  addEff=addEff[,i],
                  intercept=mean[i]-intercept)
      private$.addTrait(trait,tmp$output$varA,tmp$output$varG)
    }
    invisible(self)
  }
)

#' @title Add additive and dominance traits
#' 
#' @description 
#' Randomly assigns eligble QTLs for one or more traits with dominance. 
#' If simulating more than one trait, all traits will be pleiotrophic 
#' with correlated additive effects.
#' 
#' @section Usage: SP$addTraitAD(nQtlPerChr, mean = 0, var = 1, meanDD = 0, 
#' varDD = 0, corA = NULL, corDD = NULL, useVarA = TRUE, gamma = FALSE, 
#' shape = 1, force = FALSE)
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of desired mean genetic values for one or more traits
#' @param var a vector of desired genetic variances for one or more traits
#' @param meanDD mean dominance degree
#' @param varDD variance of dominance degree
#' @param corA a matrix of correlations between additive effects
#' @param corDD a matrix of correlations between dominance degrees
#' @param useVarA tune according to additive genetic variance if true. If 
#' FALSE, tuning is performed according to total genetic variance.
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#'  
#' @name SimParam_addTraitAD
NULL
# addTraitAD ----
SimParam$set(
  "public",
  "addTraitAD",
  function(nQtlPerChr,mean=0,var=1,meanDD=0,
           varDD=0,corA=NULL,corDD=NULL,useVarA=TRUE,
           gamma=FALSE,shape=1,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    nTraits = length(mean)
    if(length(gamma)==1) gamma = rep(gamma,nTraits)
    if(length(shape)==1) shape = rep(shape,nTraits)
    if(is.null(corA)) corA=diag(nTraits)
    if(is.null(corDD)) corDD=diag(nTraits)
    stopifnot(length(mean)==length(var),
              isSymmetric(corA),
              isSymmetric(corDD),
              length(mean)==nrow(corA))
    qtlLoci = private$.pickQtlLoci(nQtlPerChr)
    addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                        corr=corA,gamma=gamma,shape=shape)
    domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                        corDD=corDD,meanDD=meanDD,varDD=varDD)
    geno = getGeno(private$.founderPop@geno,
                   qtlLoci@lociPerChr,
                   qtlLoci@lociLoc)
    for(i in 1:nTraits){
      tmp = tuneTraitAD(geno,addEff[,i],domEff[,i],var[i],useVarA)
      intercept = tmp$output$intercept
      addEff[,i] = addEff[,i]*tmp$parameter
      domEff[,i] = domEff[,i]*tmp$parameter
      trait = new("TraitAD",
                  qtlLoci,
                  addEff=addEff[,i],
                  domEff=domEff[,i],
                  intercept=mean[i]-intercept)
      private$.addTrait(trait,tmp$output$varA,tmp$output$varG)
    }
    invisible(self)
  }
)

#' @title Add additive GxE traits
#' 
#' @description 
#' Randomly assigns eligble QTLs for one ore more additive GxE traits. 
#' If simulating more than one trait, all traits will be pleiotrophic 
#' with correlated effects.
#' 
#' @section Usage: SP$addTraitAG(nQtlPerChr, mean = 0, var = 1, varEnv = 1e-6, 
#' varGxE = 1e-6, corA = NULL, corGxE = NULL, gamma = FALSE, shape = 1)
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of desired mean genetic values for one or more traits
#' @param var a vector of desired genetic variances for one or more traits
#' @param varEnv a vector of environmental variances for one or more traits
#' @param varGxE a vector of total genotype-by-environment variances for the traits
#' @param corA a matrix of correlations between additive effects
#' @param corGxE a matrix of correlations between GxE effects
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#' 
#' @name SimParam_addTraitAG
NULL
# addTraitAG ----
SimParam$set(
  "public",
  "addTraitAG",
  function(nQtlPerChr,mean=0,var=1,varEnv=1e-6,
           varGxE=1e-6,corA=NULL,corGxE=NULL,gamma=FALSE,
           shape=1){
    private$.isRunning()
    nTraits = length(mean)
    if(length(gamma)==1) gamma = rep(gamma,nTraits)
    if(length(shape)==1) shape = rep(shape,nTraits)
    if(is.null(corA)) corA=diag(nTraits)
    if(is.null(corGxE)) corGxE=diag(nTraits)
    stopifnot(length(mean)==length(var),
              isSymmetric(corA),
              isSymmetric(corGxE),
              length(meanG)==nrow(corA),
              length(meanG)==nrow(corGxE),
              length(meanG)==length(varGxE),
              length(meanG)==length(varEnv))
    qtlLoci = private$.pickQtlLoci(nQtlPerChr)
    addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                        corr=corA,gamma=gamma,shape=shape)
    gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                        corr=corGxE,gamma=FALSE,shape=NULL)
    geno = getGeno(private$.founderPop@geno,
                   qtlLoci@lociPerChr,
                   qtlLoci@lociLoc)
    for(i in 1:nTraits){
      tmp = tuneTraitA(geno,addEff[,i],var[i])
      intercept = tmp$output$intercept
      addEff[,i] = addEff[,i]*tmp$parameter
      targetVar = varGxE[i]/varEnv[i]
      tmp = tuneTraitA(geno,gxeEff[,i],targetVar)
      gxeEff[,i] = gxeEff[,i]*tmp$parameter
      gxeInt = tmp$output$intercept
      trait = new("TraitAG",
                  qtlLoci,
                  addEff=addEff[,i],
                  intercept=mean[i]-intercept,
                  gxeEff = gxeEff[,i],
                  gxeInt = 1-gxeInt,
                  envVar = varEnv[i])
      private$.addTrait(trait,tmp$output$varA,tmp$output$varG)
    }
    invisible(self)
  }
)

#' @title Add an additive and dominance GxE trait
#' 
#' @description 
#' Randomly assigns eligble QTLs for a trait with dominance and GxE. 
#' 
#' @section Usage: SP$addTraitAG(nQtlPerChr, mean = 0, var = 1, varEnv = 1e-6, 
#' varGxE = 1e-6, meanDD = 0, varDD = 0, corA = NULL, corDD = NULL, 
#' corGxE = NULL, useVarA = TRUE, gamma = FALSE, shape = 1, force = FALSE)
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single 
#' value or nChr values.
#' @param mean a vector of desired mean genetic values for one or more traits
#' @param var a vector of desired genetic variances for one or more traits
#' @param varEnv a vector of environmental variances for one or more traits
#' @param varGxE a vector of total genotype-by-environment variances for the traits
#' @param meanDD mean dominance degree
#' @param varDD variance of dominance degree
#' @param corA a matrix of correlations between additive effects
#' @param corDD a matrix of correlations between dominance degrees
#' @param corGxE a matrix of correlations between GxE effects
#' @param useVarA tune according to additive genetic variance if true
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#'  
#' @name SimParam_addTraitADG
NULL
# addTraitADG ----
SimParam$set(
  "public",
  "addTraitADG",
  function(nQtlPerChr,mean=0,var=1,varEnv=1e-6,
           varGxE=1e-6,meanDD=0,varDD=0,corA=NULL,
           corDD=NULL,corGxE=NULL,useVarA=TRUE,gamma=FALSE,
           shape=1,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    nTraits = length(mean)
    if(length(gamma)==1) gamma = rep(gamma,nTraits)
    if(length(shape)==1) shape = rep(shape,nTraits)
    if(is.null(corA)) corA=diag(nTraits)
    if(is.null(corDD)) corDD=diag(nTraits)
    if(is.null(corGxE)) corGxE=diag(nTraits)
    stopifnot(length(mean)==length(var),
              isSymmetric(corA),
              isSymmetric(corDD),
              isSymmetric(corGxE),
              nrow(corA)==nTraits,
              nrow(corGxE)==nTraits,
              nrow(corDD)==nTraits,
              length(varGxE)==nTraits,
              length(varEnv)==nTraits,
              length(varDD)==nTraits)
    qtlLoci = private$.pickQtlLoci(nQtlPerChr)
    addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                        corr=corA,gamma=gamma,shape=shape)
    domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                        corDD=corDD,meanDD=meanDD,varDD=varDD)
    gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                        corr=corGxE,gamma=FALSE,shape=NULL)
    geno = getGeno(private$.founderPop@geno,
                   qtlLoci@lociPerChr,
                   qtlLoci@lociLoc)
    for(i in 1:nTraits){
      tmp = tuneTraitAD(geno,addEff[,i],domEff[,i],var[i],useVarA)
      intercept = tmp$output$intercept
      addEff[,i] = addEff[,i]*tmp$parameter
      domEff[,i] = domEff[,i]*tmp$parameter
      targetVar = varGxE[i]/varEnv[i]
      tmp = tuneTraitA(geno,gxeEff[,i],targetVar)
      gxeEff[,i] = gxeEff[,i]*tmp$parameter
      gxeInt = tmp$output$intercept
      trait = new("TraitADG",
                  qtlLoci,
                  addEff=addEff[,i],
                  domEff=domEff[,i],
                  intercept=mean[i]-intercept,
                  gxeEff = gxeEff[,i],
                  gxeInt = 1-gxeInt,
                  envVar = varEnv[i])
      private$.addTrait(trait,tmp$output$varA,tmp$output$varG)
    }
    invisible(self)
  }
)

#' @title Remove trait
#' 
#' @description 
#' Removes designated trait(s).
#' 
#' @section Usage: SP$removeTrait(traits, force = FALSE)
#' 
#' @param traits a vector of traits to remove
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#'  
#' @name SimParam_removeTrait
NULL
# removeTrait ----
SimParam$set(
  "public",
  "removeTrait",
  function(traits,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    traits = as.integer(traits)
    stopifnot(max(traits)<=private$.nTraits)
    private$.traits = private$.traits[-traits]
    private$.varA = private$.varA[-traits]
    private$.varG = private$.varG[-traits]
    if(is.matrix(private$.varE)){
      private$.varE = private$.varE[-traits,-traits]
      if(length(private$.varE)==0){
        private$.varE = numeric()
      }
    }else{
      private$.varE = private$.varE[-traits]
    }
    private$.nTraits = length(private$.traits)
    invisible(self)
  }
)

#' @title Switch trait
#' 
#' @description 
#' Replaces an existing trait.
#' 
#' @section Usage: SP$switchTrait(lociMap, trait, varA = NULL, varG = NULL, 
#' force = FALSE)
#' 
#' @param lociMap a new object descended from 
#' \code{\link{LociMap-class}}
#' @param trait an integer indicating which trait to replace
#' @param varA a new value for varA in the base population. 
#' If NULL, the existing value is retained.
#' @param varG a new value for varG in the base population. 
#' If NULL, the existing value is retained.
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_switchTrait
NULL
# switchTrait ----
SimParam$set(
  "public",
  "switchTrait",
  function(lociMap,trait,varA=NULL,varG=NULL,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    stopifnot(length(trait)==1,trait<=private$.nTraits)
    if(!is.null(varA)){
      private$.varA[trait] = varA
    }
    if(!is.null(varG)){
      private$.varG[trait] = varG
    }
    private$.traits[[trait]] = lociMap
    invisible(self)
  }
)

#' @title Manually add trait
#' 
#' @description 
#' Add a new trait to the simulation.
#' 
#' @section Usage: SP$manAddTrait(lociMap, varA = NULL, varG = NULL, 
#' force = FALSE)
#' 
#' @param lociMap a new object descended from 
#' \code{\link{LociMap-class}}
#' @param varA a new value for varA in the base population. 
#' @param varG a new value for varG in the base population. 
#' @param force should the check for a running simulation be 
#' ignored. Only set to TRUE if you know what you are doing.
#' 
#' @name SimParam_manAddTrait
NULL
# manAddTrait ----
SimParam$set(
  "public",
  "manAddTrait",
  function(lociMap,varA=NULL,varG=NULL,force=FALSE){
    if(!force){
      private$.isRunning()
    }
    private$.nTraits = private$.nTraits+1L
    private$.varA[private$.nTraits] = varA
    private$.varG[private$.nTraits] = varG
    private$.traits[[private$.nTraits]] = lociMap
    invisible(self)
  }
)

#' @title Rescale traits
#' 
#' @description
#' Linearly scales all traits to achieve desired 
#' values of means and variances.
#' 
#' @section Usage: SP$rescaleTraits(pop, mean = 0, var = 1, varEnv = 1e-6, 
#' varGxE = 1e-6, useVarA = TRUE)
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param mean a vector of new trait means
#' @param var a vector of new trait variances
#' @param varEnv a vector of new environmental variances
#' @param varGxE a vector of new GxE variances
#' @param useVarA tune according to additive genetic variance if true
#'
#' @note
#' You must run \code{\link{resetPop}} on existing 
#' populations to obtain the new trait values.
#' 
#' @name SimParam_rescaleTraits
NULL
# rescaleTraits ----
SimParam$set(
  "public",
  "rescaleTraits",
  function(pop,mean=0,var=1,varEnv=1e-6,
           varGxE=1e-6,useVarA=TRUE){
    isGxe = sapply(private$.traits,function(x){
      class(x)%in%c("TraitAG","TraitADG")
    })
    if(any(isGxe)){
      stopifnot(length(mean)==private$.nTraits,
                length(var)==private$.nTraits,
                length(varEnv)==private$.nTraits,
                length(varGxE)==private$.nTraits)
    }else{
      stopifnot(length(mean)==private$.nTraits,
                length(var)==private$.nTraits)
    }
    
    for(i in 1:private$.nTraits){
      trait = private$.traits[[i]]
      geno = getGeno(pop@geno,
                     trait@lociPerChr,
                     trait@lociLoc)
      if(class(trait)%in%c("TraitAD","TraitADG")){
        tmp = tuneTraitAD(geno,trait@addEff,trait@domEff,var[i],useVarA)
        trait@domEff = trait@domEff*tmp$parameter
      }else{
        tmp = tuneTraitA(geno,trait@addEff,var[i])
      }
      trait@addEff = trait@addEff*tmp$parameter
      trait@intercept = mean[i]-tmp$output$intercept
      if(class(trait)%in%c("TraitAG","TraitADG")){
        targetVar = varGxE[i]/varEnv[i]
        tmp = tuneTraitA(geno,trait@gxeEff,targetVar)
        trait@gxeEff = trait@gxeEff*tmp$parameter
        trait@gxeInt = 1-tmp$output$intercept
        trait@envVar = varEnv[i]
      }
      private$.traits[[i]] = trait
    }
    invisible(self)
  }
)

#' @title Switch founder population
#' 
#' @description
#' Switches the founder population in the founderPop 
#' field. This may be desirable if traits are to be 
#' tuned to a population derived from the original 
#' founderPop. Note that no checking is performed to verify 
#' that the genetic map and/or number of segregating sites 
#' hasn't changed. The new founderPop can be 
#' \code{\link{MapPop-class}} or \code{\link{RawPop-class}}
#' 
#' @section Usage: SP$switchFounderPop(founderPop)
#' 
#' @name SimParam_removeFounderPop
NULL
# switchFounderPop ----
SimParam$set(
  "public",
  "switchFounderPop",
  function(founderPop){
    private$.founderPop = founderPop
    invisible(self)
  }
)

#' @title Remove founder population
#' 
#' @description
#' Removes the founder population from the founderPop 
#' field. This can be ran after all traits have been 
#' added to reduce the size of the SimParam object.
#' 
#' @section Usage: SP$removeFounderPop()
#' 
#' @name SimParam_removeFounderPop
NULL
# removeFounderPop ----
SimParam$set(
  "public",
  "removeFounderPop",
  function(){
    private$.founderPop = NULL
    invisible(self)
  }
)

# addToPed ----
SimParam$set(
  "public",
  "addToPed",
  function(lastId,mother,father,isDH){
    if(!private$.isTrackPed){
      stop("isTrackPed is FALSE")
    }
    nNewInd = lastId-private$.lastId
    stopifnot(nNewInd>0)
    if(length(isDH)==1) isDH = rep(isDH,nNewInd)
    mother = as.integer(mother)
    father = as.integer(father)
    isDH = as.integer(isDH)
    stopifnot(length(mother)==nNewInd,
              length(father)==nNewInd,
              length(isDH)==nNewInd)
    if(private$.isTrackRec){
      if(length(private$.recHist)==lastId){
        #Recombination history already added
      }else if(length(private$.recHist)==private$.lastId){
        #No recombination history, assume founder individuals
        private$.recHist = c(private$.recHist,vector("list",nNewInd))
      }else{
        stop("Unexpected outcome in recombination tracking")
      }
    }
    tmp = cbind(mother,father,isDH)
    private$.pedigree = rbind(private$.pedigree,tmp)
    private$.lastId = lastId
    invisible(self)
  }
)

# addToRec ----
SimParam$set(
  "public",
  "addToRec",
  function(hist){
    stopifnot(is.list(hist))
    if(!private$.isTrackRec){
      stop("isTrackRec is FALSE")
    }
    private$.recHist = c(private$.recHist,hist)
    invisible(self)
  }
)

#' @title Switch genetic map
#' 
#' @description 
#' Replaces existing genetic map.
#' 
#' @section Usage: SP$switchGenMap(genMap)
#' 
#' @param genMap a list of length nChr containing 
#' numeric vectors for the position of each segregating 
#' site on a chromosome.
#' 
#' @name SimParam_switchGenMap
NULL
# switchGenMap ----
SimParam$set(
  "public",
  "switchGenMap",
  function(genMap){
    stopifnot(length(genMap)==private$.nChr)
    tmp = do.call("c",lapply(genMap,length))
    stopifnot(all(tmp==private$.segSites))
    private$.sepMap = FALSE
    private$.femaleMap = genMap
    private$.maleMap = NULL
    invisible(self)
  }
)

#' @title Switch female genetic map
#' 
#' @description 
#' Replaces existing female genetic map.
#' 
#' @section Usage: SP$switchFemaleMap(genMap)
#' 
#' @param genMap a list of length nChr containing 
#' numeric vectors for the position of each segregating 
#' site on a chromosome.
#' 
#' @name SimParam_switchFemaleMap
NULL
# switchFemaleMap ----
SimParam$set(
  "public",
  "switchFemaleMap",
  function(genMap){
    stopifnot(length(genMap)==private$.nChr)
    tmp = do.call("c",lapply(genMap,length))
    stopifnot(all(tmp==private$.segSites))
    if(private$.sepMap){
      private$.femaleMap = genMap
    }else{
      private$.sepMap = TRUE
      private$.maleMap = private$.femaleMap
      private$.femaleMap = genMap
    }
    invisible(self)
  }
)

#' @title Switch male genetic map
#' 
#' @description 
#' Replaces existing male genetic map.
#' 
#' @section Usage: SP$switchMaleMap(genMap)
#' 
#' @param genMap a list of length nChr containing 
#' numeric vectors for the position of each segregating 
#' site on a chromosome.
#' 
#' @name SimParam_switchMaleMap
NULL
# switchMaleMap ----
SimParam$set(
  "public",
  "switchMaleMap",
  function(genMap){
    stopifnot(length(genMap)==private$.nChr)
    tmp = do.call("c",lapply(genMap,length))
    stopifnot(all(tmp==private$.segSites))
    private$.sepMap = TRUE
    private$.maleMap = genMap
    invisible(self)
  }
)
