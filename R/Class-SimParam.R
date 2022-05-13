#' @title Simulation parameters
#' 
#' @description 
#' Container for global simulation parameters. Saving this object 
#' as SP will allow it to be accessed by function defaults.
#'
#' @export
SimParam = R6Class(
  "SimParam",
  public = list(
    #### Public ----
    
    #' @field nThreads number of threads used on platforms with OpenMP support
    nThreads = "integer",
    
    #' @field snpChips list of SNP chips
    snpChips = "list",
    
    #' @field invalidQtl list of segregating sites that aren't valid QTL
    invalidQtl = "list",
    
    #' @field invalidSnp list of segregating sites that aren't valid SNP
    invalidSnp = "list",
    
    #' @field founderPop founder population used for variance scaling
    founderPop = "MapPop",
    
    #' @field finalizePop function applied to newly created populations.
    #' Currently does nothing and should only be changed by expert users.
    finalizePop = "function",
    
    #' @field allowEmptyPop if true, population arguments with nInd=0 will 
    #' return an empty population with a warning instead of an error.
    allowEmptyPop = "logical",
    
    #' @description Starts the process of building a new simulation 
    #' by creating a new SimParam object and assigning a founder 
    #' population to the class. It is recommended that you save the 
    #' object with the name "SP", because subsequent functions will 
    #' check your global environment for an object of this name if 
    #' their simParam arguments are NULL. This allows you to call 
    #' these functions without explicitly supplying a simParam 
    #' argument with every call.
    #' 
    #' @param founderPop an object of \code{\link{MapPop-class}}
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    initialize = function(founderPop){
      stopifnot(class(founderPop)=="MapPop")
      
      # Public items
      self$nThreads = getNumThreads()
      self$v = 2.6 # Kosambi
      self$p = 0 # Single pathway gamma model
      self$quadProb = 0 # No quadrivalent pairing
      self$snpChips = list()
      self$invalidQtl = vector("list",founderPop@nChr) # All eligible
      self$invalidSnp = vector("list",founderPop@nChr) # All eligible
      self$founderPop = founderPop
      self$finalizePop = function(pop, ...){return(pop)}
      self$allowEmptyPop = FALSE # Empty populations trigger an error
      
      # Private items
      private$.restrSites = TRUE
      private$.traits = list()
      private$.segSites = founderPop@nLoci
      private$.sexes = "no"
      private$.femaleMap = lapply(founderPop@genMap, function(x) x-x[1]) # Set position 1 to 0
      private$.maleMap = NULL
      private$.sepMap = FALSE
      private$.femaleCentromere = founderPop@centromere
      private$.maleCentromere = NULL
      private$.lastId = 0L
      private$.isTrackPed = FALSE
      private$.pedigree = matrix(NA_integer_,nrow=0,ncol=3)
      private$.isTrackRec = FALSE
      private$.recHist = list()
      private$.varA = numeric()
      private$.varG = numeric()
      private$.varE = numeric()
      private$.version = packageDescription("AlphaSimR")$Version 
      private$.lastHaplo = 0L
      private$.hasHap = logical()
      private$.hap = list()
      private$.isFounder = logical()
      
      invisible(self)
    },
    
    #' @description Sets pedigree tracking for the simulation. 
    #' By default pedigree tracking is turned off. When turned on, 
    #' the pedigree of all individuals created will be tracked, 
    #' except those created by \code{\link{hybridCross}}. Turning 
    #' off pedigree tracking will turn off recombination tracking 
    #' if it is turned on.
    #' 
    #' @param isTrackPed should pedigree tracking be on.
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$setTrackPed(TRUE)
    setTrackPed = function(isTrackPed, force=FALSE){
      stopifnot(is.logical(isTrackPed))
      if(!force){
        private$.isRunning()
      }
      private$.isTrackPed = isTrackPed
      if(!isTrackPed){
        private$.isTrackRec = FALSE
      }
      invisible(self)
    },
    
    #' @description Sets recombination tracking for the simulation.
    #' By default recombination tracking is turned off. When turned
    #' on recombination tracking will also turn on pedigree tracking.
    #' Recombination tracking keeps records of all individuals created,
    #' except those created by \code{\link{hybridCross}}, because their
    #' pedigree is not tracked.
    #'
    #' @param isTrackRec should recombination tracking be on.
    #' @param force should the check for a running simulation be
    #' ignored. Only set to TRUE if you know what you are doing.
    #'
    #' @examples
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #'
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$setTrackRec(TRUE)
    setTrackRec = function(isTrackRec, force=FALSE){
      stopifnot(is.logical(isTrackRec))
      if(!force){
        private$.isRunning()
      }
      private$.isTrackRec = isTrackRec
      if(isTrackRec){
        private$.isTrackPed = TRUE
      }
      invisible(self)
    },
    
    #' @description Resets the internal lastId, the pedigree 
    #' and recombination tracking (if in use) to the 
    #' supplied lastId. Be careful using this function because 
    #' it may introduce a bug if you use individuals from
    #' the deleted portion of the pedigree.
    #' 
    #' @param lastId last ID to include in pedigree
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' 
    #' #Create population
    #' pop = newPop(founderPop, simParam=SP)
    #' pop@id # 1:10
    #' 
    #' #Create another population after reseting pedigree
    #' SP$resetPed()
    #' pop2 = newPop(founderPop, simParam=SP)
    #' pop2@id # 1:10
    resetPed =function(lastId=0L){
      private$.lastId = lastId
      private$.pedigree = private$.pedigree[0:lastId,,drop=FALSE]
      if(private$.isTrackRec){
        private$.recHist = private$.recHist[0:lastId]
      }
      invisible(self)
    },
    
    #' @description Sets restrictions on which segregating sites 
    #' can serve as SNP and/or QTL.
    #' 
    #' @param minQtlPerChr the minimum number of segSites for QTLs. 
    #' Can be a single value or a vector values for each 
    #' chromosome.
    #' @param minSnpPerChr the minimum number of segSites for SNPs. 
    #' Can be a single value or a vector values for each 
    #' chromosome.
    #' @param overlap should SNP and QTL sites be allowed to overlap.
    #' @param minSnpFreq minimum allowable frequency for SNP loci. 
    #' No minimum SNP frequency is used if value is NULL.
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$restrSegSites(minQtlPerChr=5, minSnpPerChr=5)
    restrSegSites = function(minQtlPerChr=NULL, minSnpPerChr=NULL, overlap=FALSE,
                             minSnpFreq=NULL){
      if(overlap){
        private$.restrSites = FALSE
        invisible(self)
      }else{
        # Check inputs
        if(length(minSnpPerChr)==1){
          minSnpPerChr = rep(minSnpPerChr,self$nChr)
        }
        if(length(minQtlPerChr)==1){
          minQtlPerChr = rep(minQtlPerChr,self$nChr)
        }
        stopifnot(length(minSnpPerChr)==self$nChr,
                  length(minQtlPerChr)==self$nChr)
        
        # Restrict SNPs and then QTL
        private$.restrSites = TRUE
        invisible(private$.pickLoci(minSnpPerChr, FALSE, minSnpFreq))
        invisible(private$.pickLoci(minQtlPerChr))
        invisible(self)
      }
    },
    
    #' @description Changes how sexes are determined in the simulation. 
    #' The default sexes is "no", indicating all individuals are hermaphrodites. 
    #' To add sexes to the simulation, run this function with "yes_sys" or 
    #' "yes_rand". The value "yes_sys" will systematically assign 
    #' sexes to newly created individuals as first male and then female. 
    #' Populations with an odd number of individuals will have one more male than 
    #' female. The value "yes_rand" will randomly assign a sex to each
    #' individual.
    #' 
    #' @param sexes acceptable value are "no", "yes_sys", or 
    #' "yes_rand"
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$setSexes("yes_sys")
    setSexes = function(sexes, force=FALSE){
      if(!force){
        private$.isRunning()
      }
      sexes = tolower(sexes)
      if(sexes=="no"){
        private$.sexes="no"
      }else if(sexes=="yes_sys"){
        private$.sexes="yes_sys"
      }else if(sexes=="yes_rand"){
        private$.sexes="yes_rand"
      }else{
        stop(paste0("sexes=",sexes," is not a valid option"))
      }
      invisible(self)
    },
    
    #' @description 
    #' Randomly assigns eligible SNPs to a SNP chip
    #' 
    #' @param nSnpPerChr number of SNPs per chromosome. 
    #' Can be a single value or nChr values.
    #' @param minSnpFreq minimum allowable frequency for SNP loci.
    #' If NULL, no minimum frequency is used. 
    #' @param refPop reference population for calculating SNP 
    #' frequency. If NULL, the founder population is used.
    #' @param name optional name for chip
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addSnpChip(10)
    addSnpChip = function(nSnpPerChr, minSnpFreq=NULL, refPop=NULL, name=NULL){
      if(length(nSnpPerChr)==1){
        nSnpPerChr = rep(nSnpPerChr,self$nChr)
      }
      snpChip = private$.pickLoci(nSnpPerChr, FALSE, minSnpFreq, refPop)
      if(is.null(name)){
        snpChip@name = paste0("Chip",self$nSnpChips + 1L)
      }else{
        snpChip@name = name
      }
      self$snpChips[[self$nSnpChips + 1L]] = snpChip
      invisible(self)
    },
    
    #' @description 
    #' Randomly selects the number of snps in structure and then
    #' assigns them to chips based on structure
    #' 
    #' @param nSnpPerChr number of SNPs per chromosome. 
    #' Can be a single value or nChr values.
    #' @param structure a matrix.  Rows are snp chips, columns are chips.
    #' If value is true then that snp is on that chip.
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    addStructuredSnpChip = function(nSnpPerChr,structure,force=FALSE){
      if(!force){
        private$.isRunning()
      }
      if(length(nSnpPerChr)==1){
        nSnpPerChr = rep(nSnpPerChr,self$nChr)
      }
      stopifnot(length(nSnpPerChr)==self$nChr)
      stopifnot(sapply(self$potSnp,length)>=nSnpPerChr)
      stopifnot(dim(structure)[2]==sum(nSnpPerChr))
      lociLoc = lapply(1:self$nChr,function(x){
        sort(sample(self$potSnp[[x]],nSnpPerChr[x]))
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
        self$snpChips[[self$nSnpChips+1L]] = snpChip
      }
      invisible(self)
    },
    
    ### Traits (public) ----
    
    #' @description 
    #' Randomly assigns eligible QTLs for one or more additive traits. 
    #' If simulating more than one trait, all traits will be pleiotrophic 
    #' with correlated additive effects.
    #' 
    #' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
    #' @param mean a vector of desired mean genetic values for one or more traits
    #' @param var a vector of desired genetic variances for one or more traits
    #' @param corA a matrix of correlations between additive effects
    #' @param gamma should a gamma distribution be used instead of normal
    #' @param shape the shape parameter for the gamma distribution
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    #' @param name optional name for trait(s)
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitA(10)
    addTraitA = function(nQtlPerChr,mean=0,var=1,corA=NULL,
                         gamma=FALSE,shape=1,force=FALSE,name=NULL){
      if(!force){
        private$.isRunning()
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      nTraits = length(mean)
      if(length(gamma)==1) gamma = rep(gamma,nTraits)
      if(length(shape)==1) shape = rep(shape,nTraits)
      if(is.null(corA)) corA=diag(nTraits)
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }
      stopifnot(length(mean)==length(var),
                isSymmetric(corA),
                length(mean)==nrow(corA),
                length(mean)==length(name))
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      for(i in 1:nTraits){
        trait = new("TraitA",
                    qtlLoci,
                    addEff=addEff[,i],
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, 
                           self$nThreads)
        scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        trait@addEff = trait@addEff*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        private$.addTrait(trait,var[i],var[i])
      }
      invisible(self)
    },
    
    #' @description 
    #' Randomly assigns eligible QTLs for one or more traits with dominance. 
    #' If simulating more than one trait, all traits will be pleiotrophic 
    #' with correlated effects.
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
    #' @param name optional name for trait(s)
    #'  
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitAD(10, meanDD=0.5)
    addTraitAD = function(nQtlPerChr,mean=0,var=1,meanDD=0,
                          varDD=0,corA=NULL,corDD=NULL,useVarA=TRUE,
                          gamma=FALSE,shape=1,force=FALSE,name=NULL){
      if(!force){
        private$.isRunning()
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      nTraits = length(mean)
      if(length(meanDD)==1) meanDD = rep(meanDD,nTraits)
      if(length(varDD)==1) varDD = rep(varDD,nTraits)
      if(length(gamma)==1) gamma = rep(gamma,nTraits)
      if(length(shape)==1) shape = rep(shape,nTraits)
      if(is.null(corA)) corA=diag(nTraits)
      if(is.null(corDD)) corDD=diag(nTraits)
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }
      stopifnot(length(mean)==length(var),
                isSymmetric(corA),
                isSymmetric(corDD),
                length(mean)==nrow(corA),
                length(mean)==length(name))
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                          corDD=corDD,meanDD=meanDD,varDD=varDD)
      for(i in 1:nTraits){
        trait = new("TraitAD",
                    qtlLoci,
                    addEff=addEff[,i],
                    domEff=domEff[,i],
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, 
                           self$nThreads)
        if(useVarA){
          scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        }else{
          scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
        }
        trait@addEff = trait@addEff*scale
        trait@domEff = trait@domEff*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        if(useVarA){
          private$.addTrait(trait,var[i],popVar(tmp$gv*scale)[1])
        }else{
          private$.addTrait(trait,popVar(tmp$bv*scale)[1],var[i])
        }
      }
      invisible(self)
    },
    
    #' @description 
    #' Randomly assigns eligible QTLs for one ore more additive GxE traits. 
    #' If simulating more than one trait, all traits will be pleiotrophic 
    #' with correlated effects.
    #' 
    #' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
    #' @param mean a vector of desired mean genetic values for one or more traits
    #' @param var a vector of desired genetic variances for one or more traits
    #' @param varGxE a vector of total genotype-by-environment variances for the traits
    #' @param varEnv a vector of environmental variances for one or more traits
    #' @param corA a matrix of correlations between additive effects
    #' @param corGxE a matrix of correlations between GxE effects
    #' @param gamma should a gamma distribution be used instead of normal
    #' @param shape the shape parameter for the gamma distribution
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    #' @param name optional name for trait(s)
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitAG(10, varGxE=2)
    addTraitAG = function(nQtlPerChr,mean=0,var=1,varGxE=1e-6,varEnv=0,
                          corA=NULL,corGxE=NULL,gamma=FALSE,shape=1,
                          force=FALSE,name=NULL){
      if(!force){
        private$.isRunning()
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      nTraits = length(mean)
      if(length(gamma)==1) gamma = rep(gamma,nTraits)
      if(length(shape)==1) shape = rep(shape,nTraits)
      if(length(varEnv)==1) varEnv = rep(varEnv,nTraits)
      if(is.null(corA)) corA=diag(nTraits)
      if(is.null(corGxE)) corGxE=diag(nTraits)
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }
      stopifnot(length(mean)==length(var),
                isSymmetric(corA),
                isSymmetric(corGxE),
                length(mean)==nrow(corA),
                length(mean)==nrow(corGxE),
                length(mean)==length(varGxE),
                length(mean)==length(varEnv),
                length(mean)==length(name))
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corGxE,gamma=FALSE,shape=NULL)
      for(i in 1:nTraits){
        trait = new("TraitA",
                    qtlLoci,
                    addEff=addEff[,i],
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, 
                           self$nThreads)
        scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        trait@addEff = trait@addEff*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        
        # GxE component
        traitG = new("TraitA",
                     qtlLoci,
                     addEff=gxeEff[,i],
                     intercept=0)
        tmpG = calcGenParam(traitG, self$founderPop, 
                            self$nThreads)
        if(varEnv[i]==0){
          scaleG = sqrt(varGxE[i])/sqrt(popVar(tmpG$gv)[1])
          trait = new("TraitAG",
                      trait,
                      gxeEff = gxeEff[,i]*scaleG,
                      gxeInt = 0-mean(tmpG$gv*scaleG),
                      envVar = 1)
        }else{
          scaleG = sqrt(varGxE[i]/varEnv[i])/sqrt(popVar(tmpG$gv)[1])
          trait = new("TraitAG",
                      trait,
                      gxeEff = gxeEff[,i]*scaleG,
                      gxeInt = 1-mean(tmpG$gv*scaleG),
                      envVar = varEnv[i])
        }
        
        private$.addTrait(trait,var[i],var[i])
      }
      invisible(self)
    },
    
    #' @description 
    #' Randomly assigns eligible QTLs for a trait with dominance and GxE. 
    #' 
    #' @param nQtlPerChr number of QTLs per chromosome. Can be a single 
    #' value or nChr values.
    #' @param mean a vector of desired mean genetic values for one or more traits
    #' @param var a vector of desired genetic variances for one or more traits
    #' @param varGxE a vector of total genotype-by-environment variances for the traits
    #' @param varEnv a vector of environmental variances for one or more traits
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
    #' @param name optional name for trait(s)
    #'  
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitADG(10, meanDD=0.5, varGxE=2)
    addTraitADG = function(nQtlPerChr,mean=0,var=1,varEnv=0,
                           varGxE=1e-6,meanDD=0,varDD=0,corA=NULL,
                           corDD=NULL,corGxE=NULL,useVarA=TRUE,gamma=FALSE,
                           shape=1,force=FALSE,name=NULL){
      if(!force){
        private$.isRunning()
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      nTraits = length(mean)
      if(length(meanDD)==1) meanDD = rep(meanDD,nTraits)
      if(length(varDD)==1) varDD = rep(varDD,nTraits)
      if(length(varEnv)==1) varEnv = rep(varEnv,nTraits)
      if(length(gamma)==1) gamma = rep(gamma,nTraits)
      if(length(shape)==1) shape = rep(shape,nTraits)
      if(is.null(corA)) corA=diag(nTraits)
      if(is.null(corDD)) corDD=diag(nTraits)
      if(is.null(corGxE)) corGxE=diag(nTraits)
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }
      stopifnot(length(mean)==length(var),
                isSymmetric(corA),
                isSymmetric(corDD),
                isSymmetric(corGxE),
                nrow(corA)==nTraits,
                nrow(corGxE)==nTraits,
                nrow(corDD)==nTraits,
                length(varGxE)==nTraits,
                length(varEnv)==nTraits,
                length(mean)==length(name))
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                          corDD=corDD,meanDD=meanDD,varDD=varDD)
      gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corGxE,gamma=FALSE,shape=NULL)
      for(i in 1:nTraits){
        trait = new("TraitAD",
                    qtlLoci,
                    addEff=addEff[,i],
                    domEff=domEff[,i],
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, 
                           self$nThreads)
        if(useVarA){
          scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        }else{
          scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
        }
        trait@addEff = trait@addEff*scale
        trait@domEff = trait@domEff*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        
        # GxE component
        traitG = new("TraitA",
                     qtlLoci,
                     addEff=gxeEff[,i],
                     intercept=0)
        tmpG = calcGenParam(traitG, self$founderPop, 
                            self$nThreads)
        if(varEnv[i]==0){
          scaleG = sqrt(varGxE[i])/sqrt(popVar(tmpG$gv)[1])
          trait = new("TraitADG",
                      trait,
                      gxeEff = gxeEff[,i]*scaleG,
                      gxeInt = 0-mean(tmpG$gv*scaleG),
                      envVar = 1)
        }else{
          scaleG = sqrt(varGxE[i]/varEnv[i])/sqrt(popVar(tmpG$gv)[1])
          trait = new("TraitADG",
                      trait,
                      gxeEff = gxeEff[,i]*scaleG,
                      gxeInt = 1-mean(tmpG$gv*scaleG),
                      envVar = varEnv[i])
        }
        
        if(useVarA){
          private$.addTrait(trait,var[i],popVar(tmp$gv*scale)[1])
        }else{
          private$.addTrait(trait,popVar(tmp$bv*scale)[1],var[i])
        }
      }
      invisible(self)
    },
    
    #' @description 
    #' Randomly assigns eligible QTLs for one or more additive and epistasis 
    #' traits. If simulating more than one trait, all traits will be pleiotrophic 
    #' with correlated additive effects.
    #' 
    #' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
    #' @param mean a vector of desired mean genetic values for one or more traits
    #' @param var a vector of desired genetic variances for one or more traits
    #' @param relAA the relative value of additive-by-additive variance compared 
    #' to additive variance in a diploid organism with allele frequency 0.5
    #' @param corA a matrix of correlations between additive effects
    #' @param corAA a matrix of correlations between additive-by-additive effects
    #' @param useVarA tune according to additive genetic variance if true. If 
    #' FALSE, tuning is performed according to total genetic variance.
    #' @param gamma should a gamma distribution be used instead of normal
    #' @param shape the shape parameter for the gamma distribution
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    #' @param name optional name for trait(s)
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    addTraitAE = function(nQtlPerChr,mean=0,var=1,relAA=0,corA=NULL,
                          corAA=NULL,useVarA=TRUE,gamma=FALSE,shape=1,force=FALSE,
                          name=NULL){
      if(!force){
        private$.isRunning()
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      nTraits = length(mean)
      relAA = relAA*4
      if(length(gamma)==1) gamma = rep(gamma,nTraits)
      if(length(shape)==1) shape = rep(shape,nTraits)
      if(length(relAA)==1) relAA = rep(relAA,nTraits)
      if(is.null(corA)) corA=diag(nTraits)
      if(is.null(corAA)) corAA=diag(nTraits)
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }
      stopifnot(length(mean)==length(var),
                isSymmetric(corA),
                isSymmetric(corAA),
                length(relAA)==length(mean),
                length(mean)==nrow(corA),
                (sum(nQtlPerChr)%%2L)==0L,
                length(mean)==length(name))
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape,
                          relVar=relAA)
      E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
      for(i in 1:nTraits){
        trait = new("TraitAE",
                    qtlLoci,
                    addEff=addEff[,i],
                    epiEff=cbind(E,epiEff[,i]),
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, 
                           self$nThreads)
        if(useVarA){
          scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        }else{
          scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
        }
        trait@addEff = trait@addEff*scale
        trait@epiEff[,3] = trait@epiEff[,3]*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        if(useVarA){
          private$.addTrait(trait,var[i],popVar(tmp$gv*scale)[1])
        }else{
          private$.addTrait(trait,popVar(tmp$bv*scale)[1],var[i])
        }
      }
      invisible(self)
    },
    
    #' @description 
    #' Randomly assigns eligible QTLs for one or more traits with dominance and 
    #' epistasis. If simulating more than one trait, all traits will be pleiotrophic 
    #' with correlated effects.
    #' 
    #' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
    #' @param mean a vector of desired mean genetic values for one or more traits
    #' @param var a vector of desired genetic variances for one or more traits
    #' @param meanDD mean dominance degree
    #' @param varDD variance of dominance degree
    #' @param relAA the relative value of additive-by-additive variance compared 
    #' to additive variance in a diploid organism with allele frequency 0.5
    #' @param corA a matrix of correlations between additive effects
    #' @param corDD a matrix of correlations between dominance degrees
    #' @param corAA a matrix of correlations between additive-by-additive effects
    #' @param useVarA tune according to additive genetic variance if true. If 
    #' FALSE, tuning is performed according to total genetic variance.
    #' @param gamma should a gamma distribution be used instead of normal
    #' @param shape the shape parameter for the gamma distribution
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    #' @param name optional name for trait(s)
    #'  
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitADE(10)
    addTraitADE = function(nQtlPerChr,mean=0,var=1,meanDD=0,
                           varDD=0,relAA=0,corA=NULL,corDD=NULL,corAA=NULL,
                           useVarA=TRUE,gamma=FALSE,shape=1,force=FALSE,
                           name=NULL){
      if(!force){
        private$.isRunning()
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      nTraits = length(mean)
      relAA = relAA*4
      if(length(meanDD)==1) meanDD = rep(meanDD,nTraits)
      if(length(varDD)==1) varDD = rep(varDD,nTraits)
      if(length(gamma)==1) gamma = rep(gamma,nTraits)
      if(length(shape)==1) shape = rep(shape,nTraits)
      if(length(relAA)==1) relAA = rep(relAA,nTraits)
      if(is.null(corA)) corA=diag(nTraits)
      if(is.null(corDD)) corDD=diag(nTraits)
      if(is.null(corAA)) corAA=diag(nTraits)
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }
      stopifnot(length(mean)==length(var),
                isSymmetric(corA),
                isSymmetric(corDD),
                length(mean)==nrow(corA),
                length(mean)==nrow(corAA),
                length(mean)==nrow(corDD),
                length(relAA)==length(mean),
                (sum(nQtlPerChr)%%2L)==0L,
                length(mean)==length(name))
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                          corDD=corDD,meanDD=meanDD,varDD=varDD)
      epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape,
                          relVar=relAA)
      E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
      for(i in 1:nTraits){
        trait = new("TraitADE",
                    qtlLoci,
                    addEff=addEff[,i],
                    domEff=domEff[,i],
                    epiEff=cbind(E,epiEff[,i]),
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, 
                           self$nThreads)
        if(useVarA){
          scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        }else{
          scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
        }
        trait@addEff = trait@addEff*scale
        trait@domEff = trait@domEff*scale
        trait@epiEff[,3] = trait@epiEff[,3]*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        if(useVarA){
          private$.addTrait(trait,var[i],popVar(tmp$gv*scale)[1])
        }else{
          private$.addTrait(trait,popVar(tmp$bv*scale)[1],var[i])
        }
      }
      invisible(self)
    },
    
    #' @description 
    #' Randomly assigns eligible QTLs for one or more additive and epistasis 
    #' GxE traits. If simulating more than one trait, all traits will be pleiotrophic 
    #' with correlated effects.
    #' 
    #' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
    #' @param mean a vector of desired mean genetic values for one or more traits
    #' @param var a vector of desired genetic variances for one or more traits
    #' @param relAA the relative value of additive-by-additive variance compared 
    #' to additive variance in a diploid organism with allele frequency 0.5
    #' @param varGxE a vector of total genotype-by-environment variances for the traits
    #' @param varEnv a vector of environmental variances for one or more traits
    #' @param corA a matrix of correlations between additive effects
    #' @param corAA a matrix of correlations between additive-by-additive effects
    #' @param corGxE a matrix of correlations between GxE effects
    #' @param useVarA tune according to additive genetic variance if true. If 
    #' FALSE, tuning is performed according to total genetic variance.
    #' @param gamma should a gamma distribution be used instead of normal
    #' @param shape the shape parameter for the gamma distribution
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    #' @param name optional name for trait(s)
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitAEG(10, varGxE=2)
    addTraitAEG = function(nQtlPerChr,mean=0,var=1,relAA=0,varGxE=1e-6,varEnv=0,
                           corA=NULL,corAA=NULL,corGxE=NULL,useVarA=TRUE,gamma=FALSE,
                           shape=1,force=FALSE,name=NULL){
      if(!force){
        private$.isRunning()
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      nTraits = length(mean)
      relAA = relAA*4
      if(length(gamma)==1) gamma = rep(gamma,nTraits)
      if(length(shape)==1) shape = rep(shape,nTraits)
      if(length(relAA)==1) relAA = rep(relAA,nTraits)
      if(length(varEnv)==1) varEnv = rep(varEnv,nTraits)
      if(is.null(corA)) corA=diag(nTraits)
      if(is.null(corAA)) corAA=diag(nTraits)
      if(is.null(corGxE)) corGxE=diag(nTraits)
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }
      stopifnot(length(mean)==length(var),
                length(relAA)==length(mean),
                isSymmetric(corA),
                isSymmetric(corGxE),
                isSymmetric(corAA),
                length(mean)==nrow(corA),
                length(mean)==nrow(corAA),
                length(mean)==nrow(corGxE),
                length(mean)==length(varGxE),
                length(mean)==length(varEnv),
                (sum(nQtlPerChr)%%2L)==0L,
                length(mean)==length(name))
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape,
                          relVar=relAA)
      E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
      gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corGxE,gamma=FALSE,shape=NULL)
      for(i in 1:nTraits){
        trait = new("TraitAE",
                    qtlLoci,
                    addEff=addEff[,i],
                    epiEff=cbind(E,epiEff[,i]),
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, 
                           self$nThreads)
        if(useVarA){
          scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        }else{
          scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
        }
        trait@addEff = trait@addEff*scale
        trait@epiEff[,3] = trait@epiEff[,3]*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        
        # GxE component
        traitG = new("TraitA",
                     qtlLoci,
                     addEff=gxeEff[,i],
                     intercept=0)
        tmpG = calcGenParam(traitG, self$founderPop, 
                            self$nThreads)
        if(varEnv[i]==0){
          scaleG = sqrt(varGxE[i])/sqrt(popVar(tmpG$gv)[1])
          trait = new("TraitAEG",
                      trait,
                      gxeEff = gxeEff[,i]*scaleG,
                      gxeInt = 0-mean(tmpG$gv*scaleG),
                      envVar = 1)
        }else{
          scaleG = sqrt(varGxE[i]/varEnv[i])/sqrt(popVar(tmpG$gv)[1])
          trait = new("TraitAEG",
                      trait,
                      gxeEff = gxeEff[,i]*scaleG,
                      gxeInt = 1-mean(tmpG$gv*scaleG),
                      envVar = varEnv[i])
        }
        
        if(useVarA){
          private$.addTrait(trait,var[i],popVar(tmp$gv*scale)[1])
        }else{
          private$.addTrait(trait,popVar(tmp$bv*scale)[1],var[i])
        }
      }
      invisible(self)
    },
    
    #' @description 
    #' Randomly assigns eligible QTLs for a trait with dominance, 
    #' epistasis and GxE. 
    #' 
    #' @param nQtlPerChr number of QTLs per chromosome. Can be a single 
    #' value or nChr values.
    #' @param mean a vector of desired mean genetic values for one or more traits
    #' @param var a vector of desired genetic variances for one or more traits
    #' @param varGxE a vector of total genotype-by-environment variances for the traits
    #' @param varEnv a vector of environmental variances for one or more traits
    #' @param meanDD mean dominance degree
    #' @param varDD variance of dominance degree
    #' @param relAA the relative value of additive-by-additive variance compared 
    #' to additive variance in a diploid organism with allele frequency 0.5
    #' @param corA a matrix of correlations between additive effects
    #' @param corDD a matrix of correlations between dominance degrees
    #' @param corAA a matrix of correlations between additive-by-additive effects
    #' @param corGxE a matrix of correlations between GxE effects
    #' @param useVarA tune according to additive genetic variance if true
    #' @param gamma should a gamma distribution be used instead of normal
    #' @param shape the shape parameter for the gamma distribution
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing.
    #' @param name optional name for trait(s)
    #'  
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitADEG(10, meanDD=0.5, varGxE=2)
    addTraitADEG = function(nQtlPerChr,mean=0,var=1,varEnv=0,
                            varGxE=1e-6,meanDD=0,varDD=0,relAA=0,corA=NULL,
                            corDD=NULL,corAA=NULL,corGxE=NULL,useVarA=TRUE,
                            gamma=FALSE,shape=1,force=FALSE,name=NULL){
      if(!force){
        private$.isRunning()
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      nTraits = length(mean)
      relAA = relAA*4
      if(length(meanDD)==1) meanDD = rep(meanDD,nTraits)
      if(length(varDD)==1) varDD = rep(varDD,nTraits)
      if(length(relAA)==1) relAA = rep(relAA,nTraits)
      if(length(varEnv)==1) varEnv = rep(varEnv,nTraits)
      if(length(gamma)==1) gamma = rep(gamma,nTraits)
      if(length(shape)==1) shape = rep(shape,nTraits)
      if(is.null(corA)) corA=diag(nTraits)
      if(is.null(corDD)) corDD=diag(nTraits)
      if(is.null(corGxE)) corGxE=diag(nTraits)
      if(is.null(corAA)) corAA=diag(nTraits)
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }
      stopifnot(length(mean)==length(var),
                isSymmetric(corA),
                isSymmetric(corDD),
                isSymmetric(corGxE),
                isSymmetric(corAA),
                nrow(corA)==nTraits,
                nrow(corGxE)==nTraits,
                nrow(corDD)==nTraits,
                nrow(corAA)==nTraits,
                length(varGxE)==nTraits,
                length(varEnv)==nTraits,
                length(relAA)==length(mean),
                (sum(nQtlPerChr)%%2L)==0L,
                length(mean)==length(name))
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                          corDD=corDD,meanDD=meanDD,varDD=varDD)
      epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape,
                          relVar=relAA)
      E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
      gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corGxE,gamma=FALSE,shape=NULL)
      for(i in 1:nTraits){
        trait = new("TraitADE",
                    qtlLoci,
                    addEff=addEff[,i],
                    domEff=domEff[,i],
                    epiEff=cbind(E,epiEff[,i]),
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, 
                           self$nThreads)
        if(useVarA){
          scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        }else{
          scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
        }
        trait@addEff = trait@addEff*scale
        trait@domEff = trait@domEff*scale
        trait@epiEff[,3] = trait@epiEff[,3]*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        
        # GxE component
        traitG = new("TraitA",
                     qtlLoci,
                     addEff=gxeEff[,i],
                     intercept=0)
        tmpG = calcGenParam(traitG, self$founderPop, 
                            self$nThreads)
        if(varEnv[i]==0){
          scaleG = sqrt(varGxE[i])/sqrt(popVar(tmpG$gv)[1])
          trait = new("TraitADEG",
                      trait,
                      gxeEff = gxeEff[,i]*scaleG,
                      gxeInt = 0-mean(tmpG$gv*scaleG),
                      envVar = 1)
        }else{
          scaleG = sqrt(varGxE[i]/varEnv[i])/sqrt(popVar(tmpG$gv)[1])
          trait = new("TraitADEG",
                      trait,
                      gxeEff = gxeEff[,i]*scaleG,
                      gxeInt = 1-mean(tmpG$gv*scaleG),
                      envVar = varEnv[i])
        }
        
        if(useVarA){
          private$.addTrait(trait,var[i],popVar(tmp$gv*scale)[1])
        }else{
          private$.addTrait(trait,popVar(tmp$bv*scale)[1],var[i])
        }
      }
      invisible(self)
    },
    
    #' @description 
    #' Manually add a new trait to the simulation. Trait must 
    #' be formatted as a \code{\link{LociMap-class}}. If the 
    #' trait is not already formatted, consider using importTrait.
    #' 
    #' @param lociMap a new object descended from 
    #' \code{\link{LociMap-class}}
    #' @param varE default error variance for phenotype, optional
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing
    manAddTrait = function(lociMap,varE=NA_real_,force=FALSE){
      if(!force){
        private$.isRunning()
      }
      stopifnot(is(lociMap,"LociMap"))
      tmp = calcGenParam(lociMap, self$founderPop, 
                         self$nThreads)
      varA = popVar(tmp$bv)[1]
      varG = popVar(tmp$gv)[1]
      private$.addTrait(lociMap,varA,varG,varE)
      invisible(self)
    },
    
    #' @description 
    #' Manually add a new trait(s) to the simulation. Unlike the 
    #' manAddTrait function, this function does not require 
    #' formatting the trait as a \code{\link{LociMap-class}}. 
    #' The formatting is performed automatically for the user, 
    #' with more user friendly data.frames or matrices taken as 
    #' inputs. This function only works for A and AD trait types.
    #' 
    #' @param markerNames a vector of names for the QTL
    #' @param addEff a matrix of additive effects (nLoci x nTraits). 
    #' Alternatively, a vector of length nLoci can be supplied for 
    #' a single trait.
    #' @param domEff optional dominance effects for each locus
    #' @param intercept optional intercepts for each trait
    #' @param name optional name(s) for the trait(s)
    #' @param varE default error variance for phenotype, optional
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing
    importTrait = function(markerNames, 
                           addEff, 
                           domEff=NULL,
                           intercept=NULL, 
                           name=NULL, 
                           varE=NULL,
                           force=FALSE){
      if(!force){
        private$.isRunning()
      }
      
      # Check addEff and domEff inputs
      addEff = as.matrix(addEff)
      stopifnot(length(markerNames)==nrow(addEff))
      nTraits = ncol(addEff)
      if(is.null(domEff)){
        useDom = FALSE
      }else{
        useDom = TRUE
        domEff = as.matrix(domEff)
        stopifnot(nrow(addEff)==nrow(domEff),
                  ncol(addEff)==nrow(domEff))
      }
      
      # Prepare the intercept
      if(is.null(intercept)){
        intercept = rep(0, nTraits)
      }else{
        intercept = as.numeric(intercept)
        stopifnot(length(intercept)==nTraits)
      }
      
      # Prepare varE
      if(!is.null(varE)){
        varE = as.numeric(varE)
        stopifnot(length(varE)==nTraits)
      }
      
      # Prepare trait names
      if(is.null(names)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }else{
        stopifnot(length(names)==nTraits)
      }
      
      
      # Extract genetic map and check if marker names are on the map
      genMapMarkerNames = unlist(lapply(private$.femaleMap, names))
      stopifnot(all(markerNames%in%genMapMarkerNames))
      
      # Create trait variables
      lociPerChr = integer(self$nChr)
      lociLoc = vector("list", self$nChr)
      addEffList = domEffList = vector("list", nTraits)
      for(i in 1:nTraits){
        addEffList[[i]] = domEffList[[i]] = vector("list", self$nChr)
      }
      
      # Loop through chromosomes
      for(i in 1:self$nChr){
        # Working on trait 1
        # Initialize variables
        addEffList[[1]][[i]] = domEffList[[1]][[i]] = numeric()
        lociLoc[[i]] = integer()
        
        # Find matches if they exist
        take = match(names(genMap[[i]]), markerNames)
        lociPerChr[i] = length(na.omit(take))
        
        if(lociPerChr[i]>0L){
          lociLoc[[i]] = which(!is.na(take))
          addEffList[[1]][[i]] = addEff[na.omit(take),1]
          if(useDom){
            domEffList[[1]][[i]] = domEff[na.omit(take),1]
          }
        }
        
        # Work on additional traits?
        if(nTraits>1){
          for(j in 2:nTraits){
            addEffList[[j]][[i]] = domEffList[[j]][[i]] = numeric()
            if(lociPerChr[i]>0L){
              addEffList[[j]][[i]] = addEff[na.omit(take),j]
              if(useDom){
                domEffList[[j]][[i]] = domEff[na.omit(take),j]
              }
            }
          }
        }
      }
      
      lociLoc = unlist(lociLoc)
      nLoci = sum(lociPerChr)
      
      # Create Trait(s)
      for(i in 1:nTraits){
        addEff = unlist(addEffList[[i]])
        if(useDom){
          domEff = unlist(domEffList[[i]])
          trait = new("TraitAD", 
                      addEff=addEff,
                      domEff=domEff,
                      intercept=intercept[i],
                      nLoci=nLoci,
                      lociPerChr=lociPerChr,
                      lociLoc=lociLoc,
                      name=name[i])
        }else{
          trait = new("TraitA", 
                      addEff=addEff,
                      intercept=intercept[i],
                      nLoci=nLoci,
                      lociPerChr=lociPerChr,
                      lociLoc=lociLoc,
                      name=name[i])
        }
        
        # Add trait to simParam
        self$manAddTrait(lociMap=trait, varE=varE[i], force=force)
      }
      
      invisible(self)
    },
    
    #' @description 
    #' Switch a trait in the simulation.
    #' 
    #' @param traitPos an integer indicate which trait to switch
    #' @param lociMap a new object descended from 
    #' \code{\link{LociMap-class}}
    #' @param varE default error variance for phenotype, optional
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing
    switchTrait = function(traitPos,lociMap,varE=NA_real_,force=FALSE){
      if(!force){
        private$.isRunning()
      }
      stopifnot(is(lociMap,"LociMap"),
                traitPos<=self$nTraits,
                traitPos>0)
      tmp = calcGenParam(lociMap, self$founderPop, 
                         self$nThreads)
      private$.traits[[traitPos]] = lociMap
      private$.varA[traitPos] = popVar(tmp$bv)[1]
      private$.varG[traitPos] = popVar(tmp$gv)[1]
      if(is.matrix(private$.varE)){
        private$.varE[traitPos,] = 0
        private$.varE[,traitPos] = 0
        private$.varE[traitPos,traitPos] = varE
      }else{
        private$.varE[traitPos] = varE
      }
      invisible(self)
    },
    
    #' @description
    #' Remove a trait from the simulation
    #' 
    #' @param traits an integer vector indicating which traits to remove
    #' @param force should the check for a running simulation be 
    #' ignored. Only set to TRUE if you know what you are doing
    removeTrait = function(traits,force=FALSE){
      if(!force){
        private$.isRunning()
      }
      stopifnot(max(traits)<=self$nTraits, min(traits)>0)
      private$.traits[-traits]
      private$.varA[-traits]
      private$.varG[-traits]
      if(is.matrix(private$.varE)){
        private$.varE[-traits,-traits]
      }else{
        private$.varE[-traits]
      }
      invisible(self)
    },
    
    #' @description Defines a default value for error 
    #' variances in the simulation.
    #' 
    #' @param h2 a vector of desired narrow-sense heritabilities
    #' @param H2 a vector of desired broad-sense heritabilities
    #' @param varE a vector or matrix of error variances
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitA(10)
    #' SP$setVarE(h2=0.5)
    setVarE = function(h2=NULL,H2=NULL,varE=NULL){
      if(!is.null(h2)){
        stopifnot(length(h2)==self$nTraits,
                  all(private$.varG>0),
                  all(private$.varA>0))
        varE = numeric(self$nTraits)
        for(i in 1:length(h2)){
          tmp = private$.varA[i]/h2[i]-private$.varG[i]
          if(tmp<0){
            stop(paste0("h2=",h2[i]," is not possible for trait ",i))
          }
          varE[i] = tmp
        }
        private$.varE = varE
      }else if(!is.null(H2)){
        stopifnot(length(H2)==self$nTraits)
        varE = numeric(self$nTraits)
        for(i in 1:length(H2)){
          tmp = private$.varG[i]/H2[i]-private$.varG[i]
          varE[i] = tmp
        }
        private$.varE = varE
      }else if(!is.null(varE)){
        if(is.matrix(varE)){
          stopifnot(nrow(varE)==self$nTraits,
                    ncol(varE)==self$nTraits)
        }else{
          stopifnot(length(varE)==self$nTraits)
        }
        private$.varE = varE
      }else{
        private$.varE = rep(NA_real_,self$nTraits)
      }
      invisible(self)
    },
    
    #' @description Defines a correlation structure for default 
    #' error variances. You must call \code{setVarE} first to define 
    #' the default error variances.
    #' 
    #' @param corE a correlation matrix for the error variances
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitA(10, mean=c(0,0), var=c(1,1), corA=diag(2))
    #' SP$setVarE(varE=c(1,1))
    #' E = 0.5*diag(2)+0.5 #Positively correlated error
    #' SP$setCorE(E)
    setCorE = function(corE){
      stopifnot(isSymmetric(corE),
                nrow(corE)==self$nTraits,
                length(private$.varE)==self$nTraits)
      if(is.matrix(private$.varE)){
        varE = diag(private$.varE)
      }else{
        varE = private$.varE
      }
      varE = diag(sqrt(varE),
                  nrow=self$nTraits,
                  ncol=self$nTraits)
      varE = varE%*%corE%*%varE
      private$.varE = varE
      invisible(self)
    },
    
    #' @description
    #' Linearly scales all traits to achieve desired 
    #' values of means and variances in the founder population. 
    #' 
    #' @param mean a vector of new trait means
    #' @param var a vector of new trait variances
    #' @param varEnv a vector of new environmental variances
    #' @param varGxE a vector of new GxE variances
    #' @param useVarA tune according to additive genetic variance if true
    #'
    #' @note
    #' By default the founder population is the population used to 
    #' initalize the SimParam object. This population can be changed by 
    #' replacing the population in the founderPop slot. You must run 
    #' \code{\link{resetPop}} on any existing populations to obtain the 
    #' new trait values. 
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$addTraitA(10)
    #' 
    #' #Create population
    #' pop = newPop(founderPop, simParam=SP)
    #' meanG(pop)
    #' 
    #' #Change mean to 1
    #' SP$rescaleTraits(mean=1)
    #' #Run resetPop for change to take effect
    #' pop = resetPop(pop, simParam=SP) 
    #' meanG(pop)
    rescaleTraits = function(mean=0,var=1,varEnv=0,
                             varGxE=1e-6,useVarA=TRUE){
      stopifnot(length(mean)==self$nTraits,
                length(var)==self$nTraits,
                length(varEnv)==self$nTraits,
                length(varGxE)==self$nTraits)
      for(i in 1:self$nTraits){
        trait = private$.traits[[i]]
        trait@intercept = 0
        tmp = calcGenParam(trait,
                           self$founderPop, 
                           self$nThreads)
        if(useVarA){
          scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        }else{
          scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
        }
        trait@addEff = trait@addEff*scale
        if(.hasSlot(trait,"domEff")){
          trait@domEff = trait@domEff*scale
        }
        if(.hasSlot(trait,"epiEff")){
          trait@epiEff[,3] = trait@epiEff[,3]*scale
        }
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        
        if(.hasSlot(trait,"gxeEff")){
          traitG = new("TraitA",
                       nLoci = trait@nLoci,
                       lociPerChr = trait@lociPerChr,
                       lociLoc = trait@lociLoc,
                       addEff = trait@gxeEff,
                       intercept = 0)
          tmpG = calcGenParam(traitG,
                              self$founderPop, 
                              self$nThreads)
          
          if(varEnv[i]==0){
            scaleG = sqrt(varGxE[i])/sqrt(popVar(tmpG$gv)[1])
            trait@gxeEff = trait@gxeEff*scaleG
            trait@gxeInt = 0-mean(tmpG$gv*scaleG)
            trait@envVar = 1
          }else{
            scaleG = sqrt(varGxE[i]/varEnv[i])/sqrt(popVar(tmpG$gv)[1])
            trait@gxeEff = trait@gxeEff*scaleG
            trait@gxeInt = 1-mean(tmpG$gv*scaleG)
            trait@envVar = varEnv[i]
          }
        }
        private$.varA[i] = popVar(tmp$bv*scale)[1]
        private$.varG[i] = popVar(tmp$gv*scale)[1]
        private$.traits[[i]] = trait
      }
      invisible(self)
    },
    
    #### Genetic map (public) ----
    
    #' @field v the crossover interference parameter for a gamma model of 
    #' recombination. A value of 1 indicates no crossover interference 
    #' (e.g. Haldane mapping function). A value of 2.6 approximates the 
    #' degree of crossover interference implied by the Kosambi mapping 
    #' function. (default is 2.6)
    v = "numeric",
    
    #' @field p the proportion of crossovers coming from a non-interfering
    #' pathway. (default is 0)
    p = "numeric",
    
    #' @field quadProb the probability of quadrivalent pairing in an 
    #' autopolyploid. (default is 0)
    quadProb = "numeric",
    
    #' @description Set the relative recombination rates between males 
    #' and females. This allows for sex-specific recombination rates, 
    #' under the assumption of equivalent recombination landscapes.
    #' 
    #' @param femaleRatio relative ratio of recombination in females compared to 
    #' males. A value of 2 indicate twice as much recombination in females. The 
    #' value must be greater than 0. (default is 1)
    #' 
    #' @examples 
    #' #Create founder haplotypes
    #' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    #' 
    #' #Set simulation parameters
    #' SP = SimParam$new(founderPop)
    #' SP$setRecombRatio(2) #Twice as much recombination in females
    setRecombRatio = function(femaleRatio){
      stopifnot(femaleRatio>0)
      genMap = self$genMap
      private$.sepMap = TRUE
      feSc = 2/(1/femaleRatio+1)
      maSc = 2/(femaleRatio+1)
      private$.femaleMap = lapply(genMap,
                                  function(x){
                                    feSc*x
                                  })
      private$.femaleCentromere = feSc*private$.femaleCentromere
      private$.maleMap = lapply(genMap,
                                function(x){
                                  maSc*x
                                })
      private$.maleCentromere = maSc*private$.maleCentromere
      invisible(self)
    },
    
    #' @description 
    #' Replaces existing genetic map.
    #' 
    #' @param genMap a list of length nChr containing 
    #' numeric vectors for the position of each segregating 
    #' site on a chromosome.
    #' @param centromere a numeric vector of centromere 
    #' positions. If NULL, the centromere are assumed to 
    #' be metacentric.
    switchGenMap = function(genMap, centromere=NULL){
      genMap = lapply(genMap, function(x) x-x[1]) # Set position 1 to 0
      if(is.null(centromere)){
        centromere=sapply(genMap,max)/2
      }
      stopifnot(length(genMap)==self$nChr,
                centromere<=sapply(genMap,max))
      tmp = do.call("c",lapply(genMap,length))
      stopifnot(all(tmp==private$.segSites))
      private$.sepMap = FALSE
      private$.femaleMap = genMap
      private$.maleMap = NULL
      private$.femaleCentromere = centromere
      private$.maleCentromere = NULL
      invisible(self)
    },
    
    #' @description 
    #' Replaces existing female genetic map.
    #' 
    #' @param genMap a list of length nChr containing 
    #' numeric vectors for the position of each segregating 
    #' site on a chromosome.
    #' @param centromere a numeric vector of centromere 
    #' positions. If NULL, the centromere are assumed to 
    #' be metacentric.
    switchFemaleMap = function(genMap, centromere=NULL){
      genMap = lapply(genMap, function(x) x-x[1]) # Set position 1 to 0
      if(is.null(centromere)){
        centromere=sapply(genMap,max)/2
      }
      stopifnot(length(genMap)==self$nChr,
                centromere<=sapply(genMap,max))
      tmp = do.call("c",lapply(genMap,length))
      stopifnot(all(tmp==private$.segSites))
      if(private$.sepMap){
        private$.femaleMap = genMap
        private$.femaleCentromere = centromere
      }else{
        private$.sepMap = TRUE
        private$.maleMap = private$.femaleMap
        private$.femaleMap = genMap
        private$.maleCentromere = private$.femaleCentromere
        private$.femaleCentromere = centromere
      }
      invisible(self)
    },
    
    #' @description 
    #' Replaces existing male genetic map.
    #' 
    #' @param genMap a list of length nChr containing 
    #' numeric vectors for the position of each segregating 
    #' site on a chromosome.
    #' @param centromere a numeric vector of centromere 
    #' positions. If NULL, the centromere are assumed to 
    #' be metacentric.
    switchMaleMap = function(genMap, centromere=NULL){
      genMap = lapply(genMap, function(x) x-x[1]) # Set position 1 to 0
      if(is.null(centromere)){
        centromere=sapply(genMap,max)/2
      }
      stopifnot(length(genMap)==self$nChr,
                centromere<=sapply(genMap,max))
      tmp = do.call("c",lapply(genMap,length))
      stopifnot(all(tmp==private$.segSites))
      private$.sepMap = TRUE
      private$.maleMap = genMap
      private$.maleCentromere = centromere
      invisible(self)
    },
    
    #### Internal (public) ----
    
    #' @description For internal use only.
    #' 
    #' @param lastId ID of last individual
    #' @param id the name of each individual
    #' @param mother vector of mother iids
    #' @param father vector of father iids
    #' @param isDH indicator for DH lines
    #' @param hist new recombination history
    #' @param ploidy ploidy level
    addToRec = function(lastId,id,mother,father,isDH,
                        hist,ploidy){
      nNewInd = lastId-private$.lastId
      stopifnot(nNewInd>0)
      if(length(isDH)==1) isDH = rep(isDH,nNewInd)
      mother = as.integer(mother)
      father = as.integer(father)
      isDH = as.integer(isDH)
      stopifnot(length(mother)==nNewInd,
                length(father)==nNewInd,
                length(isDH)==nNewInd)
      tmp = cbind(mother,father,isDH)
      rownames(tmp) = id
      if(is.null(hist)){
        newRecHist = vector("list",nNewInd)
        tmpLastHaplo = private$.lastHaplo
        if(all(isDH==1L)){
          for(i in 1:nNewInd){
            tmpLastHaplo = tmpLastHaplo + 1L
            newRecHist[[i]] = rep(tmpLastHaplo, ploidy)
          }
        }else{
          for(i in 1:nNewInd){
            newRecHist[[i]] = (tmpLastHaplo+1L):(tmpLastHaplo+ploidy)
            tmpLastHaplo = tmpLastHaplo + ploidy
          }
        }
        private$.hasHap = c(private$.hasHap, rep(FALSE, nNewInd))
        private$.isFounder = c(private$.isFounder, rep(TRUE, nNewInd))
        private$.recHist = c(private$.recHist, newRecHist)
        private$.lastHaplo = tmpLastHaplo
      }else{
        # Add hist to recombination history
        private$.hasHap = c(private$.hasHap, rep(FALSE, nNewInd))
        private$.isFounder = c(private$.isFounder, rep(FALSE, nNewInd))
        private$.recHist = c(private$.recHist, hist)
      }
      private$.pedigree = rbind(private$.pedigree, tmp)
      private$.lastId = lastId
      
      invisible(self)
    },
    
    #' @description For internal use only.
    #' 
    #' @param iid internal ID
    ibdHaplo = function(iid){
      if(all(private$.hasHap[iid])){
        # Return relevant haplotypes
        return(private$.hap[as.character(iid)])
      }
      
      ## Fill in missing haplotypes
      
      # Determine unique iid for needed individuals without hap data
      uid = list()
      i = 1L
      uid[[i]] = unique(iid)
      while(any(uid[[i]]!=0L)){
        i = i+1L
        uid[[i]] = unique(c(private$.pedigree[uid[[i-1]],1:2]))
      }
      uid = unique(unlist(uid))
      uid = sort(uid)[-1] # First one is always zero
      uid = uid[!private$.hasHap[uid]]
      
      # Split uid by founder and non-founder
      fuid = uid[private$.isFounder[uid]]
      nfuid = uid[!private$.isFounder[uid]]
      
      # Set hap for founders
      if(length(fuid)>0){
        nChr = length(private$.femaleMap)
        newHap = getFounderIbd(founder=private$.recHist[fuid], 
                               nChr=nChr)
        names(newHap) = as.character(fuid)
        private$.hap = c(private$.hap, newHap)
        private$.hasHap[fuid] = TRUE
      }
      
      # Set hap for non-founders
      if(length(nfuid)>0){
        for(id in nfuid){
          mother = as.character(private$.pedigree[id,1])
          father = as.character(private$.pedigree[id,2])
          private$.hap[[as.character(id)]] = 
            getNonFounderIbd(recHist=private$.recHist[[id ]],
                             mother=private$.hap[[mother]],
                             father=private$.hap[[father]])
          private$.hasHap[id] = TRUE
        }
      }
      
      # Return relevant haplotypes
      return(private$.hap[as.character(iid)])
    },
    
    #' @description For internal use only.
    #' 
    #' @param lastId last ID assigned
    updateLastId = function(lastId){
      lastId = as.integer(lastId)
      stopifnot(lastId>=private$.lastId)
      private$.lastId = lastId
      invisible(self)
    },
    
    #' @description For internal use only.
    #' 
    #' @param lastId ID of last individual
    #' @param id the name of each individual
    #' @param mother vector of mother iids
    #' @param father vector of father iids
    #' @param isDH indicator for DH lines
    addToPed = function(lastId,id,mother,father,isDH){
      nNewInd = lastId-private$.lastId
      stopifnot(nNewInd>0)
      if(length(isDH)==1) isDH = rep(isDH,nNewInd)
      mother = as.integer(mother)
      father = as.integer(father)
      isDH = as.integer(isDH)
      stopifnot(length(mother)==nNewInd,
                length(father)==nNewInd,
                length(isDH)==nNewInd)
      tmp = cbind(mother,father,isDH)
      rownames(tmp) = id
      private$.pedigree = rbind(private$.pedigree,tmp)
      private$.lastId = lastId
      invisible(self)
    }
    
  ),
  private = list(
    #### Private ----
    
    .restrSites="logical",
    .traits="list",
    .segSites="integer",
    .sexes="character",
    .femaleMap="list",
    .maleMap="list",
    .sepMap="logical",
    .femaleCentromere="numeric",
    .maleCentromere="numeric",
    .lastId="integer",
    .isTrackPed="logical",
    .pedigree="matrix",
    .isTrackRec="logical",
    .recHist="list",
    .varA="numeric",
    .varG="numeric",
    .varE="numeric",
    .version="character",
    .lastHaplo="integer",
    .hasHap="logical",
    .hap="list",
    .isFounder="logical",
    
    .isRunning = function(){
      if(private$.lastId==0L){
        invisible(self)
      }else{
        stop("lastId doesn't equal 0, you must run resetPed to proceed")
      }
    },
    
    .addTrait = function(lociMap,varA=NA_real_,varG=NA_real_,varE=NA_real_){
      stopifnot(is.numeric(varA),is.numeric(varG),is.numeric(varE),
                length(varA)==1,length(varG)==1,length(varE)==1)
      private$.traits[[self$nTraits + 1L]] = lociMap
      private$.varA = c(private$.varA,varA)
      private$.varG = c(private$.varG,varG)
      private$.varE = c(private$.varE,varE)
      invisible(self)
    },
    
    .pickLoci = function(nSitesPerChr, QTL=TRUE, minFreq=NULL, refPop=NULL){
      stopifnot(length(nSitesPerChr)==self$nChr)
      
      # Get invalid sites
      if(QTL){
        restr = self$invalidQtl
      }else{
        restr = self$invalidSnp
      }
      
      # Identify potential sites
      pot = vector('list', self$nChr)
      for(i in 1:self$nChr){
        if(is.null(restr[[i]])){
          pot[[i]] = 1:private$.segSites[i]
        }else{
          pot[[i]] = setdiff(1:private$.segSites[i], restr[[i]])
        }
      }
      
      # Filter for minimum frequency
      if(!is.null(minFreq)){
        if(is.null(refPop)){
          refPop = self$founderPop
        }
        for(chr in 1:self$nChr){
          q = calcChrFreq(refPop@geno[[chr]])
          q = 0.5-abs(q-0.5) #Convert to minor allele frequency
          tmp = which(q>=minFreq)
          pot[[chr]] = tmp[tmp%in%pot[[chr]]]
        }
      }
      stopifnot(sapply(pot,length)>=nSitesPerChr)
      
      # Sample sites
      lociLoc = lapply(1:self$nChr,function(x){
        if(nSitesPerChr[x]==0){
          return(NULL)
        }else{
          if(length(pot[[x]])==1){
            tmp = pot[[x]]
          }else{
            tmp = sort(sample(pot[[x]],nSitesPerChr[x]))
          }
          # Add site restrictions
          if(private$.restrSites){
            if(QTL){
              self$invalidSnp[[x]] = sort(union(tmp, self$invalidSnp[[x]]))
            }else{
              self$invalidQtl[[x]] = sort(union(tmp, self$invalidQtl[[x]]))
            }
          }
          return(tmp)
        }
      })
      
      # Create and return a LociMap
      lociLoc = do.call("c",lociLoc)
      loci = new("LociMap",
                 nLoci=as.integer(sum(nSitesPerChr)),
                 lociPerChr=as.integer(nSitesPerChr),
                 lociLoc=as.integer(lociLoc))
      return(loci)
    }
    
  ),
  active = list(
    #### Active ----
    
    #' @field traitNames vector of trait names
    traitNames=function(value){
      if(missing(value)){
        traitNames = sapply(private$.traits, function(x){
          x@name
        })
        return(traitNames)
      }else{
        value = as.character(value)
        if(length(value)!=self$nTraits){
          stop("length of traitNames vector must equal ",self$nTraits)
        }
        for(i in 1:self$nTraits){
          private$.traits[[i]]@name = value[i]
        }
      }
    },
    
    #' @field snpChipNames vector of chip names
    snpChipNames=function(value){
      if(missing(value)){
        snpChipNames = sapply(self$snpChips, function(x){
          x@name
        })
        return(snpChipNames)
      }else{
        value = as.character(value)
        if(length(value)!=self$nSnpChips){
          stop("length of snpChipNames vector must equal ",self$nSnpChips)
        }
        for(i in 1:self$nSnpChips){
          self$snpChips[[i]]@name = value[i]
        }
      }
    },
    
    #' @field traits list of traits
    traits=function(value){
      if(missing(value)){
        private$.traits
      }else{
        stop("`$traits` is read only, see manAddTrait",call.=FALSE)
      }
    },
    
    #' @field nChr number of chromosomes
    nChr=function(value){
      if(missing(value)){
        length(private$.segSites)
      }else{
        stop("`$nChr` is read only",call.=FALSE)
      }
    },
    
    #' @field nTraits number of traits
    nTraits=function(value){
      if(missing(value)){
        length(private$.traits)
      }else{
        stop("`$nTraits` is read only",call.=FALSE)
      }
    },
    
    #' @field nSnpChips number of SNP chips
    nSnpChips=function(value){
      if(missing(value)){
        length(self$snpChips)
      }else{
        stop("`$nSnpChips` is read only",call.=FALSE)
      }
    },
    
    #' @field segSites segregating sites per chromosome
    segSites=function(value){
      if(missing(value)){
        private$.segSites
      }else{
        stop("`$segSites` is read only",call.=FALSE)
      }
    },
    
    #' @field sexes sexes used for mating
    sexes=function(value){
      if(missing(value)){
        private$.sexes
      }else{
        stop("`$sexes` is read only",call.=FALSE)
      }
    },
    
    #' @field sepMap are there seperate genetic maps for 
    #' males and females
    sepMap=function(value){
      if(missing(value)){
        private$.sepMap
      }else{
        stop("`$sepMap` is read only",call.=FALSE)
      }
    },
    
    #' @field genMap "matrix" of chromosome genetic maps
    genMap=function(value){
      if(missing(value)){
        if(private$.sepMap){
          genMap = vector("list",self$nChr)
          for(i in 1:self$nChr){
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
    
    #' @field femaleMap "matrix" of chromosome genetic maps for 
    #' females
    femaleMap=function(value){
      if(missing(value)){
        private$.femaleMap
      }else{
        stop("`$femaleMap` is read only",call.=FALSE)
      }
    },
    
    #' @field maleMap "matrix" of chromosome genetic maps for 
    #' males
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
    
    #' @field centromere position of centromeres genetic map
    centromere=function(value){
      if(missing(value)){
        if(private$.sepMap){
          (private$.femaleCentromere+private$.maleCentromere)/2
        }else{
          private$.femaleCentromere
        }
      }else{
        stop("`$centromere` is read only",call.=FALSE)
      }
    },
    
    #' @field femaleCentromere position of centromeres on female 
    #' genetic map
    femaleCentromere=function(value){
      if(missing(value)){
        private$.femaleCentromere
      }else{
        stop("`$femaleCentromere` is read only",call.=FALSE)
      }
    },
    
    #' @field maleCentromere position of centromeres on male 
    #' genetic map
    maleCentromere=function(value){
      if(missing(value)){
        if(private$.sepMap){
          private$.maleCentromere
        }else{
          private$.femaleCentromere
        }
      }else{
        stop("`$maleCentromere` is read only",call.=FALSE)
      }
    },
    
    #' @field lastId last ID number assigned
    lastId=function(value){
      if(missing(value)){
        private$.lastId
      }else{
        stop("`$lastId` is read only",call.=FALSE)
      }
    },
    
    #' @field isTrackPed is pedigree being tracked
    isTrackPed=function(value){
      if(missing(value)){
        private$.isTrackPed
      }else{
        stop("`$isTrackPed` is read only",call.=FALSE)
      }
    },
    
    #' @field pedigree pedigree matrix for all individuals
    pedigree=function(value){
      if(missing(value)){
        private$.pedigree
      }else{
        stop("`$pedigree` is read only",call.=FALSE)
      }
    },
    
    #' @field isTrackRec is recombination being tracked
    isTrackRec=function(value){
      if(missing(value)){
        private$.isTrackRec
      }else{
        stop("`$isTrackRec` is read only",call.=FALSE)
      }
    },
    
    #' @field recHist list of historic recombination events
    recHist=function(value){
      if(missing(value)){
        private$.recHist
      }else{
        stop("`$recHist` is read only",call.=FALSE)
      }
    },
    
    #' @field haplotypes list of computed IBD haplotypes
    haplotypes=function(value){
      if(missing(value)){
        private$.hap
      }else{
        stop("`$haplotypes` is read only",call.=FALSE)
      }
    },
    
    #' @field varA additive genetic variance in founderPop
    varA=function(value){
      if(missing(value)){
        private$.varA
      }else{
        stop("`$varA` is read only",call.=FALSE)
      }
    },
    
    #' @field varG total genetic variance in founderPop
    varG=function(value){
      if(missing(value)){
        private$.varG
      }else{
        stop("`$varG` is read only",call.=FALSE)
      }
    },
    
    #' @field varE default error variance
    varE=function(value){
      if(missing(value)){
        private$.varE
      }else{
        stop("`$varE` is read only",call.=FALSE)
      }
    },
    
    #' @field version the version of AlphaSimR used to generate this object
    version=function(value){
      if(missing(value)){
        private$.version
      }else{
        stop("`$version` is read only",call.=FALSE)
      }
    }
  )
)

#### External helpers ----
sampAddEff = function(qtlLoci,nTraits,corr,gamma,shape){
  addEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%transMat(corr)
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
                  ncol=nTraits)%*%transMat(corDD)
  domEff = sweep(domEff,2,sqrt(varDD),"*")
  domEff = sweep(domEff,2,meanDD,"+")
  domEff = abs(addEff)*domEff
  return(domEff)
}

sampEpiEff = function(qtlLoci,nTraits,corr,gamma,shape,relVar){
  epiEff = matrix(rnorm(qtlLoci@nLoci*nTraits/2),
                  ncol=nTraits)%*%transMat(corr)
  if(any(gamma)){
    for(i in which(gamma)){
      x = (pnorm(epiEff[,i])-0.5)*2
      epiEff[,i] = sign(x)*qgamma(abs(x),shape=shape)
    }
  }
  epiEff = sweep(epiEff,2,sqrt(relVar),"*")
  return(epiEff)
}

# Test if object is of a SimParam class
isSimParam = function(x) {
  ret = is(x, class2 = "SimParam")
  return(ret)
}
