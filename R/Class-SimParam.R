#' @title Simulation parameters
#'
#' @description
#' Container for global simulation parameters. Saving this object
#' as SP will allow it to be accessed by function defaults.
#'
#' A simulation starts by creating a new \code{SimParam} object with
#' \code{SimParam$new()} and a founder population. It is recommended
#' that you save the object with the name \code{SP}, because subsequent
#' functions will check your global environment for an object of this
#' name if their \code{simParam} arguments are \code{NULL}. This allows
#' you to call these functions without explicitly supplying a
#' \code{simParam} argument with every call.
#'
#' The object is then changed through its methods, such as
#' \code{SP$addTraitA()}, each of which has its own help page.
#'
#' @evalRd .simParamUsage("initialize", call="SimParam$new")
#'
#' @param founderPop an object of \code{\link{MapPop-class}}
#'
#' @details
#' \subsection{Methods by task}{
#' \describe{
#'   \item{Starting a simulation}{
#'     \code{\link[=SimParam-setSexes]{setSexes()}},
#'     \code{\link[=SimParam-setTrackPed]{setTrackPed()}},
#'     \code{\link[=SimParam-setTrackRec]{setTrackRec()}},
#'     \code{\link[=SimParam-resetPed]{resetPed()}}}
#'   \item{Adding traits}{
#'     \code{\link[=SimParam-addTraitA]{addTraitA()}},
#'     \code{\link[=SimParam-addTraitAD]{addTraitAD()}},
#'     \code{\link[=SimParam-addTraitAE]{addTraitAE()}},
#'     \code{\link[=SimParam-addTraitADE]{addTraitADE()}},
#'     \code{\link[=SimParam-addTraitAG]{addTraitAG()}},
#'     \code{\link[=SimParam-addTraitADG]{addTraitADG()}},
#'     \code{\link[=SimParam-addTraitAEG]{addTraitAEG()}},
#'     \code{\link[=SimParam-addTraitADEG]{addTraitADEG()}},
#'     \code{\link[=SimParam-altAddTraitAD]{altAddTraitAD()}},
#'     \code{\link[=SimParam-manAddTrait]{manAddTrait()}},
#'     \code{\link[=SimParam-importTrait]{importTrait()}}.
#'     Each \code{addTrait} method models additive effects plus the
#'     effects named in its suffix: D for dominance, E for epistasis and
#'     G for genotype-by-environment interaction.}
#'   \item{Phenotypes}{
#'     \code{\link[=SimParam-setVarE]{setVarE()}},
#'     \code{\link[=SimParam-setCorE]{setCorE()}} (deprecated)}
#'   \item{Choosing loci and SNP chips}{
#'     \code{\link[=SimParam-restrSegSites]{restrSegSites()}},
#'     \code{\link[=SimParam-addSnpChip]{addSnpChip()}},
#'     \code{\link[=SimParam-addSnpChipByName]{addSnpChipByName()}},
#'     \code{\link[=SimParam-addStructuredSnpChip]{addStructuredSnpChip()}},
#'     \code{\link[=SimParam-setFounderHap]{setFounderHap()}}}
#'   \item{Genetic maps and recombination}{
#'     \code{\link[=SimParam-setRecombRatio]{setRecombRatio()}},
#'     \code{\link[=SimParam-switchGenMap]{switchGenMap()}},
#'     \code{\link[=SimParam-switchFemaleMap]{switchFemaleMap()}},
#'     \code{\link[=SimParam-switchMaleMap]{switchMaleMap()}}}
#'   \item{Changing traits}{
#'     \code{\link[=SimParam-switchTrait]{switchTrait()}},
#'     \code{\link[=SimParam-removeTrait]{removeTrait()}},
#'     \code{\link[=SimParam-rescaleTraits]{rescaleTraits()}}}
#' }
#' The methods \code{addToRec()}, \code{ibdHaplo()},
#' \code{updateLastId()} and \code{addToPed()} are for use by
#' AlphaSimR itself and are not documented.
#' }
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Settings for how individuals are created and recorded. Make these
#' #before creating any populations.
#' SP$setTrackPed(TRUE) #Record the pedigree
#' SP$setTrackRec(TRUE) #Also record recombination
#' SP$setSexes("yes_sys") #Assign sexes systematically
#' SP$setRecombRatio(2) #Twice as much recombination in females
#'
#' #Keep SNPs and QTL on different segregating sites
#' SP$restrSegSites(minQtlPerChr=10, minSnpPerChr=10)
#'
#' #Add SNP chips, chosen at random or by marker name
#' SP$addSnpChip(5)
#' SP$addSnpChipByName(getSnpMap(snpChip=1, simParam=SP)$id[1:2])
#'
#' #Add traits, one for each trait model
#' SP$addTraitA(10)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$addTraitAE(10, relAA=0.1)
#' SP$addTraitADE(10)
#' SP$addTraitAG(10, varGxE=2)
#' SP$addTraitADG(10, meanDD=0.5, varGxE=2)
#' SP$addTraitAEG(10, varGxE=2)
#' SP$addTraitADEG(10, meanDD=0.5, varGxE=2)
#'
#' #Set default error variances from heritabilities, used by setPheno
#' SP$setVarE(h2=rep(0.5, SP$nTraits))
#'
#' #Create a population
#' pop = newPop(founderPop, simParam=SP)
#' meanP(pop)
#'
#' @export
SimParam = R6Class(
  "SimParam",
  public = list(
    #### Public ----

    #' @field snpChips list of SNP chips
    snpChips = "list",

    #' @field invalidQtl list of segregating sites that aren't valid QTL
    invalidQtl = "list",

    #' @field invalidSnp list of segregating sites that aren't valid SNP
    invalidSnp = "list",

    #' @field founderPop founder population used for variance scaling
    founderPop = "MapPop",

    #' @field finalizePop function applied to newly created populations.
    #' It is run as the last step of creating a new population and
    #' provides a way to automatically modify populations.
    #' The function must satisfy these four requirements:
    #' 1) The first argument is a \code{\link{Pop-class}} object
    #'    which will be used in this function;
    #' 2) The second argument is a SimParam object
    #'    to provide access to global parameters;
    #' 3) The remaining arguments are passed via \code{...}; and
    #' 4) The return is a \code{\link{Pop-class}} object.
    #' See \code{\link{asCategorical}} for an example.
    finalizePop = "function",

    #' @field allowEmptyPop if true, population arguments with nInd=0 will
    #' return an empty population with a warning instead of an error.
    allowEmptyPop = "logical",

    #' @field finalizePheno function applied to newly generated phenotype values.
    #' It is run as the last step of generating new phenotype values and
    #' provides a way to automatically modify phenotype values.
    #' While \code{finalizePop} function can be used to a similar effect,
    #' \code{finalizePheno} function provides more granular support over
    #' finalizing phenotype values in new and existing populations, say,
    #' when we repeatedly call \code{\link{setPheno}} on the existing population
    #' (which would not trigger use of \code{finalizePop}).
    #' The function must satisfy these five requirements:
    #' 1) The first argument is a matrix of phenotype values,
    #'    which will be used in this function and
    #'    should have dimensions as \code{pop@pheno};
    #' 2) The second argument is a \code{\link{Pop-class}} object
    #'    to provide access to population specific values;
    #' 3) The third argument must be a SimParam object
    #'    to provide access to global parameters;
    #' 4) The remaining arguments must be passed via \code{...}; and
    #' 5) The return is a matrix of phenotype values that will form the new
    #'    \code{pop@pheno}, so the function should not change column or
    #'    row order.
    #' See \code{\link{asCategorical}} for an example.
    finalizePheno = "function",

    initialize = function(founderPop){
      stopifnot("founderPop must be a MapPop, such as one made by runMacs, quickHaplo or importHaplo"=
                  is(founderPop, "MapPop"))

      # Public items
      self$nThreads = getNumThreads()
      self$v = 2.6 # Kosambi
      self$p = 0 # Single pathway gamma model
      self$quadProb = 0 # No quadrivalent pairing
      self$snpChips = list()
      self$invalidQtl = rep(list(integer()), founderPop@nChr) # All eligible
      self$invalidSnp = rep(list(integer()), founderPop@nChr) # All eligible
      self$founderPop = founderPop
      self$finalizePop = function(pop, ...){return(pop)}
      self$allowEmptyPop = FALSE # Empty populations trigger an error
      self$finalizePheno = function(pheno, ...){return(pheno)}

      # Private items
      private$.restrSites = TRUE
      private$.traits = list()
      private$.segSites = founderPop@nLoci
      private$.sexes = "no"
      .checkGenMap(founderPop@genMap)
      private$.femaleMap = founderPop@genMap
      private$.maleMap = NULL
      private$.sepMap = FALSE
      private$.femaleCentromere = founderPop@centromere
      private$.maleCentromere = NULL
      private$.femaleGenLen = .getGenLen(founderPop)
      private$.maleGenLen = NULL
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
      private$.activeQtl = NULL
      private$.qtlIndex = list()

      invisible(self)
    },

    setTrackPed = function(isTrackPed, force=FALSE){
      stopifnot("isTrackPed must be TRUE or FALSE"=is.logical(isTrackPed))
      if(!force){
        private$.isRunning()
      }
      private$.isTrackPed = isTrackPed
      if(!isTrackPed){
        private$.isTrackRec = FALSE
      }
      invisible(self)
    },

    setTrackRec = function(isTrackRec, force=FALSE){
      stopifnot("isTrackRec must be TRUE or FALSE"=is.logical(isTrackRec))
      if(!force){
        private$.isRunning()
      }
      private$.isTrackRec = isTrackRec
      if(isTrackRec){
        private$.isTrackPed = TRUE
      }
      invisible(self)
    },

    resetPed =function(lastId=0L){
      private$.lastId = lastId
      private$.pedigree = private$.pedigree[0:lastId,,drop=FALSE]
      if(private$.isTrackRec){
        private$.recHist = private$.recHist[0:lastId]
      }
      # .hasHap, .isFounder and .hap are indexed by iid, so they have to be
      # cut back with the pedigree. Left in place, the next individual to
      # take a reused iid inherits the old individual's haplotype record.
      private$.hasHap = private$.hasHap[seq_len(lastId)]
      private$.isFounder = private$.isFounder[seq_len(lastId)]
      if(length(private$.hap)>lastId){
        private$.hap = private$.hap[seq_len(lastId)]
      }
      if(lastId==0L){
        private$.lastHaplo = 0L
      }
      invisible(self)
    },

    restrSegSites = function(minQtlPerChr=NULL, minSnpPerChr=NULL, excludeQtl=NULL,
                             excludeSnp=NULL, overlap=FALSE, minSnpFreq=NULL){
      # Handle any named QTL exclusions
      if(!is.null(excludeQtl)){
        matchList = private$.findNamedLoci(excludeQtl)

        # Make exclusions
        restr = self$invalidQtl
        for(i in seq_len(self$nChr)){
          restr[[i]] = sort(union(restr[[i]], matchList[[i]]))
        }
        self$invalidQtl = restr
      }

      # Handle any named SNP exclusions
      if(!is.null(excludeSnp)){
        # Check if the SNP list matches the QTL list
        # If so, save time by using the previously determined exclusions
        findMatch = TRUE
        if(!is.null(excludeQtl)){
          if(all(excludeSnp==excludeQtl)){
            findMatch = FALSE
          }
        }

        if(findMatch){
          matchList = private$.findNamedLoci(excludeSnp)
        }

        # Make exclusions
        restr = self$invalidSnp
        for(i in seq_len(self$nChr)){
          restr[[i]] = sort(union(restr[[i]], matchList[[i]]))
        }
        self$invalidSnp = restr
      }

      if(overlap){
        # Not setting any restrictions if overlap is allow
        # Existing restrictions will be left in place
        private$.restrSites = FALSE
        invisible(self)
      }else{
        # Check validity of inputs
        if(!is.null(minSnpPerChr)){
          if(length(minSnpPerChr)==1){
            minSnpPerChr = rep(minSnpPerChr,self$nChr)
          }else{
            stopifnot("minSnpPerChr must have length 1 or nChr"=
                        length(minSnpPerChr)==self$nChr)
          }
        }
        if(!is.null(minQtlPerChr)){
          if(length(minQtlPerChr)==1){
            minQtlPerChr = rep(minQtlPerChr,self$nChr)
          }else{
            stopifnot("minQtlPerChr must have length 1 or nChr"=
                        length(minQtlPerChr)==self$nChr)
          }
        }

        # Restrict SNPs and then QTL
        # SNPs are done first due to  potentially fewer viable loci
        private$.restrSites = TRUE
        if(!is.null(minSnpPerChr)){
          invisible(private$.pickLoci(minSnpPerChr, FALSE, minSnpFreq))
        }
        if(!is.null(minQtlPerChr)){
          invisible(private$.pickLoci(minQtlPerChr))
        }
        invisible(self)
      }
    },

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

    setFounderHap = function(hapMap){
      n = length(hapMap)
      if(is.null(names(hapMap))){
        iid = seq_len(n)
      }else{
        iid = suppressWarnings(as.integer(names(hapMap)))
        if(anyNA(iid) || anyDuplicated(iid) || any(iid<1L | iid>n)){
          stop("names of hapMap must be the internal IDs 1 to ", n)
        }
      }
      # The cache is indexed by iid and held packed
      hap = vector("list", n)
      hap[iid] = packHist(hapMap)
      private$.hap = hap
      private$.hasHap = rep(TRUE, n)
      private$.isFounder = rep(FALSE, n)
      invisible(self)
    },

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

    addSnpChipByName = function(markers, name=NULL){
      genMap = private$.femaleMap

      # Check that the markers are present on the map
      genMapMarkerNames = unlist(lapply(genMap, names))
      missing = setdiff(markers, genMapMarkerNames)
      if(length(missing)>0L){
        stop("markers not found on the genetic map: ",
             paste(head(missing, 5L), collapse=", "),
             if(length(missing)>5L) paste0(" and ", length(missing)-5L, " more"))
      }

      # Create lociPerChr and lociLoc
      lociPerChr = integer(length(genMap))
      lociLoc = vector("list", length(genMap))

      # Loop through chromosomes
      for(i in seq_len(length(genMap))){

        # Initialize lociLoc
        lociLoc[[i]] = integer()

        # Find matches if they exist
        take = match(names(genMap[[i]]), markers)
        lociPerChr[i] = length(na.omit(take))
        if(lociPerChr[i]>0L){
          lociLoc[[i]] = which(!is.na(take))
        }
      }
      lociLoc = unlist(lociLoc)

      snpChip = new("LociMap",
                    nLoci=sum(lociPerChr),
                    lociPerChr=lociPerChr,
                    lociLoc=lociLoc)

      if(is.null(name)){
        snpChip@name = paste0("Chip",self$nSnpChips + 1L)
      }else{
        snpChip@name = name
      }

      self$snpChips[[self$nSnpChips + 1L]] = snpChip
      invisible(self)
    },

    addStructuredSnpChip = function(nSnpPerChr,structure,force=FALSE){
      if(!force){
        private$.isRunning()
      }
      if(length(nSnpPerChr)==1){
        nSnpPerChr = rep(nSnpPerChr,self$nChr)
      }
      stopifnot("nSnpPerChr must have length 1 or nChr"=
                  length(nSnpPerChr)==self$nChr)
      # There is no potSnp field, so the eligible sites are worked out
      # here the same way restrSegSites does it: every segregating site on
      # the chromosome that has not been marked an invalid SNP.
      potSnp = lapply(seq_len(self$nChr), function(i){
        setdiff(1:private$.segSites[i], self$invalidSnp[[i]])
      })
      stopifnot("nSnpPerChr is more than the sites available for SNPs on a chromosome"=
                  vapply(potSnp,length,integer(1))>=nSnpPerChr)
      stopifnot("structure must have one column per SNP, sum(nSnpPerChr)"=
                  dim(structure)[2]==sum(nSnpPerChr))
      lociLoc = lapply(1:self$nChr,function(x){
        sort(sample(potSnp[[x]],nSnpPerChr[x]))
      })
      lociLoc = do.call("c",lociLoc)

      for (i in seq_len(nrow(structure))){
        snps = lociLoc[structure[i,]]
        start = 1
        numChr = numeric(length(nSnpPerChr))
        for (j in seq_len(length(nSnpPerChr))){
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

    addTraitA = function(nQtlPerChr,mean=0,var=1,corA=NULL,
                         gamma=FALSE,shape=1,force=FALSE,name=NULL,
                         nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
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
      stopifnot("mean and var must have the same length"=
                  length(mean)==length(var),
                "name must have one value per trait"=
                  length(mean)==length(name))
      checkCorMat(corA, nTraits)
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      for(i in seq_len(nTraits)){
        trait = new("TraitA",
                    qtlLoci,
                    addEff=addEff[,i],
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, nThreads)
        scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        trait@addEff = trait@addEff*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)
        private$.addTrait(trait,var[i],var[i])
      }
      invisible(self)
    },

    addTraitAD = function(nQtlPerChr,mean=0,var=1,meanDD=0,
                          varDD=0,corA=NULL,corDD=NULL,useVarA=TRUE,
                          gamma=FALSE,shape=1,force=FALSE,name=NULL,
                          nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
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
      stopifnot("mean and var must have the same length"=
                  length(mean)==length(var),
                "name must have one value per trait"=
                  length(mean)==length(name))
      checkCorMat(corA, nTraits)
      checkCorMat(corDD, nTraits)
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                          corDD=corDD,meanDD=meanDD,varDD=varDD)
      for(i in seq_len(nTraits)){
        trait = new("TraitAD",
                    qtlLoci,
                    addEff=addEff[,i],
                    domEff=domEff[,i],
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, nThreads)
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

    altAddTraitAD = function(nQtlPerChr,mean=0,varA=1,varD=0,inbrDepr=0,
                             limMeanDD=c(0,1.5),limVarDD=c(0,0.5),
                             silent=FALSE,force=FALSE,name=NULL,
                             nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
      }
      if(length(nQtlPerChr)==1){
        nQtlPerChr = rep(nQtlPerChr,self$nChr)
      }
      if(is.null(name)){
        name = paste0("Trait",self$nTraits+1)
      }

      # Pick QTL
      qtlLoci = private$.pickLoci(nQtlPerChr)

      # Create list of arguments for optimization
      argsList = argAltAD(LociMap = qtlLoci,
                          Pop = self$founderPop,
                          mean = mean,
                          varA = varA,
                          varD = varD,
                          inbrDepr = inbrDepr,
                          nThreads = nThreads)

      # Run optim to optimize meanDD and varDD
      optOut = optim(par = c(mean(limMeanDD), mean(sqrt(limVarDD))),
                     fn = objAltAD,
                     gr = NULL,
                     method = "L-BFGS-B",
                     lower = c(limMeanDD[1], sqrt(limVarDD[1])),
                     upper = c(limMeanDD[2], sqrt(limVarDD[2])),
                     args = argsList)

      # Finalize creation of trait
      output = finAltAD(input = optOut$par, args = argsList)
      trait = new("TraitAD",
                  qtlLoci,
                  addEff=c(output$a),
                  domEff=c(output$d),
                  intercept=c(output$intercept),
                  name=name)
      private$.addTrait(trait,varA,output$varG)

      # Report trait details
      if(!silent){
        cat("A new trait called", name, "was added. \n")
        cat("   varD =", output$varD, "\n")
        cat("   inbrDepr =", output$inbrDepr, "\n")
        cat("   meanDD =", output$meanDD, "\n")
        cat("   varDD =", output$varDD, "\n")
      }

      invisible(self)
    },


    addTraitAG = function(nQtlPerChr,mean=0,var=1,varGxE=1e-6,varEnv=0,
                          corA=NULL,corGxE=NULL,gamma=FALSE,shape=1,
                          force=FALSE,name=NULL,nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
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
      stopifnot("mean and var must have the same length"=
                  length(mean)==length(var),
                "varGxE must have one value per trait"=
                  length(mean)==length(varGxE),
                "varEnv must have one value per trait"=
                  length(mean)==length(varEnv),
                "name must have one value per trait"=
                  length(mean)==length(name))
      checkCorMat(corA, nTraits)
      checkCorMat(corGxE, nTraits)
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corGxE,gamma=FALSE,shape=NULL)
      for(i in seq_len(nTraits)){
        trait = new("TraitA",
                    qtlLoci,
                    addEff=addEff[,i],
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, nThreads)
        scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        trait@addEff = trait@addEff*scale
        trait@intercept = mean[i]-mean(tmp$gv*scale)

        # GxE component
        traitG = new("TraitA",
                     qtlLoci,
                     addEff=gxeEff[,i],
                     intercept=0)
        tmpG = calcGenParam(traitG, self$founderPop, nThreads)
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

    addTraitADG = function(nQtlPerChr,mean=0,var=1,varEnv=0,
                           varGxE=1e-6,meanDD=0,varDD=0,corA=NULL,
                           corDD=NULL,corGxE=NULL,useVarA=TRUE,gamma=FALSE,
                           shape=1,force=FALSE,name=NULL,nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
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
      stopifnot("mean and var must have the same length"=
                  length(mean)==length(var),
                "varGxE must have one value per trait"=
                  length(varGxE)==nTraits,
                "varEnv must have one value per trait"=
                  length(varEnv)==nTraits,
                "name must have one value per trait"=
                  length(mean)==length(name))
      checkCorMat(corA, nTraits)
      checkCorMat(corDD, nTraits)
      checkCorMat(corGxE, nTraits)
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                          corDD=corDD,meanDD=meanDD,varDD=varDD)
      gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corGxE,gamma=FALSE,shape=NULL)
      for(i in seq_len(nTraits)){
        trait = new("TraitAD",
                    qtlLoci,
                    addEff=addEff[,i],
                    domEff=domEff[,i],
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, nThreads)
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
        tmpG = calcGenParam(traitG, self$founderPop, nThreads)
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

    addTraitAE = function(nQtlPerChr,mean=0,var=1,relAA=0,corA=NULL,
                          corAA=NULL,useVarA=TRUE,gamma=FALSE,shape=1,force=FALSE,
                          name=NULL,nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
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
      stopifnot("mean and var must have the same length"=
                  length(mean)==length(var),
                "relAA must have one value per trait"=
                  length(relAA)==length(mean),
                "the total number of QTL must be even, because epistatic QTL are paired"=
                  (sum(nQtlPerChr)%%2L)==0L,
                "name must have one value per trait"=
                  length(mean)==length(name))
      checkCorMat(corA, nTraits)
      checkCorMat(corAA, nTraits)
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corAA,gamma=gamma,shape=shape,
                          relVar=relAA)
      E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
      for(i in seq_len(nTraits)){
        trait = new("TraitAE",
                    qtlLoci,
                    addEff=addEff[,i],
                    epiEff=cbind(E,epiEff[,i]),
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, nThreads)
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

    addTraitADE = function(nQtlPerChr,mean=0,var=1,meanDD=0,
                           varDD=0,relAA=0,corA=NULL,corDD=NULL,corAA=NULL,
                           useVarA=TRUE,gamma=FALSE,shape=1,force=FALSE,
                           name=NULL,nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
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
      stopifnot("mean and var must have the same length"=
                  length(mean)==length(var),
                "relAA must have one value per trait"=
                  length(relAA)==length(mean),
                "the total number of QTL must be even, because epistatic QTL are paired"=
                  (sum(nQtlPerChr)%%2L)==0L,
                "name must have one value per trait"=
                  length(mean)==length(name))
      checkCorMat(corA, nTraits)
      checkCorMat(corDD, nTraits)
      checkCorMat(corAA, nTraits)
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                          corDD=corDD,meanDD=meanDD,varDD=varDD)
      epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corAA,gamma=gamma,shape=shape,
                          relVar=relAA)
      E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
      for(i in seq_len(nTraits)){
        trait = new("TraitADE",
                    qtlLoci,
                    addEff=addEff[,i],
                    domEff=domEff[,i],
                    epiEff=cbind(E,epiEff[,i]),
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, nThreads)
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

    addTraitAEG = function(nQtlPerChr,mean=0,var=1,relAA=0,varGxE=1e-6,varEnv=0,
                           corA=NULL,corAA=NULL,corGxE=NULL,useVarA=TRUE,gamma=FALSE,
                           shape=1,force=FALSE,name=NULL,nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
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
      stopifnot("mean and var must have the same length"=
                  length(mean)==length(var),
                "relAA must have one value per trait"=
                  length(relAA)==length(mean),
                "varGxE must have one value per trait"=
                  length(mean)==length(varGxE),
                "varEnv must have one value per trait"=
                  length(mean)==length(varEnv),
                "the total number of QTL must be even, because epistatic QTL are paired"=
                  (sum(nQtlPerChr)%%2L)==0L,
                "name must have one value per trait"=
                  length(mean)==length(name))
      checkCorMat(corA, nTraits)
      checkCorMat(corAA, nTraits)
      checkCorMat(corGxE, nTraits)
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corAA,gamma=gamma,shape=shape,
                          relVar=relAA)
      E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
      gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corGxE,gamma=FALSE,shape=NULL)
      for(i in seq_len(nTraits)){
        trait = new("TraitAE",
                    qtlLoci,
                    addEff=addEff[,i],
                    epiEff=cbind(E,epiEff[,i]),
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, nThreads)
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
        tmpG = calcGenParam(traitG, self$founderPop, nThreads)
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

    addTraitADEG = function(nQtlPerChr,mean=0,var=1,varEnv=0,
                            varGxE=1e-6,meanDD=0,varDD=0,relAA=0,corA=NULL,
                            corDD=NULL,corAA=NULL,corGxE=NULL,useVarA=TRUE,
                            gamma=FALSE,shape=1,force=FALSE,name=NULL,
                            nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
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
      stopifnot("mean and var must have the same length"=
                  length(mean)==length(var),
                "varGxE must have one value per trait"=
                  length(varGxE)==nTraits,
                "varEnv must have one value per trait"=
                  length(varEnv)==nTraits,
                "relAA must have one value per trait"=
                  length(relAA)==length(mean),
                "the total number of QTL must be even, because epistatic QTL are paired"=
                  (sum(nQtlPerChr)%%2L)==0L,
                "name must have one value per trait"=
                  length(mean)==length(name))
      checkCorMat(corA, nTraits)
      checkCorMat(corDD, nTraits)
      checkCorMat(corAA, nTraits)
      checkCorMat(corGxE, nTraits)
      qtlLoci = private$.pickLoci(nQtlPerChr)
      addEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corA,gamma=gamma,shape=shape)
      domEff = sampDomEff(qtlLoci=qtlLoci,nTraits=nTraits,addEff=addEff,
                          corDD=corDD,meanDD=meanDD,varDD=varDD)
      epiEff = sampEpiEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corAA,gamma=gamma,shape=shape,
                          relVar=relAA)
      E = matrix(sample.int(sum(nQtlPerChr),sum(nQtlPerChr)),ncol=2)
      gxeEff = sampAddEff(qtlLoci=qtlLoci,nTraits=nTraits,
                          corr=corGxE,gamma=FALSE,shape=NULL)
      for(i in seq_len(nTraits)){
        trait = new("TraitADE",
                    qtlLoci,
                    addEff=addEff[,i],
                    domEff=domEff[,i],
                    epiEff=cbind(E,epiEff[,i]),
                    intercept=0,
                    name=name[i])
        tmp = calcGenParam(trait, self$founderPop, nThreads)
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
        tmpG = calcGenParam(traitG, self$founderPop, nThreads)
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

    manAddTrait = function(lociMap,varE=NA_real_,force=FALSE,nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
      }
      stopifnot("lociMap must be a LociMap, such as a TraitA object"=
                  is(lociMap,"LociMap"))
      tmp = calcGenParam(lociMap, self$founderPop, nThreads)
      varA = popVar(tmp$bv)[1]
      varG = popVar(tmp$gv)[1]
      private$.addTrait(lociMap,varA,varG,varE)
      invisible(self)
    },

    importTrait = function(markerNames,
                           addEff,
                           domEff=NULL,
                           intercept=NULL,
                           name=NULL,
                           varE=NULL,
                           force=FALSE,
                           nThreads=NULL,
                           epiPairs=NULL,
                           epiEff=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
      }

      # Duplicated names would be silently dropped by the match against
      # the genetic map below, leaving their effects unused. Factors are
      # converted so that names, not level codes, are matched to pairs.
      markerNames = as.character(markerNames)
      if(anyDuplicated(markerNames)){
        stop("markerNames must be unique",call.=FALSE)
      }

      # Check addEff and domEff inputs
      addEff = as.matrix(addEff)
      stopifnot("addEff must have one row per marker in markerNames"=
                  length(markerNames)==nrow(addEff))
      nTraits = ncol(addEff)
      if(is.null(domEff)){
        useDom = FALSE
      }else{
        useDom = TRUE
        domEff = as.matrix(domEff)
        stopifnot("domEff must have the same number of rows as addEff"=
                    nrow(addEff)==nrow(domEff),
                  "domEff must have the same number of traits as addEff"=
                    ncol(addEff)==ncol(domEff))
      }

      # Check epiPairs and epiEff inputs
      if(is.null(epiPairs)!=is.null(epiEff)){
        stop("epiPairs and epiEff must be supplied together",call.=FALSE)
      }
      useEpi = !is.null(epiPairs)
      if(useEpi){
        if(is.data.frame(epiPairs)){
          epiPairs = as.matrix(epiPairs)
        }
        if(!is.matrix(epiPairs) || ncol(epiPairs)!=2L ||
           !is.character(epiPairs)){
          stop("epiPairs must be a two-column matrix or data.frame of marker names",
               call.=FALSE)
        }
        pairedNames = as.vector(epiPairs)
        if(anyNA(pairedNames) || !all(pairedNames%in%markerNames)){
          stop("every marker in epiPairs must be in markerNames",call.=FALSE)
        }
        # The epistatic model requires each locus to be in exactly one
        # pair. Since markerNames is unique and contains every paired name,
        # equal lengths and no repeats make the pairing a partition.
        if(anyDuplicated(pairedNames) ||
           length(pairedNames)!=length(markerNames)){
          stop("each marker in markerNames must appear in exactly one pair of epiPairs",
               call.=FALSE)
        }
        epiEff = as.matrix(epiEff)
        if(nrow(epiEff)!=nrow(epiPairs) || ncol(epiEff)!=nTraits){
          stop("epiEff must have one row for each pair in epiPairs and one column for each trait",
               call.=FALSE)
        }
      }

      # Prepare the intercept
      if(is.null(intercept)){
        intercept = rep(0, nTraits)
      }else{
        intercept = as.numeric(intercept)
        stopifnot("intercept must have one value per trait"=
                    length(intercept)==nTraits)
      }

      # Prepare varE
      if(!is.null(varE)){
        varE = as.numeric(varE)
        stopifnot("varE must have one value per trait"=
                    length(varE)==nTraits)
      }else{
        varE = rep(NA_real_, nTraits)
      }

      # Prepare trait names
      if(is.null(name)){
        name = paste0("Trait",1:nTraits+self$nTraits)
      }else{
        stopifnot("name must have one value per trait"=
                    length(name)==nTraits)
      }

      # Extract genetic map and check if marker names are on the map
      genMapMarkerNames = unlist(lapply(private$.femaleMap, names))
      missing = setdiff(markerNames, genMapMarkerNames)
      if(length(missing)>0L){
        stop("markerNames not found on the genetic map: ",
             paste(head(missing, 5L), collapse=", "),
             if(length(missing)>5L) paste0(" and ", length(missing)-5L, " more"))
      }

      # Create trait variables
      lociPerChr = integer(self$nChr)
      lociLoc = lociNames = vector("list", self$nChr)
      addEffList = domEffList = vector("list", nTraits)
      for(i in seq_len(nTraits)){
        addEffList[[i]] = domEffList[[i]] = vector("list", self$nChr)
      }

      # Loop through chromosomes
      for(i in seq_len(self$nChr)){
        # Working on trait 1
        # Initialize variables
        addEffList[[1]][[i]] = domEffList[[1]][[i]] = numeric()
        lociLoc[[i]] = integer()
        lociNames[[i]] = character()

        # Find matches if they exist
        take = match(names(private$.femaleMap[[i]]), markerNames)
        lociPerChr[i] = length(na.omit(take))

        if(lociPerChr[i]>0L){
          lociLoc[[i]] = which(!is.na(take))
          lociNames[[i]] = markerNames[na.omit(take)]
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

      # Pair members index the trait's loci, which are now in map order
      # rather than the order of markerNames, so match them by name.
      if(useEpi){
        lociNames = unlist(lociNames)
        epiLoci = cbind(match(epiPairs[,1], lociNames),
                        match(epiPairs[,2], lociNames))
      }

      # Create Trait(s)
      for(i in seq_len(nTraits)){
        addEff = unlist(addEffList[[i]])
        if(useEpi){
          E = cbind(epiLoci, unname(epiEff[,i]))
          dimnames(E) = NULL
        }
        if(useDom && useEpi){
          domEff = unlist(domEffList[[i]])
          trait = new("TraitADE",
                      addEff=addEff,
                      domEff=domEff,
                      epiEff=E,
                      intercept=intercept[i],
                      nLoci=nLoci,
                      lociPerChr=lociPerChr,
                      lociLoc=lociLoc,
                      name=name[i])
        }else if(useEpi){
          trait = new("TraitAE",
                      addEff=addEff,
                      epiEff=E,
                      intercept=intercept[i],
                      nLoci=nLoci,
                      lociPerChr=lociPerChr,
                      lociLoc=lociLoc,
                      name=name[i])
        }else if(useDom){
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
        self$manAddTrait(lociMap=trait, varE=varE[i], force=force,
                         nThreads=nThreads)
      }

      invisible(self)
    },

    switchTrait = function(traitPos,lociMap,varE=NA_real_,force=FALSE,
                           nThreads=NULL){
      if(!force){
        private$.isRunning()
      }
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
      }
      traitPos = convertTraitsToIndex(traitPos, self)
      stopifnot("lociMap must be a LociMap, such as a TraitA object"=
                  is(lociMap,"LociMap"),
                "traitPos must be a single trait"=
                  length(traitPos)==1L,
                "traitPos must be between 1 and nTraits"=
                  traitPos<=self$nTraits && traitPos>0)
      tmp = calcGenParam(lociMap, self$founderPop, nThreads)
      private$.traits[[traitPos]] = lociMap
      private$.varA[traitPos] = popVar(tmp$bv)[1]
      private$.varG[traitPos] = popVar(tmp$gv)[1]
      if(is.matrix(private$.varE)){
        warning("Error correlations have been removed, use setVarE to reinstate")
        private$.varE = diag(private$.varE)
      }
      private$.varE[traitPos] = varE

      # Update activeQtl
      private$.setActiveQtl()

      invisible(self)
    },

    removeTrait = function(traits,force=FALSE){
      if(!force){
        private$.isRunning()
      }
      traits = convertTraitsToIndex(traits, self)
      stopifnot("traits must name at least one trait"=length(traits)>0L,
                "traits must be between 1 and nTraits"=
                  max(traits)<=self$nTraits && min(traits)>0)
      private$.traits = private$.traits[-traits]
      private$.varA = private$.varA[-traits]
      private$.varG = private$.varG[-traits]
      if(is.matrix(private$.varE)){
        private$.varE = private$.varE[-traits,-traits]
      }else{
        private$.varE = private$.varE[-traits]
      }

      # Update activeQtl
      private$.setActiveQtl()

      invisible(self)
    },

    setVarE = function(h2=NULL, H2=NULL, varE=NULL, corE=NULL){
      # Check validity of corE, if supplied
      if(!is.null(corE)){
        checkCorMat(corE, self$nTraits)
      }

      # Set error variances (.varE)
      if(!is.null(h2)){
        stopifnot("h2 must have one value per trait"=
                    length(h2)==self$nTraits,
                  "h2 needs every trait to have genetic variance in the founders"=
                    all(private$.varG>0),
                  "h2 needs every trait to have additive variance in the founders; use H2 or varE instead"=
                    all(private$.varA>0))
        varE = numeric(self$nTraits)
        for(i in seq_len(length(h2))){
          tmp = private$.varA[i]/h2[i]-private$.varG[i]
          if(tmp<0){
            stop(paste0("h2=",h2[i]," is not possible for trait ",i))
          }
          varE[i] = tmp
        }
        private$.varE = varE
      }else if(!is.null(H2)){
        stopifnot("H2 must have one value per trait"=
                    length(H2)==self$nTraits)
        varE = numeric(self$nTraits)
        for(i in seq_len(length(H2))){
          tmp = private$.varG[i]/H2[i]-private$.varG[i]
          varE[i] = tmp
        }
        private$.varE = varE
      }else if(!is.null(varE)){
        if(is.matrix(varE)){
          stopifnot("a varE matrix must have one row and one column per trait"=
                      nrow(varE)==self$nTraits && ncol(varE)==self$nTraits)
        }else{
          stopifnot("varE must have one value per trait"=
                      length(varE)==self$nTraits)
        }
        private$.varE = varE
      }else{
        private$.varE = rep(NA_real_, self$nTraits)
      }

      # Set error correlations
      if(!is.null(corE)){
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
      }

      invisible(self)
    },

    setCorE = function(corE){
      warning("This function has been deprecated. Use simParam$setVarE instead.")
      checkCorMat(corE, self$nTraits)
      stopifnot("setCorE needs varE set for every trait first, with setVarE"=
                  length(private$.varE)==self$nTraits && !any(is.na(private$.varE)))
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

    rescaleTraits = function(mean=0,var=1,varEnv=0,
                             varGxE=1e-6,useVarA=TRUE,nThreads=NULL){
      if(is.null(nThreads)){
        nThreads = self$nThreads
      }else{
        nThreads = as.integer(nThreads)
      }
      stopifnot("mean must have one value per trait"=
                  length(mean)==self$nTraits,
                "var must have one value per trait"=
                  length(var)==self$nTraits,
                "varEnv must have one value per trait"=
                  length(varEnv)==self$nTraits,
                "varGxE must have one value per trait"=
                  length(varGxE)==self$nTraits)
      for(i in seq_len(self$nTraits)){
        trait = private$.traits[[i]]
        trait@intercept = 0
        tmp = calcGenParam(trait,
                           self$founderPop,
                           nThreads)
        if(useVarA){
          scale = sqrt(var[i])/sqrt(popVar(tmp$bv)[1])
        }else{
          scale = sqrt(var[i])/sqrt(popVar(tmp$gv)[1])
        }
        trait@addEff = trait@addEff*scale
        if(.hasSlot(trait,"addEffMale")){
          # Parent of origin traits carry a second set of additive
          # effects. Leaving them unscaled puts the two sets on
          # different scales and breaks the requested variance.
          trait@addEffMale = trait@addEffMale*scale
        }
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
                              nThreads)

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

    setRecombRatio = function(femaleRatio){
      stopifnot("femaleRatio must be greater than 0"=femaleRatio>0)
      # Both must be read before .sepMap is set, because the accessors
      # switch to averaging the sex-specific values once it is TRUE
      genMap = self$genMap
      centromere = self$centromere
      genLen = self$genLen
      private$.sepMap = TRUE
      feSc = 2/(1/femaleRatio+1)
      maSc = 2/(femaleRatio+1)
      private$.femaleMap = lapply(genMap,
                                  function(x){
                                    feSc*x
                                  })
      private$.femaleCentromere = feSc*centromere
      private$.maleMap = lapply(genMap,
                                function(x){
                                  maSc*x
                                })
      private$.maleCentromere = maSc*centromere
      private$.femaleGenLen = feSc*genLen
      private$.maleGenLen = maSc*genLen
      invisible(self)
    },

    switchGenMap = function(genMap, centromere=NULL, genLen=NULL){
      if(is.data.frame(genMap)){
        genMap = importGenMap(genMap)
      }else{
        .checkGenMap(genMap)
      }

      genLen = .checkGenLen(genLen, genMap)
      if(is.null(centromere)){
        centromere=genLen/2
      }
      stopifnot("genMap must have one element per chromosome"=
                  length(genMap)==self$nChr,
                "centromere must not be negative"=
                  all(centromere>=0),
                "centromere must not be past the end of its chromosome"=
                  all(centromere<=genLen))

      tmp = do.call("c",lapply(genMap,length))
      stopifnot("genMap must have the same number of loci on each chromosome as the founder population"=
                  all(tmp==private$.segSites))

      private$.sepMap = FALSE
      private$.femaleMap = genMap
      private$.maleMap = NULL
      private$.femaleCentromere = centromere
      private$.maleCentromere = NULL
      private$.femaleGenLen = genLen
      private$.maleGenLen = NULL
      invisible(self)
    },

    switchFemaleMap = function(genMap, centromere=NULL, genLen=NULL){
      if(is.data.frame(genMap)){
        genMap = importGenMap(genMap)
      }else{
        .checkGenMap(genMap)
      }

      genLen = .checkGenLen(genLen, genMap)
      if(is.null(centromere)){
        centromere=genLen/2
      }
      stopifnot("genMap must have one element per chromosome"=
                  length(genMap)==self$nChr,
                "centromere must not be negative"=
                  all(centromere>=0),
                "centromere must not be past the end of its chromosome"=
                  all(centromere<=genLen))

      tmp = do.call("c",lapply(genMap,length))
      stopifnot("genMap must have the same number of loci on each chromosome as the founder population"=
                  all(tmp==private$.segSites))

      if(private$.sepMap){
        private$.femaleMap = genMap
        private$.femaleCentromere = centromere
        private$.femaleGenLen = genLen
      }else{
        private$.sepMap = TRUE
        private$.maleMap = private$.femaleMap
        private$.femaleMap = genMap
        private$.maleCentromere = private$.femaleCentromere
        private$.femaleCentromere = centromere
        private$.maleGenLen = private$.femaleGenLen
        private$.femaleGenLen = genLen
      }
      invisible(self)
    },

    switchMaleMap = function(genMap, centromere=NULL, genLen=NULL){
      if(is.data.frame(genMap)){
        genMap = importGenMap(genMap)
      }else{
        .checkGenMap(genMap)
      }

      genLen = .checkGenLen(genLen, genMap)
      if(is.null(centromere)){
        centromere=genLen/2
      }
      stopifnot("genMap must have one element per chromosome"=
                  length(genMap)==self$nChr,
                "centromere must not be negative"=
                  all(centromere>=0),
                "centromere must not be past the end of its chromosome"=
                  all(centromere<=genLen))

      tmp = do.call("c",lapply(genMap,length))
      stopifnot("genMap must have the same number of loci on each chromosome as the founder population"=
                  all(tmp==private$.segSites))

      private$.sepMap = TRUE
      private$.maleMap = genMap
      private$.maleCentromere = centromere
      private$.maleGenLen = genLen
      invisible(self)
    },

    #### Internal (public) ----

    # For internal use only. Records new individuals and their
    # recombination history; called by .newPop.
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
          for(i in seq_len(nNewInd)){
            tmpLastHaplo = tmpLastHaplo + 1L
            newRecHist[[i]] = rep(tmpLastHaplo, ploidy)
          }
        }else{
          for(i in seq_len(nNewInd)){
            newRecHist[[i]] = (tmpLastHaplo+1L):(tmpLastHaplo+ploidy)
            tmpLastHaplo = tmpLastHaplo + ploidy
          }
        }
        private$.hasHap = c(private$.hasHap, rep(FALSE, nNewInd))
        private$.isFounder = c(private$.isFounder, rep(TRUE, nNewInd))
        private$.recHist = c(private$.recHist, newRecHist)
        private$.lastHaplo = tmpLastHaplo
      }else{
        # Histories are stored packed, one integer vector per individual,
        # which is what meiosis returns. A history built in R as nested
        # lists, as doubleGenome and mergeGenome do, is packed here. The
        # IDs are not stored with it, because the pedigree's row names
        # already hold them.
        stopifnot(length(hist)==nNewInd)
        if(is.list(hist[[1]])){
          hist = packHist(hist)
        }
        private$.hasHap = c(private$.hasHap, rep(FALSE, nNewInd))
        private$.isFounder = c(private$.isFounder, rep(FALSE, nNewInd))
        private$.recHist = c(private$.recHist, unname(hist))
      }
      private$.pedigree = rbind(private$.pedigree, tmp)
      private$.lastId = lastId

      invisible(self)
    },

    # For internal use only. Returns IBD haplotypes by iid; called by
    # pullIbdHaplo.
    ibdHaplo = function(iid){
      # The cache is a list indexed by iid. It was once keyed by name, and
      # every lookup and insertion then searched or copied the whole list,
      # which made filling it quadratic in the size of the pedigree.
      if(length(private$.hap)<private$.lastId){
        length(private$.hap) = private$.lastId
      }

      if(!all(private$.hasHap[iid])){
        # Find the individuals lacking IBD haplotypes that are needed:
        # those requested and, recursively, their parents. The ancestors of
        # an individual that already has them are not needed.
        need = logical(private$.lastId)
        cur = unique(iid)
        while(length(cur)>0L){
          cur = cur[cur>0L]
          cur = cur[!need[cur] & !(private$.hasHap[cur] %in% TRUE)]
          need[cur] = TRUE
          cur = unique(c(private$.pedigree[cur,1:2]))
        }
        uid = which(need)

        # Parents always have smaller iids than their progeny, so a single
        # pass in iid order computes every parent before it is used
        private$.hap[uid] = computeIbd(uid=uid,
                                       mother=private$.pedigree[uid,1],
                                       father=private$.pedigree[uid,2],
                                       isFounder=private$.isFounder[uid],
                                       recHist=private$.recHist[uid],
                                       hap=private$.hap,
                                       nChr=length(private$.femaleMap))
        private$.hasHap[uid] = TRUE
      }

      # Return relevant haplotypes, packed as by packHist
      return(private$.hap[iid])
    },

    # For internal use only. Moves the ID counter on when the pedigree
    # is not tracked; called by .newPop.
    updateLastId = function(lastId){
      lastId = as.integer(lastId)
      stopifnot(lastId>=private$.lastId)
      private$.lastId = lastId
      invisible(self)
    },

    # For internal use only. Records new individuals in the pedigree;
    # called by .newPop.
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

    .nThreads="integer",
    .restrSites="logical",
    .traits="list",
    .segSites="integer",
    .sexes="character",
    .femaleMap="list",
    .maleMap="list",
    .sepMap="logical",
    .femaleCentromere="numeric",
    .maleCentromere="numeric",
    .femaleGenLen="numeric",
    .maleGenLen="numeric",
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
    .activeQtl = "LociMap",
    .qtlIndex = "list",

    # Determines whether not a simulation has started using lastId as an indicator
    .isRunning = function(){
      if(private$.lastId==0L){
        invisible(self)
      }else{
        stop("lastId doesn't equal 0, you must run resetPed to proceed")
      }
    },

    # Adds a trait to simulation and ensures all fields are propagated
    .addTrait = function(lociMap,varA=NA_real_,varG=NA_real_,varE=NA_real_){
      stopifnot(is.numeric(varA),is.numeric(varG),is.numeric(varE),
                length(varA)==1,length(varG)==1,length(varE)==1)

      if(self$nTraits == 0L){
        # Initializing .activeQtl
        private$.activeQtl = new("LociMap",
                                 nLoci = lociMap@nLoci,
                                 lociPerChr = lociMap@lociPerChr,
                                 lociLoc = lociMap@lociLoc)
        private$.qtlIndex[[1L]] = seq_len(private$.activeQtl@nLoci)
      }else{
        # Check if exiting activeQtl is a superset
        testSet = findLociMapSuperset(private$.activeQtl,
                                      lociMap)
        if(is.null(testSet)){
          # Existing activeQtl is a superset, just add qtlIndex
          private$.qtlIndex[[self$nTraits + 1L]] = findQtlIndex(private$.activeQtl, lociMap)
        }else{
          # New activeQtl object required, recompute all values for qtlIndex
          private$.activeQtl = testSet

          for(i in seq_len(self$nTraits)){
            private$.qtlIndex[[i]] = findQtlIndex(private$.activeQtl, self$traits[[i]])
          }

          private$.qtlIndex[[self$nTraits + 1L]] =
            findQtlIndex(private$.activeQtl, lociMap)
        }
      }

      private$.traits[[self$nTraits + 1L]] = lociMap
      private$.varA = c(private$.varA,varA)
      private$.varG = c(private$.varG,varG)
      # .varE holds a covariance matrix when setVarE was given one, and
      # c() on a matrix flattens it into a long vector. Reduce it to the
      # per trait variances before appending the new trait.
      if(is.matrix(private$.varE)){
        private$.varE = diag(private$.varE)
      }
      private$.varE = c(private$.varE,varE)
      invisible(self)
    },

    # Samples eligible loci for traits or SNP chips and ensures that they
    # are added to the exclusion list when applicable
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
      for(i in seq_len(self$nChr)){
        pot[[i]] = setdiff(1:private$.segSites[i], restr[[i]])
      }

      # Filter for minimum frequency
      if(!is.null(minFreq)){
        if(is.null(refPop)){
          refPop = self$founderPop
        }
        for(chr in seq_len(self$nChr)){
          q = calcChrFreq(refPop@geno[[chr]])
          q = 0.5-abs(q-0.5) #Convert to minor allele frequency
          tmp = which(q>=minFreq)
          pot[[chr]] = tmp[tmp%in%pot[[chr]]]
        }
      }

      stopifnot("Not enough eligible sites, see ?SimParam and method restrSegSites()" = sapply(pot,length)>=nSitesPerChr)

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
    },

    # Returns physical positions of named loci in a list format
    # Input order is not preserved. This function is intended as
    # a helper for restrSegSites
    .findNamedLoci = function(lociNames){
      # Loci names
      id = unlist(unname(lapply(private$.femaleMap, names)))
      take = match(lociNames, id)
      if(any(is.na(take))){
        stop("One or more loci are not on the genetic map. Beware of case sensitivity.")
      }

      # Find positions using an interval search strategy on the cumulative sum
      take = unique(take)
      posList = rep(list(integer()), self$nChr)
      cumSumSegSite = cumsum(private$.segSites)
      for(i in take){
        # Identify chromosome
        chr = findInterval(i, cumSumSegSite, left.open = TRUE) + 1L

        # Identify position
        if(chr>1L){
          pos = i - cumSumSegSite[chr-1L]
        }else{
          pos = i
        }

        # Add site to list
        posList[[chr]] = c(posList[[chr]], pos)
      }

      return(posList)
    },

    # Sets activeQtl and qtlIndex slots
    # Run after removing or modifying existing traits
    .setActiveQtl = function(){

      private$.activeQtl = NULL
      private$.qtlIndex = list()

      if(self$nTraits > 0L){
        # Set activeQtl using trait 1
        private$.activeQtl = new("LociMap",
                                 nLoci = self$traits[[1L]]@nLoci,
                                 lociPerChr = self$traits[[1L]]@lociPerChr,
                                 lociLoc = self$traits[[1L]]@lociLoc)

        private$.qtlIndex[[1L]] = seq_len(private$.activeQtl@nLoci)
      }

      if(self$nTraits > 1L){
        # Update activeQtl, if needed
        for(i in 2:self$nTraits){
          testSet = findLociMapSuperset(private$.activeQtl,
                                        self$traits[[i]])
          if(!is.null(testSet)){
            private$.activeQtl = testSet
          }
        }

        # Update qtlIndex
        for(i in seq_len(self$nTraits)){
          private$.qtlIndex[[i]] = findQtlIndex(private$.activeQtl,
                                                self$traits[[i]])
        }
      }

      invisible(self)
    }

  ),
  active = list(
    #### Active ----

    #' @field nThreads number of threads used with OpenMP (when available).
    #' Assign \code{NULL} to reset to \code{getNumThreads()}.
    #' See \code{vignette("parallelization", package="AlphaSimR")}
    #' for setup details.
    nThreads=function(value){
      if(missing(value)){
        private$.nThreads
      }else{
        if(is.null(value)){
          value = getNumThreads()
        }
        valueInt = as.integer(value)
        if(!is.numeric(value) || length(value)!=1L || is.na(value) ||
           !is.finite(value) || value < 1L || is.na(valueInt) ||
           value != valueInt){
          stop("`$nThreads` must be a single positive integer or NULL to reset",
               call.=FALSE)
        }
        private$.nThreads = valueInt
      }
    },

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
        for(i in seq_len(self$nTraits)){
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
        for(i in seq_len(self$nSnpChips)){
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

    #' @field sepMap are there separate genetic maps for
    #' males and females
    sepMap=function(value){
      if(missing(value)){
        private$.sepMap
      }else{
        stop("`$sepMap` is read only",call.=FALSE)
      }
    },

    #' @field genMap list of chromosome genetic maps
    genMap=function(value){
      if(missing(value)){
        if(private$.sepMap){
          genMap = vector("list",self$nChr)
          for(i in seq_len(self$nChr)){
            genMap[[i]] = (private$.femaleMap[[i]]+
                             private$.maleMap[[i]])/2
          }
          # Chromosome names are needed by getGenMap
          names(genMap) = names(private$.femaleMap)
          genMap
        }else{
          private$.femaleMap
        }
      }else{
        stop("`$genMap` is read only",call.=FALSE)
      }
    },

    #' @field femaleMap list of chromosome genetic maps for
    #' females
    femaleMap=function(value){
      if(missing(value)){
        private$.femaleMap
      }else{
        stop("`$femaleMap` is read only",call.=FALSE)
      }
    },

    #' @field maleMap list of chromosome genetic maps for
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

    #' @field centromere position of centromeres on genetic map
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

    #' @field genLen genetic length of each chromosome. A chromosome
    #' runs from position 0 to its genetic length, which can lie past
    #' its last locus.
    genLen=function(value){
      if(missing(value)){
        if(private$.sepMap){
          (private$.femaleGenLen+private$.maleGenLen)/2
        }else{
          private$.femaleGenLen
        }
      }else{
        stop("`$genLen` is read only",call.=FALSE)
      }
    },

    #' @field femaleGenLen genetic length of each chromosome on the
    #' female genetic map
    femaleGenLen=function(value){
      if(missing(value)){
        private$.femaleGenLen
      }else{
        stop("`$femaleGenLen` is read only",call.=FALSE)
      }
    },

    #' @field maleGenLen genetic length of each chromosome on the
    #' male genetic map
    maleGenLen=function(value){
      if(missing(value)){
        if(private$.sepMap){
          private$.maleGenLen
        }else{
          private$.femaleGenLen
        }
      }else{
        stop("`$maleGenLen` is read only",call.=FALSE)
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
        # Histories are stored packed and are expanded here into the
        # nested lists, named by ID, that this field has always returned
        output = unpackHist(private$.recHist)
        if(length(output)>0L){
          names(output) = rownames(private$.pedigree)[seq_along(output)]
        }
        output
      }else{
        stop("`$recHist` is read only",call.=FALSE)
      }
    },

    #' @field haplotypes list of computed IBD haplotypes, named by
    #' internal ID
    haplotypes=function(value){
      if(missing(value)){
        # Stored packed and indexed by iid, so expanded and named here
        take = which(private$.hasHap)
        take = take[take<=length(private$.hap)]
        output = unpackHist(private$.hap[take])
        if(length(output)>0L){
          names(output) = as.character(take)
        }
        output
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
    },

    #' @field activeQtl a LociMap representing all active QTL in simulation
    activeQtl=function(value){
      if(missing(value)){
        private$.activeQtl
      }else{
        stop("`$activeQtl` is read only",call.=FALSE)
      }
    },

    #' @field qtlIndex a list of vectors giving trait specific QTL indices
    #' relative to all active QTL
    qtlIndex=function(value){
      if(missing(value)){
        private$.qtlIndex
      }else{
        stop("`$qtlIndex` is read only",call.=FALSE)
      }
    }
  )
)

#### External helpers ----

#' @title Check a correlation matrix
#'
#' @description Stops with an informative error if a trait correlation
#' argument is not a valid correlation matrix of the right size
#'
#' @param R the matrix supplied by the user
#' @param nTraits number of traits
#'
#' @details
#' The matrix must be numeric, nTraits by nTraits, free of missing values,
#' symmetric, with ones on the diagonal and every entry in [-1, 1].
#' Positive semi-definiteness is not required here: \code{\link{transMat}}
#' smooths a matrix that fails it and warns.
#'
#' @returns invisibly TRUE
#'
#' @keywords internal
checkCorMat = function(R, nTraits){
  nameR = deparse(substitute(R))
  if(!is.matrix(R) || !is.numeric(R)){
    stop(nameR, " must be a numeric matrix")
  }
  if(nrow(R)!=nTraits || ncol(R)!=nTraits){
    stop(nameR, " must be a ", nTraits, " by ", nTraits,
         " matrix, with one row and column per trait")
  }
  if(anyNA(R)){
    stop(nameR, " contains missing values")
  }
  if(!isSymmetric(R)){
    stop(nameR, " is not a symmetric matrix")
  }
  # A covariance matrix passes every check above, so the diagonal is what
  # tells the two apart. The tolerance allows for a matrix that was itself
  # computed, for example by cov2cor.
  if(!isTRUE(all.equal(unname(diag(R)), rep(1, nTraits)))){
    stop(nameR, " must have ones on its diagonal to be a correlation matrix")
  }
  if(any(abs(R) > 1+sqrt(.Machine$double.eps))){
    stop(nameR, " has entries outside [-1, 1]")
  }
  invisible(TRUE)
}

#' @title Symmetric gamma transform
#'
#' @description Maps standard normal deviates to gamma deviates with a
#' random sign, preserving their order
#'
#' @param z standard normal deviates
#' @param shape gamma distribution shape parameter
#'
#' @returns a vector the same length as z
#'
#' @keywords internal
gammaTransform = function(z, shape){
  # |z| maps to the uniform deviate 2*pnorm(|z|)-1 and that to a gamma
  # quantile, while the sign of z is kept. Working with the upper tail on the
  # log scale keeps precision for large |z|, where 2*pnorm(|z|)-1 rounds to
  # one and qgamma would return Inf.
  sign(z)*qgamma(log(2)+pnorm(-abs(z),log.p=TRUE),shape=shape,
                 lower.tail=FALSE,log.p=TRUE)
}

#' @title Hermite coefficients of the symmetric gamma transform
#'
#' @description Computes the coefficients of \code{\link{gammaTransform}}
#' in the basis of orthonormal Hermite polynomials
#'
#' @param shape gamma distribution shape parameter
#' @param nCoef number of coefficients, starting from order one
#'
#' @returns a vector of length nCoef whose k-th element is the coefficient
#' of order k
#'
#' @keywords internal
gammaHermCoef = function(shape, nCoef=51L){
  # The transform is odd, so the coefficients of even order are zero and
  # only odd orders are integrated. The product of the transform and an odd
  # Hermite polynomial is even, so each integral is twice its value over
  # the positive half line. The normal density underflows well before 40.
  coef = numeric(nCoef)
  for(k in seq(1L, nCoef, by=2L)){
    integrand = function(z){
      # Orthonormal Hermite polynomials from the three-term recurrence,
      # which stays stable at high order where the explicit form does not
      h0 = rep(1, length(z))
      h1 = z
      m = 1L
      while(m<k){
        h2 = (z*h1-sqrt(m)*h0)/sqrt(m+1)
        h0 = h1
        h1 = h2
        m = m+1L
      }
      2*gammaTransform(z,shape)*h1*dnorm(z)
    }
    coef[k] = integrate(integrand, 0, 40, subdivisions=2000L,
                        rel.tol=1e-10)$value
  }
  return(coef)
}

#' @title Latent correlation for gamma distributed effects
#'
#' @description Finds the correlation matrix for the normal deviates that
#' gives the requested correlation after gamma distributed traits pass
#' through \code{\link{gammaTransform}}
#'
#' @param R the requested correlation matrix
#' @param gamma indicator of whether trait should use a gamma distribution
#' @param shape gamma distribution shape parameter
#' @param nCoef number of Hermite coefficients used for each trait
#'
#' @details
#' By Mehler's formula, when two standard normal deviates have correlation
#' rho, transformed deviates have correlation
#' sum_k b_ik b_jk rho^k / (sd_i sd_j), where b_ik are the Hermite
#' coefficients of trait i's transform. That series is inverted for each
#' pair of traits with \code{uniroot}. Pairs of normal traits and pairs with
#' a requested correlation of zero are left unchanged.
#'
#' Traits whose effects have different distributions cannot reach every
#' correlation. When a requested correlation is beyond the attainable
#' maximum, the maximum is used and a warning is given. The returned matrix
#' may not be positive semi-definite, which \code{\link{transMat}} handles.
#'
#' @returns a correlation matrix with the same dimensions as R
#'
#' @keywords internal
gammaLatentCorr = function(R, gamma, shape, nCoef=51L){
  nTraits = ncol(R)
  # A normal trait's transform is the identity, which is the Hermite
  # polynomial of order one alone. Traits sharing a shape share their
  # coefficients, so each shape is integrated once.
  coef = matrix(0, nrow=nCoef, ncol=nTraits)
  coef[1,] = 1
  for(s in unique(shape[gamma])){
    take = gamma & shape==s
    coef[,take] = gammaHermCoef(s, nCoef)
  }
  # Normalising by the variance the truncated series captures, rather than
  # the exact variance, makes two traits with the same shape reach exactly
  # one at rho=1. A requested correlation of one then stays attainable.
  coef = sweep(coef, 2, sqrt(colSums(coef^2)), "/")
  Rz = R
  for(i in seq_len(nTraits-1L)){
    for(j in (i+1L):nTraits){
      if(!(gamma[i] || gamma[j]) || R[i,j]==0){
        next
      }
      w = coef[,i]*coef[,j]
      corFun = function(rho) sum(w*rho^seq_len(nCoef))
      # Only odd orders are nonzero, so the attainable range is symmetric
      maxCor = corFun(1)
      if(abs(R[i,j])>=maxCor){
        if(abs(R[i,j])>maxCor+1e-8){
          warning(sprintf(paste("A correlation of %.3f between traits %d and",
                                "%d cannot be reached with these effect",
                                "distributions. The attainable limit of %.3f",
                                "is used instead."),
                          R[i,j], i, j, sign(R[i,j])*maxCor))
        }
        Rz[i,j] = sign(R[i,j])
      }else{
        Rz[i,j] = uniroot(function(rho) corFun(rho)-R[i,j],
                          c(-1,1), tol=1e-10)$root
      }
      Rz[j,i] = Rz[i,j]
    }
  }
  return(Rz)
}

#' @title Sample additive genotypic effects
#'
#' @description Samples deviates from a normal distribution or gamma distribution
#' with a random sign
#'
#' @param qtlLoci total number of loci
#' @param nTraits number of traits
#' @param corr correlation between traits
#' @param gamma indicator of whether trait should use a gamma distribution
#' @param shape gamma distribution shape parameter
#'
#' @details See \code{vignette("traits", package = "AlphaSimR")} for the
#' implemented trait model and genotypic effects.
#'
#' @returns a matrix with dimensions qtlLoci by nTraits
#'
#' @keywords internal
sampAddEff = function(qtlLoci,nTraits,corr,gamma,shape){
  # The gamma transform is nonlinear, which pulls correlations towards zero.
  # Sampling the normal deviates at an adjusted correlation undoes this, so
  # the transformed effects have the requested correlation.
  if(any(gamma)){
    corr = gammaLatentCorr(corr,gamma,shape)
  }
  addEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%transMat(corr)
  if(any(gamma)){
    for(i in which(gamma)){
      # shape holds one value per trait, so trait i takes shape[i].
      # Passing the whole vector recycles it down the loci instead.
      addEff[,i] = gammaTransform(addEff[,i],shape[i])
    }
  }
  return(addEff)
}

#' @title Sample dominance genotypic effects
#'
#' @description Samples dominance deviation effects from a normal distribution
#' and uses previously sampled additive genotypic effects to form dominance
#' effects
#'
#' @param qtlLoci total number of loci
#' @param nTraits number of traits
#' @param addEff previously sampled additive genotypic effects
#' @param corDD correlation between dominance degrees
#' @param meanDD mean value of dominance degrees
#' @param varDD variance of dominance degrees
#'
#' @details See \code{vignette("traits", package = "AlphaSimR")} for the
#' implemented trait model and genotypic effects.
#'
#' @returns a matrix with dimensions qtlLoci by nTraits
#'
#' @keywords internal
sampDomEff = function(qtlLoci,nTraits,addEff,corDD,
                      meanDD,varDD){
  domEff = matrix(rnorm(qtlLoci@nLoci*nTraits),
                  ncol=nTraits)%*%transMat(corDD)
  domEff = sweep(domEff,2,sqrt(varDD),"*")
  domEff = sweep(domEff,2,meanDD,"+")
  domEff = abs(addEff)*domEff
  return(domEff)
}

#' @title Sample epistatic effects
#'
#' @description Samples epistatic effects from a normal distribution or gamma distribution
#' with a variance relative to the variance of the additive genotypic effects
#'
#' @param qtlLoci total number of loci
#' @param nTraits number of traits
#' @param corr correlation between epistatic effects
#' @param gamma indicator of whether trait should use a gamma distribution
#' @param shape gamma distribution shape parameter
#' @param relVar desired variance for epistatic effects
#'
#' @details See \code{vignette("traits", package = "AlphaSimR")} for the
#' implemented trait model and genotypic effects.
#'
#' @returns a matrix with dimensions qtlLoci by nTraits
#'
#' @keywords internal
sampEpiEff = function(qtlLoci,nTraits,corr,gamma,shape,relVar){
  # As in sampAddEff, the correlation is adjusted for the gamma transform
  if(any(gamma)){
    corr = gammaLatentCorr(corr,gamma,shape)
  }
  epiEff = matrix(rnorm(qtlLoci@nLoci*nTraits/2),
                  ncol=nTraits)%*%transMat(corr)
  if(any(gamma)){
    for(i in which(gamma)){
      # As in sampAddEff, shape is per trait
      epiEff[,i] = gammaTransform(epiEff[,i],shape[i])
    }
  }
  epiEff = sweep(epiEff,2,sqrt(relVar),"*")
  return(epiEff)
}

# Test if object is of SimParam class
isSimParam = function(x) {
  ret = is(x, class2 = "SimParam")
  return(ret)
}

#' @title Find LociMap superset
#'
#' @description Compares two \code{\link{LociMap-class}} objects to determine if
#' the first one is a superset of the second. If it is, the function returns NULL.
#' If it is not, the function returns a \code{\link{LociMap-class}} object that is
#' a superset of both \code{\link{LociMap-class}} objects.
#'
#' @param lociMap1 a \code{\link{LociMap-class}} that is tested to determine if
#' it is a superset
#' @param lociMap2 a second \code{\link{LociMap-class}} that is tested
#'
#' @returns NULL if lociMap1 is a superset, or a \code{\link{LociMap-class}} if it
#' is not
#'
#' @keywords internal
findLociMapSuperset = function(lociMap1, lociMap2){

  loc1 = paste(rep(x = seq_len(length(lociMap1@lociPerChr)),
                   times = lociMap1@lociPerChr),
               lociMap1@lociLoc, sep="_")

  loc2 = paste(rep(x = seq_len(length(lociMap2@lociPerChr)),
                   times = lociMap2@lociPerChr),
               lociMap2@lociLoc, sep="_")

  inSet1 = loc2 %in% loc1

  if(all(inSet1)){
    # lociMap1 is a superset of lociMap2, no action needed
    return(NULL)
  }

  # Create an ordered matrix with columns chr and lociLoc for loci superset
  loc = unique(c(loc1, loc2))
  tab = sapply(loc, function(x){
    y = strsplit(x, split="_", fixed=TRUE)
    return(as.integer(unlist(y)))
  })
  tab = t(tab)
  tab = tab[order(tab[,1],tab[,2]),]

  # Return LociMap for superset of loci
  return(new("LociMap",
             nLoci = nrow(tab),
             lociPerChr = tabulate(tab[,1], nbins=length(lociMap1@lociPerChr)),
             lociLoc = unname(tab[,2])))
}


#' @title Find trait QTL index
#'
#' @description Finds the positions of a trait's QTL within the set of all
#' active QTL. The QTL of the trait are assumed to be a subset of the active
#' QTL, so every trait QTL is expected to be matched.
#'
#' @param activeQtl a \code{\link{LociMap-class}} representing all active QTL
#' @param traitQtl a \code{\link{LociMap-class}} representing QTL for a trait of
#' interest
#'
#' @returns an integer vector for relative positions of traitQtl on activeQtl
#'
#' @keywords internal
findQtlIndex = function(activeQtl, traitQtl){

  activeLoci = paste(rep(x = seq_len(length(activeQtl@lociPerChr)),
                         times = activeQtl@lociPerChr),
                     activeQtl@lociLoc, sep="_")

  traitLoci = paste(rep(x = seq_len(length(traitQtl@lociPerChr)),
                        times = traitQtl@lociPerChr),
                    traitQtl@lociLoc, sep="_")

  return(match(traitLoci, activeLoci))
}

# The genetic length of each chromosome on the female or male map of a
# SimParam. A SimParam saved before genLen existed has no such field, and
# reading it gives NULL, so it falls back to the position of each
# chromosome's last locus, which is where those chromosomes ended.
.simGenLen = function(simParam, female=TRUE){
  if(female){
    genLen = simParam$femaleGenLen
    if(is.null(genLen)){
      genLen = .lastLocus(simParam$femaleMap)
    }
  }else{
    genLen = simParam$maleGenLen
    if(is.null(genLen)){
      genLen = .lastLocus(simParam$maleMap)
    }
  }
  return(genLen)
}
