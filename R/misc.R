# Converts a matrix to integer type
# Intended for genotype matrices of raw type
convToImat = function(X){
  return(matrix(as.integer(X), nrow=nrow(X) ,ncol=ncol(X)))
}

#' @rdname isFemale
#' @title Test if individuals of a population are female or male
#'
#' @description Test if individuals of a population are female or male
#'
#' @param x \code{\link{Pop-class}}
#'
#' @return logical
#'
#' @examples
#' founderGenomes <- quickHaplo(nInd = 3, nChr = 1, segSites = 100)
#' SP <- SimParam$new(founderGenomes)
#' SP$setSexes(sexes = "yes_sys")
#' pop <- newPop(founderGenomes)
#' 
#' isFemale(pop)
#' isMale(pop)
#' 
#' pop[isFemale(pop)]
#' pop[isFemale(pop)]@sex
#' 
#' @export
isFemale <- function(x) {
  if (isPop(x)) {
    ret <- x@sex == "F"
  } else {
    stop("Argument x must be a Pop class object!")
  }
  return(ret)
}

##Test

#' @describeIn isFemale Test if individuals of a population are female or male
#' @export
isMale <- function(x) {
  if (isPop(x)) {
    ret <- x@sex == "M"
  } else {
    stop("Argument x must be a Pop class object!")
  }
  return(ret)
}


#' @rdname setMisc
#' @title Set miscelaneous information in a population
#'
#' @description Set miscelaneous information in a population
#'
#' @param x \code{\link{Pop-class}}
#' @param node character, name of the node to set within the \code{x@misc} slot
#' @param value, value to be saved into \code{x@misc[[*]][[node]]}; length of
#'   \code{value} should be equal to \code{nInd(x)}; if its length is 1, then
#'   it is repeated using \code{rep} (see examples)
#'
#' @details A \code{NULL} in \code{value} is ignored
#' 
#' @return \code{\link{Pop-class}} with \code{x@misc[[*]][[node]]} set
#' basePop <- newPop(founderGenomes)
#'
#' basePop <- setMisc(basePop, node = "info", value = 1)
#' basePop@misc
#' getMisc(x = basePop, node = "info")
#'
#' basePop <- setMisc(basePop, node = "info2", value = c("A", "B", "C"))
#' basePop@misc
#' getMisc(x = basePop, node = "info2")
#' 
#' n <- nInd(basePop)
#' location <- vector(mode = "list", length = n)
#' for (ind in seq_len(n)) {
#'   location[[ind]] <- runif(n = 2, min = 0, max = 100)
#' }
#' location
#' basePop <- setMisc(basePop, node = "location", value = location)
#' basePop@misc
#' getMisc(x = basePop, node = "location")
#' 
#' n <- nInd(basePop)
#' location <- vector(mode = "list", length = n)
#' for (ind in c(1, 3)) {
#'   location[[ind]] <- runif(n = 2, min = 0, max = 100)
#' }
#' location
#' basePop <- setMisc(basePop, node = "location", value = location)
#' basePop@misc
#' getMisc(x = basePop, node = "location")
#' 
#' getMisc(x = basePop)
#'
#' @export
setMisc <- function(x, node = NULL, value = NULL) {
  if (isPop(x)) {
    if (is.null(node)) {
      stop("Argument node must be provided!")
    }
    if (is.null(value)) {
      stop("Argument value must be provided!")
    }
    n <- nInd(x)
    if (length(value) == 1 && n > 1) {
      value <- rep(x = value, times = n)
    }
    if (length(value) != n) {
      stop("Argument value must be of length 1 or nInd(x)!")
    }
    for (ind in seq_len(n)) {
      if (!is.null(value[ind])) {
        x@misc[[ind]][node] <- value[ind]
      }
    }
  } else {
    stop("Argument x must be a Pop class object!")
  }
  return(x)
}

#' @rdname getMisc
#' @title Get miscelaneous information in a population
#'
#' @description Get miscelaneous information in a population
#'
#' @param x \code{\link{Pop-class}}
#' @param node character, name of the node to get from the \code{x@misc} slot;
#'   if \code{NULL} the whole \code{x@misc} slot is returned
#'
#' @return The \code{x@misc} slot or its nodes \code{x@misc[[*]][[node]]}
#'
#' @examples
#' founderGenomes <- quickHaplo(nInd = 3, nChr = 1, segSites = 100)
#' SP <- SimParam$new(founderGenomes)
#' basePop <- newPop(founderGenomes)
#'
#' basePop <- setMisc(basePop, node = "info", value = 1)
#' basePop@misc
#' getMisc(x = basePop, node = "info")
#'
#' basePop <- setMisc(basePop, node = "info2", value = c("A", "B", "C"))
#' basePop@misc
#' getMisc(x = basePop, node = "info2")
#' 
#' n <- nInd(basePop)
#' location <- vector(mode = "list", length = n)
#' for (ind in seq_len(n)) {
#'   location[[ind]] <- runif(n = 2, min = 0, max = 100)
#' }
#' location
#' basePop <- setMisc(basePop, node = "location", value = location)
#' basePop@misc
#' getMisc(x = basePop, node = "location")
#' 
#' n <- nInd(basePop)
#' location <- vector(mode = "list", length = n)
#' for (ind in c(1, 3)) {
#'   location[[ind]] <- runif(n = 2, min = 0, max = 100)
#' }
#' location
#' basePop <- setMisc(basePop, node = "location", value = location)
#' basePop@misc
#' getMisc(x = basePop, node = "location")
#' 
#' getMisc(x = basePop)
#' 
#' @export
getMisc <- function(x, node = NULL) {
  if (isPop(x)) {
    if (is.null(node)) {
      ret <- x@misc
    } else {
      nInd <- nInd(x)
      ret <- vector(mode = "list", length = nInd)
      for (ind in seq_len(nInd)) {
        if (!is.null(x@misc[[ind]][[node]])) {
          ret[ind] <- x@misc[[ind]][node]
        }
      }
    }
  } else {
    stop("Argument x must be a Pop class object!")
  }
  return(ret)
}

#' @title Get pedigree
#' 
#' @description 
#' Returns the population's pedigree as stored in the 
#' id, mother and father slots. NULL is returned if the 
#' input population lacks the required.
#' 
#' @param pop a population
#' 
#' @examples 
#' # Create a founder population
#' founderPop = quickHaplo(2,1,2)
#' 
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' # Create a population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' # Get the pedigree
#' getPed(pop)
#' 
#' # Returns NULL when a population lacks a pedigree
#' getPed(founderPop)
#' 
#' @export
getPed = function(pop){
  if(.hasSlot(pop, "id") & 
     .hasSlot(pop, "mother") & 
     .hasSlot(pop, "father")){
    df = data.frame(id = pop@id,
                    mother = pop@mother,
                    father = pop@father)
    return(df)
  }else{
    return(NULL)
  }
}


#' @title Selection intensity
#' 
#' @description 
#' Calculates the standardized selection intensity
#' 
#' @param p the proportion of individuals selected
#' 
#' @examples 
#' selInt(0.1)
#' 
#' @export
selInt = function(p){
  return(dnorm(qnorm(1-p))/p)
}

#' @title Calculate Smith-Hazel weights
#' 
#' @description
#' Calculates weights for Smith-Hazel index given economice weights 
#' and phenotypic and genotypic variance-covariance matrices.
#' 
#' @param econWt vector of economic weights
#' @param varG the genetic variance-covariance matrix
#' @param varP the phenotypic variance-covariance matrix
#' 
#' @return a vector of weight for calculating index values
#' 
#' @examples
#' G = 1.5*diag(2)-0.5
#' E = diag(2)
#' P = G+E
#' wt = c(1,1)
#' smithHazel(wt, G, P)
#' 
#' @export
smithHazel = function(econWt,varG,varP){
  return(solve(varP)%*%varG%*%econWt)
}

#' @title Selection index
#' 
#' @description
#' Calculates values of a selection index given trait values and 
#' weights. This function is intended to be used in combination with 
#' selection functions working on populations such as 
#' \code{\link{selectInd}}.
#' 
#' @param Y a matrix of trait values
#' @param b a vector of weights
#' @param scale should Y be scaled and centered
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' #Model two genetically correlated traits
#' G = 1.5*diag(2)-0.5 #Genetic correlation matrix
#' SP$addTraitA(10, mean=c(0,0), var=c(1,1), corA=G)
#' SP$setVarE(h2=c(0.5,0.5))
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Calculate Smith-Hazel weights
#' econWt = c(1, 1)
#' b = smithHazel(econWt, varG(pop), varP(pop))
#' 
#' #Selection 2 best individuals using Smith-Hazel index
#' #selIndex is used as a trait
#' pop2 = selectInd(pop, nInd=2, trait=selIndex, 
#'                  simParam=SP, b=b)
#' 
#' @export
selIndex = function(Y,b,scale=FALSE){
  if(scale){
    return(scale(Y)%*%b)
  }
  return(Y%*%b)
}

#' @title Edit genome
#' 
#' @description
#' Edits selected loci of selected individuals to a homozygous 
#' state for either the 1 or 0 allele. The gv slot is recalculated to 
#' reflect the any changes due to editing, but other slots remain the same.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param ind a vector of individuals to edit
#' @param chr a vector of chromosomes to edit. 
#' Length must match length of segSites.
#' @param segSites a vector of segregating sites to edit. Length must 
#' match length of chr.
#' @param allele either 0 or 1 for desired allele
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Change individual 1 to homozygous for the 1 allele 
#' #at locus 1, chromosome 1
#' pop2 = editGenome(pop, ind=1, chr=1, segSites=1, 
#'                   allele=1, simParam=SP)
#' 
#' @export
editGenome = function (pop, ind, chr, segSites, allele, simParam = NULL) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  ind = unique(as.integer(ind))
  stopifnot(all(ind %in% (1:pop@nInd)))
  chr = as.integer(chr)
  segSites = as.integer(segSites)
  stopifnot(length(chr) == length(segSites))
  allele = as.integer(allele)
  stopifnot(all(allele == 0L | allele == 1L))
  allele = as.raw(allele)
  if(length(allele) == 1L){
    allele = rep(allele, length(segSites))
  }
  stopifnot(length(allele) == length(segSites))
  for (selChr in unique(chr)) {
    sel = which(chr == selChr)
    for (i in sel) {
      BYTE = (segSites[i] - 1L)%/%8L + 1L
      BIT = (segSites[i] - 1L)%%8L + 1L
      for (selInd in ind) {
        for (j in 1:pop@ploidy) {
          TMP = pop@geno[[selChr]][BYTE, j, selInd]
          TMP = rawToBits(TMP)
          TMP[BIT] = allele[i]
          TMP = packBits(TMP)
          pop@geno[[selChr]][BYTE, j, selInd] = TMP
        }
      }
    }
  }
  PHENO = pop@pheno
  EBV = pop@ebv
  pop = resetPop(pop = pop, simParam = simParam)
  pop@pheno = PHENO
  pop@ebv = EBV
  return(pop)
}

#' @title Edit genome - the top QTL
#' 
#' @description
#' Edits the top QTL (with the largest additive effect) to a homozygous 
#' state for the allele increasing. Only nonfixed QTL are edited The gv slot is
#' recalculated to reflect the any changes due to editing, but other slots remain the same.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param ind a vector of individuals to edit
#' @param nQtl number of QTL to edit
#' @param trait which trait effects should guide selection of the top QTL
#' @param increase should the trait value be increased or decreased
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Change up to 10 loci for individual 1 
#' pop2 = editGenomeTopQtl(pop, ind=1, nQtl=10, simParam=SP)
#'                   
#' @export
editGenomeTopQtl = function(pop, ind, nQtl, trait = 1, increase = TRUE, simParam = NULL) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  ind = unique(as.integer(ind))
  stopifnot(all(ind %in% (1:pop@nInd)))
  nQtl = as.integer(nQtl)
  stopifnot(nQtl > 0 & nQtl <= simParam$traits[[trait]]@nLoci)
  
  findTopQtl = function(pop, ind, nQtl, trait, increase, simParam) {
    # @title Find the top non fixed QTL for use in editGenome()
    # @param pop an object of \code{\link{Pop-class}}
    # @param ind a vector of individuals to edit
    # @param nQtl number of QTL to edit
    # @param trait which trait effects should guide selection of the top QTL
    # @param increase should the trait value be increased or decreased
    # @param simParam an object of \code{\link{SimParam}}
    # @return: a list of four vectors with the:
    #         first  indicating which QTL (of all genome QTL) are the top,
    #         second indicating which segsite (of all segsites within a chromosome) are the top,
    #         third  indicating chromosome of the QTL
    #         fourth indicates which allele we want to fix (edit to)
    QtlGeno = pullQtlGeno(pop=pop[ind],trait=trait,simParam=simParam)
    
    QtlEff = simParam$traits[[trait]]@addEff
    ret = vector(mode = "list", length = 4)
    ret[[1]] = ret[[2]] = ret[[3]] = ret[[4]] = rep(NA, times = nQtl)
    QtlEffRank = order(abs(QtlEff), decreasing = TRUE)
    nQtlInd = 0
    Qtl = 0
    
    while (nQtlInd < nQtl) {
      Qtl = Qtl + 1
      if(Qtl>ncol(QtlGeno)){
        ret[[1]] = ret[[1]][1:nQtlInd]
        ret[[2]] = ret[[2]][1:nQtlInd]
        ret[[3]] = ret[[3]][1:nQtlInd]
        ret[[4]] = ret[[4]][1:nQtlInd]
        nQtl = nQtlInd
        break()
      }
      QtlGenoLoc = QtlGeno[QtlEffRank[Qtl]]
      if (QtlEff[QtlEffRank[Qtl]] > 0) {
        if (QtlGenoLoc < 2) {
          nQtlInd = nQtlInd + 1
          ret[[1]][nQtlInd] = QtlEffRank[Qtl]
          ret[[2]][nQtlInd] = simParam$traits[[trait]]@lociLoc[QtlEffRank[Qtl]]
          if (increase) {
            ret[[4]][nQtlInd] = 1
          } else {
            ret[[4]][nQtlInd] = 0
          }
        }
      } else {
        if (QtlGenoLoc > 0) {
          nQtlInd = nQtlInd + 1
          ret[[1]][nQtlInd] = QtlEffRank[Qtl]
          ret[[2]][nQtlInd] = simParam$traits[[trait]]@lociLoc[QtlEffRank[Qtl]]
          if (increase) {
            ret[[4]][nQtlInd] = 0
          } else {
            ret[[4]][nQtlInd] = 1
          }
        }
      }
    }
    
    # Locate QTL segsite to chromosomes
    tmp = cumsum(simParam$traits[[trait]]@lociPerChr)
    for (Qtl in 1:nQtl) {
      ret[[3]][Qtl] = which(ret[[1]][Qtl] <= tmp)[1]
    }
    ret
  }
  
  for (ind2 in ind) {
    targetQtl = findTopQtl(pop = pop, ind = ind2, nQtl = nQtl, trait = trait,
                           increase = increase, simParam = simParam)
    pop = editGenome(pop = pop,
                     ind = ind2,
                     chr = targetQtl[[3]],
                     segSites = targetQtl[[2]],
                     allele = targetQtl[[4]],
                     simParam = simParam)
  }
  return(pop)
}

#' @title Usefulness criterion
#' 
#' @description Calculates the usefulness criterion
#' 
#' @param pop and object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values (\code{gv}, default), estimated
#' breeding values (\code{ebv}), breeding values (\code{bv}), 
#' or phenotypes (\code{pheno})
#' @param p the proportion of individuals selected
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns a numeric value
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Determine usefulness of population 
#' usefulness(pop, simParam=SP)
#' 
#' #Should be equivalent to GV of best individual
#' max(gv(pop))
#' 
#' @export
usefulness = function(pop,trait=1,use="gv",p=0.1,
                      selectTop=TRUE,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  response = getResponse(pop=pop, trait=trait, use=use,
                         simParam=simParam, ...)
  response = sort(response, decreasing=selectTop)
  response = response[1:ceiling(p*length(response))]
  return(mean(response))
}

#' @title Linear transformation matrix
#' 
#' @description 
#' Creates an m by m linear transformation matrix that 
#' can be applied to n by m uncorrelated deviates 
#' sampled from a standard normal distribution to produce
#' correlated deviates with an arbitrary correlation 
#' of R. If R is not positive semi-definite, the function 
#' returns smoothing and returns a warning (see details).
#' 
#' @param R a correlation matrix
#' 
#' @details 
#' An eigendecomposition is applied to the correlation 
#' matrix and used to test if it is positive semi-definite. 
#' If the matrix is not positive semi-definite, it is not a 
#' valid correlation matrix. In this case, smoothing is 
#' applied to the matrix (as described in the 'cor.smooth' of 
#' the 'psych' library) to obtain a valid correlation matrix. 
#' The resulting deviates will thus not exactly match the 
#' desired correlation, but will hopefully be close if the 
#' input matrix wasn't too far removed from a valid 
#' correlation matrix.
#' 
#' @examples 
#' # Create an 2x2 correlation matrix
#' R = 0.5*diag(2) + 0.5
#' 
#' # Sample 1000 uncorrelated deviates from a 
#' # bivariate standard normal distribution
#' X = matrix(rnorm(2*1000), ncol=2)
#' 
#' # Compute the transformation matrix
#' T = transMat(R)
#' 
#' # Apply the transformation to the deviates
#' Y = X%*%T
#' 
#' # Measure the sample correlation
#' cor(Y)
#' 
#' @export
transMat = function(R){
  # Check if matrix is symmetric 
  # Stop if it is not
  nameR = deparse(substitute(R))
  if(!isSymmetric(R)){
    stop(nameR, " is not a symmetric matrix")
  }
  
  # Check if matrix is positive semi-definite
  # Provide a warning if it is not
  eig = eigen(R, symmetric=TRUE)
  
  if(min(eig$values)<.Machine$double.eps){
    warning("Matrix is not positive semi-definite, see ?transMat for details")
    # Performing correlation matrix smoothing
    eig$values[eig$values<.Machine$double.eps] = 100*.Machine$double.eps
    m = ncol(R)
    totVar = sum(eig$values)
    eig$values = eig$values * m/totVar
    newR = eig$vectors%*%diag(eig$values)%*%t(eig$vectors)
    newR = cov2cor(newR)
    eig = eigen(newR, symmetric=TRUE)
  }
  
  return(
    t(eig$vectors %*% 
        (t(eig$vectors)*sqrt(pmax(eig$values, 0)))
    )
  )
}

#' @title Add Random Mutations
#' 
#' @description
#' Adds random mutations to individuals in a 
#' population. Note that any existing phenotypes 
#' or EBVs are kept. Thus, the user will need to run 
#' \code{\link{setPheno}} and/or \code{\link{setEBV}} 
#' to generate new phenotypes or EBVs that reflect 
#' changes introduced by the new mutations.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param mutRate rate of new mutations
#' @param returnPos should the positions of mutations be returned
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return an object of \code{\link{Pop-class}} if 
#' returnPos=FALSE or a list containing a 
#' \code{\link{Pop-class}} and a data.frame containing the 
#' postions of mutations if returnPos=TRUE
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Introduce mutations
#' pop = mutate(pop, simParam=SP)
#' 
#' @export
mutate = function(pop, mutRate=2.5e-8, returnPos=FALSE, simParam=NULL){

  # Mutation history variable
  IND=NULL; CHR=NULL; HAP=NULL; SITE=NULL
  
  # Number of haplotypes per chromosome
  nHap = pop@nInd*pop@ploidy
  
  # Number of total sites
  s = sum(pop@nLoci)
  
  # Number of mutations per haplotype
  nMut = rbinom(nHap, s, mutRate)
  
  if(any(nMut>0L)){
    for(take in which(nMut>0L)){
      # Determine haplotype and individual
      ind = (take-1L)%/%pop@ploidy + 1L
      hap = (take-1L)%%pop@ploidy + 1L
      
      # Sample mutation sites
      sites = sampleInt(nMut[take], s) + 1L
      
      # Resolve all mutations
      chr = 1L
      for(i in sites){
        # Find chromosome
        repeat{
          if(i > sum(pop@nLoci[1L:chr])){
            chr = chr + 1L
          }else{
            break
          }
        }
        
        # Find site
        if(chr>1L){
          site = i - sum(pop@nLoci[1L:(chr-1L)]) 
        }else{
          site = i
        }
        
        # Create mutation
        BYTE = (site-1L)%/%8L + 1L
        BIT = (site-1L)%%8L + 1L
        TMP = pop@geno[[chr]][BYTE,hap,ind]
        TMP = rawToBits(TMP)
        TMP[BIT] = ifelse(TMP[BIT], as.raw(0L), as.raw(1L))
        TMP = packBits(TMP)
        pop@geno[[chr]][BYTE,hap,ind] = TMP
        
        # Record results
        if(returnPos){
          IND = c(IND, ind)
          CHR = c(CHR, chr)
          HAP = c(HAP, hap)
          SITE = c(SITE, site)
        }
      }
    }
    
    # Reset population
    PHENO = pop@pheno
    EBV = pop@ebv
    pop = resetPop(pop=pop, simParam=simParam)
    pop@pheno = PHENO
    pop@ebv = EBV
  }
  
  # Return results
  if(returnPos){
    return(list(pop,data.frame(individual=IND,chromosome=CHR,haplotype=HAP,site=SITE)))
  }else{
    return(pop)
  }
}

#' @title Lose individuals at random
#' 
#' @description
#' Samples individuals at random to remove from the population. 
#' The user supplies a probability for the individuals to be 
#' removed from the population.
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param p the expected proportion of individuals that will
#' be lost to attrition.
#'
#' @return an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=100, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Lose an expected 5% of individuals
#' pop = attrition(pop, p=0.05)
#' 
#' @export
attrition = function(pop, p){
  take = as.logical(rbinom(pop@nInd, size=1, prob=1-p))
  return(pop[take])
}

# Sample deviates from a standard normal distribution
# n is the number of deviates
# u is a deviate from a uniform distribution [0,1]
# Seed is generated from u
rnormWithSeed = function(n, u){
  glbEnv = globalenv()
  origSeed = glbEnv$.Random.seed
  on.exit({
    if(is.null(origSeed)){
      rm(list =".Random.seed", envir=glbEnv)
    }else{
      assign(".Random.seed", value=origSeed, 
             envir=glbEnv)
    }
  })
  set.seed(as.integer((u-0.5)*2*2147483647))
  rnorm(n)
}
