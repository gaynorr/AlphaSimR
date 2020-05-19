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
#' @param chr a vector of chromosomes to edit. Length must match 
#' length of segSites.
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
editGenome = function(pop,ind,chr,segSites,allele,
                      simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  ind = unique(as.integer(ind))
  stopifnot(all(ind%in%(1:pop@nInd)))
  chr = as.integer(chr)
  segSites = as.integer(segSites)
  stopifnot(length(chr)==length(segSites))
  allele = as.integer(allele)
  stopifnot(allele==0L | allele==1L)
  allele = as.raw(allele)
  for(selChr in unique(chr)){
    sel = chr==selChr
    selSegSites = segSites[sel]
    selAllele = allele[sel]
    for(selSite in selSegSites){
      BYTE = (selSite-1L)%/%8L + 1L
      BIT = (selSite-1L)%%8L + 1L
      for(selInd in ind){
        for(i in 1:pop@ploidy){
          TMP = pop@geno[[selChr]][BYTE,i,selInd]
          TMP = rawToBits(TMP)
          TMP[BIT] = selAllele
          TMP = packBits(TMP)
          pop@geno[[selChr]][BYTE,i,selInd] = TMP
        }
      }
    }
  }
  # Reset population
  PHENO = pop@pheno
  EBV = pop@ebv
  pop = resetPop(pop=pop, simParam=simParam)
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

#' @title Writes a Pop-class as PLINK files
#' 
#' @description
#' Writes a Pop-class as PLINK PED and MAP files
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param baseName a character. Basename of PED and MAP files.
#' @param trait an integer. Which phenotype trait should be used.
#' @param snpChip an integer. Which SNP array should be used.
#' @param simParam an object of \code{\link{SimParam}}
#' @param chromLength an integer. The size of chromosomes in base
#' pairs; assuming all chromosomes are of the same size.
#'
#' @examples 
#' \dontrun{
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$setGender(gender = "yes_rand")
#' SP$addTraitA(nQtlPerChr = 10)
#' SP$addSnpChip(nSnpPerChr = 5)
#' 
#' #Create population
#' pop = newPop(rawPop = founderPop)
#' pop = setPheno(pop, varE = SP$varA)
#' writePlink(pop, baseName="test")
#' 
#' #Test
#' test = read.table(file = "test.ped")
#' #...gender
#' if (!identical(x = c("M", "F")[test[[5]]], y = pop@gender)) { stop() }
#' #...pheno (issues with rounding)
#' # if (!identical(x = test[[6]], y = pop@pheno[, 1])) { stop() }
#' #...genotypes
#' x = test[, -(1:6)]  - 1
#' x[, 1] = x[, 1] + x[, 2]
#' x[, 2] = x[, 3] + x[, 4]
#' x[, 3] = x[, 5] + x[, 6]
#' x[, 4] = x[, 7] + x[, 8]
#' x[, 5] = x[, 9] + x[, 10]
#' y = pullSnpGeno(pop)
#' if (sum(x[, 1:5] - y) != 0) { stop() }
#' }
#' @export
writePlink = function(pop, baseName, trait = 1L, snpChip = 1L, simParam = NULL,
                      chromLength = 10L^8) {
  if (is.null(simParam)) {
    simParam = get(x = "SP", envir = .GlobalEnv)
  }
  if (pop@ploidy != 2L) {
    stop(paste0("writePlink() will write ", pop@ploidy, " alleles for each locus!"))
  }  
  
  # ---- Map ----
  
  # This assumes equal number of markers per chromosome!
  map = data.frame(chr = rep(x = 1L:simParam$nChr,
                             each = simParam$snpChips[[snpChip]]@lociPerChr[1L]),
                   loc = paste0("SNP_", 1L:simParam$snpChips[[snpChip]]@nLoci),
                   posGenetic = 0L,
                   pos = simParam$snpChips[[snpChip]]@lociLoc)
  for (chr in 1L:simParam$nChr) {
    # chr = 1
    sel = map$chr == chr
    map$posGenetic[sel] = simParam$genMap[[chr]][map$pos[sel]]
  }
  map$pos = round(map$posGenetic * chromLength)
  write.table(x = map, file = paste0(baseName, ".map"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # ---- Ped ----
  
  # First the FAM format, which covers the first 6 columns of the PED format
  fam = data.frame(family = rep(x = 1L, times = pop@nInd),
                   id     = as.integer(pop@id),
                   father = as.integer(pop@father),
                   mother = as.integer(pop@mother),
                   gender = 0L,
                   pheno  = 0)
  if (!any(pop@gender == "")) {
    fam$gender = (pop@gender == "F") + 1L
  }
  if (!anyNA(pop@pheno[, trait])) {
    fam$pheno = pop@pheno[, trait]
  }
  
  # Select loci on the SNP array
  tmp = selectLoci(chr          = 1L:simParam$nChr,
                   inLociPerChr = simParam$snpChips[[snpChip]]@lociPerChr,
                   inLociLoc    = simParam$snpChips[[snpChip]]@lociLoc)
  # Add loci alleles to fam and write to file directly from C++
  writePlinkPed(fam    = fam,
                haplo  = getHaplo(geno       = pop@geno,
                                  lociPerChr = tmp$lociPerChr,
                                  lociLoc    = tmp$lociLoc,
                                  nThreads   = simParam$nThreads),
                nInd   = nrow(fam),
                ploidy = pop@ploidy,
                nLoc   = sum(tmp$lociPerChr),
                file   = paste0(baseName, ".ped"))
}

#Create rotation matrix for sampling random deviates
#Uses SVD method for stability
rotMat = function(X){
  ans = svd(X)
  u = t(ans$u)*sqrt(pmax(ans$d,0))
  return(t(ans$v%*%u))
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
