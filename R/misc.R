#' @title Selection intensity
#' 
#' @description 
#' Calculates the standardized selection intensity
#' 
#' @param p the proportion of individuals selected
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
    pop@geno[[selChr]][selSegSites,,ind] = selAllele
  }
  pop@gxe = vector("list",simParam$nTraits)
  pop@gv = matrix(NA_real_,nrow=pop@nInd,
                  ncol=simParam$nTraits)
  if(simParam$nTraits>=1){
    for(i in 1:simParam$nTraits){
      tmp = getGv(simParam$traits[[i]],pop)
      pop@gv[,i] = tmp[[1]]
      if(length(tmp)>1){
        pop@gxe[[i]] = tmp[[2]]
      }
    }
  }
  validObject(pop)
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
    ret = vector(mode = "list", length = 2)
    ret[[1]] = ret[[2]] = ret[[3]] = ret[[4]] = rep(NA, times = nQtl)
    QtlEffRank = order(abs(QtlEff), decreasing = TRUE)
    nQtlInd = 0
    Qtl = 0
    
    while (nQtlInd < nQtl) {
      Qtl = Qtl + 1
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
  validObject(pop)
  return(pop)
}

#' @title Correlated vector
#' 
#' @description
#' Creates a correlated vector by adding random error. 
#'
#' @param x a numeric vector
#' @param rho desired correlation. Must be greater than 
#' 0 and less than or equal to 1.
#' @param shrink should the output vector be shrunken 
#' to have similar variance as the input vector
#'
#' @return a numeric vector
#'
#' @export
corVec = function(x,rho,shrink=FALSE){
  stopifnot(rho>0, rho<=1)
  x = as.vector(x)
  varX = var(x)
  varE = varX/(rho^2)-varX
  y = x+rnorm(length(x),sd=sqrt(varE))
  if(shrink){
    meanX = mean(x)
    y = (y-meanX)/sqrt(varX+varE)
    y = y*sqrt(varX)+meanX
  }
  return(y)
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
#' @export
usefulness = function(pop,trait=1,use="gv",p=0.1,
                      selectTop=TRUE,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  response = getResponse(pop=pop,trait=trait,use=use,
                         simParam=simParam,...)
  response = sort(response,decreasing=selectTop)
  response = response[1:ceiling(p*length(response))]
  return(mean(response))
}

#' @title Variance to correlation
#' 
#' @description
#' Converts a variance-covariance matrix to a 
#' correlation matrix.
#'
#' @param var a variance-covariance matrix 
#'
#' @return a numeric matrix
#'
#' @export
var2cor = function(var){
  tmp = diag(1/sqrt(diag(var)))
  return(tmp%*%var%*%tmp)
}

#' @title Converts a Pop-class to PLINK data
#' 
#' @description
#' Converts a Pop-class to PLINK FAM and MAP data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param trait an integer. Which phenotype trait should be used.
#' @param snpChip an integer. Which SNP array should be used.
#' @param simParam an object of \code{\link{SimParam}}
#' @param chromLength an integer. The size of chromosomes in base
#' pairs; assuming all chromosomes are of the same size.
#' 
#' @return a list with the fam, ped, and map data.frames
#'
#' @export
pop2Plink = function(pop, trait = 1L, snpChip = 1L, simParam = NULL, chromLength = 10L^8) {
  # ---- Setup ----
  
  if (is.null(simParam)) {
    simParam = get(x = "SP", envir = .GlobalEnv)
  }
  Ret = vector(mode = "list", length = 3L)
  names(Ret) = c("fam", "ped", "map")
  
  # ---- Fam ----
  
  Ret$fam = data.frame(Family = rep(x = 1L, times = pop@nInd),
                       IId  = pop@id,
                       FId  = pop@father,
                       MId  = pop@mother,
                       Sex  = 0L,
                       Phen = 0)
  if (!any(pop@gender == "")) {
    Ret$fam$Sex = (pop@gender == "F") + 1L
  }
  if (!anyNA(pop@pheno[, trait])) {
    Ret$fam$Phen = pop@pheno[, trait]
  }
  
  # ---- Ped ----
  
  nLoc = simParam$snpChips[[snpChip]]@nLoci
  Ret$ped = as.data.frame(matrix(data = 0L, nrow = pop@nInd, ncol = 6L + 2L * nLoc))
  Ret$ped[, 1L:6L] = Ret$fam
  Tmp = pullSnpHaplo(pop = pop, snpChip = snpChip)
  Sel1 = seq(from = 1L, to = 2L * pop@nInd,   by = 2L)
  Sel2 = seq(from = 2L, to = 2L * pop@nInd, by = 2L)
  k = 7L
  for (Loc in 1L:nLoc) {
    # Loc = 1L
    Ret$ped[, k] = Tmp[Sel1, Loc] + 1L
    k = k + 1L
    Ret$ped[, k] = Tmp[Sel2, Loc] + 1L
    k = k + 1L
  }
  
  # ---- Map ----
  
  # Assuming equal number of markers per chromosome!
  Ret$map = data.frame(Chr = rep(x = 1L:simParam$nChr, each = simParam$snpChips[[snpChip]]@lociPerChr[1L]),
                       Loc = colnames(Tmp),
                       PosGenetic = 0L,
                       Pos = simParam$snpChips[[snpChip]]@lociLoc)
  for (Chr in 1L:simParam$nChr) {
    Sel = Ret$map$Chr == Chr
    Ret$map$PosGenetic[Sel] = simParam$genMap[[Chr]][Ret$map$Pos[Sel]]
  }
  Ret$map$Pos = round(Ret$map$PosGenetic * chromLength)
  
  # ---- Return ----
  
  Ret
}

#' @title Writes a Pop-class to PLINK data
#' 
#' @description
#' Writes a Pop-class to PLINK FAM and MAP data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param baseName a character. Output file basename.
#' @param trait an integer. Which phenotype trait should be used.
#' @param snpChip an integer. Which SNP array should be used.
#' @param simParam an object of \code{\link{SimParam}}
#' @param chromLength an integer. The size of chromosomes in base
#' pairs; assuming all chromosomes are of the same size.
#' 
#' @return a list with the fam, ped, and map data.frames
#'
#' @export
writePlink = function(pop, baseName, trait = 1L, snpChip = 1L,
                      simParam = NULL, chromLength = 10L^8) {
  x = pop2Plink(pop = pop, trait = trait, snpChip = snpChip,
                simParam = simParam, chromLength = chromLength)
  write.table(x = x$ped, file = paste0(baseName, ".ped"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(x = x$map, file = paste0(baseName, ".map"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
}