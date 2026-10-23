#' Find loci on specific chromosomes
#' 
#' This function alters the lociPerChr and lociLoc vectors to reflect 
#' only loci on specific chromosomes.
#'
#' @param chr chromosomes to select
#' @param inLociPerChr original lociPerChr vector
#' @param inLociLoc original lociLoc vector
#' 
#' @return a list with lociPerChr and lociLoc
#'
#' @keywords internal
selectLoci = function(chr, inLociPerChr, inLociLoc){
  if(is.null(chr)){
    return(list(lociPerChr=inLociPerChr,
                lociLoc=inLociLoc))
  }
  nChr = length(inLociPerChr)
  if(!any(chr%in%(1:nChr))){
    stop("chr does not name any chromosome in the population")
  }
  if(max(chr)>nChr){
    stop("chr exceeds the number of chromosomes in the population")
  }
  outLociPerChr = numeric(nChr)
  outLociPerChr[chr] = inLociPerChr[chr]
  outLociLoc = numeric(sum(outLociPerChr))
  inStart = outStart = inEnd = outEnd = 0L
  for(i in seq_len(nChr)){
    inStart = inStart + 1L
    inEnd = inEnd + inLociPerChr[i]
    if(outLociPerChr[i]>0){
      outStart = outStart + 1L
      outEnd = outEnd + inLociPerChr[i]
      outLociLoc[outStart:outEnd] = inLociLoc[inStart:inEnd]
      outStart = outEnd
    }
    inStart = inEnd
  }
  return(list(lociPerChr=outLociPerChr,
              lociLoc=outLociLoc))
}

#' Check a haplotype selector
#'
#' Checks the haplo argument of the pull*Haplo functions. A single
#' haplotype is read through a raw pointer in C++, so a value outside
#' 1 to ploidy would read memory outside the genotypes.
#'
#' @param haplo either "all" or a single haplotype number
#' @param ploidy ploidy level of the population
#'
#' @return "all", or haplo as an integer
#'
#' @keywords internal
checkHaplo = function(haplo, ploidy){
  if(identical(haplo, "all")){
    return(haplo)
  }
  # The range is checked before the whole number test so that a very large
  # value never reaches as.integer, which would turn it into NA
  if((length(haplo)!=1L) || !is.numeric(haplo) || is.na(haplo) ||
     (haplo<1) || (haplo>ploidy) || (haplo!=round(haplo))){
    stop(paste0("haplo must be \"all\" or a whole number from 1 to ",
                ploidy))
  }
  return(haplo)
}

#' Retrieves marker names from genMap
#'
#' @param lociPerChr number of loci per chromosome
#' @param lociLoc index of loci within each chromosome's segregating sites
#' @param genMap internal AlphaSimR genetic map with names
#'
#' @keywords internal
getLociNames = function(lociPerChr, lociLoc, genMap){
  lociNames = character(length(lociLoc))
  start = end = 0L
  for(chr in seq_len(length(lociPerChr))){
    if(lociPerChr[chr]>0){
      start = end + 1L
      end = end + lociPerChr[chr]
      take = lociLoc[start:end]
      lociNames[start:end] = names(genMap[[chr]])[take]
    }
  }
  return(lociNames)
}

#' Finds positions of loci by marker name
#' 
#' @description Used to generate lociPerChr and lociLoc 
#' objects for a set of markers. These objects can be passed to 
#' other functions for pulling genotypes or haplotypes.
#' 
#' @param markers a vector of marker names
#' @param genMap a genetic map in AlphaSimR's internal 
#' genetic map format
#'
#' @return A list containing lociPerChr and lociLoc 
#' 
#' @keywords internal
mapLoci = function(markers, genMap){
  # Check that the markers are present on the map
  genMapMarkerNames = unlist(lapply(genMap, names))
  if(!all(markers%in%genMapMarkerNames)){
    stop("markers contains names that are not in the genetic map")
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
  
  return(list(lociPerChr=lociPerChr,
              lociLoc=lociLoc))
}

#' @title Get genetic map
#' 
#' @description Retrieves the genetic map for all loci.
#' 
#' @param object where to retrieve the genetic map. 
#' Can be an object of \code{\link{SimParam}} or 
#' \code{\link{MapPop-class}}. If NULL, the function will 
#' look for a SimParam object called "SP" in your 
#' global environment.
#' @param sex determines which sex specific map 
#' is returned. Options are "A" for average map, "F" 
#' for female map, and "M" for male map. All options are 
#' equivalent if not using sex specific maps or using 
#' pulling from a MapPop.
#'
#' @return Returns a data.frame with:
#' \describe{
#'   \item{id}{Unique identifier for locus}
#'   \item{chr}{Chromosome containing the locus}
#'   \item{pos}{Genetic map position}
#' }
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' getGenMap(founderPop)
#' 
#' @export
getGenMap = function(object=NULL, sex="A"){
  
  if(is.null(object)){
    object = get("SP",envir=.GlobalEnv)
  }
  
  if(is(object, "SimParam")){
    # Extract the sex specific map
    if(sex=="A"){
      genMap = object$genMap
    }else if(sex=="F"){
      genMap = object$femaleMap
    }else if(sex=="M"){
      genMap = object$maleMap
    }else{
      stop(paste("sex =",sex,"is not a valid option"))
    }
  }else if(is(object, "MapPop")){
    # No sex specific maps
    genMap = object@genMap
  }else{
    stop("object is class ",
         class(object), 
         ", which is not a valid source for the genetic map")
  }
  
  # Loci names
  id = unlist(unname(lapply(genMap, names)))
  
  # Chromosome names
  nLoci = sapply(genMap, length)
  chr = rep(names(genMap), nLoci)
  
  # Map positions
  pos = unname(do.call("c", genMap))
  
  # Output data.frame
  return(data.frame(id=id, chr=chr, pos=pos))
}

#' @title Get SNP genetic map
#' 
#' @description Retrieves the genetic map for a 
#' given SNP chip.
#' 
#' @param snpChip an integer. Indicates which SNP
#' chip's map to retrieve.
#' @param sex determines which sex specific map 
#' is returned. Options are "A" for average map, "F" 
#' for female map, and "M" for male map. All options are 
#' equivalent if not using sex specific maps.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#'
#' @return Returns a data.frame with:
#' \describe{
#'   \item{id}{Unique identifier for the SNP}
#'   \item{chr}{Chromosome containing the SNP}
#'   \item{site}{Segregating site on the chromosome}
#'   \item{pos}{Genetic map position}
#' }
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addSnpChip(5)
#' 
#' #Pull SNP map
#' getSnpMap(snpChip=1, simParam=SP)
#' 
#' @export
getSnpMap = function(snpChip=1, sex="A", simParam=NULL){
  
  simParam = getSimParam(simParam)
  
  #Extract genetic map and SNP positions
  if(sex=="A"){
    genMap = simParam$genMap
  }else if(sex=="F"){
    genMap = simParam$femaleMap
  }else if(sex=="M"){
    genMap = simParam$maleMap
  }else{
    stop(paste("sex =",sex,"is not a valid option"))
  }
  
  if(is.character(snpChip)){
    # Suspect snpChip is a name
    chipNames = simParam$snpChipNames
    take = match(snpChip, chipNames)
    if(is.na(take)){
      stop("'",snpChip,"' did not match any SNP chip names")
    }
    snpChip = take
  }
  
  snp = simParam$snpChips[[snpChip]] #SNP positions
  
  #Create a list of SNP positions on the genetic map
  #Each list element corresponds to a chromosome
  snpMap = lapply(1:simParam$nChr, function(x){
    if(snp@lociPerChr[x]==0){
      #No SNPs on chromosome
      return(NULL)
    }else{
      if(x==1){
        #First chromosome, start at position 1
        take = 1:snp@lociPerChr[x]
      }else{
        #All other chromosomes
        take = (sum(snp@lociPerChr[1:(x-1)])+1):sum(snp@lociPerChr[1:x])
      }
      return(genMap[[x]][snp@lociLoc[take]])
    }
  })
  
  #Create a data.frame with SNP positions on genetic map
  output = data.frame(id=getLociNames(snp@lociPerChr, snp@lociLoc, genMap),
                      chr=rep(names(genMap),snp@lociPerChr),
                      site=snp@lociLoc,
                      pos=do.call("c",snpMap))
  return(output)
}

#' @title Get QTL genetic map
#' 
#' @description Retrieves the genetic map for the 
#' QTL of a given trait.
#' 
#' @param trait an integer or a trait name indicating which trait's QTL
#' to return
#' @param sex determines which sex specific map 
#' is returned. Options are "A" for average map, "F" 
#' for female map, and "M" for male map. All options are 
#' equivalent if not using sex specific maps.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#'
#' @return Returns a data.frame with:
#' \describe{
#'   \item{id}{Unique identifier for the QTL}
#'   \item{chr}{Chromosome containing the QTL}
#'   \item{site}{Segregating site on the chromosome}
#'   \item{pos}{Genetic map position}
#' }
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(5)
#' 
#' #Pull SNP map
#' getQtlMap(trait=1, simParam=SP)
#' 
#' @export
getQtlMap = function(trait=1, sex="A", simParam=NULL){
  
  simParam = getSimParam(simParam)
  
  if(is.character(trait)){
    # Suspect trait is a name
    traitNames = simParam$traitNames
    take = match(trait, traitNames)
    if(is.na(take)){
      stop("'",trait,"' did not match any trait names")
    }
    trait = take
  }
  
  #Extract genetic map and SNP positions
  if(sex=="A"){
    genMap = simParam$genMap
  }else if(sex=="F"){
    genMap = simParam$femaleMap
  }else if(sex=="M"){
    genMap = simParam$maleMap
  }else{
    stop(paste("sex =",sex,"is not a valid option"))
  }
  qtl = simParam$traits[[trait]] #QTL positions
  
  #Create a list of QTL positions on the genetic map
  #Each list element corresponds to a chromosome
  qtlMap = lapply(1:simParam$nChr, function(x){
    if(qtl@lociPerChr[x]==0){
      #No QTL on chromosome
      return(NULL)
    }else{
      if(x==1){
        #First chromosome, start at position 1
        take = 1:qtl@lociPerChr[x]
      }else{
        #All other chromosomes
        take = (sum(qtl@lociPerChr[1:(x-1)])+1):sum(qtl@lociPerChr[1:x])
      }
      return(genMap[[x]][qtl@lociLoc[take]])
    }
  })
  
  #Create a data.frame with QTL positions on genetic map
  output = data.frame(id=getLociNames(qtl@lociPerChr, qtl@lociLoc, genMap),
                      chr=rep(names(genMap),qtl@lociPerChr),
                      site=qtl@lociLoc,
                      pos=do.call("c",qtlMap))
  return(output)
}

#' @name pullGeno
#' @rdname pullGeno
#'
#' @title Pull genotypes
#'
#' @description
#' Retrieves genotype data from:    \cr
#' - all segregating sites (\code{\link{pullSegSiteGeno}})    \cr
#' - quantitative trait loci (QTL) for a trait (\code{\link{pullQtlGeno}})    \cr
#' - SNP chip markers (\code{\link{pullSnpGeno}})    \cr
#' - a user-defined set of markers by name (\code{\link{pullMarkerGeno}})    \cr
#'
#' Supports \code{\link{Pop-class}}, \code{\link{HybridPop-class}}, and
#' \code{\link{MultiPop-class}} inputs. For \code{MultiPop} objects, output can
#' optionally be simplified to a chosen nesting level.
#'
#' \code{\link{MapPop-class}} objects are supported by
#' \code{\link{pullSegSiteGeno}} and \code{\link{pullMarkerGeno}}.
#'
#' @param pop An object of \code{\link{Pop-class}},
#'   \code{\link{HybridPop-class}}, or \code{\link{MultiPop-class}}.
#'   \code{\link{MapPop-class}} objects are also supported by
#'   \code{\link{pullSegSiteGeno}} and \code{\link{pullMarkerGeno}}.
#' @param chr Integer vector of chromosomes to retrieve. If \code{NULL},
#'   all chromosomes are retrieved.
#' @param trait An integer or a trait name. Indicates which trait's
#'   QTL genotypes to retrieve.
#' @param snpChip SNP chip whose genotypes are retrieved. Either an integer
#'   index or a SNP chip name.
#' @param markers Character vector of marker names to retrieve.
#' @param asRaw Logical; if \code{TRUE}, returns genotypes in raw (byte) format.
#' @param simParam An object of class \code{\link{SimParam}} for non-
#'   \code{\link{MapPop-class}} inputs. If \code{NULL}, the object named
#'   \code{SP} in the global environment is used. Ignored when \code{pop} is a
#'   \code{\link{MapPop-class}}.
#' @param nThreads Number of threads to use if OpenMP is available. If
#'   \code{NULL}, uses \code{\link{getNumThreads}} for
#'   \code{\link{MapPop-class}} objects and \code{simParam$nThreads} otherwise.
#' @param simplify Logical. Used only when \code{pop} is a
#'   \code{\link{MultiPop-class}}. If \code{TRUE}, flattens \code{pop} to
#'   \code{level} with \code{\link{flattenMultiPop}}, combines terminal results
#'   with \code{\link{rbind}}, and attaches a \code{"source"} attribute
#'   indicating row origins.
#' @param level Integer scalar \eqn{\ge 1}. Used only when \code{pop} is a
#'   \code{\link{MultiPop-class}} and \code{simplify=TRUE}. Number of
#'   \code{MultiPop} levels to preserve. Passed to
#'   \code{\link{flattenMultiPop}}. Ignored when \code{simplify=FALSE}.
#'
#' @details
#' With \code{simplify=FALSE}, \code{MultiPop} structure is preserved and output
#' follows the same nesting as \code{pop}.
#'
#' With \code{simplify=TRUE}, output is simplified to the requested
#' \code{level}. If \code{level} exceeds the nesting depth of \code{pop}, the
#' original structure is returned unchanged.
#'
#' Simplified output includes a \code{"source"} attribute: a data frame with
#' columns \code{level1}, \code{level2}, etc., giving the source population
#' name (or index) for each row at each retained level.
#' 
#' @seealso \code{\link{pullHaplo}}, \code{\link{pullAf}}
#'
#' @return
#' For a non-\code{MultiPop} input, a genotype matrix of allele dosages, 
#' (\emph{\eqn{0=} homozygous ancestral genotype, rows = individuals, columns =
#' loci}) with row names set to individual IDs when available.
#'
#' For a \code{\link{MultiPop-class}} input:    \cr
#' - if \code{simplify=FALSE}, a nested list matching the \code{MultiPop}
#'   structure, with one genotype matrix per terminal population.    \cr
#' - if \code{simplify=TRUE}, a simplified object where terminal genotype
#'   matrices are combined (within each retained branch) using
#'   \code{\link{rbind}}. A \code{"source"} attribute records row origins.
#'   When \code{level=1}, the result is a single genotype matrix.
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#'
#' # Retrieve genotypes from all segregating loci in founder haplotypes
#' pullSegSiteGeno(founderPop)
#'
#' # Retrieve genotypes from a specific set of markers in founder haplotypes
#' pullMarkerGeno(founderPop, markers=c("1_1","1_2"))
#'
#' # Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads=1L}
#'
#' # Add an additive trait with 10 QTL
#' SP$addTraitA(10, name="Trait1")
#'
#' # Add a SNP chip
#' SP$addSnpChip(nSnpPerChr=5, name="5SNPchip")
#'
#' # Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' # Retrieve genotypes from all segregating loci
#' pullSegSiteGeno(pop, simParam=SP)
#'
#' # Retrieve genotypes from QTL for a specific trait
#' pullQtlGeno(pop, trait="Trait1", simParam=SP)
#'
#' # Retrieve genotypes from a SNP chip
#' pullSnpGeno(pop, snpChip="5SNPchip", simParam=SP)
#'
#' # Retrieve genotypes from a specific set of markers
#' pullMarkerGeno(pop, markers=c("1_1","1_2"), simParam=SP)
#'
#' # Create a nested MultiPop
#' pop2 = randCross(pop, nCrosses=3, nProgeny=4)
#' mp1 = splitPop(
#'   pop2,
#'   by = list(
#'     function(x) rep(LETTERS[1:2], length.out=length(x)),
#'     function(x) paste(x@mother, x@father, sep="_")
#'   )
#' )
#'
#' # Retrieve QTL genotypes at different simplification levels
#' pullQtlGeno(mp1, simplify=FALSE)
#' pullQtlGeno(mp1, simplify=TRUE, level=2)
#' pullQtlGeno(mp1, simplify=TRUE, level=1)
NULL

#' @rdname pullGeno
#' @export
pullSegSiteGeno = function(
  pop,
  chr = NULL,
  asRaw = FALSE,
  simParam = NULL,
  nThreads = NULL,
  simplify = FALSE,
  level = 1L
) {
  calcPopValue(
    x = pop,
    FUN = function(
      pop,
      chr = chr,
      asRaw = asRaw,
      simParam = simParam,
      nThreads = nThreads
    ) {
      if (is(pop, "MapPop")) {
        allLoci = unlist(c(sapply(pop@nLoci, function(x) 1:x)))
        lociTot = pop@nLoci
        if (is.null(nThreads)) {
          nThreads = getNumThreads()
        } else {
          nThreads = as.integer(nThreads)
        }
        map = pop@genMap
      } else {
        simParam = getSimParam(simParam)

        allLoci = unlist(c(sapply(simParam$segSites, function(x) 1:x)))
        lociTot = simParam$segSites
        if (is.null(nThreads)) {
          nThreads = simParam$nThreads
        } else {
          nThreads = as.integer(nThreads)
        }
        map = simParam$genMap
      }

      tmp = selectLoci(chr, lociTot, allLoci)

      output = getGeno(pop@geno, tmp$lociPerChr, tmp$lociLoc, nThreads)

      if (!asRaw) {
        output = convToImat(output)
      }

      if (is(pop, "Pop") | is(pop, "NamedMapPop")) {
        rownames(output) = pop@id
      } else {
        rownames(output) = as.character(1:pop@nInd)
      }

      colnames(output) = getLociNames(tmp$lociPerChr, tmp$lociLoc, map)

      return(output)
    },
    simplify = simplify,
    level = level,
    chr = chr,
    asRaw = asRaw,
    simParam = simParam,
    nThreads = nThreads
  )
}

#' @rdname pullGeno
#' @export
pullQtlGeno = function(
  pop,
  trait = 1,
  chr = NULL,
  asRaw = FALSE,
  simParam = NULL,
  nThreads = NULL,
  simplify = FALSE,
  level = 1L
) {
  calcPopValue(
    x = pop,
    FUN = function(
      pop,
      trait = trait,
      chr = chr,
      asRaw = asRaw,
      simParam = simParam,
      nThreads = nThreads,
      simplify = simplify,
      level = level
    ) {
      simParam = getSimParam(simParam)
      if (is.null(nThreads)) {
        nThreads = simParam$nThreads
      } else {
        nThreads = as.integer(nThreads)
      }

      if (is.character(trait)) {
        # Suspect trait is a name
        traitNames = simParam$traitNames
        take = match(trait, traitNames)
        if (is.na(take)) {
          stop("'", trait, "' did not match any trait names")
        }
        trait = take
      }

      tmp = selectLoci(
        chr,
        simParam$traits[[trait]]@lociPerChr,
        simParam$traits[[trait]]@lociLoc
      )

      output = getGeno(pop@geno, tmp$lociPerChr, tmp$lociLoc, nThreads)

      if (!asRaw) {
        output = convToImat(output)
      }

      if (is(pop, "Pop")) {
        rownames(output) = pop@id
      } else {
        rownames(output) = as.character(1:pop@nInd)
      }

      colnames(output) = getLociNames(
        tmp$lociPerChr,
        tmp$lociLoc,
        simParam$genMap
      )

      return(output)
    },
    simplify = simplify,
    level = level,
    trait = trait,
    chr = chr,
    asRaw = asRaw,
    simParam = simParam,
    nThreads = nThreads
  )
}

#' @rdname pullGeno
#' @export
pullSnpGeno = function(
  pop,
  snpChip = 1,
  chr = NULL,
  asRaw = FALSE,
  simParam = NULL,
  nThreads = NULL,
  simplify = FALSE,
  level = 1L
) {
  calcPopValue(
    x = pop,
    FUN = function(
      pop,
      snpChip = snpChip,
      chr = chr,
      asRaw = asRaw,
      simParam = simParam,
      nThreads = nThreads
    ) {
      simParam = getSimParam(simParam)
      if (is.null(nThreads)) {
        nThreads = simParam$nThreads
      } else {
        nThreads = as.integer(nThreads)
      }

      if (is.character(snpChip)) {
        # Suspect snpChip is a name
        chipNames = simParam$snpChipNames
        take = match(snpChip, chipNames)
        if (is.na(take)) {
          stop("'", snpChip, "' did not match any SNP chip names")
        }
        snpChip = take
      }

      tmp = selectLoci(
        chr,
        simParam$snpChips[[snpChip]]@lociPerChr,
        simParam$snpChips[[snpChip]]@lociLoc
      )

      output = getGeno(pop@geno, tmp$lociPerChr, tmp$lociLoc, nThreads)

      if (!asRaw) {
        output = convToImat(output)
      }

      if (is(pop, "Pop")) {
        rownames(output) = pop@id
      } else {
        rownames(output) = as.character(1:pop@nInd)
      }

      colnames(output) = getLociNames(
        tmp$lociPerChr,
        tmp$lociLoc,
        simParam$genMap
      )

      return(output)
    },
    snpChip = snpChip,
    chr = chr,
    asRaw = asRaw,
    simParam = simParam,
    nThreads = nThreads,
    simplify = simplify,
    level = level
  )
}

#' @rdname pullGeno
#' @export
pullMarkerGeno = function(
  pop,
  markers,
  asRaw = FALSE,
  simParam = NULL,
  nThreads = NULL,
  simplify = FALSE,
  level = 1L
) {
  calcPopValue(
    x = pop,
    FUN = function(
      pop,
      markers,
      asRaw = asRaw,
      simParam = simParam,
      nThreads = nThreads
    ) {
      # Get genetic map and nThreads
      if (is(pop, "MapPop")) {
        if (is.null(nThreads)) {
          nThreads = getNumThreads()
        } else {
          nThreads = as.integer(nThreads)
        }
        genMap = pop@genMap
      } else {
        simParam = getSimParam(simParam)
        if (is.null(nThreads)) {
          nThreads = simParam$nThreads
        } else {
          nThreads = as.integer(nThreads)
        }
        genMap = simParam$genMap
      }

      # Map markers to genetic map
      lociMap = mapLoci(markers, genMap)

      # Get genotypes
      output = getGeno(pop@geno, lociMap$lociPerChr, lociMap$lociLoc, nThreads)

      if (!asRaw) {
        output = convToImat(output)
      }

      if (is(pop, "Pop") | is(pop, "NamedMapPop")) {
        rownames(output) = pop@id
      } else {
        rownames(output) = as.character(1:pop@nInd)
      }

      colnames(output) = getLociNames(
        lociMap$lociPerChr,
        lociMap$lociLoc,
        genMap
      )

      output = output[, match(markers, colnames(output)), drop = FALSE]

      return(output)
    },
    markers = markers,
    asRaw = asRaw,
    simParam = simParam,
    nThreads = nThreads,
    simplify = simplify,
    level = level
  )
}

#' @title Pull SNP haplotypes
#'
#' @description Retrieves SNP haplotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param snpChip an integer. Indicates which SNP
#' chip's haplotypes to retrieve.
#' @param haplo either "all" for all haplotypes or an integer
#' for a single set of haplotypes. Use a value of 1 for female
#' haplotypes and a value of 2 for male haplotypes in diploids.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosomes are retrieved.
#' @param asRaw return in raw (byte) format
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns a matrix of SNP haplotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullSnpHaplo(pop, simParam=SP)
#' 
#' @export
pullSnpHaplo = function(pop, snpChip=1, haplo="all",
                        chr=NULL, asRaw=FALSE, simParam=NULL,
                        nThreads=NULL){
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  if(is.character(snpChip)){
    # Suspect snpChip is a name
    chipNames = simParam$snpChipNames
    take = match(snpChip, chipNames)
    if(is.na(take)){
      stop("'",snpChip,"' did not match any SNP chip names")
    }
    snpChip = take
  }
  
  tmp = selectLoci(chr,
                   simParam$snpChips[[snpChip]]@lociPerChr,
                   simParam$snpChips[[snpChip]]@lociLoc)
  
  lociPerChr = tmp$lociPerChr
  lociLoc = tmp$lociLoc
  
  haplo = checkHaplo(haplo, pop@ploidy)
  if(haplo=="all"){
    output = getHaplo(pop@geno,lociPerChr,lociLoc,nThreads)
    
    if(!asRaw){
      output = convToImat(output)
    }
    
    if(is(pop,"Pop")){
      rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(rep(1:pop@nInd,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }
  }else{
    output = getOneHaplo(pop@geno,lociPerChr,lociLoc,
                         haplo,nThreads)
    
    if(!asRaw){
      output = convToImat(output)
    }
    
    if(is(pop,"Pop")){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  
  colnames(output) = getLociNames(tmp$lociPerChr, tmp$lociLoc, simParam$genMap)
  
  return(output)
}

#' @title Pull QTL haplotypes
#'
#' @description Retrieves QTL haplotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param trait an integer or a trait name. Indicates which trait's
#' QTL haplotypes to retrieve.
#' @param haplo either "all" for all haplotypes or an integer
#' for a single set of haplotypes. Use a value of 1 for female
#' haplotypes and a value of 2 for male haplotypes in diploids.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosomes are retrieved.
#' @param asRaw return in raw (byte) format
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns a matrix of QTL haplotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullQtlHaplo(pop, simParam=SP)
#' 
#' @export
pullQtlHaplo = function(pop, trait=1, haplo="all",
                        chr=NULL, asRaw=FALSE, simParam=NULL,
                        nThreads=NULL){
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  if(is.character(trait)){
    # Suspect trait is a name
    traitNames = simParam$traitNames
    take = match(trait, traitNames)
    if(is.na(take)){
      stop("'",trait,"' did not match any trait names")
    }
    trait = take
  }
  
  tmp = selectLoci(chr,
                   simParam$traits[[trait]]@lociPerChr,
                   simParam$traits[[trait]]@lociLoc)
  
  lociPerChr = tmp$lociPerChr
  
  lociLoc = tmp$lociLoc
  
  haplo = checkHaplo(haplo, pop@ploidy)
  if(haplo=="all"){
    output = getHaplo(pop@geno,lociPerChr,lociLoc,nThreads)
    
    if(!asRaw){
      output = convToImat(output)
    }
    
    if(is(pop,"Pop")){
      rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(rep(1:pop@nInd,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }
  }else{
    output = getOneHaplo(pop@geno,lociPerChr,lociLoc,
                         haplo,nThreads)
    
    if(!asRaw){
      output = convToImat(output)
    }
    
    if(is(pop,"Pop")){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  
  colnames(output) = getLociNames(tmp$lociPerChr, tmp$lociLoc, simParam$genMap)
  
  return(output)
}

#' @title Pull seg site haplotypes
#'
#' @description
#' Retrieves haplotype data for all segregating sites
#'
#' @param pop an object of \code{\link{RawPop-class}} or
#' \code{\link{MapPop-class}}
#' @param haplo either "all" for all haplotypes or an integer
#' for a single set of haplotypes. Use a value of 1 for female
#' haplotypes and a value of 2 for male haplotypes in diploids.
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosomes are retrieved.
#' @param asRaw return in raw (byte) format
#' @param simParam an object of class \code{\link{SimParam}} for non-
#' \code{\link{MapPop-class}} inputs. If \code{NULL}, the function uses the
#' object named \code{SP} from the global environment. Ignored when
#' \code{pop} is a \code{\link{MapPop-class}}.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from
#' \code{getNumThreads()} for \code{\link{MapPop-class}} objects and
#' \code{simParam$nThreads} for \code{\link{RawPop-class}}.
#'
#' @return Returns a matrix of haplotypes
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullSegSiteHaplo(pop, simParam=SP)
#' 
#' @export
pullSegSiteHaplo = function(pop, haplo="all",
                            chr=NULL, asRaw=FALSE, simParam=NULL,
                            nThreads=NULL){
  if(is(pop,"MapPop")){
    allLoci = unlist(c(sapply(pop@nLoci, function(x) 1:x)))
    lociTot = pop@nLoci
    if(is.null(nThreads)){
      nThreads = getNumThreads()
    }else{
      nThreads = as.integer(nThreads)
    }
    map = pop@genMap
  }else{
    simParam = getSimParam(simParam)
    allLoci = unlist(c(sapply(simParam$segSites, function(x) 1:x)))
    lociTot = simParam$segSites
    if(is.null(nThreads)){
      nThreads = simParam$nThreads
    }else{
      nThreads = as.integer(nThreads)
    }
    map = simParam$genMap
  }
  if(!is.null(chr)){
    tmp = selectLoci(chr,lociTot,allLoci)
    lociTot = tmp$lociPerChr
    allLoci = tmp$lociLoc
  }else{
    chr = 1:pop@nChr
  }
  haplo = checkHaplo(haplo, pop@ploidy)
  if(haplo=="all"){
    output = getHaplo(pop@geno,
                      lociTot,
                      allLoci,
                      nThreads)
    if(!asRaw){
      output = convToImat(output)
    }
    
    if(is(pop,"Pop") | is(pop,"NamedMapPop")){
      rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(rep(1:pop@nInd,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }
    
  }else{
    output = getOneHaplo(pop@geno,
                         lociTot,
                         allLoci,
                         haplo,
                         nThreads)
    
    if(!asRaw){
      output = convToImat(output)
    }
    
    if(is(pop,"Pop") | is(pop,"NamedMapPop")){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  
  colnames(output) = getLociNames(lociTot, allLoci, map)
  
  return(output)
}


#' @title Pull IBD haplotypes
#'
#' @description Retrieves IBD haplotype data
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param chr a vector of chromosomes to retrieve. If NULL,
#' all chromosomes are retrieved.
#' @param snpChip an integer indicating which SNP array loci 
#' are to be retrieved. If NULL, all sites are retrieved.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns a matrix of IBD haplotypes.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' SP$setTrackRec(TRUE)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pullIbdHaplo(pop, simParam=SP)
#' 
#' @export
pullIbdHaplo = function(pop, chr=NULL, snpChip=NULL, simParam=NULL,
                        nThreads=NULL){
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  if(is.character(snpChip)){
    # Suspect snpChip is a name
    chipNames = simParam$snpChipNames
    take = match(snpChip, chipNames)
    if(is.na(take)){
      stop("'",snpChip,"' did not match any SNP chip names")
    }
    snpChip = take
  }
  
  if(!simParam$isTrackRec){
    stop(
      "pullIbdHaplo can only be used with the trackRec option, see trackRec in ?SimParam"
      )
  }
  
  if(is.null(chr)){
    chr = 1:pop@nChr
  }
  
  # Retrieve IBD data
  ibd = simParam$ibdHaplo(pop@iid)
  
  # Fill in output matrix
  output = createIbdMat(ibd=ibd, chr=chr, 
                        nLoci=pop@nLoci, 
                        ploidy=pop@ploidy,
                        nThreads=nThreads)
  
  rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                           rep(1:pop@ploidy,pop@nInd),sep="_")
  
  colnames(output) = unlist(lapply(simParam$genMap[chr], names))
  
  if(!is.null(snpChip)){
    nLoci = pop@nLoci[chr]
    tmp = getSnpMap(snpChip=snpChip,simParam=simParam)
    # getSnpMap reports chromosomes by name while chr is an index, so the
    # names are looked up before the two are matched. Comparing an index
    # to a name silently returned no markers for any named genetic map.
    chrNames = names(simParam$genMap)[chr]
    tmp = tmp[tmp$chr%in%chrNames,]
    if(length(chr)>1){
      for(i in 2:length(chr)){
        j = chrNames[i]
        tmp[tmp$chr==j,"site"] =
          tmp[tmp$chr==j,"site"] + sum(nLoci[1:(i-1)])
      }
    }
    output = output[,tmp$site,drop=FALSE]
  }
  
  return(output)
}

#' @title Pull marker haplotypes
#'
#' @description Retrieves haplotype data for user
#' specified loci
#'
#' @param pop an object of \code{\link{RawPop-class}} or
#' \code{\link{MapPop-class}}
#' @param markers a character vector. Indicates the
#' names of the loci to be retrieved
#' @param haplo either "all" for all haplotypes or an integer
#' for a single set of haplotypes. Use a value of 1 for female
#' haplotypes and a value of 2 for male haplotypes in diploids.
#' @param asRaw return in raw (byte) format
#' @param simParam an object of class \code{\link{SimParam}} for non-
#' \code{\link{MapPop-class}} inputs. If \code{NULL}, the function uses the
#' object named \code{SP} from the global environment. Ignored when
#' \code{pop} is a \code{\link{MapPop-class}}.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from
#' \code{getNumThreads()} for \code{\link{MapPop-class}} objects and
#' \code{simParam$nThreads} for \code{\link{RawPop-class}}.
#'
#' @return Returns a matrix of genotypes.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$addSnpChip(5)
#' SP$setTrackRec(TRUE)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Pull haplotype data for first two markers on chromosome one.
#' #Marker name is consistent with default naming in AlphaSimR.
#' pullMarkerHaplo(pop, markers=c("1_1","1_2"), simParam=SP)
#'
#' @export
pullMarkerHaplo = function(pop, markers, haplo="all", asRaw=FALSE,
                           simParam=NULL, nThreads=NULL){
  # Get genetic map and nThreads
  if(is(pop,"MapPop")){
    if(is.null(nThreads)){
      nThreads = getNumThreads()
    }else{
      nThreads = as.integer(nThreads)
    }
    genMap = pop@genMap
  }else{
    simParam = getSimParam(simParam)
    if(is.null(nThreads)){
      nThreads = simParam$nThreads
    }else{
      nThreads = as.integer(nThreads)
    }
    genMap = simParam$genMap
  }
  
  # Map markers to genetic map
  lociMap = mapLoci(markers, genMap)
  
  haplo = checkHaplo(haplo, pop@ploidy)
  if(haplo=="all"){
    output = getHaplo(pop@geno, lociMap$lociPerChr, lociMap$lociLoc, nThreads)
    
    if(!asRaw){
      output = convToImat(output)
    }
    
    if(is(pop,"Pop")){
      rownames(output) = paste(rep(pop@id,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(rep(1:pop@nInd,each=pop@ploidy),
                               rep(1:pop@ploidy,pop@nInd),sep="_")
    }
  }else{
    output = getOneHaplo(pop@geno, lociMap$lociPerChr, lociMap$lociLoc,
                         haplo, nThreads)
    
    if(!asRaw){
      output = convToImat(output)
    }
    
    if(is(pop,"Pop")){
      rownames(output) = paste(pop@id,rep(haplo,pop@nInd),sep="_")
    }else{
      rownames(output) = paste(1:pop@nInd,rep(haplo,pop@nInd),sep="_")
    }
  }
  
  colnames(output) = getLociNames(lociMap$lociPerChr, 
                                  lociMap$lociLoc, 
                                  genMap)
  
  output = output[,match(markers, colnames(output)),drop=FALSE]
  
  return(output)
}

#' @title Set marker haplotypes 
#'
#' @description Manually sets the haplotypes in a population 
#' for all individuals at one or more loci.
#'
#' @param pop an object of \code{\link{RawPop-class}} or
#' \code{\link{MapPop-class}}
#' @param haplo a matrix of haplotypes, see details
#' @param simParam an object of class \code{\link{SimParam}} for non-
#' \code{\link{MapPop-class}} inputs. If \code{NULL}, the function uses the
#' object named \code{SP} from the global environment. Ignored when
#' \code{pop} is a \code{\link{MapPop-class}}.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from
#' \code{getNumThreads()} for \code{\link{MapPop-class}} objects and
#' \code{simParam$nThreads} for \code{\link{RawPop-class}}.
#' 
#' @details The format of the haplotype matrix should match 
#' the format of the output from \code{\link{pullMarkerHaplo}}
#' with the option haplo="all". Thus, it is recommended that this 
#' function is first used to extract the haplotypes and that any 
#' desired changes be made to the output of pullMarkerHaplo before 
#' passing the matrix to setMarkerHaplo. Any changes made to QTL 
#' may potentially result in changes to an individual's genetic 
#' value. These changes will be reflected in the gv and/or gxe slot. 
#' All other slots will remain unchanged, so the ebv and pheno slots 
#' will not reflect the new genotypes.
#'
#' @return an object of the same class as the "pop" input
#'
#' @examples 
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' # Extract haplotypes for marker "1_1"
#' H = pullMarkerHaplo(founderPop, markers="1_1", nThreads=1L)
#' H[1,1]
#' 
#' # Set the first haplotype to 1
#' H[1,1] = 1L
#' 
#' # Set marker haplotypes
#' founderPop = setMarkerHaplo(founderPop, haplo=H, nThreads=1L)
#' pullMarkerHaplo(founderPop, markers="1_1", nThreads=1L)[1,1]
#' 
#' @export
setMarkerHaplo = function(pop, haplo, simParam=NULL, nThreads=NULL){
  # Check validity of rows
  if(nrow(haplo)!=(pop@nInd*pop@ploidy)){
    stop("nrow(haplo) must equal nInd(pop)*ploidy")
  }
  
  # Get genetic map
  if(is(pop,"MapPop")){
    genMap = pop@genMap
    if(is.null(nThreads)){
      nThreads = getNumThreads()
    }else{
      nThreads = as.integer(nThreads)
    }
  }else{
    simParam = getSimParam(simParam)
    genMap = simParam$genMap
    if(is.null(nThreads)){
      nThreads = simParam$nThreads
    }else{
      nThreads = as.integer(nThreads)
    }
  }
  
  # Map markers to the genetic map
  markers = colnames(haplo)
  lociMap = mapLoci(markers, genMap)
  
  # Order haplotype data
  orderedMapNames = getLociNames(lociMap$lociPerChr, 
                                 lociMap$lociLoc, 
                                 genMap)
  # The columns have to be put into the order setHaplo expects, which
  # means indexing haplo by where each ordered map name sits in markers.
  # Matching the other way round produces the inverse permutation.
  haplo = haplo[,match(orderedMapNames, markers), drop=FALSE]
  
  # Set haplotypes
  geno = setHaplo(pop@geno, haplo, lociMap$lociPerChr, 
                  lociMap$lociLoc, nThreads)
  dim(geno) = NULL # Account for matrix bug in RcppArmadillo
  pop@geno = geno
  
  if(is(pop, "Pop")){
    PHENO = pop@pheno
    EBV = pop@ebv
    pop = resetPop(pop=pop, simParam=simParam, nThreads=nThreads)
    pop@pheno = PHENO
    pop@ebv = EBV
  }
  
  return(pop)
}
