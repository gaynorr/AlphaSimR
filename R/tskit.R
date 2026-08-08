.recordedIbdHaplo = function(samplePops, chr, simParam){
  pedigree = simParam$pedigree
  recHist = simParam$recHist
  nIndividuals = length(recHist)
  needed = rep(FALSE, nIndividuals)
  frontier = unique(as.integer(unlist(
    lapply(samplePops, function(x) x@iid),
    use.names=FALSE
  )))
  while(length(frontier) > 0L){
    frontier = frontier[!needed[frontier]]
    if(length(frontier) == 0L){
      break
    }
    needed[frontier] = TRUE
    parents = unique(as.integer(pedigree[frontier, 1:2, drop=FALSE]))
    frontier = parents[parents > 0L]
  }

  neededIid = which(needed)
  founderIid = neededIid[
    pedigree[neededIid, 1L] == 0L & pedigree[neededIid, 2L] == 0L
  ]
  ibd = vector("list", nIndividuals)
  founderIbd = getFounderIbd(recHist[founderIid], samplePops[[1L]]@nChr)
  for(i in seq_along(founderIid)){
    ibd[[founderIid[i]]] = founderIbd[[i]]
  }
  for(iid in neededIid){
    if(pedigree[iid, 1L] == 0L && pedigree[iid, 2L] == 0L){
      next
    }
    mother = pedigree[iid, 1L]
    father = pedigree[iid, 2L]
    if(mother < 1L || father < 1L ||
       is.null(ibd[[mother]]) || is.null(ibd[[father]])){
      stop("Recorded pedigree is incomplete for the sampled population")
    }
    ibd[[iid]] = getNonFounderIbd(
      recHist=recHist[[iid]],
      mother=ibd[[mother]],
      father=ibd[[father]]
    )
  }

  do.call("rbind", lapply(samplePops, function(x){
    createIbdMat(
      ibd=ibd[x@iid],
      chr=chr,
      nLoci=x@nLoci,
      ploidy=x@ploidy,
      nThreads=simParam$nThreads
    )
  }))
}

#' Convert recorded ancestry to tree sequences
#'
#' Converts the pedigree and recombination history recorded by
#' \code{SimParam$setTrackRec(TRUE)} into one succinct tree sequence per
#' chromosome. AlphaSimR continues to use its native genotype representation;
#' conversion occurs only when this function is called.
#'
#' @param pop a \code{\link{Pop-class}} or \code{\link{MultiPop-class}}
#'   whose haplotypes should be marked as samples. MultiPop components may
#'   have different ploidies but must share chromosome and locus structure.
#' @param chr chromosomes to export. The default exports all chromosomes.
#' @param includeVariants include simulated sites and founder alleles. Set to
#'   \code{FALSE} to export ancestry only.
#' @param simplify simplify each tree sequence to the sampled haplotypes.
#' @param simParam the \code{\link{SimParam}} object that recorded the
#'   population.
#'
#' @return A named list of \code{RcppTskit::TreeSequence} objects, with class
#'   \code{AlphaSimRTreeSequence}.
#'
#' @details Coordinates are locus based. For a chromosome with \code{n} loci,
#' the sequence spans \code{[0, n)} and locus \code{i} is stored at
#' \code{i - 0.5}. When variants are requested, the current population is
#' checked against the founder alleles and recorded inheritance. Direct genome
#' changes made by functions such as \code{\link{mutate}} and
#' \code{\link{editGenome}} cannot currently be represented and produce an
#' informative error. Inbred founders retain their shared haplotype origins,
#' and samples may span generations and even ploidies in a MultiPop.
#' \code{HybridPop} is not supported because that lightweight class stores
#' genetic values rather than haplotypes and recombination histories.
#'
#' @examples
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$setTrackRec(TRUE)
#' pop = newPop(founderPop, simParam=SP)
#' progeny = randCross(pop, nCrosses=2, simParam=SP)
#' ts = asTreeSequence(progeny, simParam=SP)
#'
#' @export
asTreeSequence = function(pop, chr=NULL, includeVariants=TRUE,
                          simplify=TRUE, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP", envir=.GlobalEnv)
  }
  if(is(pop, "HybridPop")){
    stop(paste0(
      "HybridPop does not store haplotypes or recombination history; ",
      "create it with returnHybridPop=FALSE"
    ))
  }
  flattenPops = function(x){
    if(is(x, "Pop")){
      return(list(x))
    }
    if(is(x, "MultiPop")){
      return(unlist(lapply(x@pops, flattenPops), recursive=FALSE))
    }
    stop("pop must be a Pop or MultiPop object")
  }
  samplePops = flattenPops(pop)
  if(length(samplePops) == 0L){
    stop("pop must contain at least one population")
  }
  nChr = vapply(samplePops, function(x) x@nChr, integer(1))
  if(any(nChr != nChr[1L])){
    stop("MultiPop components must have the same chromosomes")
  }
  nChr = nChr[1L]
  nLoci = lapply(samplePops, function(x) x@nLoci)
  if(any(vapply(
    nLoci,
    function(x) !identical(x, nLoci[[1L]]),
    logical(1)
  ))){
    stop("MultiPop components must have the same loci per chromosome")
  }
  nLoci = nLoci[[1L]]
  sampleIid = as.integer(unlist(
    lapply(samplePops, function(x) x@iid),
    use.names=FALSE
  ))
  if(length(sampleIid) == 0L){
    stop("pop must contain at least one individual")
  }
  if(anyDuplicated(sampleIid)){
    stop("pop must contain unique individual IDs")
  }
  if(!simParam$isTrackRec){
    stop("asTreeSequence requires SP$setTrackRec(TRUE) before creating populations")
  }
  stopifnot(length(includeVariants)==1L, !is.na(includeVariants),
            is.logical(includeVariants),
            length(simplify)==1L, !is.na(simplify), is.logical(simplify))

  if(is.null(chr)){
    chr = seq_len(nChr)
  }
  chr = as.integer(chr)
  if(length(chr) == 0L || anyNA(chr) || anyDuplicated(chr) ||
     any(chr < 1L | chr > nChr)){
    stop("chr contains invalid or duplicated chromosome numbers")
  }

  recHist = simParam$recHist
  pedigree = simParam$pedigree
  if(length(recHist) != nrow(pedigree) ||
     any(sampleIid < 1L | sampleIid > length(recHist))){
    stop("pop and the recorded recombination history do not match")
  }
  for(history in recHist){
    if(is.list(history)){
      for(chromosome in chr){
        chromosomeHistory = history[[chromosome]]
        if(any(vapply(
          chromosomeHistory,
          function(segments) any(segments[,1L] >= 100L),
          logical(1)
        ))){
          stop(paste0(
            "Recombination history contains unresolved quadrivalent ",
            "homolog labels and must be regenerated"
          ))
        }
      }
    }
  }

  founderIid = as.integer(which(
    pedigree[,1L] == 0L & pedigree[,2L] == 0L
  ))
  nFounder = simParam$founderPop@nInd
  expectedFounderIid = seq_len(nFounder)
  if(includeVariants && !identical(founderIid, expectedFounderIid)){
    stop(paste0(
      "Variant export currently requires the original SimParam founder ",
      "population to be the only recorded founders; use includeVariants=FALSE"
    ))
  }

  if(includeVariants){
    ibd = .recordedIbdHaplo(samplePops, chr, simParam)
    founderHap = pullSegSiteHaplo(simParam$founderPop, chr=chr)
    currentHap = do.call("rbind", lapply(
      samplePops,
      pullSegSiteHaplo,
      chr=chr,
      simParam=simParam
    ))
    unrecordedGenomeMessage = paste0(
      "Current genomes contain changes not represented by recombination ",
      "history (for example mutate() or editGenome()); use ",
      "includeVariants=FALSE for ancestry-only export"
    )
    if(!identical(dim(ibd), dim(currentHap)) ||
       ncol(founderHap) != ncol(currentHap)){
      stop(unrecordedGenomeMessage)
    }
    founderOrigins = as.integer(unlist(
      recHist[expectedFounderIid],
      use.names=FALSE
    ))
    if(length(founderOrigins) != nrow(founderHap) ||
       anyNA(founderOrigins) || any(founderOrigins < 1L)){
      stop("Recorded founder origins do not match SimParam$founderPop")
    }
    uniqueOrigins = unique(founderOrigins)
    originRows = match(uniqueOrigins, founderOrigins)
    originHap = founderHap[originRows,,drop=FALSE]
    copyOriginRows = match(founderOrigins, uniqueOrigins)
    for(row in seq_len(nrow(founderHap))){
      if(!identical(
        unname(founderHap[row,]),
        unname(originHap[copyOriginRows[row],])
      )){
        stop("Copies of a founder haplotype origin contain different alleles")
      }
    }
    ibdOriginRows = matrix(
      match(as.integer(ibd), uniqueOrigins),
      nrow=nrow(ibd),
      ncol=ncol(ibd)
    )
    if(anyNA(ibdOriginRows)){
      stop("Recorded founder origins do not match SimParam$founderPop")
    }
    expectedHap = matrix(0L, nrow=nrow(ibd), ncol=ncol(ibd))
    for(j in seq_len(ncol(ibd))){
      expectedHap[,j] = originHap[ibdOriginRows[,j],j]
    }
    if(!identical(unname(currentHap), expectedHap)){
      stop(unrecordedGenomeMessage)
    }
  }

  individualId = names(recHist)
  if(is.null(individualId) || length(individualId) != length(recHist)){
    individualId = as.character(seq_along(recHist))
  }
  version = as.character(utils::packageVersion("AlphaSimR"))
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")

  output = vector("list", length(chr))
  names(output) = paste0("chr", chr)
  for(i in seq_along(chr)){
    if(includeVariants){
      founderCopyChr = pullSegSiteHaplo(
        simParam$founderPop,
        chr=chr[i]
      )
      founderChr = founderCopyChr[originRows,,drop=FALSE]
    }else{
      founderChr = matrix(integer(), nrow=0L, ncol=0L)
    }
    xptr = buildTreeSequenceCpp(
      recHist=recHist,
      pedigree=pedigree,
      sampleIid=sampleIid,
      individualId=individualId,
      chromosome=chr[i]-1L,
      nLoci=nLoci[chr[i]],
      founderIid=founderIid,
      founderHaplotypes=founderChr,
      includeVariants=includeVariants,
      simplify=simplify,
      version=version,
      timestamp=timestamp
    )
    if(utils::packageVersion("RcppTskit") >= "0.3.0"){
      tables = RcppTskit::TableCollection$new(xptr=xptr)
    }else{
      tables = RcppTskit::TableCollection$new(pointer=xptr)
    }
    output[[i]] = tables$tree_sequence()
  }
  attr(output, "chromosome") = chr
  class(output) = c("AlphaSimRTreeSequence", "list")
  return(output)
}

#' Write AlphaSimR tree sequences
#'
#' @param x an object returned by \code{\link{asTreeSequence}}.
#' @param file output path. For multiple chromosomes, chromosome labels are
#'   inserted before the \code{.trees} extension.
#' @param overwrite overwrite existing files.
#'
#' @return The written paths, invisibly.
#'
#' @export
writeTreeSequence = function(x, file, overwrite=FALSE){
  if(!inherits(x, "AlphaSimRTreeSequence") || length(x) < 1L){
    stop("x must be a non-empty result from asTreeSequence()")
  }
  if(!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)){
    stop("file must be one non-empty path")
  }
  if(!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)){
    stop("overwrite must be TRUE or FALSE")
  }

  chromosomes = attr(x, "chromosome")
  if(length(x) == 1L){
    paths = file
  }else{
    stem = sub("\\.trees$", "", file, ignore.case=TRUE)
    paths = paste0(stem, "_chr", chromosomes, ".trees")
  }
  directories = unique(dirname(paths))
  if(any(!dir.exists(directories))){
    stop("Output directory does not exist")
  }
  if(!overwrite && any(file.exists(paths))){
    stop("Output file already exists; use overwrite=TRUE")
  }
  for(i in seq_along(x)){
    x[[i]]$dump(paths[i])
  }
  invisible(normalizePath(paths, mustWork=TRUE))
}
