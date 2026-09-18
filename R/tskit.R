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
#' checked against the founder alleles and recorded inheritance. Changes made
#' by functions such as \code{\link{mutate}} and \code{\link{editGenome}} are
#' represented as mutations on terminal sample nodes. Alleles from
#' additional founders are inferred from sampled descendants where possible.
#' Known founder alleles are used only when sampled tracked founders verify
#' their correspondence to \code{SimParam$founderPop}; otherwise founder
#' alleles are inferred conservatively from the requested samples. Inbred
#' founders retain their shared haplotype origins, and samples may span
#' generations and even ploidies in a MultiPop. Sample and recorded-individual
#' node times are pedigree depths measured backwards from the deepest recorded
#' individual. Synthetic founder-origin roots and sample ancestry proxies use
#' the explicitly documented offsets stored in the top-level metadata.
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
  if(!inherits(simParam, "SimParam")){
    stop("simParam must be a SimParam object")
  }
  if(!is.logical(includeVariants) || length(includeVariants) != 1L ||
     is.na(includeVariants)){
    stop("includeVariants must be TRUE or FALSE")
  }
  if(!is.logical(simplify) || length(simplify) != 1L || is.na(simplify)){
    stop("simplify must be TRUE or FALSE")
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
  if(nChr != simParam$nChr ||
     !identical(nLoci, simParam$founderPop@nLoci)){
    stop("pop and simParam chromosome/locus structures do not match")
  }
  sampleIid = as.integer(unlist(
    lapply(samplePops, function(x) x@iid),
    use.names=FALSE
  ))
  sampleId = as.character(unlist(
    lapply(samplePops, function(x) x@id),
    use.names=FALSE
  ))
  if(length(sampleIid) == 0L){
    stop("pop must contain at least one individual")
  }
  if(anyNA(sampleIid)){
    stop("pop contains missing individual IDs")
  }
  if(anyDuplicated(sampleIid)){
    stop("pop must contain unique individual IDs")
  }
  if(length(sampleId) != length(sampleIid) || anyNA(sampleId)){
    stop("pop must contain one non-missing ID per individual")
  }
  if(!simParam$isTrackRec){
    stop("asTreeSequence requires SP$setTrackRec(TRUE) before creating populations")
  }
  if(is.null(chr)){
    chr = seq_len(nChr)
  }
  if(!is.numeric(chr) || anyNA(chr) || any(chr != as.integer(chr))){
    stop("chr must contain whole chromosome numbers")
  }
  chr = as.integer(chr)
  if(length(chr) == 0L || anyNA(chr) || anyDuplicated(chr) ||
     any(chr < 1L | chr > nChr)){
    stop("chr contains invalid or duplicated chromosome numbers")
  }

  recHist = simParam$recHist
  pedigree = simParam$pedigree
  if(ncol(pedigree) < 2L || anyNA(pedigree[,1:2,drop=FALSE]) ||
     length(recHist) != nrow(pedigree) ||
     any(sampleIid < 1L | sampleIid > length(recHist))){
    stop("pop and the recorded recombination history do not match")
  }
  for(samplePop in samplePops){
    for(iid in samplePop@iid){
      history = recHist[[iid]]
      recordedPloidy = if(is.integer(history)){
        rep(length(history), length(chr))
      }else if(is.list(history)){
        vapply(chr, function(chromosome){
          length(history[[chromosome]])
        }, integer(1))
      }else{
        integer()
      }
      if(length(recordedPloidy) != length(chr) ||
         any(recordedPloidy != samplePop@ploidy)){
        stop("pop ploidy and the recorded recombination history do not match")
      }
    }
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
  if(includeVariants){
    ibd = .recordedIbdHaplo(samplePops, chr, simParam)
    sortedChr = sort(chr)
    blockEnd = cumsum(nLoci[sortedChr])
    blockStart = c(1L, utils::head(blockEnd, -1L) + 1L)
    blockOrder = match(chr, sortedChr)
    columnOrder = unlist(Map(
      seq.int,
      blockStart[blockOrder],
      blockEnd[blockOrder]
    ), use.names=FALSE)
    founderCopyHap = pullSegSiteHaplo(
      simParam$founderPop,
      chr=sortedChr
    )[,columnOrder,drop=FALSE]
    currentHap = do.call("rbind", lapply(samplePops, function(x){
      pullSegSiteHaplo(
        x,
        chr=sortedChr,
        simParam=simParam
      )[,columnOrder,drop=FALSE]
    }))
    sampleRowIid = as.integer(unlist(lapply(samplePops, function(x){
      rep(x@iid, each=x@ploidy)
    }), use.names=FALSE))
    if(!identical(dim(ibd), dim(currentHap)) ||
       ncol(founderCopyHap) != ncol(currentHap)){
      stop("Current genomes and recorded ancestry have different dimensions")
    }
    if(anyNA(currentHap) || any(currentHap != 0L & currentHap != 1L)){
      stop("Current haplotypes must contain only 0 and 1")
    }

    founderOrigins = as.integer(unlist(
      recHist[founderIid],
      use.names=FALSE
    ))
    if(anyNA(founderOrigins) || any(founderOrigins < 1L)){
      stop("Recorded founder origins must be positive integers")
    }
    uniqueOrigins = unique(founderOrigins)
    originHap = matrix(
      0L,
      nrow=length(uniqueOrigins),
      ncol=ncol(currentHap)
    )
    originKnown = rep(FALSE, nrow(originHap))

    originalFoundersRecorded =
      length(expectedFounderIid) <= nrow(pedigree) &&
      all(pedigree[expectedFounderIid, 1L] == 0L) &&
      all(pedigree[expectedFounderIid, 2L] == 0L) &&
      all(vapply(recHist[expectedFounderIid], is.integer, logical(1))) &&
      all(expectedFounderIid %in% sampleIid)
    if(originalFoundersRecorded){
      founderSampleRows = as.integer(unlist(lapply(
        expectedFounderIid,
        function(iid) which(sampleRowIid == iid)
      ), use.names=FALSE))
      originalFoundersRecorded =
        length(founderSampleRows) == nrow(founderCopyHap) &&
        identical(
          unname(currentHap[founderSampleRows,,drop=FALSE]),
          unname(founderCopyHap)
        )
    }
    if(originalFoundersRecorded){
      originalOrigins = as.integer(unlist(
        recHist[expectedFounderIid],
        use.names=FALSE
      ))
      if(length(originalOrigins) != nrow(founderCopyHap)){
        stop("Recorded original founders do not match SimParam$founderPop")
      }
      for(row in seq_len(nrow(founderCopyHap))){
        originRow = match(originalOrigins[row], uniqueOrigins)
        if(is.na(originRow)){
          stop("Recorded founder origins are internally inconsistent")
        }
        if(originKnown[originRow] &&
           !identical(
             unname(founderCopyHap[row,]),
             unname(originHap[originRow,])
           )){
          stop("Copies of a founder haplotype origin contain different alleles")
        }
        originHap[originRow,] = founderCopyHap[row,]
        originKnown[originRow] = TRUE
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

    variantEncoding = resolveVariantEncodingCpp(
      originRows=ibdOriginRows,
      currentHaplotypes=currentHap,
      originHaplotypes=originHap,
      knownOrigins=originKnown
    )
    originHap = variantEncoding$originHaplotypes
    sampleAlleleOverrides = variantEncoding$sampleAlleleOverrides
  }

  individualId = names(recHist)
  if(is.null(individualId) || length(individualId) != length(recHist)){
    individualId = as.character(seq_along(recHist))
  }else{
    missingId = is.na(individualId)
    individualId[missingId] = as.character(which(missingId))
  }
  individualId[sampleIid] = sampleId
  individualId = enc2utf8(individualId)
  if(anyNA(iconv(individualId, from="UTF-8", to="UTF-8"))){
    stop("Individual IDs must be valid UTF-8 strings")
  }
  version = as.character(utils::packageVersion("AlphaSimR"))
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC")

  output = vector("list", length(chr))
  names(output) = paste0("chr", chr)
  columnOffset = 0L
  for(i in seq_along(chr)){
    if(includeVariants){
      columns = columnOffset + seq_len(nLoci[chr[i]])
      originChr = originHap[,columns,drop=FALSE]
      overridesChr = sampleAlleleOverrides[,columns,drop=FALSE]
      columnOffset = columnOffset + nLoci[chr[i]]
    }else{
      originChr = matrix(integer(), nrow=0L, ncol=0L)
      overridesChr = matrix(integer(), nrow=0L, ncol=0L)
    }
    xptr = buildTreeSequenceCpp(
      recHist=recHist,
      pedigree=pedigree,
      sampleIid=sampleIid,
      individualId=individualId,
      chromosome=chr[i]-1L,
      nLoci=nLoci[chr[i]],
      founderIid=founderIid,
      originHaplotypes=originChr,
      sampleAlleleOverrides=overridesChr,
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
  attr(output, "coordinate_system") = "locus_index"
  attr(output, "sample_iid") = sampleIid
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
