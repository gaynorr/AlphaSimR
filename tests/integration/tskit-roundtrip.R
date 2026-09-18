library(AlphaSimR)

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 1L){
  stop("Usage: Rscript tests/integration/tskit-roundtrip.R OUTPUT_DIR")
}
outputDir = args[[1L]]
dir.create(outputDir, recursive=TRUE, showWarnings=FALSE)
caseCount = 0L

newSP = function(founderPop){
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  SP
}

writeExpectedStructure = function(label, samplePops, SP, chr){
  expectedOrigins = AlphaSimR:::.recordedIbdHaplo(
    samplePops,
    chr,
    SP
  )
  write.table(
    expectedOrigins,
    file.path(outputDir, paste0(label, ".origins.tsv")),
    sep="\t",
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE
  )
  pedigree = SP$pedigree[,1:2,drop=FALSE]
  depth = integer(nrow(pedigree))
  for(iid in seq_len(nrow(pedigree))){
    parents = pedigree[iid,]
    if(all(parents > 0L)){
      depth[iid] = max(depth[parents]) + 1L
    }
  }
  maxDepth = max(depth)
  expectedSamples = do.call(rbind, lapply(samplePops, function(pop){
    cbind(
      iid=rep(pop@iid, each=pop@ploidy),
      homolog=rep(seq_len(pop@ploidy), times=pop@nInd),
      time=rep(maxDepth-depth[pop@iid], each=pop@ploidy)
    )
  }))
  write.table(
    expectedSamples,
    file.path(outputDir, paste0(label, ".samples.tsv")),
    sep="\t",
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE
  )
  expectedIds = as.character(unlist(lapply(samplePops, function(pop){
    rep(pop@id, each=pop@ploidy)
  }), use.names=FALSE))
  idHex = vapply(expectedIds, function(id){
    paste0(
      "x",
      paste(
        sprintf("%02x", as.integer(charToRaw(enc2utf8(id)))),
        collapse=""
      )
    )
  }, character(1))
  writeLines(
    idHex,
    file.path(outputDir, paste0(label, ".ids.txt")),
    useBytes=TRUE
  )
  write.table(
    pedigree,
    file.path(outputDir, paste0(label, ".pedigree.tsv")),
    sep="\t",
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE
  )
}

emit = function(name, pop, SP, chromosomes=seq_len(pop@nChr),
                 simplifyValues=c(TRUE, FALSE)){
  for(chr in chromosomes){
    expected = pullSegSiteHaplo(pop, chr=chr, simParam=SP)
    for(simplify in simplifyValues){
      label = paste0(
        name,
        "_chr", chr,
        if(simplify) "_simple" else "_full"
      )
      ts = asTreeSequence(
        pop,
        chr=chr,
        simplify=simplify,
        simParam=SP
      )
      stopifnot(
        ts[[1]]$num_samples() == nrow(expected),
        ts[[1]]$num_sites() == ncol(expected),
        ts[[1]]$sequence_length() == pop@nLoci[chr]
      )
      writeTreeSequence(
        ts,
        file.path(outputDir, paste0(label, ".trees")),
        overwrite=TRUE
      )
      write.table(
        expected,
        file.path(outputDir, paste0(label, ".tsv")),
        sep="\t",
        row.names=FALSE,
        col.names=FALSE,
        quote=FALSE
      )
      writeExpectedStructure(label, list(pop), SP, chr)
      caseCount <<- caseCount + 1L
    }
  }
}

emitTogether = function(name, pop, SP, chromosomes){
  for(simplify in c(TRUE, FALSE)){
    suffix = if(simplify) "simple" else "full"
    ts = asTreeSequence(
      pop,
      chr=chromosomes,
      simplify=simplify,
      simParam=SP
    )
    stopifnot(
      identical(names(ts), paste0("chr", chromosomes)),
      identical(attr(ts, "chromosome"), as.integer(chromosomes))
    )
    paths = writeTreeSequence(
      ts,
      file.path(outputDir, paste0(name, "_", suffix, ".trees")),
      overwrite=TRUE
    )
    stopifnot(length(paths) == length(chromosomes))
    for(i in seq_along(chromosomes)){
      chr = chromosomes[i]
      label = paste0(name, "_", suffix, "_chr", chr)
      expected = pullSegSiteHaplo(pop, chr=chr, simParam=SP)
      stopifnot(
        ts[[i]]$num_samples() == nrow(expected),
        ts[[i]]$num_sites() == ncol(expected),
        ts[[i]]$sequence_length() == pop@nLoci[chr]
      )
      write.table(
        expected,
        file.path(outputDir, paste0(label, ".tsv")),
        sep="\t",
        row.names=FALSE,
        col.names=FALSE,
        quote=FALSE
      )
      writeExpectedStructure(label, list(pop), SP, chr)
      caseCount <<- caseCount + 1L
    }
  }
}

emitMulti = function(name, pops, SP, chromosomes){
  multi = do.call(newMultiPop, pops)
  flattenPops = function(pop){
    if(is(pop, "Pop")){
      return(list(pop))
    }
    unlist(lapply(pop@pops, flattenPops), recursive=FALSE)
  }
  samplePops = unlist(lapply(pops, flattenPops), recursive=FALSE)
  for(chr in chromosomes){
    expected = do.call(rbind, lapply(
      samplePops,
      pullSegSiteHaplo,
      chr=chr,
      simParam=SP
    ))
    for(simplify in c(TRUE, FALSE)){
      label = paste0(
        name,
        "_chr", chr,
        if(simplify) "_simple" else "_full"
      )
      ts = asTreeSequence(
        multi,
        chr=chr,
        simplify=simplify,
        simParam=SP
      )
      writeTreeSequence(
        ts,
        file.path(outputDir, paste0(label, ".trees")),
        overwrite=TRUE
      )
      write.table(
        expected,
        file.path(outputDir, paste0(label, ".tsv")),
        sep="\t",
        row.names=FALSE,
        col.names=FALSE,
        quote=FALSE
      )
      writeExpectedStructure(label, samplePops, SP, chr)
      caseCount <<- caseCount + 1L
    }
  }
}

set.seed(20260808)

# Variable chromosome lengths, sex-specific maps, generations, and ordering.
founder = quickHaplo(8, 3, c(7, 19, 31))
SP = newSP(founder)
SP$setSexes("yes_sys")
SP$setRecombRatio(2)
pop = newPop(founder, simParam=SP)
generation1 = randCross(pop, nCrosses=8, ignoreSexes=FALSE, simParam=SP)
generation2 = randCross(generation1, nCrosses=6, simParam=SP)
emit("diploid_multigeneration", generation2, SP, c(3, 1, 2))
emitTogether("permuted_chromosomes", generation2, SP, c(3, 1, 2))
emit("reordered_samples", generation2[c(6, 1, 4, 2)], SP, c(1, 3))
emit(
  "mixed_generation",
  mergePops(list(pop[1:2], generation1[1:2], generation2[1:2])),
  SP,
  2
)

# Separate male/female maps, random sexes, crossover interference, and a
# non-interfering crossover pathway.
founder = quickHaplo(12, 2, c(23, 37))
SP = newSP(founder)
SP$setSexes("yes_rand")
femaleMap = list(
  seq(0, 4, length.out=23),
  seq(0, 2.5, length.out=37)
)
maleMap = list(
  seq(0, 0.2, length.out=23),
  seq(0, 0.4, length.out=37)
)
for(chr in seq_along(femaleMap)){
  names(femaleMap[[chr]]) = names(SP$femaleMap[[chr]])
  names(maleMap[[chr]]) = names(SP$femaleMap[[chr]])
}
SP$switchFemaleMap(femaleMap)
SP$switchMaleMap(maleMap)
SP$v = 1
SP$p = 0.35
pop = newPop(founder, simParam=SP)
if(!all(c("F", "M") %in% pop@sex)){
  pop@sex[1:2] = c("F", "M")
}
sexSpecific = randCross(
  pop,
  nCrosses=12,
  nProgeny=2,
  balance=FALSE,
  simParam=SP
)
emit("sex_specific_maps_interference", sexSpecific, SP)
emit(
  "sex_specific_second_generation",
  randCross(
    sexSpecific,
    nCrosses=8,
    ignoreSexes=TRUE,
    simParam=SP
  ),
  SP
)

founder = quickHaplo(6, 1, 10)
SP = newSP(founder)
highRecombinationMap = list(seq(0, 20, length.out=10))
names(highRecombinationMap[[1]]) = names(SP$genMap[[1]])
SP$switchGenMap(highRecombinationMap)
pop = newPop(founder, simParam=SP)
emit(
  "dense_crossovers_between_loci",
  randCross(pop, nCrosses=6, simParam=SP),
  SP
)

# Shared inbred origins, selfing, doubled haploids, and hybrid Pop output.
founder = quickHaplo(5, 2, c(13, 21), inbred=TRUE)
SP = newSP(founder)
pop = newPop(founder, simParam=SP)
emit("inbred_founders", pop, SP)
emit("inbred_hybrid", hybridCross(pop[1:2], pop[3:5], simParam=SP), SP)
emit("selfed", self(pop, nProgeny=2, simParam=SP), SP)
emit("dh_female", makeDH(pop, nDH=2, useFemale=TRUE, simParam=SP), SP)
emit("dh_male", makeDH(pop, nDH=2, useFemale=FALSE, simParam=SP), SP)

# Minimal chromosomes, user IDs requiring JSON escaping, and a founder subset
# that does not initially include the full SimParam founder population.
founder = quickHaplo(4, 1, 1)
SP = newSP(founder)
pop = newPop(
  founder,
  id=c("quote\"id", "back\\slash", "line\nbreak", "unicode-α"),
  simParam=SP
)
emit("single_locus_special_ids", randCross(
  pop,
  nCrosses=3,
  ignoreSexes=TRUE,
  simParam=SP
), SP)

founder = quickHaplo(8, 2, c(9, 16), inbred=TRUE)
SP = newSP(founder)
founderSubset = newPop(founder[c(8, 3, 5)], simParam=SP)
emit("initial_inbred_founder_subset", founderSubset, SP)
emit("initial_inbred_founder_subset_descendants", randCross(
  founderSubset,
  nCrosses=5,
  simParam=SP
), SP)

# Even autopolyploids and forced/probabilistic quadrivalent pairing.
for(ploidy in c(4L, 6L, 8L)){
  founder = quickHaplo(6, 2, c(17, 25), ploidy=ploidy)
  SP = newSP(founder)
  SP$quadProb = if(ploidy == 4L) 1 else 0.5
  pop = newPop(founder, simParam=SP)
  polyploidProgeny = randCross(pop, nCrosses=5, simParam=SP)
  emit(
    paste0("autopolyploid_", ploidy),
    polyploidProgeny,
    SP
  )
  emit(
    paste0("autopolyploid_self_", ploidy),
    self(
      polyploidProgeny[1:2],
      nProgeny=2,
      keepParents=FALSE,
      simParam=SP
    ),
    SP,
    1
  )
  emit(
    paste0("autopolyploid_second_generation_", ploidy),
    randCross(polyploidProgeny, nCrosses=4, simParam=SP),
    SP,
    2
  )
}

# Genome doubling, reduction, merging, and unequal-ploidy offspring.
founder = quickHaplo(5, 2, c(15, 29))
SP = newSP(founder)
diploid = newPop(founder, simParam=SP)
tetraploid = doubleGenome(diploid, simParam=SP)
octoploid = doubleGenome(tetraploid, simParam=SP)
emit("double_tetraploid", tetraploid, SP)
emit("double_octoploid", octoploid, SP)
emit("reduce_recombining", reduceGenome(
  tetraploid,
  nProgeny=2,
  simRecomb=TRUE,
  simParam=SP
), SP)
emit("reduce_nonrecombining", reduceGenome(
  tetraploid,
  nProgeny=2,
  simRecomb=FALSE,
  simParam=SP
), SP)

haploid = reduceGenome(
  diploid,
  nProgeny=2,
  useFemale=FALSE,
  keepParents=FALSE,
  simRecomb=TRUE,
  simParam=SP
)
emit("diploid_to_haploid_male_map", haploid, SP)
emit(
  "haploid_genome_doubling",
  doubleGenome(haploid, keepParents=FALSE, simParam=SP),
  SP
)

crossPlan = matrix(c(1, 2, 2, 3), ncol=2)
mergedTetraploid = mergeGenome(diploid, diploid, crossPlan, simParam=SP)
mergedHexaploid = mergeGenome(diploid, tetraploid, crossPlan, simParam=SP)
triploidA = makeCross2(diploid, tetraploid, crossPlan, simParam=SP)
triploidB = makeCross2(tetraploid, diploid, crossPlan, simParam=SP)
emit("merge_diploid_diploid", mergedTetraploid, SP)
emit("merge_diploid_tetraploid", mergedHexaploid, SP)
emit("triploid_maternal_diploid", triploidA, SP)
emit("triploid_maternal_tetraploid", triploidB, SP)
emit(
  "triploid_genome_doubling",
  doubleGenome(triploidA, keepParents=FALSE, simParam=SP),
  SP
)
emit(
  "merge_triploid_diploid",
  mergeGenome(triploidA, diploid, crossPlan, simParam=SP),
  SP
)
emit(
  "cross_tetraploid_hexaploid",
  makeCross2(
    tetraploid,
    mergedHexaploid,
    matrix(c(1L, 1L, 2L, 2L), ncol=2),
    simParam=SP
  ),
  SP
)
emitMulti(
  "multipop_mixed_ploidy_generation",
  list(diploid[c(5, 1)], tetraploid[2:3], triploidA[1]),
  SP,
  c(2, 1)
)
emitMulti(
  "nested_multipop",
  list(newMultiPop(
    diploid[5],
    newMultiPop(tetraploid[2], triploidA[1])
  )),
  SP,
  1
)

# Current-genome overlays and founders introduced after the original founders.
knownHap = pullSegSiteHaplo(diploid, chr=1, simParam=SP)
knownOne = which(colSums(knownHap[1:2,,drop=FALSE]) > 0L)[1L]
edited = editGenome(
  diploid,
  ind=1,
  chr=1,
  segSites=knownOne,
  allele=0L,
  simParam=SP
)
emit("edited_back_mutation", edited, SP, 1)
emit("mutated_all_sites", mutate(diploid, mutRate=1, simParam=SP), SP, 1)

# Retaining every unchanged original founder lets the exporter verify the
# founder-origin mapping and place a descendant's back mutation on its shared
# ancestral state rather than treating the mapping as known by convention.
backMutationChild = randCross(diploid, nCrosses=1, simParam=SP)
childHap = pullSegSiteHaplo(backMutationChild, chr=1, simParam=SP)
childOne = which(colSums(childHap) > 0L)[1L]
stopifnot(!is.na(childOne))
editedBackMutationChild = editGenome(
  backMutationChild,
  ind=1,
  chr=1,
  segSites=childOne,
  allele=0L,
  simParam=SP
)
emitMulti(
  "verified_founders_back_mutation",
  list(diploid, editedBackMutationChild),
  SP,
  1
)

# A historical sample edited after its descendants were generated must not
# alter those descendants in the exported tree sequence.
homozygousLocus = which(knownHap[1,] == knownHap[2,])[1L]
stopifnot(!is.na(homozygousLocus))
preEditDescendants = makeCross(
  diploid,
  crossPlan=matrix(c(1L, 2L), ncol=2),
  nProgeny=4,
  simParam=SP
)
postHocEditedAncestor = editGenome(
  diploid,
  ind=1,
  chr=1,
  segSites=homozygousLocus,
  allele=1L-knownHap[1,homozygousLocus],
  simParam=SP
)
emitMulti(
  "posthoc_edited_ancestor_and_descendants",
  list(postHocEditedAncestor[1], preEditDescendants),
  SP,
  1
)

ancestryOnly = asTreeSequence(
  edited,
  chr=1,
  includeVariants=FALSE,
  simParam=SP
)
writeTreeSequence(
  ancestryOnly,
  file.path(outputDir, "ancestry_only.trees"),
  overwrite=TRUE
)
writeExpectedStructure("ancestry_only", list(edited), SP, 1)
caseCount = caseCount + 1L
additionalFounders = newPop(founder, simParam=SP)
additionalProgeny = randCross(additionalFounders, nCrosses=4, simParam=SP)
emit("additional_founders", additionalFounders, SP)
emit("additional_founder_descendants", additionalProgeny, SP)

# User-supplied pedigree and reset pedigree/cache state.
founder = quickHaplo(6, 2, c(18, 27))
SP = newSP(founder)
pop = newPop(founder, simParam=SP)
pedigree = pedigreeCross(
  pop,
  id=c("A", "B", "C", "D", "E"),
  mother=c("0", "0", "A", "C", "D"),
  father=c("0", "0", "B", "C", "C"),
  DH=c(FALSE, FALSE, FALSE, FALSE, TRUE),
  nSelf=c(0, 0, 0, 1, 0),
  simParam=SP
)
emit("pedigree_cross", pedigree, SP)

founder = quickHaplo(5, 2, c(12, 22))
SP = newSP(founder)
namedFounders = newPop(
  founder,
  id=paste0("F", seq_len(founder@nInd)),
  simParam=SP
)
matchedPedigree = pedigreeCross(
  namedFounders,
  id=c("F1", "F2", "F3", "C1", "C2"),
  mother=c("0", "0", "0", "F1", "C1"),
  father=c("0", "0", "0", "F2", "F3"),
  matchID=TRUE,
  nSelf=c(0, 0, 0, 1, 0),
  simParam=SP
)
emit("pedigree_cross_matched_founders", matchedPedigree, SP)

founder = quickHaplo(4, 1, 11)
SP = newSP(founder)
oldFounders = newPop(founder, simParam=SP)
invisible(randCross(oldFounders, nCrosses=4, simParam=SP))
SP$resetPed()
emit("reset_pedigree", newPop(founder, simParam=SP), SP)

founder = quickHaplo(5, 1, 19)
SP = newSP(founder)
emit(
  "permuted_initial_founders",
  newPop(founder[c(5, 2, 4, 1, 3)], simParam=SP),
  SP
)
SP$resetPed()
emit(
  "permuted_reset_founders",
  newPop(founder[c(3, 5, 1, 4, 2)], simParam=SP),
  SP
)

founder = quickHaplo(6, 2, c(14, 24))
SP = newSP(founder)
founders = newPop(founder, simParam=SP)
generation1 = randCross(founders, nCrosses=6, simParam=SP)
invisible(randCross(generation1, nCrosses=5, simParam=SP))
SP$resetPed(max(generation1@iid))
replacementGeneration2 = randCross(
  generation1,
  nCrosses=5,
  simParam=SP
)
emit("partial_pedigree_reset", replacementGeneration2, SP)

# Randomized multi-generation breeding chains exercise recombination tracking
# with multiple threads, interference settings, quadrivalent probabilities,
# current-genome mutation overlays, and samples spanning generations.
for(caseSeed in 1:36){
  set.seed(91000 + caseSeed)
  ploidy = c(2L, 4L, 6L, 8L)[(caseSeed - 1L) %% 4L + 1L]
  loci = sample(10:30, 2)
  founder = quickHaplo(8, 2, loci, ploidy=ploidy)
  SP = newSP(founder)
  SP$nThreads = if(caseSeed %% 2L) 1L else 2L
  SP$quadProb = if(ploidy == 2L) 0 else c(0, 0.5, 1)[
    (caseSeed - 1L) %% 3L + 1L
  ]
  SP$v = c(0.7, 1, 2.6)[(caseSeed - 1L) %/% 12L + 1L]
  SP$p = c(0, 0.2)[caseSeed %% 2L + 1L]
  recombinationRatio = c(0.5, 1, 3)[
    ((caseSeed - 1L) %/% 3L) %% 3L + 1L
  ]
  SP$setRecombRatio(recombinationRatio)
  founders = newPop(founder, simParam=SP)
  generation1 = randCross(founders, nCrosses=7, simParam=SP)
  generation2 = self(
    generation1[1:4],
    nProgeny=2,
    keepParents=caseSeed %% 2L == 0L,
    simParam=SP
  )
  generation3 = randCross(generation2, nCrosses=6, simParam=SP)
  mutatedGeneration3 = mutate(
    generation3,
    mutRate=0.05,
    simParam=SP
  )
  selectedChr = caseSeed %% 2L + 1L
  emitMulti(
    paste0("randomized_breeding_chain_", caseSeed, "_ploidy_", ploidy),
    list(
      founders[1],
      generation1[2],
      generation2[3],
      mutatedGeneration3[1:2]
    ),
    SP,
    selectedChr
  )
}

writeLines(as.character(caseCount), file.path(outputDir, "case-count.txt"))
cat(caseCount, "tree-sequence configurations written\n")
