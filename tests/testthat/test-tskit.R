context("tskit")

test_that("recorded diploid ancestry converts and writes", {
  set.seed(101)
  founderPop = quickHaplo(nInd=4, nChr=2, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  founders = newPop(founderPop, simParam=SP)
  progeny = randCross(founders, nCrosses=6, simParam=SP)

  ts = asTreeSequence(progeny, simplify=FALSE, simParam=SP)
  expect_s3_class(ts, "AlphaSimRTreeSequence")
  expect_named(ts, c("chr1", "chr2"))
  expect_equal(as.numeric(ts[[1]]$num_samples()), 12)
  expect_equal(as.numeric(ts[[1]]$num_individuals()), 10)
  expect_equal(as.numeric(ts[[1]]$num_nodes()), 28 + 12)
  expect_equal(as.numeric(ts[[1]]$num_sites()), 20)
  expect_equal(ts[[1]]$sequence_length(), 20)

  prefix = tempfile(fileext=".trees")
  paths = writeTreeSequence(ts, prefix)
  expect_length(paths, 2)
  expect_true(all(file.exists(paths)))
  expect_match(paths, "_chr[12]\\.trees$")
  loaded = RcppTskit::ts_load(paths[1])
  expect_equal(as.numeric(loaded$num_samples()), 12)
  expect_equal(as.numeric(loaded$num_sites()), 20)
})

test_that("simplification retains requested samples and variants", {
  set.seed(102)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=25)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  generation1 = randCross(pop, nCrosses=8, simParam=SP)
  generation2 = randCross(generation1, nCrosses=5, simParam=SP)

  full = asTreeSequence(generation2, simplify=FALSE, simParam=SP)
  small = asTreeSequence(generation2, simplify=TRUE, simParam=SP)
  expect_equal(as.numeric(small[[1]]$num_samples()), 10)
  expect_equal(as.numeric(small[[1]]$num_sites()), 25)
  expect_lte(as.numeric(small[[1]]$num_nodes()),
             as.numeric(full[[1]]$num_nodes()))
  expect_lte(as.numeric(small[[1]]$num_edges()),
             as.numeric(full[[1]]$num_edges()))
})

test_that("chromosome requests retain their requested order", {
  set.seed(109)
  founderPop = quickHaplo(nInd=5, nChr=3, segSites=c(7, 11, 19))
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  founders = newPop(founderPop, simParam=SP)
  progeny = randCross(founders, nCrosses=4, simParam=SP)

  reordered = asTreeSequence(progeny, chr=c(3, 1), simParam=SP)
  expect_named(reordered, c("chr3", "chr1"))
  expect_equal(attr(reordered, "chromosome"), c(3L, 1L))
  expect_equal(
    unname(vapply(reordered, function(x) x$sequence_length(), numeric(1))),
    c(19, 7)
  )
  expect_equal(
    unname(vapply(
      reordered,
      function(x) as.numeric(x$num_sites()),
      numeric(1)
    )),
    c(19, 7)
  )
})

test_that("doubled haploid ancestry is supported", {
  set.seed(103)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=30)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  dh = makeDH(pop, nDH=2, simParam=SP)

  ts = asTreeSequence(dh, simParam=SP)
  expect_equal(as.numeric(ts[[1]]$num_samples()), dh@nInd * dh@ploidy)
  expect_equal(as.numeric(ts[[1]]$num_sites()), dh@nLoci[1])
})

test_that("inbred founders preserve shared haplotype origins", {
  set.seed(106)
  founderPop = quickHaplo(
    nInd=4,
    nChr=2,
    segSites=20,
    inbred=TRUE
  )
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  founders = newPop(founderPop, simParam=SP)

  founderTs = asTreeSequence(founders, simplify=FALSE, simParam=SP)
  expect_equal(as.numeric(founderTs[[1]]$num_samples()), 8)
  expect_equal(as.numeric(founderTs[[1]]$num_nodes()), 12 + 8)

  progeny = hybridCross(
    founders[1:2],
    founders[3:4],
    simParam=SP
  )
  progenyTs = asTreeSequence(progeny, simParam=SP)
  expect_equal(as.numeric(progenyTs[[1]]$num_samples()), 8)
  expect_equal(as.numeric(progenyTs[[1]]$num_sites()), 20)
})

test_that("genome edits and mutations are represented on samples", {
  set.seed(104)
  founderPop = quickHaplo(nInd=3, nChr=1, segSites=15)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  originalAllele = pullSegSiteHaplo(pop, simParam=SP)[1,1]
  edited = editGenome(pop, ind=1, chr=1, segSites=1,
                      allele=1L-originalAllele,
                      simParam=SP)

  originalTs = asTreeSequence(pop, simParam=SP)
  editedTs = asTreeSequence(edited, simParam=SP)
  expect_gt(as.numeric(editedTs[[1]]$num_mutations()),
            as.numeric(originalTs[[1]]$num_mutations()))

  mutated = mutate(pop, mutRate=1, simParam=SP)
  mutatedTs = asTreeSequence(mutated, simParam=SP)
  expect_equal(
    as.numeric(mutatedTs[[1]]$num_sites()),
    mutated@nLoci[1]
  )
  expect_gt(as.numeric(mutatedTs[[1]]$num_mutations()), 0)

  ancestry = asTreeSequence(edited, includeVariants=FALSE, simParam=SP)
  expect_equal(as.numeric(ancestry[[1]]$num_sites()), 0)
})

test_that("input and overwrite safeguards are enforced", {
  founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  pop = newPop(founderPop, simParam=SP)
  expect_error(asTreeSequence(pop, simParam=SP), "setTrackRec")

  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  expect_error(asTreeSequence(pop[c(1,1)], simParam=SP), "unique")
  missingId = pop
  missingId@id[1] = NA_character_
  expect_error(asTreeSequence(missingId, simParam=SP), "non-missing")
  invalidId = pop
  invalidId@id[1] = rawToChar(as.raw(0xff))
  expect_error(asTreeSequence(invalidId, simParam=SP), "valid UTF-8")
  duplicatePublicId = pop
  duplicatePublicId@id[] = "shared"
  expect_s3_class(
    asTreeSequence(duplicatePublicId, simParam=SP),
    "AlphaSimRTreeSequence"
  )
  otherFounderPop = quickHaplo(nInd=2, nChr=1, segSites=11)
  otherSP = SimParam$new(otherFounderPop)
  otherSP$nThreads = 1L
  otherSP$setTrackRec(TRUE)
  invisible(newPop(otherFounderPop, simParam=otherSP))
  expect_error(
    asTreeSequence(pop, includeVariants=FALSE, simParam=otherSP),
    "chromosome/locus structures"
  )
  expect_error(asTreeSequence(pop, chr=integer(), simParam=SP), "invalid")
  expect_error(asTreeSequence(pop, chr=1.5, simParam=SP), "whole")
  expect_error(asTreeSequence(pop, simplify=NA, simParam=SP), "simplify")
  expect_error(asTreeSequence(pop, includeVariants=1, simParam=SP),
               "includeVariants")
  ts = asTreeSequence(pop, simParam=SP)
  path = tempfile(fileext=".trees")
  writeTreeSequence(ts, path)
  expect_error(writeTreeSequence(ts, path), "already exists")
  expect_silent(writeTreeSequence(ts, path, overwrite=TRUE))
})

test_that("polyploid histories and additional ancestry-only founders work", {
  set.seed(105)
  founderPop = quickHaplo(nInd=4, nChr=2, segSites=24)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  tetraploid = doubleGenome(pop, keepParents=FALSE, simParam=SP)
  SP$quadProb = 1
  progeny = randCross(tetraploid, nCrosses=4, simParam=SP)

  ts = asTreeSequence(progeny, chr=2, simParam=SP)
  expect_named(ts, "chr2")
  expect_equal(as.numeric(ts[[1]]$num_samples()),
               progeny@nInd * progeny@ploidy)
  expect_equal(ts[[1]]$sequence_length(), progeny@nLoci[2])

  additionalFounders = newPop(founderPop, simParam=SP)
  variants = asTreeSequence(additionalFounders, simParam=SP)
  expect_equal(as.numeric(variants[[1]]$num_samples()),
               additionalFounders@nInd * additionalFounders@ploidy)
  expect_equal(as.numeric(variants[[1]]$num_sites()),
               additionalFounders@nLoci[1])

  changedFounder = editGenome(
    additionalFounders,
    ind=1,
    chr=1,
    segSites=1,
    allele=1L-pullSegSiteHaplo(
      additionalFounders,
      simParam=SP
    )[1,1],
    simParam=SP
  )
  expect_s3_class(
    asTreeSequence(changedFounder, simParam=SP),
    "AlphaSimRTreeSequence"
  )
})

test_that("MultiPop samples can span generations and ploidies", {
  set.seed(107)
  founderPop = quickHaplo(nInd=5, nChr=2, segSites=18)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  founders = newPop(founderPop, simParam=SP)
  progeny = randCross(founders, nCrosses=4, simParam=SP)
  tetraploid = doubleGenome(progeny[1:2], simParam=SP)
  samples = newMultiPop(progeny[3:4], tetraploid)

  ts = asTreeSequence(samples, simplify=FALSE, simParam=SP)
  expect_named(ts, c("chr1", "chr2"))
  expect_equal(as.numeric(ts[[1]]$num_samples()), 12)
  expect_equal(as.numeric(ts[[1]]$num_sites()), 18)

  nested = newMultiPop(newMultiPop(progeny[1], tetraploid[1]), progeny[2])
  nestedTs = asTreeSequence(nested, simParam=SP)
  expect_equal(as.numeric(nestedTs[[1]]$num_samples()), 8)
})

test_that("HybridPop reports its representation limitation", {
  founderPop = quickHaplo(nInd=3, nChr=1, segSites=12, inbred=TRUE)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  founders = newPop(founderPop, simParam=SP)
  hybrid = hybridCross(
    founders,
    founders[1],
    returnHybridPop=TRUE,
    simParam=SP
  )

  expect_error(
    asTreeSequence(hybrid, simParam=SP),
    "does not store haplotypes or recombination history"
  )
})

test_that("pedigree resets do not reuse stale IBD ancestry", {
  set.seed(108)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=14)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  founders = newPop(founderPop, simParam=SP)
  invisible(randCross(founders, nCrosses=4, simParam=SP))

  SP$resetPed()
  replacementFounders = newPop(founderPop, simParam=SP)
  ts = asTreeSequence(replacementFounders, simParam=SP)
  expect_equal(as.numeric(ts[[1]]$num_samples()), 8)
  expect_equal(as.numeric(ts[[1]]$num_sites()), 14)
})

test_that("randomized ancestry exports remain structurally valid", {
  for(seed in 201:206){
    set.seed(seed)
    nLoci = sample(8:20, 2)
    founderPop = quickHaplo(
      nInd=6,
      nChr=2,
      segSites=nLoci,
      ploidy=if(seed %% 2L) 2L else 4L
    )
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$quadProb = if(seed %% 2L) 0 else (seed %% 3L) / 2
    SP$setTrackRec(TRUE)
    founders = newPop(founderPop, simParam=SP)
    generation1 = randCross(founders, nCrosses=5, simParam=SP)
    generation2 = randCross(generation1, nCrosses=4, simParam=SP)
    samples = generation2[sample.int(generation2@nInd)]

    for(simplify in c(FALSE, TRUE)){
      ts = asTreeSequence(samples, simplify=simplify, simParam=SP)
      expect_equal(as.numeric(ts[[1]]$num_samples()),
                   samples@nInd * samples@ploidy)
      expect_equal(as.numeric(ts[[2]]$num_sites()), nLoci[2])
      expect_equal(ts[[1]]$sequence_length(), nLoci[1])
      expect_equal(attr(ts, "sample_iid"), samples@iid)
      expect_identical(attr(ts, "coordinate_system"), "locus_index")
    }
  }
})

test_that("variant encoding handles inferred, recurrent, and back mutations", {
  encoding = AlphaSimR:::resolveVariantEncodingCpp(
    originRows=matrix(c(1L, 1L, 2L, 2L), nrow=2),
    currentHaplotypes=matrix(c(0L, 1L, 0L, 1L), nrow=2),
    originHaplotypes=matrix(c(0L, 1L, 0L, 1L), nrow=2),
    knownOrigins=c(FALSE, TRUE)
  )

  expect_equal(encoding$originHaplotypes[1,], c(0L, 0L))
  expect_equal(encoding$sampleAlleleOverrides[,1], c(-1L, 1L))
  expect_equal(encoding$sampleAlleleOverrides[,2], c(0L, -1L))
})
