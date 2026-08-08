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
  expect_equal(as.numeric(ts[[1]]$num_nodes()), 28)
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
  expect_equal(as.numeric(founderTs[[1]]$num_nodes()), 12)

  progeny = hybridCross(
    founders[1:2],
    founders[3:4],
    simParam=SP
  )
  progenyTs = asTreeSequence(progeny, simParam=SP)
  expect_equal(as.numeric(progenyTs[[1]]$num_samples()), 8)
  expect_equal(as.numeric(progenyTs[[1]]$num_sites()), 20)
})

test_that("unrecorded genome changes are explicit", {
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

  expect_error(
    asTreeSequence(edited, simParam=SP),
    "changes not represented by recombination history"
  )
  expect_s3_class(
    asTreeSequence(edited, includeVariants=FALSE, simParam=SP),
    "AlphaSimRTreeSequence"
  )
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
  expect_error(asTreeSequence(pop, chr=integer(), simParam=SP), "invalid")
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
  expect_error(
    asTreeSequence(additionalFounders, simParam=SP),
    "only recorded founders"
  )
  ancestry = asTreeSequence(
    additionalFounders,
    includeVariants=FALSE,
    simParam=SP
  )
  expect_equal(as.numeric(ancestry[[1]]$num_samples()),
               additionalFounders@nInd * additionalFounders@ploidy)
  expect_equal(as.numeric(ancestry[[1]]$num_sites()), 0)
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
