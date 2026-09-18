context("vectorNProgeny")

# Issue #256 asked for vectorised family sizes. nProgeny was vectorised for
# the crossing functions; nDH in makeDH and nProgeny in reduceGenome are the
# remaining cases, and they are the awkward ones because the C++ behind them
# assumed a single value and wrote output at a fixed stride.
#
# The properties worth pinning down are the ones a stride bug breaks without
# failing loudly: the right number of progeny per parent, each progeny
# attributed to the correct parent, and a result that does not depend on how
# the work was divided between threads.

vecPop = function(nInd=5, nChr=2, segSites=20, ploidy=2L, seed=10301,
                  trackRec=FALSE, nThreads=1L){
  set.seed(seed)
  founderPop = quickHaplo(nInd=nInd, nChr=nChr, segSites=segSites,
                          ploidy=ploidy)
  SP = SimParam$new(founderPop)
  SP$nThreads = nThreads
  if(trackRec){
    SP$setTrackRec(TRUE)
  }
  SP$addTraitA(nQtlPerChr=5)
  SP$setVarE(h2=0.5)
  pop = newPop(founderPop, simParam=SP)
  return(list(pop=pop, SP=SP))
}

# Every progeny should name the parent it actually came from, in parent
# order, each parent appearing as many times as it was asked for
expect_parent_mapping = function(progeny, parent, counts, label=""){
  expect_equal(nInd(progeny), as.integer(sum(counts)), info=label)
  expect_equal(progeny@mother, rep(parent@id, times=counts), info=label)
  expect_equal(progeny@father, rep(parent@id, times=counts), info=label)
}

test_that("makeDH accepts one nDH per individual", {
  d = vecPop(nInd=5)
  nDH = c(3L, 1L, 4L, 1L, 2L)

  dh = makeDH(d$pop, nDH=nDH, keepParents=FALSE, simParam=d$SP)

  expect_parent_mapping(dh, d$pop, nDH)
  expect_equal(dh@ploidy, 2L)
  expect_equal(dh@nLoci, d$pop@nLoci)

  # Still doubled haploids: the two haplotypes of each are identical
  h = pullSegSiteHaplo(dh, simParam=d$SP)
  odd = seq(1, nrow(h), by=2)
  expect_equal(unname(h[odd,,drop=FALSE]), unname(h[odd+1,,drop=FALSE]))

  # And each DH's genotype is one the parent could have produced, so every
  # locus dosage is 0 or 2 and the parent is not homozygous for the other
  # allele there
  geno = pullSegSiteGeno(dh, simParam=d$SP)
  parGeno = pullSegSiteGeno(d$pop, simParam=d$SP)
  parentOf = rep(seq_len(nInd(d$pop)), times=nDH)
  expect_true(all(geno %in% c(0, 2)))
  for(i in seq_len(nrow(geno))){
    p = parGeno[parentOf[i],]
    # A DH carrying the 1 allele needs a parent that had one to give
    expect_true(all(p[geno[i,]==2] >= 1),
                info=paste("DH", i, "has an allele its parent lacks"))
    expect_true(all(p[geno[i,]==0] <= 1),
                info=paste("DH", i, "lacks an allele its parent was fixed for"))
  }
})

test_that("a vector of one repeated value matches the single value", {
  nDH = 3L
  a = vecPop(nInd=4, seed=10311)
  set.seed(555)
  fromScalar = makeDH(a$pop, nDH=nDH, simParam=a$SP)

  b = vecPop(nInd=4, seed=10311)
  set.seed(555)
  fromVector = makeDH(b$pop, nDH=rep(nDH, 4), simParam=b$SP)

  expect_equal(nInd(fromScalar), nInd(fromVector))
  expect_equal(unname(pullSegSiteHaplo(fromScalar, simParam=a$SP)),
               unname(pullSegSiteHaplo(fromVector, simParam=b$SP)))
})

test_that("an nDH of zero gives that individual no DH lines", {
  d = vecPop(nInd=4, seed=10321)
  nDH = c(2L, 0L, 0L, 3L)

  dh = makeDH(d$pop, nDH=nDH, keepParents=FALSE, simParam=d$SP)

  expect_equal(nInd(dh), 5L)
  expect_parent_mapping(dh, d$pop, nDH)
  # The skipped parents appear nowhere
  expect_false(any(dh@mother %in% d$pop@id[2:3]))
})

test_that("makeDH returns an empty population when no DH are requested", {
  d = vecPop(nInd=3, seed=10331)

  fromVector = makeDH(d$pop, nDH=c(0L, 0L, 0L), simParam=d$SP)
  expect_true(isPop(fromVector))
  expect_equal(nInd(fromVector), 0L)
  expect_equal(fromVector@nChr, d$pop@nChr)
  expect_equal(fromVector@ploidy, d$pop@ploidy)
  expect_equal(nrow(gv(fromVector)), 0L)
  expect_equal(ncol(gv(fromVector)), d$SP$nTraits)
  expect_true(isTRUE(validObject(fromVector, test=TRUE)))

  # A single zero does the same thing
  fromScalar = makeDH(d$pop, nDH=0, simParam=d$SP)
  expect_equal(nInd(fromScalar), 0L)

  # An empty result still merges with a real population
  some = makeDH(d$pop, nDH=1, simParam=d$SP)
  expect_equal(nInd(mergePops(list(fromVector, some))), nInd(some))
})

test_that("makeDH checks the length and the values of nDH", {
  d = vecPop(nInd=4, seed=10341)

  expect_error(makeDH(d$pop, nDH=c(1L,2L), simParam=d$SP),
               "Length of nDH")
  expect_error(makeDH(d$pop, nDH=rep(1L,5), simParam=d$SP),
               "Length of nDH")
  expect_error(makeDH(d$pop, nDH=c(1L,-1L,1L,1L), simParam=d$SP),
               "non-negative")
  expect_error(makeDH(d$pop, nDH=c(1L,NA_integer_,1L,1L), simParam=d$SP),
               "non-negative")
})

test_that("makeDH refuses a vector nDH for a MultiPop", {
  # Each population in a MultiPop has its own size, so one vector cannot
  # be right for all of them
  d = vecPop(nInd=4, seed=10351)
  multi = newMultiPop(d$pop, d$pop[1:2])

  expect_error(makeDH(multi, nDH=c(1L,2L,1L,1L), simParam=d$SP),
               "single value")

  # A single value is still fine and reaches every population
  out = makeDH(multi, nDH=2, simParam=d$SP)
  expect_equal(nInd(out@pops[[1]]), 8L)
  expect_equal(nInd(out@pops[[2]]), 4L)
})

test_that("a vector nDH gives the same answer whatever the thread count", {
  # The output slice for an individual is found from a running total, and
  # the point of computing it up front is that it cannot depend on how the
  # work was split
  skip_on_cran()
  nDH = c(4L, 1L, 0L, 3L, 2L, 5L)

  one = vecPop(nInd=6, segSites=50, seed=10361, nThreads=1L)
  set.seed(778)
  a = makeDH(one$pop, nDH=nDH, simParam=one$SP)

  two = vecPop(nInd=6, segSites=50, seed=10361, nThreads=4L)
  set.seed(778)
  b = makeDH(two$pop, nDH=nDH, simParam=two$SP)

  expect_equal(nInd(a), nInd(b))
  expect_equal(unname(pullSegSiteHaplo(a, simParam=one$SP)),
               unname(pullSegSiteHaplo(b, simParam=two$SP)))
})

test_that("recombination tracking lines up with a vector nDH", {
  # The recombination history is written at the same slice index as the
  # genotype, so a stride bug here shows up as a mismatch between the two
  skip_on_cran()
  d = vecPop(nInd=4, nChr=2, segSites=30, seed=10371, trackRec=TRUE)
  nDH = c(2L, 0L, 3L, 1L)

  dh = makeDH(d$pop, nDH=nDH, keepParents=FALSE, simParam=d$SP)
  expect_equal(nInd(dh), 6L)

  ibd = pullIbdHaplo(dh, simParam=d$SP)
  expect_equal(nrow(ibd), nInd(dh)*dh@ploidy)

  # A DH may only carry founder haplotypes belonging to its own parent
  parentIbd = pullIbdHaplo(d$pop, simParam=d$SP)
  originsOf = function(m, i, ploidy){
    return(unique(c(m[((i-1L)*ploidy + 1L):(i*ploidy), ])))
  }
  parentOf = rep(seq_len(nInd(d$pop)), times=nDH)
  for(i in seq_len(nInd(dh))){
    allowed = originsOf(parentIbd, parentOf[i], d$pop@ploidy)
    expect_true(all(originsOf(ibd, i, dh@ploidy) %in% allowed),
                info=paste("DH", i, "carries an unrelated origin"))
  }
})

test_that("reduceGenome accepts one nProgeny per individual", {
  d = vecPop(nInd=5, ploidy=4L, seed=10381)
  nProgeny = c(1L, 3L, 2L, 1L, 2L)

  red = reduceGenome(d$pop, nProgeny=nProgeny, keepParents=FALSE,
                     simParam=d$SP)

  expect_parent_mapping(red, d$pop, nProgeny)
  expect_equal(red@ploidy, 2L)
  expect_equal(red@nLoci, d$pop@nLoci)

  geno = pullSegSiteGeno(red, simParam=d$SP)
  expect_true(all(geno >= 0))
  expect_true(all(geno <= red@ploidy))
  expect_true(isTRUE(validObject(red, test=TRUE)))
})

test_that("reduceGenome handles zeros and an empty result", {
  d = vecPop(nInd=4, ploidy=4L, seed=10391)

  some = reduceGenome(d$pop, nProgeny=c(0L,2L,0L,1L), keepParents=FALSE,
                      simParam=d$SP)
  expect_equal(nInd(some), 3L)
  expect_parent_mapping(some, d$pop, c(0L,2L,0L,1L))

  none = reduceGenome(d$pop, nProgeny=0, simParam=d$SP)
  expect_true(isPop(none))
  expect_equal(nInd(none), 0L)
  # Reducing halves the ploidy, and the empty population has to say so
  expect_equal(none@ploidy, 2L)
  expect_true(isTRUE(validObject(none, test=TRUE)))
})

test_that("reduceGenome checks the length and the values of nProgeny", {
  d = vecPop(nInd=4, ploidy=4L, seed=10401)

  expect_error(reduceGenome(d$pop, nProgeny=c(1L,2L), simParam=d$SP),
               "Length of nProgeny")
  expect_error(reduceGenome(d$pop, nProgeny=c(1L,-2L,1L,1L), simParam=d$SP),
               "non-negative")
})

test_that("a vector nProgeny matches the single value and the thread count", {
  # reduceGenome divides its work over progeny rather than individuals, so
  # each output slice looks its parent up in a table built up front
  skip_on_cran()
  a = vecPop(nInd=4, ploidy=4L, segSites=40, seed=10411, nThreads=1L)
  set.seed(991)
  fromScalar = reduceGenome(a$pop, nProgeny=2, simParam=a$SP)

  b = vecPop(nInd=4, ploidy=4L, segSites=40, seed=10411, nThreads=1L)
  set.seed(991)
  fromVector = reduceGenome(b$pop, nProgeny=rep(2L,4), simParam=b$SP)

  expect_equal(unname(pullSegSiteHaplo(fromScalar, simParam=a$SP)),
               unname(pullSegSiteHaplo(fromVector, simParam=b$SP)))

  cc = vecPop(nInd=4, ploidy=4L, segSites=40, seed=10411, nThreads=4L)
  set.seed(991)
  manyThreads = reduceGenome(cc$pop, nProgeny=c(3L,0L,1L,4L), simParam=cc$SP)
  dd = vecPop(nInd=4, ploidy=4L, segSites=40, seed=10411, nThreads=1L)
  set.seed(991)
  oneThread = reduceGenome(dd$pop, nProgeny=c(3L,0L,1L,4L), simParam=dd$SP)

  expect_equal(nInd(manyThreads), 8L)
  expect_equal(unname(pullSegSiteHaplo(manyThreads, simParam=cc$SP)),
               unname(pullSegSiteHaplo(oneThread, simParam=dd$SP)))
})
