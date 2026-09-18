context("polyploids")

# reduceGenome, doubleGenome and mergeGenome change the ploidy of a
# population. There was no test file for any of them, and a bug in
# reduceGenome's choice of centromeres has already shipped once.
#
# These check the properties that hold whatever the ploidy: the result has the
# ploidy it should, carries the same loci, and its dosages stay inside the
# range the new ploidy allows.

polyPop = function(ploidy, nInd=6, nChr=2, segSites=12, nQtl=6, seed=9001){
  set.seed(seed)
  founderPop = quickHaplo(nInd=nInd, nChr=nChr, segSites=segSites,
                          ploidy=ploidy)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=nQtl)
  SP$setVarE(h2=0.5)
  pop = newPop(founderPop, simParam=SP)
  return(list(pop=pop, SP=SP))
}

# Dosages must lie between zero and the ploidy, and the number of columns must
# match the loci the population claims to have
expect_sane_geno = function(pop, SP, label=""){
  geno = pullSegSiteGeno(pop, simParam=SP)
  expect_equal(nrow(geno), nInd(pop), info=label)
  expect_equal(ncol(geno), sum(pop@nLoci), info=label)
  expect_true(all(geno >= 0), info=label)
  expect_true(all(geno <= pop@ploidy), info=label)
  expect_true(all(geno == round(geno)), info=label)
}

test_that("doubleGenome doubles the ploidy and keeps every dosage doubled", {
  for(ploidy in c(2L, 4L)){
    d = polyPop(ploidy)
    before = pullSegSiteGeno(d$pop, simParam=d$SP)

    dbl = doubleGenome(d$pop, simParam=d$SP)

    expect_equal(dbl@ploidy, 2L*ploidy)
    expect_equal(nInd(dbl), nInd(d$pop))
    expect_equal(dbl@nLoci, d$pop@nLoci)
    expect_sane_geno(dbl, d$SP, label=paste("doubled from", ploidy))

    # Every haplotype is copied, so every dosage is exactly doubled. This is
    # the whole of what doubling means and it holds with no randomness.
    after = pullSegSiteGeno(dbl, simParam=d$SP)
    expect_equal(unname(after), unname(2L*before))

    # A doubled individual is fully homozygous in the sense that its
    # haplotypes come in identical pairs
    h = pullSegSiteHaplo(dbl, simParam=d$SP)
    odd = seq(1, nrow(h), by=2)
    expect_equal(unname(h[odd,,drop=FALSE]), unname(h[odd+1,,drop=FALSE]))
  }
})

test_that("reduceGenome halves the ploidy", {
  # Meiosis at ploidy 8 is expensive relative to what it adds over ploidy 4
  skip_on_cran()
  for(ploidy in c(4L, 8L)){
    d = polyPop(ploidy)
    red = reduceGenome(d$pop, simParam=d$SP)

    expect_equal(red@ploidy, ploidy %/% 2L)
    expect_equal(nInd(red), nInd(d$pop))
    expect_equal(red@nLoci, d$pop@nLoci)
    expect_sane_geno(red, d$SP, label=paste("reduced from", ploidy))
  }
})

test_that("reduceGenome makes the requested number of progeny", {
  d = polyPop(4L, nInd=5)
  red = reduceGenome(d$pop, nProgeny=3, simParam=d$SP)
  expect_equal(nInd(red), 15L)
  expect_equal(red@ploidy, 2L)
})

test_that("reduceGenome refuses an odd ploidy", {
  d = polyPop(4L)
  odd = d$pop
  odd@ploidy = 3L
  expect_error(reduceGenome(odd, simParam=d$SP), "odd ploidy")
})

test_that("reduceGenome without recombination copies whole haplotypes", {
  # With recombination switched off every gamete is one of the parent's
  # haplotypes taken whole, so each locus of the gamete must match the parent
  # at some single haplotype
  d = polyPop(4L, nInd=4, nChr=1, segSites=20)
  red = reduceGenome(d$pop, simRecomb=FALSE, simParam=d$SP)
  expect_equal(red@ploidy, 2L)

  parent = pullSegSiteHaplo(d$pop, simParam=d$SP)
  gamete = pullSegSiteHaplo(red, simParam=d$SP)
  expect_equal(nrow(gamete), nInd(red) * red@ploidy)

  for(i in seq_len(nInd(red))){
    candidates = parent[((i-1L)*4L+1L):(i*4L),,drop=FALSE]
    for(k in seq_len(red@ploidy)){
      g = gamete[(i-1L)*red@ploidy + k,]
      matches = apply(candidates, 1, function(h) all(h == g))
      expect_true(any(matches),
                  info=paste("gamete", i, "haplotype", k,
                             "is not a whole parent haplotype"))
    }
  }
})

test_that("mergeGenome adds the two ploidies together", {
  d2 = polyPop(2L, nInd=4, seed=9101)

  # Same ploidy on both sides
  crossPlan = cbind(1:4, 4:1)
  merged = mergeGenome(d2$pop, d2$pop, crossPlan=crossPlan, simParam=d2$SP)
  expect_equal(merged@ploidy, 4L)
  expect_equal(nInd(merged), 4L)
  expect_sane_geno(merged, d2$SP, label="2 + 2")

  # The merged dosage is the sum of the two parents' dosages, because no
  # meiosis happens: whole genomes are placed side by side
  parents = pullSegSiteGeno(d2$pop, simParam=d2$SP)
  got = pullSegSiteGeno(merged, simParam=d2$SP)
  expect_equal(unname(got),
               unname(parents[crossPlan[,1],,drop=FALSE] +
                      parents[crossPlan[,2],,drop=FALSE]))

  # Different ploidies on the two sides, both descended from one SimParam
  d4 = polyPop(4L, nInd=4, seed=9103)
  reduced = reduceGenome(d4$pop, simParam=d4$SP)
  expect_equal(reduced@ploidy, 2L)
  mixed = mergeGenome(reduced, d4$pop, crossPlan=crossPlan, simParam=d4$SP)
  expect_equal(mixed@ploidy, 6L)
  expect_equal(nInd(mixed), 4L)
  expect_sane_geno(mixed, d4$SP, label="2 + 4")
})

test_that("mergeGenome refuses a crossPlan that points outside the parents", {
  d = polyPop(2L, nInd=4)
  expect_error(mergeGenome(d$pop, d$pop, crossPlan=cbind(1:4, c(1,2,3,99)),
                           simParam=d$SP),
               "Invalid crossPlan")
  expect_error(mergeGenome(d$pop, d$pop, crossPlan=cbind(1:4, c(0,2,3,4)),
                           simParam=d$SP),
               "Invalid crossPlan")
})

test_that("mergeGenome can match parents by id", {
  d = polyPop(2L, nInd=4)
  byIndex = mergeGenome(d$pop, d$pop, crossPlan=cbind(1:2, 3:4),
                        simParam=d$SP)
  byId = mergeGenome(d$pop, d$pop,
                     crossPlan=cbind(d$pop@id[1:2], d$pop@id[3:4]),
                     simParam=d$SP)
  expect_equal(unname(pullSegSiteGeno(byIndex, simParam=d$SP)),
               unname(pullSegSiteGeno(byId, simParam=d$SP)))

  expect_error(mergeGenome(d$pop, d$pop,
                           crossPlan=cbind(c("notAnId","x"), d$pop@id[1:2]),
                           simParam=d$SP),
               "Failed to match supplied IDs")
})

test_that("doubling then reducing returns to the starting ploidy", {
  # Doubles to ploidy 8 before reducing back
  skip_on_cran()
  for(ploidy in c(2L, 4L)){
    d = polyPop(ploidy)
    roundTrip = reduceGenome(doubleGenome(d$pop, simParam=d$SP),
                             simParam=d$SP)
    expect_equal(roundTrip@ploidy, ploidy)
    expect_sane_geno(roundTrip, d$SP, label=paste("round trip at", ploidy))
  }
})

test_that("changing ploidy leaves a usable population", {
  # Genetic values have to be recalculated for the new ploidy, and the result
  # has to behave like any other population
  d = polyPop(4L, nInd=8)
  red = reduceGenome(d$pop, simParam=d$SP)

  expect_true(isTRUE(validObject(red, test=TRUE)))
  expect_equal(nrow(gv(red)), nInd(red))
  expect_true(all(is.finite(c(gv(red)))))

  phenotyped = setPheno(red, varE=1, simParam=d$SP)
  expect_true(all(is.finite(c(pheno(phenotyped)))))

  selected = selectInd(phenotyped, nInd=2, simParam=d$SP)
  expect_equal(nInd(selected), 2L)
  expect_equal(selected@ploidy, red@ploidy)
})

test_that("reduceGenome follows the map it is told to use", {
  # A hundred sites and two sex specific maps, for one qualitative comparison
  skip_on_cran()
  # useFemale chooses which of the sex specific maps supplies both the
  # genetic positions and the centromeres. With a recombination ratio set,
  # the two maps differ, so the two gametes should differ as well.
  set.seed(9201)
  founderPop = quickHaplo(nInd=20, nChr=1, segSites=100, ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setRecombRatio(4)
  SP$addTraitA(nQtlPerChr=10)
  pop = newPop(founderPop, simParam=SP)

  expect_false(isTRUE(all.equal(SP$femaleMap, SP$maleMap)))

  set.seed(9202)
  fem = reduceGenome(pop, useFemale=TRUE, simParam=SP)
  set.seed(9202)
  mal = reduceGenome(pop, useFemale=FALSE, simParam=SP)

  expect_equal(fem@ploidy, 2L)
  expect_equal(mal@ploidy, 2L)
  expect_false(isTRUE(all.equal(pullSegSiteHaplo(fem, simParam=SP),
                                pullSegSiteHaplo(mal, simParam=SP))))
})
