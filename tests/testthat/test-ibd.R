context("ibd")

# Recombination tracking was never switched on anywhere in the test suite, so
# setTrackRec and pullIbdHaplo had no coverage at all.
#
# IBD haplotypes name, for every locus of every haplotype, which founder
# haplotype it descends from. That gives properties that hold exactly and need
# no reference implementation: founders descend from themselves, progeny carry
# only their parents' founder haplotypes, and the number of switches along a
# chromosome is the number of crossovers that happened.

ibdPop = function(nInd=10, nChr=2, segSites=30, seed=9701){
  set.seed(seed)
  founderPop = quickHaplo(nInd=nInd, nChr=nChr, segSites=segSites)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  SP$addTraitA(nQtlPerChr=5)
  SP$addSnpChip(nSnpPerChr=5)
  SP$setVarE(h2=0.5)
  pop = newPop(founderPop, simParam=SP)
  return(list(pop=pop, SP=SP))
}

test_that("pullIbdHaplo needs recombination tracking", {
  set.seed(9702)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(5)
  pop = newPop(founderPop, simParam=SP)

  expect_false(SP$isTrackRec)
  expect_error(pullIbdHaplo(pop, simParam=SP), "trackRec")
})

test_that("setTrackRec also turns on pedigree tracking", {
  set.seed(9703)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(FALSE)
  SP$setTrackRec(TRUE)

  expect_true(SP$isTrackRec)
  expect_true(SP$isTrackPed)
})

test_that("founders descend from their own haplotypes", {
  d = ibdPop(nInd=6, nChr=2, segSites=20)
  ibd = pullIbdHaplo(d$pop, simParam=d$SP)

  expect_equal(nrow(ibd), nInd(d$pop) * d$pop@ploidy)
  expect_equal(ncol(ibd), sum(d$pop@nLoci))

  # Each founder haplotype is one founder haplotype the whole way along, and
  # they are numbered one per row in order
  for(i in seq_len(nrow(ibd))){
    expect_equal(length(unique(ibd[i,])), 1L,
                 info=paste("founder haplotype", i, "is not of one origin"))
  }
  # Distinct founders have distinct origins
  expect_equal(length(unique(ibd[,1])), nrow(ibd))

  # Every value names a haplotype that exists
  expect_true(all(ibd >= 1))
  expect_true(all(ibd <= nInd(d$pop) * d$pop@ploidy))
})

test_that("progeny carry only their parents' founder haplotypes", {
  d = ibdPop(nInd=10, nChr=2, segSites=30)
  progeny = randCross(d$pop, nCrosses=8, simParam=d$SP)

  parentIbd = pullIbdHaplo(d$pop, simParam=d$SP)
  childIbd = pullIbdHaplo(progeny, simParam=d$SP)

  expect_equal(nrow(childIbd), nInd(progeny) * progeny@ploidy)
  expect_equal(ncol(childIbd), ncol(parentIbd))

  # The rows of an IBD matrix follow the individuals in order, two rows each,
  # which the row names state, so a parent's origins can be read off directly
  originsOf = function(ibd, i){
    return(unique(c(ibd[((i-1L)*2L + 1L):(i*2L), ])))
  }
  parentOf = function(id) which(d$pop@id == id)

  for(i in seq_len(nInd(progeny))){
    allowed = c(originsOf(parentIbd, parentOf(progeny@mother[i])),
                originsOf(parentIbd, parentOf(progeny@father[i])))
    expect_true(all(originsOf(childIbd, i) %in% allowed),
                info=paste("individual", i, "carries an unrelated origin"))
  }
})

test_that("a doubled haploid is identical in both its haplotypes", {
  d = ibdPop(nInd=6, nChr=2, segSites=30)
  f1 = randCross(d$pop, nCrosses=4, simParam=d$SP)
  dh = makeDH(f1, nDH=1, simParam=d$SP)

  ibd = pullIbdHaplo(dh, simParam=d$SP)
  odd = seq(1, nrow(ibd), by=2)
  expect_equal(unname(ibd[odd,,drop=FALSE]), unname(ibd[odd+1,,drop=FALSE]))
})

test_that("a selfed individual of one founder stays that founder", {
  # A single founder selfed to itself can only pass on its own two
  # haplotypes, however many crossovers happen
  d = ibdPop(nInd=4, nChr=2, segSites=30)
  one = d$pop[1]
  allowed = unique(c(pullIbdHaplo(one, simParam=d$SP)))
  expect_equal(length(allowed), 2L)

  kids = self(one, nProgeny=5, simParam=d$SP)
  ibd = pullIbdHaplo(kids, simParam=d$SP)
  expect_true(all(ibd %in% allowed))
})

test_that("switches along a chromosome count the crossovers", {
  # Two hundred sites across forty gametes, to count crossovers
  skip_on_cran()
  # Without recombination a gamete is one founder haplotype from end to end,
  # so a chromosome of a doubled haploid made without crossovers has no
  # switches. With recombination there should be some.
  d = ibdPop(nInd=20, nChr=1, segSites=200, seed=9801)
  progeny = randCross(d$pop, nCrosses=20, simParam=d$SP)
  ibd = pullIbdHaplo(progeny, simParam=d$SP)

  switches = apply(ibd, 1, function(r) sum(r[-1] != r[-length(r)]))

  # Over a whole Morgan and forty gametes, at least one crossover is certain
  # enough to assert
  expect_gt(sum(switches), 0)

  # And no gamete should look like noise. A one Morgan chromosome averages
  # one crossover, so a gamete with dozens of switches would mean the IBD
  # positions are not in map order.
  expect_lt(max(switches), 20)
})

test_that("pullIbdHaplo can be narrowed to chromosomes and to a chip", {
  d = ibdPop(nInd=6, nChr=3, segSites=20)
  progeny = randCross(d$pop, nCrosses=4, simParam=d$SP)

  full = pullIbdHaplo(progeny, simParam=d$SP)
  expect_equal(ncol(full), sum(progeny@nLoci))

  one = pullIbdHaplo(progeny, chr=2, simParam=d$SP)
  expect_equal(ncol(one), progeny@nLoci[2])

  two = pullIbdHaplo(progeny, chr=c(1,3), simParam=d$SP)
  expect_equal(ncol(two), sum(progeny@nLoci[c(1,3)]))

  snpMap = getSnpMap(snpChip=1, simParam=d$SP)
  chip = pullIbdHaplo(progeny, snpChip=1, simParam=d$SP)
  expect_equal(ncol(chip), nrow(snpMap))

  # The chip columns are a subset of the full set, taken at the chip's loci
  expect_true(all(colnames(chip) %in% colnames(full)))
})

test_that("IBD haplotypes are named after the individuals they belong to", {
  d = ibdPop(nInd=5, nChr=1, segSites=20)
  ibd = pullIbdHaplo(d$pop, simParam=d$SP)

  expect_equal(rownames(ibd),
               paste(rep(d$pop@id, each=2), rep(1:2, 5), sep="_"))
  expect_equal(unname(colnames(ibd)), names(d$SP$genMap[[1]]))
})

test_that("tracking recombination does not change the genotypes produced", {
  # The IBD bookkeeping rides alongside the simulation and must not alter it
  set.seed(9901)
  founderPop = quickHaplo(nInd=10, nChr=2, segSites=30)

  makeSP = function(track){
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$setTrackRec(track)
    SP$addTraitA(nQtlPerChr=5)
    SP$setVarE(h2=0.5)
    return(SP)
  }

  SPoff = makeSP(FALSE)
  set.seed(9902)
  popOff = randCross(newPop(founderPop, simParam=SPoff), nCrosses=6,
                     simParam=SPoff)

  SPon = makeSP(TRUE)
  set.seed(9902)
  popOn = randCross(newPop(founderPop, simParam=SPon), nCrosses=6,
                    simParam=SPon)

  expect_equal(unname(pullSegSiteHaplo(popOff, simParam=SPoff)),
               unname(pullSegSiteHaplo(popOn, simParam=SPon)))
})
