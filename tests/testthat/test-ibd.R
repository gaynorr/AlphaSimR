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

test_that("IBD agrees with genotypes in tetraploids with quadrivalents", {
  # Founder haplotypes are numbered one per row in order, so the allele at
  # each locus of a progeny haplotype must equal the allele of the founder
  # haplotype its IBD names. Quadrivalent gametes are drawn from all four
  # chromosomes, and a history that names the wrong chromosome shows up here
  # as a mismatch or as an error from pullIbdHaplo.
  for(quadProb in c(0, 1)){
    set.seed(9911)
    founderPop = quickHaplo(nInd=10, nChr=2, segSites=50, ploidy=4L)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$setTrackRec(TRUE)
    SP$quadProb = quadProb
    pop = newPop(founderPop, simParam=SP)
    founderHaplo = pullSegSiteHaplo(pop, simParam=SP)

    progeny = list(cross=randCross(pop, nCrosses=10, simParam=SP),
                   reduced=reduceGenome(pop, nProgeny=2, simParam=SP))
    progeny$cross2 = randCross(progeny$cross, nCrosses=10, simParam=SP)

    for(type in names(progeny)){
      ibd = pullIbdHaplo(progeny[[type]], simParam=SP)
      haplo = pullSegSiteHaplo(progeny[[type]], simParam=SP)
      expect_equal(dim(ibd), dim(haplo))
      expected = founderHaplo[cbind(c(ibd), rep(seq_len(ncol(ibd)),
                                                each=nrow(ibd)))]
      expect_equal(unname(c(haplo)), unname(expected),
                   info=paste(type, "with quadProb", quadProb))
    }
  }
})


# Checks every IBD label against the genotypes. Founder haplotypes are
# numbered one per row in order, so the allele at each locus of a haplotype
# must equal the allele of the founder haplotype its IBD names.
expectIbdMatchesGeno = function(pop, founderHaplo, SP, label){
  ibd = pullIbdHaplo(pop, simParam=SP)
  haplo = pullSegSiteHaplo(pop, simParam=SP)
  expect_equal(dim(ibd), dim(haplo), info=label)
  expected = founderHaplo[cbind(c(ibd), rep(seq_len(ncol(ibd)),
                                            each=nrow(ibd)))]
  expect_equal(unname(c(haplo)), unname(expected), info=label)
}

test_that("IBD agrees with genotypes through every kind of meiosis", {
  set.seed(9921)
  founderPop = quickHaplo(nInd=10, nChr=2, segSites=50)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  founderHaplo = pullSegSiteHaplo(pop, simParam=SP)

  f1 = randCross(pop, nCrosses=10, simParam=SP)
  progeny = list(f1=f1,
                 dh=makeDH(f1, nDH=2, simParam=SP),
                 self=self(f1, nProgeny=2, simParam=SP),
                 haploid=reduceGenome(f1, nProgeny=2, simParam=SP))
  # doubleGenome and mergeGenome build their histories in R rather than
  # in meiosis
  progeny$doubled = doubleGenome(progeny$haploid, simParam=SP)
  progeny$merged = mergeGenome(progeny$haploid[1:5], progeny$haploid[6:10],
                               crossPlan=cbind(1:5, 5:1), simParam=SP)
  progeny$f2 = randCross(c(progeny$dh, progeny$doubled, progeny$merged),
                         nCrosses=10, simParam=SP)

  for(type in names(progeny)){
    expectIbdMatchesGeno(progeny[[type]], founderHaplo, SP, type)
  }
})

test_that("IBD pulled in pieces matches IBD pulled at once", {
  # Haplotypes are cached as they are computed, so a later call reuses the
  # ancestors computed by an earlier one
  run = function(piecemeal){
    set.seed(9931)
    founderPop = quickHaplo(nInd=8, nChr=2, segSites=30)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$setTrackRec(TRUE)
    pop = newPop(founderPop, simParam=SP)
    f1 = randCross(pop, nCrosses=8, simParam=SP)
    f2 = randCross(f1, nCrosses=8, simParam=SP)
    f3 = randCross(f2, nCrosses=8, simParam=SP)
    if(piecemeal){
      pullIbdHaplo(f1[1:3], simParam=SP)
      pullIbdHaplo(f2, simParam=SP)
    }
    return(list(ibd=pullIbdHaplo(f3, simParam=SP), SP=SP))
  }
  once = run(FALSE)
  pieces = run(TRUE)
  expect_identical(pieces$ibd, once$ibd)

  # The cache is named by internal ID. The piecemeal run also holds the
  # individuals pulled earlier that are not ancestors of f3, and the two
  # agree wherever both hold an individual.
  onceHap = once$SP$haplotypes
  piecesHap = pieces$SP$haplotypes
  expect_identical(names(onceHap),
                   as.character(sort(as.integer(names(onceHap)))))
  expect_true(all(names(onceHap) %in% names(piecesHap)))
  expect_true(all(as.character(25:32) %in% names(onceHap)))
  expect_identical(piecesHap[names(onceHap)], onceHap)
})

test_that("recHist is returned as nested lists named by ID", {
  d = ibdPop(nInd=4, nChr=2, segSites=20)
  progeny = randCross(d$pop, nCrosses=3, simParam=d$SP)
  recHist = d$SP$recHist

  expect_identical(names(recHist), c(d$pop@id, progeny@id))
  # Founders list their founder haplotypes
  expect_identical(recHist[[d$pop@id[2]]], 3:4)
  # Progeny hold a two column matrix per copy of each chromosome
  for(id in progeny@id){
    expect_equal(length(recHist[[id]]), 2L)
    for(chr in recHist[[id]]){
      expect_equal(length(chr), 2L)
      for(m in chr){
        expect_true(is.integer(m) && is.matrix(m) && ncol(m)==2L)
        expect_true(all(m[,1] %in% 1:2))
        expect_equal(m[1,2], 1L)
      }
    }
  }
})

test_that("packed recombination records round trip", {
  x = list(1:2,
           list(list(matrix(c(1L,1L), 1, 2), matrix(c(2L,1L,1L,7L), 2, 2)),
                list(matrix(c(2L,1L), 1, 2), matrix(c(1L,2L,1L,1L,4L,9L),
                                                    3, 2))))
  packed = AlphaSimR:::packHist(x)
  expect_identical(packed[[1]], 1:2)
  expect_true(is.integer(packed[[2]]) && packed[[2]][1] == -2L)
  unpacked = AlphaSimR:::unpackHist(packed)
  expect_identical(unpacked[[1]], 1:2)
  for(chr in 1:2){
    for(k in 1:2){
      expect_identical(unpacked[[2]][[chr]][[k]], x[[2]][[chr]][[k]])
    }
  }
})

test_that("resetPed discards the cached IBD of removed individuals", {
  set.seed(9941)
  founderPop = quickHaplo(nInd=6, nChr=1, segSites=30)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  founderHaplo = pullSegSiteHaplo(pop, simParam=SP)
  f1 = randCross(pop, nCrosses=6, simParam=SP)
  pullIbdHaplo(f1, simParam=SP)

  # The new individuals take the iids of the ones removed
  SP$resetPed(6L)
  f1b = randCross(pop, nCrosses=6, simParam=SP)
  expect_identical(f1b@iid, f1@iid)
  expectIbdMatchesGeno(f1b, founderHaplo, SP, "after resetPed")
})

test_that("founder haplotypes can be supplied with setFounderHap", {
  set.seed(9951)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=20)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  hapMap = lapply(1:4, function(i)
    list(list(matrix(c(100L+i, 1L), 1, 2), matrix(c(200L+i, 1L), 1, 2))))
  names(hapMap) = as.character(4:1)
  expect_error(SP$setFounderHap(setNames(hapMap, c("1","2","3","9"))),
               "internal IDs")
  SP$setFounderHap(hapMap)

  ibd = pullIbdHaplo(pop, simParam=SP)
  expect_equal(unname(ibd[,1]), c(104L, 204L, 103L, 203L,
                                  102L, 202L, 101L, 201L))
  progeny = randCross(pop, nCrosses=6, simParam=SP)
  expect_true(all(pullIbdHaplo(progeny, simParam=SP) %in% c(ibd)))
})
