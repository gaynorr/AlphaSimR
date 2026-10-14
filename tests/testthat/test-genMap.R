context("genMap")

# Genetic map positions are absolute: they are measured from the start of the
# chromosome, which is position 0, and the first locus need not sit there.
# Meiosis samples chiasmata from position 0, so some crossovers fall before
# the first locus. Those change which parental chromosome a gamete starts on
# without adding a crossover between loci, and a recombination history still
# begins at site 1 with strictly increasing sites after it.
#
# The structural check used below holds exactly: when every locus sits at the
# same position, all chiasmata fall before the first locus, so each gamete
# has to be one parental haplotype, never a mixture of two.

# A MapPop with every locus at the same position, away from position 0
stackedMapPop = function(ploidy, nInd=6, nLoci=40, pos=2){
  H = matrix(sample(0:1, nInd*ploidy*nLoci, replace=TRUE), ncol=nLoci)
  return(newMapPop(genMap=list(rep(pos, nLoci)), haplotypes=list(H),
                   ploidy=ploidy))
}

# every row of child must equal some row of parent
expect_parental_haplotypes = function(child, parent, label=""){
  key = function(H) apply(H, 1, paste, collapse="")
  expect_true(all(key(child) %in% key(parent)),
              info=paste(label, "a gamete mixes two parental haplotypes"))
}

# Collects every recombination history matrix. Founders are recorded as
# plain integer vectors, so they are skipped.
recHistMatrices = function(x){
  if(is.numeric(x)){
    if(is.matrix(x)){
      return(list(x))
    }
    return(list())
  }
  return(do.call(c, lapply(x, recHistMatrices)))
}

# ---------------------------------------------------------------------------
# PART 1  Map positions are kept as supplied
# ---------------------------------------------------------------------------

test_that("KEEP runMacs returns absolute positions within genLen", {
  set.seed(6101)
  pop = runMacs(nInd=6, nChr=2, segSites=15, nThreads=1)
  for(chr in seq_along(pop@genMap)){
    pos = pop@genMap[[chr]]
    expect_true(all(pos>0))
    expect_true(all(pos<1))
  }
})

test_that("CENTER runMacs and runMacs2 put the centromere at half of genLen", {
  set.seed(6106)
  pop = runMacs(nInd=4, nChr=2, segSites=10, nThreads=1)
  expect_equal(pop@centromere, c(0.5, 0.5))
  SP = SimParam$new(pop)
  expect_equal(unname(SP$centromere), c(0.5, 0.5))

  # A different genLen moves the centromere with it, not with the sites
  pop = runMacs2(nInd=4, nChr=2, segSites=10, genLen=1.5, nThreads=1)
  expect_equal(pop@centromere, c(0.75, 0.75))
  pop = runMacs(nInd=4, nChr=2, segSites=10, manualGenLen=c(1, 2),
                nThreads=1)
  expect_equal(pop@centromere, c(0.5, 1))
})

test_that("KEEP importGenMap, importHaplo and importInbredGeno keep positions", {
  genMap = data.frame(markerName=letters[1:5],
                      chromosome=c(1,1,1,2,2),
                      position=c(0.2,0.5,1,0.15,0.4))
  map = importGenMap(genMap)
  expect_equal(unname(map[["1"]]), c(0.2,0.5,1))
  expect_equal(unname(map[["2"]]), c(0.15,0.4))

  haplo = rbind(c(1,1,0,1,0),
                c(1,1,0,1,0),
                c(0,1,1,0,0),
                c(0,1,1,0,0))
  colnames(haplo) = letters[1:5]
  pop = importHaplo(haplo=haplo, genMap=genMap, ploidy=2L)
  expect_equal(pop@genMap, map)
  pop = importInbredGeno(geno=2*haplo, genMap=genMap)
  expect_equal(pop@genMap, map)

  # Dropping a marker with no genotype leaves the others where they were
  pop = importHaplo(haplo=haplo[,-1,drop=FALSE], genMap=genMap, ploidy=2L)
  expect_equal(unname(pop@genMap[["1"]]), c(0.5,1))
  pop = importInbredGeno(geno=2*haplo[,-1,drop=FALSE], genMap=genMap)
  expect_equal(unname(pop@genMap[["1"]]), c(0.5,1))
})

test_that("KEEP SimParam and switchGenMap keep positions", {
  map = list(seq(0.3, 1.3, length.out=10))
  founderPop = newMapPop(genMap=map,
                         haplotypes=list(matrix(0L, nrow=4, ncol=10)))
  SP = SimParam$new(founderPop)
  expect_equal(unname(SP$genMap[[1]]), map[[1]])

  newMap = list(seq(2, 3, length.out=10))
  SP$switchGenMap(newMap)
  expect_equal(unname(SP$genMap[[1]]), newMap[[1]])
  SP$switchFemaleMap(newMap)
  expect_equal(unname(SP$femaleMap[[1]]), newMap[[1]])
  SP$switchMaleMap(map)
  expect_equal(unname(SP$maleMap[[1]]), map[[1]])
})

test_that("KEEP addSegSite before the first locus leaves the map in place", {
  founderPop = newMapPop(genMap=list(c(0.5,1)),
                         haplotypes=list(matrix(0L, nrow=4, ncol=2)))
  out = addSegSite(founderPop, siteName="x", chr=1, mapPos=0.2,
                   haplo=c(0,1,0,1), nThreads=1L)
  expect_equal(unname(out@genMap[[1]]), c(0.2,0.5,1))
  expect_error(addSegSite(founderPop, siteName="y", chr=1, mapPos=-0.1,
                          haplo=c(0,1,0,1), nThreads=1L),
               "non-negative")
})

# ---------------------------------------------------------------------------
# PART 2  Invalid maps are rejected rather than shifted
# ---------------------------------------------------------------------------

test_that("REJECT negative, unsorted or non-finite positions", {
  H = list(matrix(0L, nrow=4, ncol=3))
  expect_error(newMapPop(genMap=list(c(-0.1,0.5,1)), haplotypes=H),
               "non-negative")
  expect_error(newMapPop(genMap=list(c(0.5,0.1,1)), haplotypes=H),
               "non-decreasing")
  expect_error(newMapPop(genMap=list(c(0,NA,1)), haplotypes=H),
               "finite")
  expect_error(importGenMap(data.frame(markerName=letters[1:2],
                                       chromosome=c(1,1),
                                       position=c(-1,0.5))),
               "non-negative")

  founderPop = newMapPop(genMap=list(c(0,0.5,1)), haplotypes=H)
  SP = SimParam$new(founderPop)
  expect_error(SP$switchGenMap(list(c(-0.1,0.5,1))), "non-negative")
  expect_error(SP$switchFemaleMap(list(c(1,0.5,0))), "non-decreasing")
  expect_error(SP$switchMaleMap(list(c(0,0.5,1)), centromere=-0.1))
})

# ---------------------------------------------------------------------------
# PART 3  Meiosis with crossovers before the first locus
# ---------------------------------------------------------------------------

test_that("RESOLVE crossovers before the first locus in bivalents", {
  set.seed(6102)
  founderPop = stackedMapPop(ploidy=2L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  H = pullSegSiteHaplo(pop, simParam=SP)

  kids = randCross(pop, nCrosses=200, simParam=SP)
  kidH = pullSegSiteHaplo(kids, simParam=SP)
  expect_parental_haplotypes(kidH, H, "makeCross")
  # Both chromosomes of a parent are passed on, not only the first
  expect_gt(length(unique(apply(kidH, 1, paste, collapse=""))), 6L)

  dh = makeDH(pop, nDH=20, simParam=SP)
  expect_parental_haplotypes(pullSegSiteHaplo(dh, simParam=SP), H, "makeDH")

  # With no crossover between loci every history is one row at site 1
  M = recHistMatrices(SP$recHist)
  expect_gt(length(M), 0L)
  expect_true(all(sapply(M, function(m) nrow(m)==1L && m[1,2]==1L)))

  # The recorded history must name the chromosome actually transmitted,
  # which is not always the first one passed to meiosis
  ibd = pullIbdHaplo(kids, simParam=SP)
  expect_equal(dim(ibd), c(400L, 40L))
  expect_true(all(ibd==ibd[,1]))
  expect_equal(unname(kidH), unname(H[ibd[,1],]))
})

test_that("RESOLVE crossovers before the first locus in quadrivalents", {
  set.seed(6103)
  founderPop = stackedMapPop(ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$quadProb = 1
  pop = newPop(founderPop, simParam=SP)
  H = pullSegSiteHaplo(pop, simParam=SP)

  kids = randCross(pop, nCrosses=100, simParam=SP)
  expect_parental_haplotypes(pullSegSiteHaplo(kids, simParam=SP), H,
                             "quadrivalent makeCross")
  gametes = reduceGenome(pop, nProgeny=20, simParam=SP)
  expect_parental_haplotypes(pullSegSiteHaplo(gametes, simParam=SP), H,
                             "quadrivalent reduceGenome")
})

test_that("RESOLVE quadrivalents with the centromere beyond the last locus", {
  # runMacs places the centromere at half of genLen, which can lie past the
  # last site returned when there are few of them
  set.seed(6107)
  founderPop = newMapPop(
    genMap=list(seq(0.1, 0.4, length.out=20)),
    haplotypes=list(matrix(sample(0:1, 24*20, replace=TRUE), nrow=24)),
    ploidy=4L)
  founderPop@centromere = 0.5
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$quadProb = 1
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  kids = randCross(pop, nCrosses=100, simParam=SP)
  gametes = reduceGenome(pop, nProgeny=20, simParam=SP)
  expect_equal(nInd(kids), 100L)
  expect_equal(nInd(gametes), 120L)
  M = recHistMatrices(SP$recHist)
  expect_true(all(sapply(M, function(m){
    m[1,2]==1L && !is.unsorted(m[,2], strictly=TRUE)
  })))
})

test_that("RESOLVE recombination histories keep their form on an offset map", {
  set.seed(6104)
  founderPop = newMapPop(
    genMap=list(seq(0.5, 1.5, length.out=100), seq(1, 3, length.out=100)),
    haplotypes=list(matrix(sample(0:1, 2000, replace=TRUE), nrow=20),
                    matrix(sample(0:1, 2000, replace=TRUE), nrow=20)))
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  for(i in 1:3){
    pop = randCross(pop, nCrosses=50, simParam=SP)
  }
  M = recHistMatrices(SP$recHist)
  expect_true(any(sapply(M, nrow)>1L))
  expect_true(all(sapply(M, function(m){
    m[1,2]==1L && !is.unsorted(m[,2], strictly=TRUE) &&
      all(m[-1,1]!=m[-nrow(m),1])
  })))
  expect_silent(pullIbdHaplo(pop, simParam=SP))
})

test_that("RESOLVE recombination between loci does not depend on an offset", {
  skip_on_cran() # Monte Carlo comparison
  rf = function(offset, n=20000){
    founderPop = newMapPop(genMap=list(c(offset, offset+0.1)),
                           haplotypes=list(rbind(c(0,0), c(1,1))))
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    pop = newPop(founderPop, simParam=SP)
    dh = makeDH(pop, nDH=n, simParam=SP)
    H = pullSegSiteHaplo(dh, haplo=1, simParam=SP)
    return(mean(H[,1]!=H[,2]))
  }
  set.seed(6105)
  # Each estimate has a standard error of about 0.002
  expect_equal(rf(0), rf(3), tolerance=0.012, scale=1)
})
