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

# ---------------------------------------------------------------------------
# PART 4  Chromosomes extending past the last locus
# ---------------------------------------------------------------------------

# A chromosome runs from position 0 to its genetic length, genLen, which can
# lie past the last locus. Crossovers are sampled over the whole chromosome,
# and those after the last locus affect no locus. With every locus at
# position 0 and a longer chromosome, every chiasma falls after the last
# locus, so as in PART 3 each gamete has to be one parental haplotype.

test_that("STORE genLen defaults to the last locus and is kept when given", {
  H = list(matrix(0L, nrow=4, ncol=3), matrix(0L, nrow=4, ncol=2))
  map = list(c(0.1,0.5,0.8), c(0.2,0.4))
  expect_equal(AlphaSimR:::.getGenLen(newMapPop(map, H)), c(0.8,0.4))
  pop = newMapPop(map, H, genLen=c(1,2))
  expect_equal(pop@genLen, c(1,2))
  expect_equal(pop@centromere, c(0.5,1))
  expect_equal(newMapPop(map, H, genLen=1.5)@genLen, c(1.5,1.5))

  # cChr joins the lengths, and sampleHaplo keeps them
  expect_equal(sampleHaplo(pop, nInd=3)@genLen, c(1,2))
  expect_equal(cChr(pop, newMapPop(map[1], H[1]))@genLen, c(1,2,0.8))

  # quickHaplo's chromosomes are genLen long even with a single site
  expect_equal(quickHaplo(nInd=2, nChr=2, segSites=1, genLen=c(1,2))@genLen,
               c(1,2))

  # A SimParam takes the lengths from its founders
  SP = SimParam$new(pop)
  expect_equal(SP$genLen, c(1,2))
  expect_equal(SP$femaleGenLen, c(1,2))
  expect_equal(SP$maleGenLen, c(1,2))
})

test_that("STORE runMacs keeps the requested genetic length", {
  set.seed(6201)
  pop = runMacs(nInd=4, nChr=2, segSites=10, nThreads=1)
  expect_equal(pop@genLen, c(1,1))
  pop = runMacs2(nInd=4, nChr=2, segSites=10, genLen=1.5, nThreads=1)
  expect_equal(pop@genLen, c(1.5,1.5))
  pop = runMacs(nInd=4, nChr=2, segSites=10, manualGenLen=c(1,2),
                nThreads=1)
  expect_equal(pop@genLen, c(1,2))
  expect_true(isTRUE(validObject(pop, test=TRUE)))
})

test_that("STORE importHaplo and importInbredGeno take genLen", {
  genMap = data.frame(markerName=letters[1:5],
                      chromosome=c(1,1,1,2,2),
                      position=c(0.2,0.5,1,0.15,0.4))
  haplo = rbind(c(1,1,0,1,0),
                c(1,1,0,1,0),
                c(0,1,1,0,0),
                c(0,1,1,0,0))
  colnames(haplo) = letters[1:5]
  expect_equal(importHaplo(haplo, genMap)@genLen, c(1,0.4))
  expect_equal(importHaplo(haplo, genMap, genLen=c(1.2,0.9))@genLen,
               c(1.2,0.9))
  expect_equal(importInbredGeno(2*haplo[c(1,3),], genMap, genLen=2)@genLen,
               c(2,2))
  expect_error(importHaplo(haplo, genMap, genLen=c(0.5,0.9)),
               "last locus")
})

test_that("STORE crossing a MapPop keeps its genLen", {
  set.seed(6202)
  inbreds = newMapPop(list(seq(0,0.5,length.out=20)),
                      list(matrix(sample(0:1, 5*20, replace=TRUE), nrow=5)),
                      inbred=TRUE, genLen=1)
  out = hybridCross(inbreds[1:3], inbreds[4:5], nThreads=1L)
  expect_equal(out@genLen, 1)

  out = pedigreeCross(inbreds[1:2], id=as.character(1:4),
                      mother=c(NA, NA, "1", "3"), father=c(NA, NA, "2", "3"))
  expect_equal(out@genLen, 1)
})

test_that("REJECT a genLen shorter than the map or mismatched", {
  H = list(matrix(0L, nrow=4, ncol=3))
  map = list(c(0.1,0.5,0.8))
  expect_error(newMapPop(map, H, genLen=0.5), "last locus")
  expect_error(newMapPop(map, H, genLen=c(1,2)), "one per chromosome")
  expect_error(newMapPop(map, H, genLen=NA_real_), "finite")

  pop = newMapPop(map, H, genLen=1)
  expect_error(c(pop, newMapPop(map, H, genLen=2)), "genLen does not match")

  SP = SimParam$new(pop)
  expect_error(SP$switchGenMap(map, genLen=0.5), "last locus")
  expect_error(SP$switchGenMap(map, centromere=1.5, genLen=1))
})

test_that("SWITCH map replacement sets genLen and setRecombRatio scales it", {
  H = list(matrix(0L, nrow=4, ncol=3))
  pop = newMapPop(list(c(0.1,0.5,0.8)), H, genLen=2)
  SP = SimParam$new(pop)

  # Without genLen a new map ends at its last locus
  SP$switchGenMap(list(c(0.2,0.4,0.6)))
  expect_equal(SP$genLen, 0.6)
  expect_equal(SP$centromere, 0.3)
  SP$switchGenMap(list(c(0.2,0.4,0.6)), genLen=3)
  expect_equal(SP$genLen, 3)
  expect_equal(SP$centromere, 1.5)

  SP$switchFemaleMap(list(c(0.2,0.4,0.6)), genLen=4)
  expect_equal(SP$femaleGenLen, 4)
  expect_equal(SP$maleGenLen, 3)
  expect_equal(SP$genLen, 3.5)
  SP$switchMaleMap(list(c(0.2,0.4,0.6)), genLen=2)
  expect_equal(SP$maleGenLen, 2)

  SP$switchGenMap(list(c(0.2,0.4,0.6)), genLen=3)
  SP$setRecombRatio(2)
  expect_equal(SP$femaleGenLen, 3*2/(1/2+1))
  expect_equal(SP$maleGenLen, 3*2/(2+1))
})

test_that("ADDSITE a site past the chromosome needs extendGenLen", {
  founderPop = newMapPop(genMap=list(c(0.5,1)),
                         haplotypes=list(matrix(0L, nrow=4, ncol=2)),
                         genLen=1.5)
  out = addSegSite(founderPop, siteName="x", chr=1, mapPos=1.2,
                   haplo=c(0,1,0,1), nThreads=1L)
  expect_equal(out@genLen, 1.5)
  expect_error(addSegSite(founderPop, siteName="y", chr=1, mapPos=2,
                          haplo=c(0,1,0,1), nThreads=1L),
               "extendGenLen")
  out = addSegSite(founderPop, siteName="y", chr=1, mapPos=2,
                   haplo=c(0,1,0,1), nThreads=1L, extendGenLen=TRUE)
  expect_equal(out@genLen, 2)
  expect_equal(unname(out@genMap[[1]]), c(0.5,1,2))
  expect_equal(out@centromere, founderPop@centromere)
})

test_that("FALLBACK objects without genLen end at their last locus", {
  H = list(matrix(0L, nrow=4, ncol=3))
  pop = newMapPop(list(c(0.1,0.5,0.8)), H, genLen=2)

  # A MapPop saved before genLen existed has no slot at all
  old = pop
  attr(old, "genLen") = NULL
  expect_false(methods::.hasSlot(old, "genLen"))
  expect_equal(AlphaSimR:::.getGenLen(old), 0.8)
  expect_equal(SimParam$new(old)$genLen, 0.8)

  # One built directly with new() can leave it empty
  built = new("MapPop", nInd=pop@nInd, nChr=pop@nChr, ploidy=pop@ploidy,
              nLoci=pop@nLoci, geno=pop@geno, genMap=pop@genMap,
              centromere=pop@centromere, inbred=pop@inbred)
  expect_equal(AlphaSimR:::.getGenLen(built), 0.8)

  # A SimParam saved before genLen existed gives NULL for the field
  oldSP = list(femaleMap=pop@genMap, maleMap=pop@genMap)
  expect_equal(AlphaSimR:::.simGenLen(oldSP, female=TRUE), 0.8)
  expect_equal(AlphaSimR:::.simGenLen(oldSP, female=FALSE), 0.8)
})

test_that("RESOLVE crossovers after the last locus in bivalents", {
  set.seed(6203)
  nLoci = 40
  H = matrix(sample(0:1, 12*nLoci, replace=TRUE), ncol=nLoci)
  founderPop = newMapPop(genMap=list(rep(0, nLoci)), haplotypes=list(H),
                         genLen=2)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackRec(TRUE)
  pop = newPop(founderPop, simParam=SP)
  H = pullSegSiteHaplo(pop, simParam=SP)

  kids = randCross(pop, nCrosses=200, simParam=SP)
  kidH = pullSegSiteHaplo(kids, simParam=SP)
  expect_parental_haplotypes(kidH, H, "makeCross")
  dh = makeDH(pop, nDH=20, simParam=SP)
  expect_parental_haplotypes(pullSegSiteHaplo(dh, simParam=SP), H, "makeDH")

  # Every history is one row at site 1, and it names the chromosome that
  # was transmitted
  M = recHistMatrices(SP$recHist)
  expect_true(all(sapply(M, function(m) nrow(m)==1L && m[1,2]==1L)))
  ibd = pullIbdHaplo(kids, simParam=SP)
  expect_true(all(ibd==ibd[,1]))
  expect_equal(unname(kidH), unname(H[ibd[,1],]))
})

test_that("RESOLVE crossovers after the last locus in quadrivalents", {
  set.seed(6204)
  nLoci = 40
  H = matrix(sample(0:1, 24*nLoci, replace=TRUE), ncol=nLoci)
  founderPop = newMapPop(genMap=list(rep(0, nLoci)), haplotypes=list(H),
                         ploidy=4L, genLen=2)
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

test_that("REACH genLen is what meiosis samples crossovers over", {
  # Extra length after the last locus leaves bivalent recombination between
  # loci unchanged in distribution, so the only sign that it reaches
  # meiosis is that more chiasmata are drawn. From the same seed, changing
  # nothing but genLen must then change the progeny.
  set.seed(6207)
  map = list(seq(0, 1, length.out=50))
  H = list(matrix(sample(0:1, 16*50, replace=TRUE), nrow=16))
  run = function(genLen, ploidy){
    founderPop = newMapPop(map, H, ploidy=ploidy, genLen=genLen)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    pop = newPop(founderPop, simParam=SP)
    set.seed(1)
    if(ploidy==2L){
      return(list(cross=pullSegSiteHaplo(randCross(pop, 20, simParam=SP),
                                         simParam=SP),
                  dh=pullSegSiteHaplo(makeDH(pop, 5, simParam=SP),
                                      simParam=SP)))
    }
    return(list(gamete=pullSegSiteHaplo(reduceGenome(pop, 5, simParam=SP),
                                        simParam=SP)))
  }
  short = run(1, 2L)
  long = run(4, 2L)
  expect_false(identical(short$cross, long$cross))
  expect_false(identical(short$dh, long$dh))
  expect_false(identical(run(1, 4L)$gamete, run(4, 4L)$gamete))
})

test_that("RESOLVE histories keep their form with length at both ends", {
  set.seed(6205)
  founderPop = newMapPop(
    genMap=list(seq(0.5, 1.5, length.out=100), seq(1, 2, length.out=100)),
    haplotypes=list(matrix(sample(0:1, 2000, replace=TRUE), nrow=20),
                    matrix(sample(0:1, 2000, replace=TRUE), nrow=20)),
    genLen=c(2.5, 3))
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
      all(m[,2]<=100L) && all(m[-1,1]!=m[-nrow(m),1])
  })))
  expect_silent(pullIbdHaplo(pop, simParam=SP))
})

test_that("RESOLVE recombination between loci does not depend on genLen", {
  skip_on_cran() # Monte Carlo comparison
  rf = function(genLen, n=20000){
    founderPop = newMapPop(genMap=list(c(0.5, 0.6)),
                           haplotypes=list(rbind(c(0,0), c(1,1))),
                           genLen=genLen)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    pop = newPop(founderPop, simParam=SP)
    dh = makeDH(pop, nDH=n, simParam=SP)
    H = pullSegSiteHaplo(dh, haplo=1, simParam=SP)
    return(mean(H[,1]!=H[,2]))
  }
  set.seed(6206)
  # Each estimate has a standard error of about 0.002
  expect_equal(rf(0.6), rf(3), tolerance=0.012, scale=1)
})
