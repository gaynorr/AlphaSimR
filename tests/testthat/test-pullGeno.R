context("pullGeno")

# The genotype and haplotype accessors gather bits out of the packed cube in
# tiles of 64 individuals by 256 loci. The counts below are chosen so that
# neither the individuals nor the loci divide evenly into a tile, and so that
# the loci do not divide evenly into the 8 loci held in a byte.

makeTestMap = function(nPerChr, prefix = "M") {
  data.frame(
    markerName = paste0(prefix, seq_len(sum(nPerChr))),
    chromosome = rep(seq_along(nPerChr), nPerChr),
    position = unlist(lapply(nPerChr, function(n) seq(0, 1, length.out = n))),
    stringsAsFactors = FALSE
  )
}

# Builds a population from a haplotype matrix that is known in advance, so the
# expected values do not come from AlphaSimR itself
knownPop = function(nInd, nPerChr, ploidy = 2L, prefix = "M") {
  genMap = makeTestMap(nPerChr, prefix)
  nMarker = sum(nPerChr)
  haplo = matrix(sample(0:1, nInd * ploidy * nMarker, replace = TRUE),
                 nrow = nInd * ploidy, ncol = nMarker)
  # importHaplo takes marker names from the column names
  colnames(haplo) = genMap$markerName
  founderPop = importHaplo(haplo = haplo, genMap = genMap, ploidy = ploidy)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop, simParam = SP)
  # Order the reference to match the markers that come back
  ord = match(colnames(pullSegSiteHaplo(pop, simParam = SP)), colnames(haplo))
  list(pop = pop, SP = SP, haplo = haplo[, ord, drop = FALSE],
       nInd = nInd, ploidy = ploidy, nMarker = nMarker)
}

test_that("genotypes and haplotypes match a known haplotype matrix", {
  set.seed(5001)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))

  H = pullSegSiteHaplo(x$pop, simParam = x$SP)
  expect_equal(dim(H), dim(x$haplo))
  expect_true(all(H == x$haplo))

  odd = seq(1, nrow(x$haplo), by = 2)
  even = seq(2, nrow(x$haplo), by = 2)

  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  expect_equal(dim(G), c(x$nInd, x$nMarker))
  expect_true(all(G == x$haplo[odd, ] + x$haplo[even, ]))

  # getMaternalGeno and getPaternalGeno
  expect_true(all(pullSegSiteHaplo(x$pop, haplo = 1, simParam = x$SP) ==
                    x$haplo[odd, ]))
  expect_true(all(pullSegSiteHaplo(x$pop, haplo = 2, simParam = x$SP) ==
                    x$haplo[even, ]))
})

test_that("tetraploid genotypes and haplotypes match a known matrix", {
  set.seed(5002)
  x = knownPop(nInd = 20L, nPerChr = c(90L, 70L), ploidy = 4L, prefix = "T")

  H = pullSegSiteHaplo(x$pop, simParam = x$SP)
  expect_true(all(H == x$haplo))

  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  expected = t(sapply(seq_len(x$nInd), function(i)
    colSums(x$haplo[((i - 1L) * 4L + 1L):(i * 4L), , drop = FALSE])))
  expect_equal(dim(G), dim(expected))
  expect_true(all(G == expected))

  # A single haplotype of a polyploid
  expect_true(all(pullSegSiteHaplo(x$pop, haplo = 1, simParam = x$SP) ==
                    x$haplo[seq(1, nrow(x$haplo), by = 4L), ]))
})

test_that("gathering is correct at the edges of a tile", {
  set.seed(5003)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))

  # 64 individuals is one tile
  for (n in c(1L, 2L, 63L, 64L, 65L)) {
    sub = x$pop[seq_len(n)]
    rows = seq_len(n * 2L)
    G = pullSegSiteGeno(sub, simParam = x$SP)
    H = pullSegSiteHaplo(sub, simParam = x$SP)
    expect_equal(dim(G), c(n, x$nMarker))
    expect_true(all(H == x$haplo[rows, , drop = FALSE]))
    expect_true(all(G == x$haplo[rows[c(TRUE, FALSE)], , drop = FALSE] +
                      x$haplo[rows[c(FALSE, TRUE)], , drop = FALSE]))
  }

  # A single locus
  one = pullMarkerHaplo(x$pop, markers = colnames(x$haplo)[7],
                        simParam = x$SP)
  expect_true(all(as.vector(one) == x$haplo[, 7]))
})

test_that("a SNP chip returns the same values as the full site set", {
  set.seed(5004)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))
  x$SP$addSnpChip(nSnpPerChr = c(90L, 70L))

  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  H = pullSegSiteHaplo(x$pop, simParam = x$SP)
  snpG = pullSnpGeno(x$pop, simParam = x$SP)
  snpH = pullSnpHaplo(x$pop, simParam = x$SP)

  keep = match(colnames(snpG), colnames(G))
  expect_false(anyNA(keep))
  expect_true(all(snpG == G[, keep, drop = FALSE]))
  expect_true(all(snpH == H[, keep, drop = FALSE]))
})

test_that("the raw genotype path agrees with the integer path", {
  set.seed(5005)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))
  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  gRaw = pullSegSiteGeno(x$pop, asRaw = TRUE, simParam = x$SP)
  expect_true(all(as.integer(gRaw) == as.integer(G)))
})

test_that("gathering does not depend on the number of threads", {
  skip_on_cran()
  nThreads = getNumThreads()
  skip_if_not(nThreads > 1L, "only one thread available")

  set.seed(5006)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))

  expect_identical(pullSegSiteGeno(x$pop, nThreads = 1L, simParam = x$SP),
                   pullSegSiteGeno(x$pop, nThreads = nThreads, simParam = x$SP))
  expect_identical(pullSegSiteHaplo(x$pop, nThreads = 1L, simParam = x$SP),
                   pullSegSiteHaplo(x$pop, nThreads = nThreads, simParam = x$SP))
  expect_identical(pullSegSiteHaplo(x$pop, haplo = 1, nThreads = 1L,
                                    simParam = x$SP),
                   pullSegSiteHaplo(x$pop, haplo = 1, nThreads = nThreads,
                                    simParam = x$SP))
})

test_that("genotypes of a crossed population equal the sum of its haplotypes", {
  set.seed(5007)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))
  prog = randCross(x$pop, nCrosses = 40L, nProgeny = 2L, simParam = x$SP)

  G = pullSegSiteGeno(prog, simParam = x$SP)
  H = pullSegSiteHaplo(prog, simParam = x$SP)
  expect_true(all(G == H[seq(1, nrow(H), by = 2), ] +
                    H[seq(2, nrow(H), by = 2), ]))
  expect_true(all(G >= 0L) && all(G <= 2L))
})

# A single haplotype and the loci are read through raw pointers in C++, so a
# value outside the population would read memory beyond its genotypes. Both
# the R wrappers and the C++ code refuse them.
test_that("REFUSE a haplo outside 1 to ploidy", {
  set.seed(5008)
  x = knownPop(nInd = 5L, nPerChr = c(20L, 12L))
  # Traits have to be added before any individuals are made
  x$SP$resetPed()
  x$SP$addTraitA(nQtlPerChr = 2L)
  x$pop = newPop(x$SP$founderPop, simParam = x$SP)
  markers = colnames(pullSegSiteHaplo(x$pop, simParam = x$SP))
  for (bad in list(0, -1, 3, NA, 1.5, c(1, 2), "1")) {
    expect_error(pullSegSiteHaplo(x$pop, haplo = bad, simParam = x$SP),
                 "haplo must be")
    expect_error(pullQtlHaplo(x$pop, haplo = bad, simParam = x$SP),
                 "haplo must be")
    expect_error(pullMarkerHaplo(x$pop, markers = markers, haplo = bad,
                                 simParam = x$SP),
                 "haplo must be")
  }
  # Valid selectors still work, and as a double as well as an integer
  expect_equal(pullSegSiteHaplo(x$pop, haplo = 2, simParam = x$SP),
               pullSegSiteHaplo(x$pop, haplo = 2L, simParam = x$SP))
  expect_equal(nrow(pullSegSiteHaplo(x$pop, haplo = "all", simParam = x$SP)),
               2L * x$pop@nInd)

  geno = x$pop@geno
  lociPerChr = x$pop@nLoci
  lociLoc = c(seq_len(20L), seq_len(12L))
  for (bad in c(0L, -1L, 3L, NA_integer_)) {
    expect_error(AlphaSimR:::getOneHaplo(geno, lociPerChr, lociLoc, bad, 1L),
                 "haplo must be between")
  }
})

test_that("REFUSE loci outside their chromosome", {
  set.seed(5009)
  x = knownPop(nInd = 5L, nPerChr = c(20L, 12L))
  geno = x$pop@geno
  lociPerChr = x$pop@nLoci
  lociLoc = c(seq_len(20L), seq_len(12L))
  # The last chromosome holds 12 loci in 2 bytes, so 17 is past its end
  for (bad in c(0L, -1L, 17L, NA_integer_)) {
    badLoc = lociLoc
    badLoc[32L] = bad
    expect_error(AlphaSimR:::getGeno(geno, lociPerChr, badLoc, 1L),
                 "outside its chromosome")
    expect_error(AlphaSimR:::getHaplo(geno, lociPerChr, badLoc, 1L),
                 "outside its chromosome")
    expect_error(AlphaSimR:::getMaternalGeno(geno, lociPerChr, badLoc, 1L),
                 "outside its chromosome")
    expect_error(AlphaSimR:::getPaternalGeno(geno, lociPerChr, badLoc, 1L),
                 "outside its chromosome")
    expect_error(AlphaSimR:::getOneHaplo(geno, lociPerChr, badLoc, 1L, 1L),
                 "outside its chromosome")
  }
  # lociPerChr has to account for every locus and every chromosome
  expect_error(AlphaSimR:::getGeno(geno, c(20L, 13L), lociLoc, 1L),
               "greater than length")
  expect_error(AlphaSimR:::getGeno(geno, c(20L, 11L), lociLoc, 1L),
               "less than length")
  expect_error(AlphaSimR:::getGeno(geno, 32L, lociLoc, 1L),
               "number of chromosomes")
  expect_error(AlphaSimR:::getGeno(geno, c(-1L, 33L), lociLoc, 1L),
               "negative")
  # setHaplo writes rather than reads, and checks the haplotypes it is given
  haplo = AlphaSimR:::getHaplo(geno, lociPerChr, lociLoc, 1L)
  badLoc = lociLoc
  badLoc[1L] = 0L
  expect_error(AlphaSimR:::setHaplo(geno, haplo, lociPerChr, badLoc, 1L),
               "outside its chromosome")
  expect_error(AlphaSimR:::setHaplo(geno, haplo[-1L, ], lociPerChr,
                                    lociLoc, 1L),
               "haplo must have")
  # The loci that are valid still give the genotypes
  expect_equal(AlphaSimR:::getGeno(geno, lociPerChr, lociLoc, 1L),
               unname(pullSegSiteGeno(x$pop, asRaw = TRUE, simParam = x$SP)))
})

test_that("REFUSE a manual trait with loci outside their chromosome", {
  set.seed(5010)
  x = knownPop(nInd = 5L, nPerChr = c(20L, 12L))
  # Traits have to be added before any individuals are made
  x$SP$resetPed()
  x$SP$addTraitA(nQtlPerChr = 2L)
  trait = x$SP$traits[[1L]]
  # The LociMap validity method only checks counts, so these are valid
  # objects that do not fit the genotypes
  for (bad in c(0L, 100L)) {
    badTrait = trait
    badTrait@lociLoc[1L] = bad
    expect_true(validObject(badTrait))
    expect_error(x$SP$manAddTrait(badTrait), "outside its chromosome")
  }
  expect_equal(x$SP$nTraits, 1L)
})
