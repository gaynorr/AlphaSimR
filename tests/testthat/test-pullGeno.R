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

test_that("genotypes, haplotypes and allele frequencies match a known haplotype matrix", {
  set.seed(5001)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))

  H = pullSegSiteHaplo(x$pop, simParam = x$SP)
  expect_equal(dim(H), dim(x$haplo))
  expect_true(all(H == x$haplo))

  odd = seq(1, nrow(x$haplo), by = 2)
  even = seq(2, nrow(x$haplo), by = 2)

  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  expectedG = x$haplo[odd, ] + x$haplo[even, ]
  expect_equal(dim(G), c(x$nInd, x$nMarker))
  expect_true(all(G == expectedG))

  # getMaternalGeno and getPaternalGeno
  expect_true(all(pullSegSiteHaplo(x$pop, haplo = 1, simParam = x$SP) ==
                    x$haplo[odd, ]))
  expect_true(all(pullSegSiteHaplo(x$pop, haplo = 2, simParam = x$SP) ==
                    x$haplo[even, ]))

  Af = pullSegSiteAf(x$pop, simParam = x$SP)
  expect_equal(dim(Af), c(1, x$nMarker))
  expect_true(all(Af == colMeans(expectedG) / x$ploidy))
})

test_that("tetraploid genotypes, haplotypes and allele frequencies match a known matrix", {
  set.seed(5002)
  x = knownPop(nInd = 20L, nPerChr = c(90L, 70L), ploidy = 4L, prefix = "T")

  H = pullSegSiteHaplo(x$pop, simParam = x$SP)
  expect_true(all(H == x$haplo))

  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  expectedG = t(sapply(seq_len(x$nInd), function(i)
    colSums(x$haplo[((i - 1L) * 4L + 1L):(i * 4L), , drop = FALSE])))
  expect_equal(dim(G), dim(expectedG))
  expect_true(all(G == expectedG))

  # A single haplotype of a polyploid
  expect_true(all(pullSegSiteHaplo(x$pop, haplo = 1, simParam = x$SP) ==
                    x$haplo[seq(1, nrow(x$haplo), by = 4L), ]))

  Af = pullSegSiteAf(x$pop, simParam = x$SP)
  expect_equal(dim(Af), c(1, x$nMarker))
  expect_true(all(Af == colMeans(expectedG) / x$ploidy))
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
    expectedG = x$haplo[rows[c(TRUE, FALSE)], , drop = FALSE] +
                 x$haplo[rows[c(FALSE, TRUE)], , drop = FALSE]
    expect_equal(dim(G), c(n, x$nMarker))
    expect_true(all(H == x$haplo[rows, , drop = FALSE]))
    expect_true(all(G == expectedG))

    Af = pullSegSiteAf(sub, simParam = x$SP)
    expect_true(all(Af == colMeans(expectedG) / x$ploidy))
  }

  # A single locus
  one = pullMarkerHaplo(x$pop, markers = colnames(x$haplo)[7],
                        simParam = x$SP)
  expect_true(all(as.vector(one) == x$haplo[, 7]))

  markerAf = pullMarkerAf(x$pop, markers = colnames(x$haplo)[7],
                          simParam = x$SP)
  expect_equal(as.numeric(markerAf), mean(x$haplo[, 7]))
})

test_that("a SNP chip returns the same values as the full site set", {
  set.seed(5004)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))
  x$SP$addSnpChip(nSnpPerChr = c(90L, 70L), name = "160SNPchip")

  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  H = pullSegSiteHaplo(x$pop, simParam = x$SP)
  snpG = pullSnpGeno(x$pop, simParam = x$SP)
  expect_equal(snpG, pullSnpGeno(x$pop, snpChip = "160SNPchip", simParam = x$SP))
  snpH = pullSnpHaplo(x$pop, simParam = x$SP)
  expect_equal(snpH, pullSnpHaplo(x$pop, snpChip = "160SNPchip",
                                  haplo = "all", simParam = x$SP))

  keep = match(colnames(snpG), colnames(G))
  expect_false(anyNA(keep))
  expect_true(all(snpG == G[, keep, drop = FALSE]))
  expect_true(all(snpH == H[, keep, drop = FALSE]))

  segAf = pullSegSiteAf(x$pop, simParam = x$SP)
  snpAf = pullSnpAf(x$pop, simParam = x$SP)
  expect_equal(snpAf, pullSnpAf(x$pop, snpChip = "160SNPchip",
                                nThreads = 1L, simParam = x$SP))
  expect_true(all(snpAf == segAf[, keep, drop = FALSE]))
  expect_error(
    pullSnpGeno(x$pop, snpChip = "SNPchip", simParam = x$SP),
    paste0("'SNPchip' did not match any SNP chip names")
  )
  expect_error(
    pullSnpHaplo(x$pop, snpChip = "SNPchip", simParam = x$SP),
    paste0("'SNPchip' did not match any SNP chip names")
  )
  expect_error(
    pullSnpAf(x$pop, snpChip = "SNPchip", simParam = x$SP),
    paste0("'SNPchip' did not match any SNP chip names")
  )
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
  x$SP$resetPed()
  x$SP$addTraitA(10, name="Trait1")
  x$SP$restrSegSites(10, 10, overlap = FALSE)
  x$SP$addSnpChip(nSnpPerChr = 10, name = "20SNPchip")

  expect_identical(
    pullSegSiteGeno(x$SP$founderPop, nThreads = 1L, simParam = x$SP),
    pullSegSiteGeno(x$SP$founderPop, nThreads = nThreads, simParam = x$SP)
  )
  expect_identical(
    pullSegSiteHaplo(x$SP$founderPop, nThreads = 1L, simParam = x$SP),
    pullSegSiteHaplo(x$SP$founderPop, nThreads = nThreads, simParam = x$SP)
  )
  expect_identical(
    pullSegSiteHaplo(
      x$SP$founderPop,
      haplo = 1,
      nThreads = 1L,
      simParam = x$SP
    ),
    pullSegSiteHaplo(
      x$SP$founderPop,
      haplo = 1,
      nThreads = nThreads,
      simParam = x$SP
    )
  )
  expect_identical(
    pullSegSiteAf(x$SP$founderPop, nThreads = 1L, simParam = x$SP),
    pullSegSiteAf(x$SP$founderPop, nThreads = nThreads, simParam = x$SP)
  )

  expect_identical(pullSegSiteGeno(x$pop, nThreads = 1L, simParam = x$SP),
                   pullSegSiteGeno(x$pop, nThreads = nThreads, simParam = x$SP))
  expect_identical(pullSegSiteHaplo(x$pop, nThreads = 1L, simParam = x$SP),
                   pullSegSiteHaplo(x$pop, nThreads = nThreads, simParam = x$SP))
  expect_identical(pullSegSiteHaplo(x$pop, haplo = 1, nThreads = 1L,
                                    simParam = x$SP),
                   pullSegSiteHaplo(x$pop, haplo = 1, nThreads = nThreads,
                                    simParam = x$SP))
  expect_identical(
    pullSegSiteAf(x$pop, nThreads = 1L, simParam = x$SP),
    pullSegSiteAf(x$pop, nThreads = nThreads, simParam = x$SP)
  )

  expect_identical(
    pullQtlGeno(x$pop, nThreads = 1L, simParam = x$SP),
    pullQtlGeno(x$pop, nThreads = nThreads, simParam = x$SP)
  )
  expect_identical(
    pullQtlHaplo(x$pop, nThreads = 1L, simParam = x$SP),
    pullQtlHaplo(x$pop, nThreads = nThreads, simParam = x$SP)
  )
  expect_identical(
    pullQtlHaplo(x$pop, haplo = 1, nThreads = 1L, simParam = x$SP),
    pullQtlHaplo(x$pop, haplo = 1, nThreads = nThreads, simParam = x$SP)
  )
  expect_identical(
    pullQtlAf(x$pop, nThreads = 1L, simParam = x$SP),
    pullQtlAf(x$pop, nThreads = nThreads, simParam = x$SP)
  )

  expect_identical(
    pullSnpGeno(x$pop, nThreads = 1L, simParam = x$SP),
    pullSnpGeno(x$pop, nThreads = nThreads, simParam = x$SP)
  )
  expect_identical(
    pullSnpHaplo(x$pop, nThreads = 1L, simParam = x$SP),
    pullSnpHaplo(x$pop, nThreads = nThreads, simParam = x$SP)
  )
  expect_identical(
    pullSnpHaplo(x$pop, haplo = 1, nThreads = 1L, simParam = x$SP),
    pullSnpHaplo(x$pop, haplo = 1, nThreads = nThreads, simParam = x$SP)
  )
  expect_identical(
    pullSnpAf(x$pop, nThreads = 1L, simParam = x$SP),
    pullSnpAf(x$pop, nThreads = nThreads, simParam = x$SP)
  )

  markers <- paste0('M', sample(300, 10))

  expect_identical(
    pullMarkerGeno(
      x$SP$founderPop,
      markers = markers,
      nThreads = 1L,
      simParam = x$SP
    ),
    pullMarkerGeno(
      x$SP$founderPop,
      markers = markers,
      nThreads = nThreads,
      simParam = x$SP
    )
  )
  expect_identical(
    pullMarkerHaplo(
      x$SP$founderPop,
      markers = markers,
      nThreads = 1L,
      simParam = x$SP
    ),
    pullMarkerHaplo(
      x$SP$founderPop,
      markers = markers,
      nThreads = nThreads,
      simParam = x$SP
    )
  )
  expect_identical(
    pullMarkerHaplo(
      x$SP$founderPop,
      markers = markers,
      haplo = 1,
      nThreads = 1L,
      simParam = x$SP
    ),
    pullMarkerHaplo(
      x$SP$founderPop,
      markers = markers,
      haplo = 1,
      nThreads = nThreads,
      simParam = x$SP
    )
  )
  expect_identical(
    pullMarkerAf(
      x$SP$founderPop,
      markers = markers,
      nThreads = 1L,
      simParam = x$SP
    ),
    pullMarkerAf(
      x$SP$founderPop,
      markers = markers,
      nThreads = nThreads,
      simParam = x$SP
    )
  )

  expect_identical(
    pullMarkerGeno(x$pop, markers = markers, nThreads = 1L, simParam = x$SP),
    pullMarkerGeno(
      x$pop,
      markers = markers,
      nThreads = nThreads,
      simParam = x$SP
    )
  )
  expect_identical(
    pullMarkerHaplo(x$pop, markers = markers, nThreads = 1L, simParam = x$SP),
    pullMarkerHaplo(
      x$pop,
      markers = markers,
      nThreads = nThreads,
      simParam = x$SP
    )
  )
  expect_identical(
    pullMarkerHaplo(
      x$pop,
      markers = markers,
      haplo = 1,
      nThreads = 1L,
      simParam = x$SP
    ),
    pullMarkerHaplo(
      x$pop,
      markers = markers,
      haplo = 1,
      nThreads = nThreads,
      simParam = x$SP
    )
  )
  expect_identical(
    pullMarkerAf(x$pop, markers = markers, nThreads = 1L, simParam = x$SP),
    pullMarkerAf(x$pop, markers = markers, nThreads = nThreads, simParam = x$SP)
  )
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

  Af = pullSegSiteAf(prog, simParam = x$SP)
  expect_true(all(Af == colMeans(G) / x$ploidy))
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

# The accessors read the genotypes in place rather than through a copy, so
# nothing they do may write to the population. A deep copy taken first is the
# only reference that would show such a write: an ordinary assignment in R
# shares memory with the original until one of them is modified.
deepCopy = function(x) unserialize(serialize(x, NULL))

test_that("KEEP the genotypes unchanged by every accessor", {
  set.seed(5011)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L), ploidy = 4L)
  before = deepCopy(x$pop@geno)
  lociPerChr = x$pop@nLoci
  lociLoc = c(seq_len(170L), seq_len(130L))
  invisible(pullSegSiteGeno(x$pop, simParam = x$SP))
  invisible(pullSegSiteHaplo(x$pop, simParam = x$SP))
  invisible(pullSegSiteHaplo(x$pop, haplo = 3, simParam = x$SP))
  invisible(pullSegSiteAf(x$pop, simParam = x$SP))
  invisible(AlphaSimR:::getMaternalGeno(x$pop@geno, lociPerChr, lociLoc, 2L))
  invisible(AlphaSimR:::getPaternalGeno(x$pop@geno, lociPerChr, lociLoc, 2L))
  invisible(AlphaSimR:::calcGenoFreq(x$pop@geno, lociPerChr, lociLoc, 2L))
  invisible(AlphaSimR:::calcChrFreq(x$pop@geno[[2L]]))
  expect_identical(x$pop@geno, before)
})

test_that("SPLIT the haplotypes between the maternal and paternal dosages", {
  set.seed(5012)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L), ploidy = 4L)
  lociPerChr = x$pop@nLoci
  lociLoc = c(seq_len(170L), seq_len(130L))
  # Rows of the known haplotypes for haplotype h of every individual
  hapRows = function(h) seq(h, nrow(x$haplo), by = x$ploidy)
  mat = AlphaSimR:::getMaternalGeno(x$pop@geno, lociPerChr, lociLoc, 1L)
  pat = AlphaSimR:::getPaternalGeno(x$pop@geno, lociPerChr, lociLoc, 1L)
  expect_true(all(mat == x$haplo[hapRows(1), ] + x$haplo[hapRows(2), ]))
  expect_true(all(pat == x$haplo[hapRows(3), ] + x$haplo[hapRows(4), ]))
  # Together they are the genotype dosage
  expect_identical(
    as.integer(mat) + as.integer(pat),
    as.integer(AlphaSimR:::getGeno(x$pop@geno, lociPerChr, lociLoc, 1L)))
})

test_that("COUNT allele frequencies exactly from the packed genotypes", {
  set.seed(5013)
  for (ploidy in c(2L, 4L)) {
    x = knownPop(nInd = 70L, nPerChr = c(170L, 130L), ploidy = ploidy)
    G = pullSegSiteGeno(x$pop, simParam = x$SP)
    lociPerChr = x$pop@nLoci
    lociLoc = c(seq_len(170L), seq_len(130L))
    # The column sums are whole numbers, so dividing them in this order
    # gives the frequencies to the last bit
    expected = unname(colSums(G))/x$nInd/ploidy
    public = pullSegSiteAf(x$pop, simParam = x$SP)
    expect_identical(c(public), expected)
    for (nThreads in c(1L, 2L)) {
      p = AlphaSimR:::calcGenoFreq(x$pop@geno, lociPerChr, lociLoc, nThreads)
      expect_identical(c(p), expected)
    }
    # calcChrFreq covers every bit of a chromosome's bytes, and the bits
    # past the last locus hold no allele
    q = c(AlphaSimR:::calcChrFreq(x$pop@geno[[2L]]))
    expect_equal(length(q), 8L*dim(x$pop@geno[[2L]])[1L])
    expect_identical(q[seq_len(130L)], expected[170L + seq_len(130L)])
    expect_true(all(q[-seq_len(130L)] == 0))
  }
  # A population without individuals has no frequencies. It gave an error
  # when the frequencies were the means of a dosage matrix, and still does.
  empty = x$pop@geno
  empty = lapply(empty, function(chr) chr[, , integer(), drop = FALSE])
  expect_error(AlphaSimR:::calcGenoFreq(empty, lociPerChr, lociLoc, 1L))
  expect_error(AlphaSimR:::calcChrFreq(empty[[1L]]))
})

test_that("REFUSE chromosomes that disagree on ploidy or individuals", {
  set.seed(5014)
  x = knownPop(nInd = 5L, nPerChr = c(20L, 12L))
  lociPerChr = x$pop@nLoci
  lociLoc = c(seq_len(20L), seq_len(12L))
  short = x$pop@geno
  short[[2L]] = short[[2L]][, , 1:3, drop = FALSE]
  expect_error(AlphaSimR:::getGeno(short, lociPerChr, lociLoc, 1L),
               "same ploidy and individuals")
  haploid = x$pop@geno
  haploid[[2L]] = haploid[[2L]][, 1L, , drop = FALSE]
  expect_error(AlphaSimR:::getHaplo(haploid, lociPerChr, lociLoc, 1L),
               "same ploidy and individuals")
  notRaw = x$pop@geno
  notRaw[[1L]] = array(as.integer(notRaw[[1L]]), dim = dim(notRaw[[1L]]))
  expect_error(AlphaSimR:::getGeno(notRaw, lociPerChr, lociLoc, 1L),
               "raw vectors")
})

test_that("KEEP the input population unchanged in setMarkerHaplo", {
  set.seed(5011)
  founder = quickHaplo(nInd = 4, nChr = 2, segSites = 16)
  SP = SimParam$new(founder)
  SP$nThreads = 1L
  pop = newPop(founder, simParam = SP)
  # A deep copy, so it cannot share memory with pop@geno
  before = unserialize(serialize(pop@geno, NULL))
  other = pop
  H = pullSegSiteHaplo(pop, simParam = SP)[, 1:3]
  H[] = 1L - H
  edited = setMarkerHaplo(pop, H, simParam = SP)
  expect_identical(pop@geno, before)
  expect_identical(other@geno, before)
  expect_equal(pullSegSiteHaplo(pop, simParam = SP)[, 1:3], 1L - H)
  # The result carries the new haplotypes and nothing else changes
  newHaplo = pullSegSiteHaplo(edited, simParam = SP)
  oldHaplo = pullSegSiteHaplo(pop, simParam = SP)
  expect_equal(newHaplo[, 1:3], H)
  expect_equal(newHaplo[, -(1:3)], oldHaplo[, -(1:3)])
  # The second chromosome has no loci set, so it is the same as before
  expect_identical(edited@geno[[2L]], before[[2L]])
})

test_that("pullQtl returns the same values as the full site set", {
  set.seed(5004)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))
  x$SP$resetPed()
  x$SP$addTraitA(10, name="Trait1")

  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  H = pullSegSiteHaplo(x$pop, simParam = x$SP)
  qtlG = pullQtlGeno(x$pop, simParam = x$SP)
  expect_equal(qtlG, pullQtlGeno(x$pop, trait = "Trait1", simParam = x$SP))
  qtlH = pullQtlHaplo(x$pop, simParam = x$SP)
  expect_equal(qtlH, pullQtlHaplo(x$pop, trait = "Trait1",
                                  haplo = "all", simParam = x$SP))

  keep = match(colnames(qtlG), colnames(G))
  expect_false(anyNA(keep))
  expect_true(all(qtlG == G[, keep, drop = FALSE]))
  expect_true(all(qtlH == H[, keep, drop = FALSE]))

  segAf = pullSegSiteAf(x$pop, simParam = x$SP)
  qtlAf = pullQtlAf(x$pop, simParam = x$SP)
  expect_equal(qtlAf, pullQtlAf(x$pop, trait = "Trait1", simParam = x$SP))
  expect_true(all(qtlAf == segAf[, keep, drop = FALSE]))

  expect_error(
    pullQtlGeno(x$pop, trait = "noTrait", simParam = x$SP),
    paste0("'noTrait' did not match any trait names")
  )
  expect_error(
    pullQtlHaplo(x$pop, trait = "noTrait", simParam = x$SP),
    paste0("'noTrait' did not match any trait names")
  )
  expect_error(
    pullQtlAf(x$pop, trait = "noTrait", simParam = x$SP),
    paste0("'noTrait' did not match any trait names")
  )
})

test_that("pullMarker returns the same values as the full site set", {
  set.seed(5004)
  x = knownPop(nInd = 70L, nPerChr = c(170L, 130L))

  G = pullSegSiteGeno(x$pop, simParam = x$SP)
  H = pullSegSiteHaplo(x$pop, simParam = x$SP)

  markers = paste0("M", sample(300, 10))
  markerG = pullMarkerGeno(x$pop, markers = markers, simParam = x$SP)
  expect_equal(markerG, pullMarkerGeno(x$pop, markers = markers,
                                       nThreads = 1L, simParam = x$SP))
  expect_identical(markers, colnames(markerG))
  markerH = pullMarkerHaplo(x$pop, markers = markers, simParam = x$SP)
  expect_equal(markerH, pullMarkerHaplo(x$pop, markers = markers, nThreads = 1L,
                                        haplo = "all", simParam = x$SP))

  expect_true(all(markerG == G[, markers, drop = FALSE]))
  expect_true(all(markerH == H[, markers, drop = FALSE]))

  segAf = pullSegSiteAf(x$pop, simParam = x$SP)
  markerAf = pullMarkerAf(x$pop, markers = markers, simParam = x$SP)
  expect_equal(markerAf, pullMarkerAf(x$pop, markers = markers,
                                      nThreads = 1L, simParam = x$SP))
  expect_true(all(markerAf == segAf[, markers, drop = FALSE]))
})

test_that("Sampling individuals for pullAf", {

  set.seed(5004)
  x = knownPop(nInd = 50L, nPerChr = c(170L, 130L))

  # MapPop object in `pullSegSite`
  expect_equal(
    pullSegSiteAf(x$SP$founderPop, simParam = x$SP),
    pullSegSiteAf(x$SP$founderPop, sample = 50L, simParam = x$SP)
  )
  expect_equal(
    pullSegSiteAf(x$SP$founderPop, simParam = x$SP),
    pullSegSiteAf(x$SP$founderPop, sample = 45L, simParam = x$SP),
    tolerance = 0.1
  )
  expect_error(
    pullSegSiteAf(x$SP$founderPop, sample = 51L, simParam = x$SP),
    "`pop` has fewer individuals than `sample`"
  )

  # MapPop object in `pullMarker`
  expect_equal(
    pullMarkerAf(x$SP$founderPop, markers = c('M10','M20'), simParam = x$SP),
    pullMarkerAf(x$SP$founderPop, markers = c('M10','M20'), sample = 50L,
                 simParam = x$SP)
  )
  expect_equal(
    pullMarkerAf(x$SP$founderPop, markers = c('M10','M20'), simParam = x$SP),
    pullMarkerAf(x$SP$founderPop, markers = c('M10','M20'), sample = 45L,
                 simParam = x$SP),
    tolerance = 0.1
  )
  expect_error(
    pullMarkerAf(x$SP$founderPop, markers = c('M10','M20'), sample = 51L,
                 simParam = x$SP),
    "`pop` has fewer individuals than `sample`"
  )

  # Pop object in pullSegSite
  expect_equal(
    pullSegSiteAf(x$pop, simParam = x$SP),
    pullSegSiteAf(x$pop, sample = 50L, simParam = x$SP)
  )
  expect_equal(
    pullSegSiteAf(x$pop, simParam = x$SP),
    pullSegSiteAf(x$pop, sample = 45L, simParam = x$SP),
    tolerance = 0.1
  )
  expect_error(
    pullSegSiteAf(x$pop, sample = 51L, simParam = x$SP),
    "`pop` has fewer individuals than `sample`"
  )

  # Pop object in pullQtl
  x$SP$resetPed()
  x$SP$addTraitA(10, name="Trait1")
  expect_equal(
    pullQtlAf(x$pop, trait = "Trait1", simParam = x$SP),
    pullQtlAf(x$pop, trait = "Trait1", sample = 50L, simParam = x$SP)
  )
  expect_equal(
    pullQtlAf(x$pop, trait = "Trait1", simParam = x$SP),
    pullQtlAf(x$pop, trait = "Trait1", sample = 45L, simParam = x$SP),
    tolerance = 0.1
  )
  expect_error(
    pullQtlAf(x$pop, trait = "Trait1", sample = 51L, simParam = x$SP),
    "`pop` has fewer individuals than `sample`"
  )

  # Pop object in pullSnp
  x$SP$addSnpChip(nSnpPerChr = c(90L, 70L), name = "160SNPchip")
  expect_equal(
    pullSnpAf(x$pop, snpChip = "160SNPchip", simParam = x$SP),
    pullSnpAf(x$pop, snpChip = "160SNPchip", sample = 50L, simParam = x$SP)
  )
  expect_equal(
    pullSnpAf(x$pop, snpChip = "160SNPchip", simParam = x$SP),
    pullSnpAf(x$pop, snpChip = "160SNPchip", sample = 45L, simParam = x$SP),
    tolerance = 0.1
  )
  expect_error(
    pullSnpAf(x$pop, snpChip = "160SNPchip", sample = 51L, simParam = x$SP),
    "`pop` has fewer individuals than `sample`"
  )

  # Pop object in pullMarker
  expect_equal(
    pullMarkerAf(x$pop, markers = c("M25","M1"), simParam = x$SP),
    pullMarkerAf(x$pop, markers = c("M25","M1"), sample = 50L, simParam = x$SP)
  )
  expect_equal(
    pullMarkerAf(x$pop, markers = c("M25","M1"), simParam = x$SP),
    pullMarkerAf(x$pop, markers = c("M25","M1"), sample = 45L, simParam = x$SP),
    tolerance = 0.1
  )
  expect_error(
    pullMarkerAf(x$pop, markers = c("M25","M1"), sample = 51L, simParam = x$SP),
    "`pop` has fewer individuals than `sample`"
  )

  # Wrong `sample`
  for (bad in list(0, -1, NA_integer_, 1.5, c(1, 2), "1", TRUE, Inf)) {
    expect_error(
        pullSegSiteAf(x$SP$founderPop, sample = bad, simParam = x$SP),
        "`sample` must be a positive non-NA integer value."
      )
    expect_error(
        pullMarkerAf(x$SP$founderPop, markers = c('M10','M20'), sample = bad,
                     simParam = x$SP),
        "`sample` must be a positive non-NA integer value."
      )
    expect_error(
        pullSegSiteAf(x$pop, sample = bad, simParam = x$SP),
        "`sample` must be a positive non-NA integer value."
      )
    expect_error(
        pullQtlAf(x$pop, trait = "Trait1", sample = bad, simParam = x$SP),
        "`sample` must be a positive non-NA integer value."
      )
    expect_error(
        pullSnpAf(x$pop, snpChip = "160SNPchip", sample = bad, simParam = x$SP),
        "`sample` must be a positive non-NA integer value."
      )
    expect_error(
        pullMarkerAf(x$pop, markers = c("M25","M1"), sample = bad,
                     simParam = x$SP),
        "`sample` must be a positive non-NA integer value."
      )
  }
})

test_that("MultiPop input for pullSegSite", {
  set.seed(5004)

  x = knownPop(nInd = 32L, nPerChr = c(130L, 170L))

  df_split = data.frame(
    level1 = rep(LETTERS[1:4], each = 8),
    level2 = rep(c('R1', 'R2'), each = 4, times = 4),
    row.names = x$pop@id
  )
  mp = splitPop(x$pop, by = df_split)

  # Haplotypes
  segSiteHaplo = pullSegSiteHaplo(x$pop, haplo = 1, simParam = x$SP)
  segSiteHaplo_l1 = split(as.data.frame(segSiteHaplo), df_split$level1)
  segSiteHaplo_l1 = lapply(segSiteHaplo_l1, as.matrix)
  segSiteHaplo_l2 = lapply(segSiteHaplo_l1, \(x) {
    split(as.data.frame(x), rep(c('R1', 'R2'), each = 4))
  })
  segSiteHaplo_l2 = lapply(segSiteHaplo_l2, \(x) lapply(x, as.matrix))

  expect_equivalent(
    pullSegSiteHaplo(mp, haplo = 1, simParam = x$SP, simplify = TRUE),
    segSiteHaplo
  )
  expect_equivalent(
    pullSegSiteHaplo(
      mp,
      haplo = 1,
      simParam = x$SP,
      simplify = TRUE,
      level = 2
    ),
    segSiteHaplo_l1
  )
  expect_equivalent(
    pullSegSiteHaplo(mp, haplo = 1, simParam = x$SP),
    segSiteHaplo_l2
  )

  # Genotypes
  segSiteGeno = pullSegSiteGeno(x$pop, simParam = x$SP)
  segSiteGeno_l1 = split(as.data.frame(segSiteGeno), df_split$level1)
  segSiteGeno_l1 = lapply(segSiteGeno_l1, as.matrix)
  segSiteGeno_l2 = lapply(segSiteGeno_l1, \(x) {
    split(as.data.frame(x), rep(c('R1', 'R2'), each = 4))
  })
  segSiteGeno_l2 = lapply(segSiteGeno_l2, \(x) lapply(x, as.matrix))

  expect_equivalent(
    pullSegSiteGeno(mp, simParam = x$SP, simplify = TRUE),
    segSiteGeno
  )
  expect_equivalent(
    pullSegSiteGeno(mp, simParam = x$SP, simplify = TRUE, level = 2),
    segSiteGeno_l1
  )
  expect_equivalent(
    pullSegSiteGeno(mp, simParam = x$SP),
    segSiteGeno_l2
  )

  # Allele frequencies
  segSiteAf = do.call(
    'rbind',
    lapply(segSiteGeno_l2, \(x) {
      do.call('rbind', lapply(x, \(y) colMeans(y) / 2))
    })
  )
  segSiteAf_l1 = lapply(segSiteGeno_l2, \(x) {
    do.call('rbind', lapply(x, \(y) colMeans(y) / 2))
  })
  segSiteAf_l2 = lapply(segSiteGeno_l2, \(x) {
    lapply(x, \(y) as.matrix(colMeans(y) / 2))
  })

  expect_equivalent(
    pullSegSiteAf(mp, simParam = x$SP, simplify = TRUE),
    segSiteAf
  )
  expect_equivalent(
    pullSegSiteAf(mp, simParam = x$SP, simplify = TRUE, level = 2),
    segSiteAf_l1
  )
  expect_equivalent(
    pullSegSiteAf(mp, simParam = x$SP),
    segSiteAf_l2
  )
})
