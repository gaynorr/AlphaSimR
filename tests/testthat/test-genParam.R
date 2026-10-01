context("genParam")

# See also test-addTrait.R & test-importData.R for some other unit tests

# Tests below are against:
# * Falconer (1961) Introduction to quantitative genetics, chapter 7, examples 7.1, 7.2 and 7.4
#   with two different allele frequencies, but assuming Hardy-Weinberg equilibrium
#   proportions, while here we test both HWE proportions and a heterozygote
#   deficit at the same allele frequencies
#   https://archive.org/details/introductiontoqu0000falc/page/112
#   (prefer the book as the primary reference: it covers the broader set of concepts0
# * Falconer (1985), A note on Fisher's average effect and average excess
#   https://doi.org/10.1017/S0016672300022825
#   Table 1; equations 3, 11, 13; discussion of breeding values on p. 346.
# * Additional trait expansions to test AlphaSimR's implementation

local({
  # ---- Trait guide ----
  #
  # The first three traits are each controlled by one different locus:
  #   Trait  Locus   p     q     Observed dosage frequencies (0, 1, 2)   F
  #     1      a   0.9   0.1              0.01, 0.18, 0.81             0
  #     2      b   0.6   0.4              0.16, 0.48, 0.36             0
  #     3      c   0.6   0.4              0.34, 0.12, 0.54             0.75
  #     4      d   0.6   0.4              0.16, 0.48, 0.36             0
  #     5      e   0.6   0.4              0.40, 0.00, 0.60             1
  #     6     f&g  0.6   0.4              HWE at each locus; LE
  #     7     h&i  0.6   0.4              HWE at each locus; partial LD
  #     8     j&k  0.6   0.4              HWE at each locus; perfect LD
  #     9     f&g  as trait 6, with AA effect 2
  #    10     h&i  as trait 7, with AA effect 2
  #    11     j&k  as trait 8, with AA effect 2
  #    12     a&b  p = (0.9, 0.6); HWE marginals, asymmetric LD; no AA
  #    13     a&b  as trait 12, with AA effect 2
  #    14     b&c  p = (0.6, 0.6); F = (0, 0.75); AD, no AA
  #    15     b&c  as trait 14, with AA effect 2
  #    16      l   0.0   1.0              1.00, 0.00, 0.00             undefined
  #    17     b&l  HWE segregating b + fixed l, with AA effect 2
  #    18   f&g,h&i  two AA pairs, combining traits 9 and 10
  #    19   m&n,o&p  unequal a = (1,2,3,4), d = (0.1,0.2,0.3,0.4), AA = (0.05,0)
  #    20     l&b  same effects as trait 17, with fixed l first in the AA pair
  #
  # All have intercept = 10; except trait 19, a = 4 and d = 2 (trait 4 has d = 0)
  # Traits 1 vs 2 isolate allele-frequency differences under HWE.
  # Traits 2 vs 3 isolate departures from HWE at the same allele frequency.
  # Trait 14 combines traits 2 and 3: one HWE and one non-HWE locus.
  # Their existing association is retained; this pair is not at LE.
  # Trait indices refer to this guide. Section headings specify
  # whether calculations use the idealized or observed reference population.

  # ---- Haplotypes and genotypes ----

  # Haplotypes (filled below to get targeted genotype frequencies and LD)
  haplo = matrix(0, nrow = 200, ncol = 3,
                 dimnames = list(NULL, c("a", "b", "c")))

  # Trait 1, Locus a - at Hardy-Weinberg genotype proportions (F = 0)
  # q = 0.1 as in Falconer (1961) examples 7.1 and 7.2
  haplo[1:2, "a"] = 0 # 1 individual with 0/0
  haplo[3:38, "a"] = rep(c(0, 1), times = 18) # 18 individuals with 0/1
  haplo[39:200, "a"] = 1 # 81 individuals with 1/1

  # Trait 2, Locus b - at Hardy-Weinberg genotype proportions (F = 0)
  # q = 0.4 as in Falconer (1961) examples 7.1 and 7.2
  haplo[1:32, "b"] = 0 # 16 individuals with 0/0
  haplo[33:128, "b"] = rep(c(0, 1), times = 48) # 48 individuals with 0/1
  haplo[129:200, "b"] = 1 # 36 individuals with 1/1

  # Trait 3, Locus c - with a heterozygote deficit relative to HWE (F = 0.75)
  # q = 0.4 as in Falconer (1961) examples 7.1 and 7.2
  haplo[1:68, "c"] = 0 # 34 individuals with 0/0
  haplo[69:92, "c"] = rep(c(0, 1), times = 12) # 12 individuals with 0/1
  haplo[93:200, "c"] = 1 # 54 individuals with 1/1

  # Traits 4-13: exact frequencies need 2,500 individuals = 5,000 haplotypes
  # Expand the original population, retaining its representative row indices, and genotype frequencies
  haplo = haplo[rep(seq_len(nrow(haplo)), times = 25), ]
  haplo = cbind(haplo, d = haplo[, "b"], e = rep(c(0, 1), times = c(2000, 3000)))

  # Traits 6&9, loci f&g: independent gametes at p = 0.6 for both loci.
  # Haplotype frequencies (00, 01, 10, 11) = (0.16, 0.24, 0.24, 0.36)
  gametes = rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
  gametesLE = gametes[rep(1:4, times = c(4, 6, 6, 9)), ]
  # cor(gametesLE)

  # Traits 7&10, loci h&i: partial LD, with the same allele frequencies.
  # D = 0.48 - 0.6^2 = 0.12; gametic correlation = 0.12/0.24 = 0.5.
  gametesPartialLD = gametes[rep(1:4, times = c(7, 3, 3, 12)), ]

  # Traits 8&11, loci j&k: perfect coupling LD, with the same allele frequencies
  gametesLD = gametes[rep(1:4, times = c(10, 0, 0, 15)), ]
  pairs = expand.grid(first = 1:25, second = 1:25)
  pairs = pairs[rep(seq_len(nrow(pairs)), times = 4), ]

  # Traits 12&13 are reusing loci a&b

  haploLE = haploPartialLD = haploLD = matrix(0, nrow = 5000, ncol = 2)
  haploLE[seq(from = 1, to = 5000, by = 2), ] = gametesLE[pairs$first, ]
  haploLE[seq(from = 2, to = 5000, by = 2), ] = gametesLE[pairs$second, ]
  haploPartialLD[seq(from = 1, to = 5000, by = 2), ] = gametesPartialLD[pairs$first, ]
  haploPartialLD[seq(from = 2, to = 5000, by = 2), ] = gametesPartialLD[pairs$second, ]
  haploLD[seq(from = 1, to = 5000, by = 2), ] = gametesLD[pairs$first, ]
  haploLD[seq(from = 2, to = 5000, by = 2), ] = gametesLD[pairs$second, ]
  # cor(haploLE)
  # cor(haploPartialLD)
  # cor(haploLD)
  colnames(haploLE) = c("f", "g")
  colnames(haploPartialLD) = c("h", "i")
  colnames(haploLD) = c("j", "k")
  haplo = cbind(haplo, haploLE, haploPartialLD, haploLD)
  # Traits 16&17: fixed locus l, alone and as an AA partner of b.
  haplo = cbind(haplo, l = 0)
  # Trait 19: reuse a,b,c,e genotypes at separately named loci with unequal effects
  haplo = cbind(haplo, m = haplo[, "a"], n = haplo[, "b"],
                o = haplo[, "c"], p = haplo[, "e"])
  # cor(haplo)

  nInd = nrow(haplo) / 2
  nLoc = ncol(haplo)
  locusNames = colnames(haplo)
  # Traits 1-5 each use one locus, in this explicit trait order
  singleLoci = c("a", "b", "c", "d", "e")

  p = colMeans(haplo)
  q = 1 - p
  P_HW = p^2
  H_HW = 2 * p * q
  Q_HW = q^2

  geno = do.call(
    rbind,
    args = lapply(seq(from = 1, to = nrow(haplo), by = 2),
                  FUN = function(i) haplo[i, ] + haplo[i + 1, ])
  )
  colnames(geno) = locusNames
  # cor(geno)
  P = apply(geno, MARGIN = 2, FUN = function(x) sum(x == 2) / length(x))
  H = apply(geno, MARGIN = 2, FUN = function(x) sum(x == 1) / length(x))
  Q = apply(geno, MARGIN = 2, FUN = function(x) sum(x == 0) / length(x))

  # F = (ExpHet - ObsHet) / ExpHet
  # Falconer (1961) Introduction to quantitative genetics
  # https://archive.org/details/introductiontoqu0000falc/page/66
  # page 66, equation 3.15: H/H_HW = 1-F
  F = (H_HW - H) / H_HW
  # F is undefined at fixed locus l: both observed and expected heterozygosity are zero

  # HWE frequencies: Falconer (1961), p. 114, Table 7.1.
  # https://archive.org/details/introductiontoqu0000falc/page/114
  # Non-HWE frequencies: Falconer (1985), p. 339, Table 1.
  # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/26DFA92B3BA3EA92CD76847BFA21C5C8/S0016672300022825a.pdf#page=3
  # Joint-frequency and correlation checks below verify our constructed fixture.
  test_that("genParam fixture has the expected allele and genotype frequencies", {
    # Loci a-e control traits 1-5; f&g control traits 6&9; h&i control traits 7&10; j&k control traits 8&11.
    # Traits 12&13 also use loci a&b, with different allele frequencies.
    # Expected vectors below follow locus order a-l; m-p copy a,b,c,e
    expect_identical(colnames(geno), locusNames)
    for (quantity in list(p, q, P_HW, H_HW, Q_HW, P, H, Q, F)) {
      expect_identical(names(quantity), locusNames)
    }
    expect_equal(unname(p[letters[1:12]]), c(0.9, rep(0.6, times = 10), 0), tolerance = 1e-6)
    expect_equal(unname(q[letters[1:12]]), c(0.1, rep(0.4, times = 10), 1), tolerance = 1e-6)
    # Idealized genotype frequencies at the observed allele frequencies.
    expect_equal(unname(P_HW[letters[1:12]]), c(0.81, rep(0.36, times = 10), 0), tolerance = 1e-6)
    expect_equal(unname(H_HW[letters[1:12]]), c(0.18, rep(0.48, times = 10), 0), tolerance = 1e-6)
    expect_equal(unname(Q_HW[letters[1:12]]), c(0.01, rep(0.16, times = 10), 1), tolerance = 1e-6)
    # Observed genotype frequencies: c has F = 0.75; e is fully homozygous.
    expect_equal(unname(P[letters[1:12]]), c(0.81, 0.36, 0.54, 0.36, 0.6, rep(0.36, times = 6), 0), tolerance = 1e-6)
    expect_equal(unname(H[letters[1:12]]), c(0.18, 0.48, 0.12, 0.48, 0, rep(0.48, times = 6), 0), tolerance = 1e-6)
    expect_equal(unname(Q[letters[1:12]]), c(0.01, 0.16, 0.34, 0.16, 0.4, rep(0.16, times = 6), 1), tolerance = 1e-6)
    expect_equal(unname(F[letters[1:12]]), c(0, 0, 0.75, 0, 1, rep(0, times = 6), NaN), tolerance = 1e-6)

    expect_equal(unname(geno[, c("m", "n", "o", "p")]),
                 unname(geno[, c("a", "b", "c", "e")]))
    for (quantity in list(p, q, P_HW, H_HW, Q_HW, P, H, Q, F)) {
      expect_equal(unname(quantity[c("m", "n", "o", "p")]),
                   unname(quantity[c("a", "b", "c", "e")]), tolerance = 1e-6)
    }

    # Traits 6&9, loci f&g: LE haplotype frequencies equal products of allele frequencies.
    expect_equal(matrix(as.vector(table(haplo[, "f"], haplo[, "g"])), nrow = 2, ncol = 2)/nrow(haplo),
                 unname(outer(X = c(q["f"], p["f"]), Y = c(q["g"], p["g"]))), tolerance = 1e-6)
    # Joint genotype frequencies also equal products of marginal frequencies.
    expect_equal(matrix(as.vector(table(geno[, "f"], geno[, "g"])), nrow = 3, ncol = 3)/nInd,
                 unname(outer(X = c(Q["f"], H["f"], P["f"]), Y = c(Q["g"], H["g"], P["g"]))), tolerance = 1e-6)
    # Traits 7&10, loci h&i: HWE marginals but partial LD, so independence fails.
    expect_failure(expect_equal(
      matrix(as.vector(table(haplo[, "h"], haplo[, "i"])), nrow = 2, ncol = 2)/nrow(haplo),
      unname(outer(X = c(q["h"], p["h"]), Y = c(q["i"], p["i"]))), tolerance = 1e-6
    ))
    expect_failure(expect_equal(
      matrix(as.vector(table(geno[, "h"], geno[, "i"])), nrow = 3, ncol = 3)/nInd,
      unname(outer(X = c(Q["h"], H["h"], P["h"]), Y = c(Q["i"], H["i"], P["i"]))), tolerance = 1e-6
    ))
    # Traits 8&11, loci j&k: independence also fails under perfect LD.
    expect_failure(expect_equal(
      matrix(as.vector(table(haplo[, "j"], haplo[, "k"])), nrow = 2, ncol = 2)/nrow(haplo),
      unname(outer(X = c(q["j"], p["j"]), Y = c(q["k"], p["k"]))), tolerance = 1e-6
    ))
    expect_failure(expect_equal(
      matrix(as.vector(table(geno[, "j"], geno[, "k"])), nrow = 3, ncol = 3)/nInd,
      unname(outer(X = c(Q["j"], H["j"], P["j"]), Y = c(Q["k"], H["k"], P["k"]))), tolerance = 1e-6
    ))
    # ... but we know it's perfect corr
    expect_equal(cor(haplo[, "j"], haplo[, "k"]), 1, tolerance = 1e-6)
    expect_equal(cor(geno[, "j"], geno[, "k"]), 1, tolerance = 1e-6)
  })

  # ---- AlphaSimR setup of founders ----

  genMap = data.frame(
    markerName = locusNames,
    chromosome = rep(1, times = nLoc),
    position = seq(from = 0, by = 0.1, length.out = nLoc)
  )

  ped = data.frame(
    id = as.character(1:nInd),
    mother = rep(0, times = nInd),
    father = rep(0, times = nInd)
  )

  founderPop = importHaplo(
    haplo = haplo,
    genMap = genMap,
    ploidy = 2L,
    ped = ped
  )

  SP = SimParam$new(founderPop = founderPop)
  SP$nThreads = 1L

  # ---- AlphaSimR setup of traits ----

  intercept = 10
  # Additive and dominance effects by locus a-p; AA effects by trait 1-20
  a = c(a=4, b=4, c=4, d=4, e=4, f=4, g=4, h=4, i=4, j=4, k=4, l=4, m=1, n=2, o=3, p=4)
  d = c(a=2, b=2, c=2, d=0, e=2, f=2, g=2, h=2, i=2, j=2, k=2, l=2, m=0.1, n=0.2, o=0.3, p=0.4)
  aa = c(0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 2, 0, 2, 0, 2, 2, 0.05, 2)

  # Trait 1: q = 0.1, observed genotype frequencies at HWE.
  SP$importTrait(
    markerNames = "a",
    addEff = a["a"],
    domEff = d["a"],
    intercept = intercept,
    name = "Falconer7.1_q=0.1_HWE"
  )
  # Trait 2: q = 0.4, observed genotype frequencies at HWE.
  SP$importTrait(
    markerNames = "b",
    addEff = a["b"],
    domEff = d["b"],
    intercept = intercept,
    name = "Falconer7.1_q=0.4_HWE"
  )
  # Trait 3: q = 0.4, observed heterozygote deficit (F = 0.75).
  SP$importTrait(
    markerNames = "c",
    addEff = a["c"],
    domEff = d["c"],
    intercept = intercept,
    name = "Falconer7.1_q=0.4_not_HWE"
  )
  # Trait 4: locus d copies b, but has no dominance.
  SP$importTrait(
    markerNames = "d",
    addEff = a["d"],
    domEff = d["d"],
    intercept = intercept,
    name = "A_HWE"
  )
  # Trait 5: locus e is fully homozygous; dominance effects remain present.
  SP$importTrait(
    markerNames = "e",
    addEff = a["e"],
    domEff = d["e"],
    intercept = intercept,
    name = "AD_inbred"
  )
  # Traits 6-8: two loci at LE, partial LD, and perfect LD.
  SP$importTrait(
    markerNames = c("f", "g"),
    addEff = a[c("f", "g")],
    domEff = d[c("f", "g")],
    intercept = intercept,
    name = "AD_LE"
  )
  SP$importTrait(
    markerNames = c("h", "i"),
    addEff = a[c("h", "i")],
    domEff = d[c("h", "i")],
    intercept = intercept,
    name = "AD_partial_LD"
  )
  SP$importTrait(
    markerNames = c("j", "k"),
    addEff = a[c("j", "k")],
    domEff = d[c("j", "k")],
    intercept = intercept,
    name = "AD_perfect_LD"
  )
  # Traits 9-11 overlay AA epistasis on the same three pairs of loci.
  # Using SP$manAddTrait(new("TraitADE", ...) because SP$importTrait() does not support epistatis
  SP$manAddTrait(new("TraitADE", nLoci = 2L, lociPerChr = 2L,
                    lociLoc = match(c("f", "g"), table = locusNames), addEff = a[c("f", "g")], domEff = d[c("f", "g")],
                    epiEff = matrix(c(1, 2, aa[9]), nrow = 1),
                    intercept = intercept, name = "ADE_LE"))
  SP$manAddTrait(new("TraitADE", nLoci = 2L, lociPerChr = 2L,
                    lociLoc = match(c("h", "i"), table = locusNames), addEff = a[c("h", "i")], domEff = d[c("h", "i")],
                    epiEff = matrix(c(1, 2, aa[10]), nrow = 1),
                    intercept = intercept, name = "ADE_partial_LD"))
  SP$manAddTrait(new("TraitADE", nLoci = 2L, lociPerChr = 2L,
                    lociLoc = match(c("j", "k"), table = locusNames), addEff = a[c("j", "k")], domEff = d[c("j", "k")],
                    epiEff = matrix(c(1, 2, aa[11]), nrow = 1),
                    intercept = intercept, name = "ADE_perfect_LD"))
  # Traits 12&13: asymmetric allele frequencies, using the existing a&b genotypes.
  SP$importTrait(
    markerNames = c("a", "b"),
    addEff = a[c("a", "b")],
    domEff = d[c("a", "b")],
    intercept = intercept,
    name = "AD_asymmetric_LD"
  )
  SP$manAddTrait(new("TraitADE", nLoci = 2L, lociPerChr = 2L,
                    lociLoc = match(c("a", "b"), table = locusNames), addEff = a[c("a", "b")], domEff = d[c("a", "b")],
                    epiEff = matrix(c(1, 2, aa[13]), nrow = 1),
                    intercept = intercept, name = "ADE_asymmetric_LD"))
  # Trait 14: combine b (HWE) and c (F = 0.75)
  SP$importTrait(
    markerNames = c("b", "c"),
    addEff = a[c("b", "c")],
    domEff = d[c("b", "c")],
    intercept = intercept,
    name = "AD_mixed_HWE"
  )
  # Trait 15: the AA counterpart of trait 14, with the same mixed HWE marginals
  SP$manAddTrait(new("TraitADE", nLoci = 2L, lociPerChr = 2L,
                    lociLoc = match(c("b", "c"), table = locusNames), addEff = a[c("b", "c")], domEff = d[c("b", "c")],
                    epiEff = matrix(c(1, 2, aa[15]), nrow = 1),
                    intercept = intercept, name = "ADE_mixed_HWE"))
  # Trait 16: fixed locus l, with nonzero input a and d
  SP$importTrait(
    markerNames = "l",
    addEff = a["l"],
    domEff = d["l"],
    intercept = intercept,
    name = "AD_fixed_p0"
  )
  # Trait 17: b plus fixed l, with l second in the AA pair
  SP$manAddTrait(new("TraitADE", nLoci = 2L, lociPerChr = 2L,
                    lociLoc = match(c("b", "l"), table = locusNames), addEff = a[c("b", "l")], domEff = d[c("b", "l")],
                    epiEff = matrix(c(1, 2, aa[17]), nrow = 1),
                    intercept = intercept, name = "ADE_fixed_partner"))
  # Trait 18: two AA pairs, reusing the LE f&g and partial-LD h&i genotypes
  SP$manAddTrait(new("TraitADE", nLoci = 4L, lociPerChr = 4L,
                    lociLoc = match(c("f", "g", "h", "i"), table = locusNames),
                    addEff = a[c("f", "g", "h", "i")], domEff = d[c("f", "g", "h", "i")],
                    epiEff = rbind(c(1, 2, aa[18]), c(3, 4, aa[18])),
                    intercept = intercept, name = "ADE_two_pairs"))
  # Trait 19: unequal locus effects, with one nonzero AA interaction
  SP$manAddTrait(new("TraitADE", nLoci = 4L, lociPerChr = 4L,
                    lociLoc = match(c("m", "n", "o", "p"), table = locusNames),
                    addEff = a[c("m", "n", "o", "p")], domEff = d[c("m", "n", "o", "p")],
                    epiEff = rbind(c(1, 2, aa[19]), c(3, 4, 0)),
                    intercept = intercept, name = "ADE_unequal_effects"))
  # Trait 20: reverse the AA pair so fixed l is first and variable b is second
  # Keep the genomic locus order b,l; pair indices specify the interaction order
  SP$manAddTrait(new("TraitADE", nLoci = 2L, lociPerChr = 2L,
                    lociLoc = match(c("b", "l"), table = locusNames),
                    addEff = a[c("b", "l")], domEff = d[c("b", "l")],
                    epiEff = matrix(c(2, 1, aa[20]), nrow = 1),
                    intercept = intercept, name = "ADE_fixed_first"))
  # TODO: Add Imprinting trait(s)
  # TODO: Add GxE trait(s)
  # SP$traits
  nTrt = length(SP$traits)
  pop = newPop(founderPop, simParam = SP)
  gp = genParam(pop, simParam = SP)

  # Reuse an ADE heterozygote and homozygote for the edge cases in each section
  singleRows = c(which(geno[, "f"] == 1 & geno[, "g"] == 1)[1],
                 which(geno[, "f"] != 1 & geno[, "g"] != 1)[1])
  singleCases = list()
  for (traits in list(9L, c(1L, 4L, 9L, 16L))) {
    singleSP = SimParam$new(founderPop)
    singleSP$nThreads = 1L
    for (trait in traits) singleSP$manAddTrait(SP$traits[[trait]])
    for (row in singleRows) {
      singlePop = newPop(founderPop[row], simParam = singleSP)
      singleCases[[length(singleCases) + 1L]] = list(
        row = row, traits = traits, SP = singleSP, pop = singlePop,
        gp = genParam(singlePop, simParam = singleSP)
      )
    }
  }

  # ---- Explicit thread arguments work ----

  test_that("genParam: explicit thread arguments work", {
    expect_identical(genParam(pop, simParam = SP, nThreads = 1), gp)
    skip_if(getNumThreads() < 2)
    expect_identical(genParam(pop, simParam = SP, nThreads = 2), gp)
    expect_identical(SP$nThreads, 1L)
  })

  # ---- Genetic value ----

  # Falconer (1961), p. 113, Fig. 7.1: midpoint-relative genotype scale
  # https://archive.org/details/introductiontoqu0000falc/page/113
  # At dosage x = 0, 1, 2, additive contributions are -a, 0, +a
  # and dominance contributions are 0, d, 0 (Falconer's -a, d, +a scale).
  # Columns are traits 1-17; all expectations use the input effects.
  myGv_a = matrix(0, nrow = nInd, ncol = nTrt)
  myGv_d = matrix(0, nrow = nInd, ncol = nTrt)
  myGv_aa = matrix(0, nrow = nInd, ncol = nTrt)

  # Single locus traits
  for (trait in 1:5) {
    locus = singleLoci[trait]
    myGv_a[, trait] = a[locus]*(geno[, locus]-1)
    myGv_d[, trait] = d[locus]*(geno[, locus]==1)
  }
  # Trait 4: d = 0. Trait 5: d = 2, but no heterozygotes
  # Both therefore have zero dominance contributions; neither has AA epistasis

  # Multiple locus traits
  # Traits 6&9 use f&g (LE), 7&10 use h&i (partial LD), and
  # 8&11 use j&k (perfect LD). Sum the contributions from both loci
  myGv_a[, 6] = a["f"]*(geno[, "f"]-1) + a["g"]*(geno[, "g"]-1)
  myGv_a[, 7] = a["h"]*(geno[, "h"]-1) + a["i"]*(geno[, "i"]-1)
  myGv_a[, 8] = a["j"]*(geno[, "j"]-1) + a["k"]*(geno[, "k"]-1)
  myGv_d[, 6] = d["f"]*(geno[, "f"]==1) + d["g"]*(geno[, "g"]==1)
  myGv_d[, 7] = d["h"]*(geno[, "h"]==1) + d["i"]*(geno[, "i"]==1)
  myGv_d[, 8] = d["j"]*(geno[, "j"]==1) + d["k"]*(geno[, "k"]==1)

  # Traits 9-11 retain those additive and dominance contributions and add AA
  myGv_a[, 9:11] = myGv_a[, 6:8]
  myGv_d[, 9:11] = myGv_d[, 6:8]
  myGv_aa[, 9] = aa[9]*(geno[, "f"]-1)*(geno[, "g"]-1)
  myGv_aa[, 10] = aa[10]*(geno[, "h"]-1)*(geno[, "i"]-1)
  myGv_aa[, 11] = aa[11]*(geno[, "j"]-1)*(geno[, "k"]-1)

  # Traits 12&13 use a&b with unequal allele frequencies: AD, then ADE
  myGv_a[, 12] = a["a"]*(geno[, "a"]-1) + a["b"]*(geno[, "b"]-1)
  myGv_d[, 12] = d["a"]*(geno[, "a"]==1) + d["b"]*(geno[, "b"]==1)
  myGv_a[, 13] = myGv_a[, 12]
  myGv_d[, 13] = myGv_d[, 12]
  myGv_aa[, 13] = aa[13]*(geno[, "a"]-1)*(geno[, "b"]-1)

  # Trait 14: sum the contributions already calculated for b and c
  # There is no AA effect, so myGv_aa[, 14] remains zero
  myGv_a[, 14] = myGv_a[, 2] + myGv_a[, 3]
  myGv_d[, 14] = myGv_d[, 2] + myGv_d[, 3]

  # Trait 15: retain trait 14's AD contributions and add AA between b and c
  myGv_a[, 15] = myGv_a[, 14]
  myGv_d[, 15] = myGv_d[, 14]
  myGv_aa[, 15] = aa[15]*(geno[, "b"]-1)*(geno[, "c"]-1)

  # Fixed single locus trait 16: constant -a, no heterozygotes or AA
  myGv_a[, 16] = -a["l"]
  # Its dominance and AA contribution columns remain zero

  # Trait 17: reuse b and l's contributions; l has additive dosage -1
  myGv_a[, 17] = myGv_a[, 2] + myGv_a[, 16]
  myGv_d[, 17] = myGv_d[, 2]
  myGv_aa[, 17] = -aa[17]*(geno[, "b"]-1)

  # Trait 18: add both existing pairs, counting the intercept only once below
  myGv_a[, 18] = myGv_a[, 9] + myGv_a[, 10]
  myGv_d[, 18] = myGv_d[, 9] + myGv_d[, 10]
  myGv_aa[, 18] = myGv_aa[, 9] + myGv_aa[, 10]

  # Trait 19: unequal effects at m-p, with AA only between m and n
  myGv_a[, 19] = (geno[, c("m", "n", "o", "p")]-1) %*% a[c("m", "n", "o", "p")]
  myGv_d[, 19] = (geno[, c("m", "n", "o", "p")]==1) %*% d[c("m", "n", "o", "p")]
  myGv_aa[, 19] = aa[19]*(geno[, "m"]-1)*(geno[, "n"]-1)

  # Trait 20: reversing the AA pair leaves every genetic contribution unchanged
  myGv_a[, 20] = myGv_a[, 17]
  myGv_d[, 20] = myGv_d[, 17]
  myGv_aa[, 20] = aa[20]*(geno[, "l"]-1)*(geno[, "b"]-1)

  # Non-additive contributions combine dominance and AA epistasis.
  # Add the intercept once per trait to obtain the genetic value.
  myGv_n = myGv_d + myGv_aa
  myGv = intercept + myGv_a + myGv_n

  test_that("genParam: Genetic value", {
    # Single locus traits
    # Falconer (1961), pp. 113-114, Example 7.1: values 6, 12, 14; midpoint 10.
    # https://archive.org/details/introductiontoqu0000falc/page/114
    expect_equal(unname(gp$gv[, 1:5]), myGv[, 1:5], tolerance = 1e-6)

    # Fixed single locus trait: gv is constant at 10-4 = 6.
    expect_equal(unname(gp$gv[, 16]), myGv[, 16], tolerance = 1e-6)

    # Multiple locus traits
    # Falconer (1961), p. 116: summing genotypic contributions across loci.
    # https://archive.org/details/introductiontoqu0000falc/page/116
    # AD: LE, partial LD, perfect LD, asymmetric a&b, and mixed HWE b&c
    for (trait in c(6, 7, 8, 12, 14)) {
      expect_equal(unname(gp$gv[, trait]), myGv[, trait], tolerance = 1e-6)
    }
    # ADE counterparts: the same loci and effects, plus AA epistasis
    for (trait in c(9, 10, 11, 13, 15, 17, 18, 19, 20)) {
      expect_equal(unname(gp$gv[, trait]), myGv[, trait], tolerance = 1e-6)
    }
    # The public accessor returns the same genetic values for all traits
    expect_equal(dim(gv(pop)), c(nInd, nTrt))
    expect_equal(unname(gv(pop)), myGv, tolerance = 1e-6)
    expect_identical(colnames(gv(pop)), SP$traitNames)
  })

  test_that("genParam: Genetic value - one individual", {
    for (case in singleCases) {
      expected = myGv[case$row, case$traits, drop = FALSE]
      colnames(expected) = case$SP$traitNames
      expect_equal(case$gp$gv, expected, tolerance = 1e-6)
      expect_equal(gv(case$pop), expected, tolerance = 1e-6)
    }
  })

  # ---- Genetic value decomposition with genotypic parameterization ----

  # Midpoint scale: Falconer (1961), p. 113, Fig. 7.1; summation: p. 116
  # https://archive.org/details/introductiontoqu0000falc/page/113
  # https://archive.org/details/introductiontoqu0000falc/page/116
  # This identity tests AlphaSimR's input-effect parameterization; it is not
  # Falconer's population-centered breeding value parameterization G = A + D + ... (p. 126)
  # https://archive.org/details/introductiontoqu0000falc/page/126
  test_that("genParam: Genetic value decomposition with genotypic parameterization", {
    for (trait in 1:nTrt) {
      expect_equal(gp$gv[, trait],
                   gp$gv_mu[trait] + gp$gv_a[, trait] + gp$gv_d[, trait] + gp$gv_aa[, trait],
                   tolerance = 1e-6)
    }
  })

  # ---- Trait intercept ----

  # Falconer (1961), p. 114, Example 7.1: the homozygote midpoint is 10
  # https://archive.org/details/introductiontoqu0000falc/page/114
  # The same intercept is assigned once per trait
  test_that("genParam: Trait intercept", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$gv_mu[trait]), intercept, tolerance = 1e-6)
    }
  })

  test_that("genParam: Trait intercept - one individual", {
    for (case in singleCases) {
      expected = setNames(rep(intercept, length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$gv_mu, expected, tolerance = 1e-6)
    }
  })

  # ---- Additive genotypic effect contributions to genetic values ----

  # Falconer (1961), p. 113, Fig. 7.1: a is half the homozygote difference
  # https://archive.org/details/introductiontoqu0000falc/page/113
  # At each locus, dosage 0, 1, 2 contributes -a, 0, +a.
  test_that("genParam: Additive genotypic effect contributions to genetic values", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$gv_a[, trait]), myGv_a[, trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Additive genotypic effect contributions to genetic values - one individual", {
    for (case in singleCases) {
      expected = myGv_a[case$row, case$traits, drop = FALSE]
      colnames(expected) = case$SP$traitNames
      expect_equal(case$gp$gv_a, expected, tolerance = 1e-6)
    }
  })

  # ---- Dominance genotypic effect contributions to genetic values ----

  # Falconer (1961), p. 113, Fig. 7.1: d is the heterozygote's midpoint offset
  # https://archive.org/details/introductiontoqu0000falc/page/113
  # Only heterozygotes contribute d; traits 4&5 have zero contributions.
  test_that("genParam: Dominance genotypic effect contributions to genetic values", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$gv_d[, trait]), myGv_d[, trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Dominance genotypic effect contributions to genetic values - one individual", {
    for (case in singleCases) {
      expected = myGv_d[case$row, case$traits, drop = FALSE]
      colnames(expected) = case$SP$traitNames
      expect_equal(case$gp$gv_d, expected, tolerance = 1e-6)
    }
  })

  # ---- Additive-by-additive epistatic genotypic effect contributions ----

  # Falconer (1961), p. 116: epistasis means non-additive combination across loci
  # https://archive.org/details/introductiontoqu0000falc/page/116
  # Only traits 9-11, 13, 15&17 have AA effects: aa*(x1-1)*(x2-1).
  test_that("genParam: Additive-by-additive epistatic genotypic effect contributions", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$gv_aa[, trait]), myGv_aa[, trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Additive-by-additive epistatic genotypic effect contributions - one individual", {
    for (case in singleCases) {
      expected = myGv_aa[case$row, case$traits, drop = FALSE]
      colnames(expected) = case$SP$traitNames
      expect_equal(case$gp$gv_aa, expected, tolerance = 1e-6)
    }
  })

  # ---- Non-additive genotypic effect contributions to genetic values ----

  # Non-additive contributions are dominance plus AA epistasis.
  test_that("genParam: Non-additive genotypic effect contributions to genetic values", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$gv_n[, trait]), myGv_n[, trait], tolerance = 1e-6)
    }
    expect_equal(gp$gv_n, gp$gv_d + gp$gv_aa, tolerance = 1e-6)
  })

  test_that("genParam: Non-additive genotypic effect contributions to genetic values - one individual", {
    for (case in singleCases) {
      expected = (myGv_d + myGv_aa)[case$row, case$traits, drop = FALSE]
      colnames(expected) = case$SP$traitNames
      expect_equal(case$gp$gv_n, expected, tolerance = 1e-6)
    }
  })

  # ---- Imprinting genotypic effect contributions to genetic values ----

  test_that("genParam: Imprinting genotypic effect contributions to genetic values", {
    skip("TODO: add tests for this section")
  })

  # ---- Genotype-by-environment effect contributions to genetic values ----

  test_that("genParam: Genotype-by-environment effect contributions to genetic values", {
    skip("TODO: add tests for this section")
  })

  # ---- Genetic value decomposition with breeding value parameterization ----

  # Falconer (1961), p. 122, equation 7.8: centered G = A + D at one locus
  # https://archive.org/details/introductiontoqu0000falc/page/122
  # Multiple-locus interaction deviations extend the decomposition on p. 126
  # https://archive.org/details/introductiontoqu0000falc/page/126
  # Here gp$gv is uncentered, so restore the observed population mean:
  # gv = mu + bv + dd + aa = mu + bv + nd for the AD/ADE traits tested here
  # mu is the population mean, not the input intercept gv_mu
  # bv, dd and aa each have zero observed mean and depend on the reference population

  test_that("genParam: Genetic value decomposition with breeding value parameterization", {
    for (trait in 1:nTrt) {
      expect_equal(gp$gv[, trait],
                   gp$mu[trait] + gp$bv[, trait] + gp$dd[, trait] + gp$aa[, trait],
                   tolerance = 1e-6)
    }
  })

  # ---- Population mean genetic value (idealized population) ----

  # All traits use their observed allele frequencies, but assume HWE and LE
  # At each locus, E(x-1) = p-q and Pr(x==1) = 2pq
  # Falconer's midpoint-relative M = a(p-q) + 2pqd
  # Falconer (1961), p. 115, equation 7.2 (Table 7.1 on p. 114).
  # https://archive.org/details/introductiontoqu0000falc/page/115
  myMeanLocus_HW = a*(p-q) + 2*p*q*d
  myMu_HW = rep(0, times = nTrt)

  # Single locus traits
  myMu_HW[1:5] = intercept + myMeanLocus_HW[c("a", "b", "c", "d", "e")]
  # Trait 5: HWE includes hypothetical heterozygotes despite full observed inbreeding
  # Its mean is 10 + 4*(0.6-0.4) + 2*0.6*0.4*2 = 11.76

  # Multiple locus traits
  # Falconer (1961), p. 116, equation 7.3: sum the single-locus means.
  # https://archive.org/details/introductiontoqu0000falc/page/116
  # Traits 6&9: f&g; 7&10: h&i; 8&11: j&k. Their idealized references all have LE
  myMu_HW[6] = intercept + sum(myMeanLocus_HW[c("f", "g")])
  myMu_HW[7] = intercept + sum(myMeanLocus_HW[c("h", "i")])
  myMu_HW[8] = intercept + sum(myMeanLocus_HW[c("j", "k")])
  # Derived here for AlphaSimR's AA product coding using independence under LE;
  # this AA extension is not part of Falconer's equation 7.3.
  # Under LE, E((x1-1)*(x2-1)) = (p1-q1)*(p2-q2)
  myMu_HW[9] = myMu_HW[6] + aa[9]*(p["f"]-q["f"])*(p["g"]-q["g"])
  myMu_HW[10] = myMu_HW[7] + aa[10]*(p["h"]-q["h"])*(p["i"]-q["i"])
  myMu_HW[11] = myMu_HW[8] + aa[11]*(p["j"]-q["j"])*(p["k"]-q["k"])
  # Traits 12&13: a&b have unequal allele frequencies; use each locus's p and q
  myMu_HW[12] = intercept + sum(myMeanLocus_HW[c("a", "b")])
  myMu_HW[13] = myMu_HW[12] + aa[13]*(p["a"]-q["a"])*(p["b"]-q["b"])

  # Trait 14: both loci use HWE here, giving 10 + 1.76 + 1.76 = 13.52.
  myMu_HW[14] = intercept + sum(myMeanLocus_HW[c("b", "c")])

  # Trait 15: under HWE and LE, AA adds 2*0.2*0.2 = 0.08, giving 13.6.
  myMu_HW[15] = myMu_HW[14] + aa[15]*(p["b"]-q["b"])*(p["c"]-q["c"])

  # Fixed single locus trait 16: HWE also assigns all weight to the fixed genotype.
  # Falconer (1961), p. 115: fixation gives midpoint-relative mean -a or +a.
  # https://archive.org/details/introductiontoqu0000falc/page/115
  myMu_HW[16] = intercept + myMeanLocus_HW["l"]
  # Trait 17: 10 + 1.76 - 4 + 2*0.2*(-1) = 7.36.
  myMu_HW[17] = intercept + sum(myMeanLocus_HW[c("b", "l")]) + aa[17]*(p["b"]-q["b"])*(p["l"]-q["l"])

  myMu_HW[18] = myMu_HW[9] + myMu_HW[10] - intercept

  myMu_HW[19] = intercept + sum(myMeanLocus_HW[c("m", "n", "o", "p")]) +
    aa[19]*(p["m"]-q["m"])*(p["n"]-q["n"])

  myMu_HW[20] = myMu_HW[17]

  test_that("genParam: Population mean genetic value (idealized population)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$mu_HW[trait]), myMu_HW[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Population mean genetic value (idealized population) - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      for (i in seq_along(case$traits)) {
        trait = case$traits[i]
        loci = if (trait == 9L) c("f", "g") else c("a", "d", "l")[match(trait, c(1L, 4L, 16L))]
        dosage = setNames(as.numeric(geno[case$row, loci]), loci)
        h = 2 * (dosage / 2) * (1 - dosage / 2)
        expected[i] = intercept + sum(a[loci] * (dosage - 1) + d[loci] * h)
        if (trait == 9L) expected[i] = expected[i] + aa[trait] * prod(dosage - 1)
      }
      expect_equal(case$gp$mu_HW, expected, tolerance = 1e-6)
    }
  })

  # ---- Population mean genetic value (observed population) ----

  # Falconer (1961), p. 114: frequency-weighted genotype means equal the
  # arithmetic mean over individuals; p. 117 illustrates joint-genotype weighting.
  # https://archive.org/details/introductiontoqu0000falc/page/114
  # https://archive.org/details/introductiontoqu0000falc/page/117
  # Average the already tested genetic values over the actual individuals
  # This retains their observed genotype frequencies and associations between loci
  # Trait 14: 11.76 + 11.04 - 10 = 12.8, counting the intercept once.
  test_that("genParam: Population mean genetic value (observed population)", {
    summary = meanG(pop)
    expect_equal(unname(summary), colMeans(myGv), tolerance = 1e-6)
    expect_identical(names(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$mu[trait]), mean(gp$gv[, trait]), tolerance = 1e-6)
    }
  })

  test_that("genParam: Population mean genetic value (observed population) - one individual", {
    for (case in singleCases) {
      expected = setNames(as.vector(myGv[case$row, case$traits]), case$SP$traitNames)
      expect_equal(case$gp$mu, expected, tolerance = 1e-6)
      expect_equal(meanG(case$pop), expected, tolerance = 1e-6)
    }
  })

  # ---- Average effect of allele substitution (idealized population) ----

  # Diploid AD locus: alpha_HW = a + d(q-p).
  # Falconer (1961), p. 119, equation 7.5; also Falconer (1985), p. 338, equation 3.
  # https://archive.org/details/introductiontoqu0000falc/page/119
  # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/26DFA92B3BA3EA92CD76847BFA21C5C8/S0016672300022825a.pdf#page=2

  # Use HWE at the observed allele frequencies and LE between loci
  myAlphaLocus_HW = a + d*(q-p)
  myAlpha_HW = vector("list", length = nTrt)

  # Single locus traits (1-5)
  for (trait in 1:5) {
    locus = singleLoci[trait]
    myAlpha_HW[[trait]] = myAlphaLocus_HW[locus]
  }
  # Traits 1-3: 2.4, 3.6, 3.6. Trait 4: d = 0, hence alpha_HW = a = 4
  # Trait 5: hypothetical HWE heterozygotes give 3.6 despite observed F = 1

  # Multiple locus traits 6-8: AD at LE, partial LD and perfect LD; all use HWE and LE here
  myAlpha_HW[[6]] = myAlphaLocus_HW[c("f", "g")]
  myAlpha_HW[[7]] = myAlphaLocus_HW[c("h", "i")]
  myAlpha_HW[[8]] = myAlphaLocus_HW[c("j", "k")]

  # Multiple locus traits 9-11: ADE at LE, partial LD and perfect LD; both slopes are 4
  # Derived here by averaging AlphaSimR's AA product over the other locus
  # Related marginal-averaging concept: Falconer (1961), pp. 127-128, Example 7.7
  # https://archive.org/details/introductiontoqu0000falc/page/127
  # AA adds aa*E(x_other-1) = aa*(p_other-q_other) to each locus's slope
  myAlpha_HW[[9]] = myAlpha_HW[[6]] + unname(aa[9]*(p[c("g", "f")]-q[c("g", "f")]))
  myAlpha_HW[[10]] = myAlpha_HW[[7]] + unname(aa[10]*(p[c("i", "h")]-q[c("i", "h")]))
  myAlpha_HW[[11]] = myAlpha_HW[[8]] + unname(aa[11]*(p[c("k", "j")]-q[c("k", "j")]))

  # Multiple locus traits 12&13: asymmetric AD/ADE. Keep locus order a,b and use the OTHER
  # locus's mean for AA: (2.4, 3.6) + 2*(0.2, 0.8) = (2.8, 5.2)
  myAlpha_HW[[12]] = myAlphaLocus_HW[c("a", "b")]
  myAlpha_HW[[13]] = myAlpha_HW[[12]] + unname(aa[13]*(p[c("b", "a")]-q[c("b", "a")]))

  # Multiple locus trait 14: b&c both have p = 0.6, giving (3.6, 3.6).
  # The idealized reference replaces c's non-HWE frequencies with HWE.
  myAlpha_HW[[14]] = myAlphaLocus_HW[c("b", "c")]

  # Trait 15: AA adds the other locus's mean dosage, giving (4, 4).
  myAlpha_HW[[15]] = myAlpha_HW[[14]] + unname(aa[15]*(p[c("c", "b")]-q[c("c", "b")]))

  # Fixed single locus trait 16: Var(x) = 0, so the regression is undefined.
  # genParam returns zero by convention (src/calcGenParam.cpp, divide-by-zero guard).
  # Do not extrapolate a+d(q-p) to fixation: that would give 6 here.
  myAlpha_HW[[16]] = c(l=0)

  # Trait 17: fixed l contributes -aa to b's slope; its own slope is zero.
  myAlpha_HW[[17]] = c(b=unname(myAlpha_HW[[2]]) - aa[17], l=0) # (1.6, 0)

  myAlpha_HW[[18]] = c(myAlpha_HW[[9]], myAlpha_HW[[10]])

  myAlpha_HW[[19]] = myAlphaLocus_HW[c("m", "n", "o", "p")]
  myAlpha_HW[[19]][c("m", "n")] = myAlpha_HW[[19]][c("m", "n")] +
    aa[19]*unname(p[c("n", "m")]-q[c("n", "m")])

  # Slopes retain genomic order b,l even when the AA pair is reversed
  myAlpha_HW[[20]] = myAlpha_HW[[17]]

  test_that("genParam: Average effect of allele substitution (idealized population)", {
    for (trait in 1:nTrt) {
      expect_equal(as.vector(gp$alpha_HW[[trait]]), unname(myAlpha_HW[[trait]]), tolerance = 1e-6)
    }
  })

  test_that("genParam: Average effect of allele substitution (idealized population) - one individual", {
    for (case in singleCases) {
      expected = setNames(vector("list", length(case$traits)), case$SP$traitNames)
      for (i in seq_along(case$traits)) {
        trait = case$traits[i]
        loci = if (trait == 9L) c("f", "g") else c("a", "d", "l")[match(trait, c(1L, 4L, 16L))]
        dosage = setNames(as.numeric(geno[case$row, loci]), loci)
        h = 2 * (dosage / 2) * (1 - dosage / 2)
        slope = a[loci] + d[loci] * (1 - dosage)
        if (trait == 9L) slope = slope + aa[trait] * rev(dosage - 1)
        slope[h == 0] = 0
        expected[[i]] = matrix(slope, ncol = 1)
      }
      expect_equal(case$gp$alpha_HW, expected, tolerance = 1e-6)
    }
  })

  # ---- Average effect of allele substitution (observed population) ----

  # Diploid AD locus: alpha = a + d(q-p)(1-F)/(1+F).
  # Falconer (1985), pp. 340-341, definition C and equations 10-11:
  # regression on dosage using observed frequencies (PDF page 5 = printed p. 341).
  # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/26DFA92B3BA3EA92CD76847BFA21C5C8/S0016672300022825a.pdf#page=5
  myAlphaLocus = a + d*(q-p)*(1-F)/(1+F)
  myAlpha = vector("list", length = nTrt)

  # Single locus traits
  for (trait in 1:5) {
    locus = singleLoci[trait]
    myAlpha[[trait]] = myAlphaLocus[locus]
  }
  # Trait 3: F = 0.75 gives 4 - 0.4/7 = 3.942857.
  # Trait 4: d = 0; trait 5: F = 1. Both have alpha = a = 4.
  # Falconer (1985), p. 342: these two limiting cases.
  # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/26DFA92B3BA3EA92CD76847BFA21C5C8/S0016672300022825a.pdf#page=6

  # Multiple locus traits
  # genParam uses marginal genotype frequencies, removing associations between
  # loci for this calculation. Regressing total observed gv on one locus under
  # LD also captures associated effects at the other locus.
  myAlpha[[6]] = myAlphaLocus[c("f", "g")]
  myAlpha[[7]] = myAlphaLocus[c("h", "i")]
  myAlpha[[8]] = myAlphaLocus[c("j", "k")]

  # Multiple locus traits 9-11: ADE at LE, partial LD and perfect LD.
  # AA still adds aa*(p_other-q_other): mean dosage is unchanged by F.
  myAlpha[[9]] = myAlpha[[6]] + unname(aa[9]*(p[c("g", "f")]-q[c("g", "f")]))
  myAlpha[[10]] = myAlpha[[7]] + unname(aa[10]*(p[c("i", "h")]-q[c("i", "h")]))
  myAlpha[[11]] = myAlpha[[8]] + unname(aa[11]*(p[c("k", "j")]-q[c("k", "j")]))

  # Multiple locus traits 12&13: asymmetric AD/ADE.
  myAlpha[[12]] = myAlphaLocus[c("a", "b")]
  myAlpha[[13]] = myAlpha[[12]] + unname(aa[13]*(p[c("b", "a")]-q[c("b", "a")]))
  # Traits 6-13 have HWE marginals, so their alpha equals alpha_HW
  # despite LE / partial LD / perfect LD / asymmetric LD.

  # Multiple locus trait 14: b&c have F = (0, 0.75), giving (3.6, 3.942857).
  # Only c's slope differs from its HWE reference; use each locus's own F.
  myAlpha[[14]] = myAlphaLocus[c("b", "c")]

  # Trait 15: AA adds 0.4 to each slope, giving (4, 4.342857).
  # Inbreeding changes the dominance correction, but not mean dosage in the AA term.
  myAlpha[[15]] = myAlpha[[14]] + unname(aa[15]*(p[c("c", "b")]-q[c("c", "b")]))

  # Fixed single locus trait 16: the same zero-variance convention applies.
  # F and cov(gv, x)/var(x) are undefined, so neither is used for these expectations.
  myAlpha[[16]] = c(l=0)

  # Trait 17: b has HWE frequencies, and fixed l shifts its slope by -aa.
  myAlpha[[17]] = c(b=unname(myAlpha[[2]]) - aa[17], l=0) # (1.6, 0)

  myAlpha[[18]] = c(myAlpha[[9]], myAlpha[[10]])

  myAlpha[[19]] = myAlphaLocus[c("m", "n", "o", "p")]
  myAlpha[[19]][c("m", "n")] = myAlpha[[19]][c("m", "n")] +
    aa[19]*unname(p[c("n", "m")]-q[c("n", "m")])

  myAlpha[[20]] = myAlpha[[17]]

  test_that("genParam: Average effect of allele substitution (observed population)", {
    for (trait in 1:nTrt) {
      expect_equal(as.vector(gp$alpha[[trait]]), unname(myAlpha[[trait]]), tolerance = 1e-6)
    }
    # Independent regression check for every single-locus trait, using input gv.
    for (trait in 1:5) {
      locus = singleLoci[trait]
      expect_equal(as.vector(gp$alpha[[trait]]),
                   cov(myGv[, trait], geno[, locus])/var(geno[, locus]), tolerance = 1e-6)
    }
    expect_equal(unname(unlist(myAlpha[1:5])), c(2.4, 3.6, 4-0.4/7, 4, 4), tolerance = 1e-6)
  })

  test_that("genParam: Average effect of allele substitution (observed population) - one individual", {
    for (case in singleCases) {
      expected = setNames(vector("list", length(case$traits)), case$SP$traitNames)
      for (i in seq_along(case$traits)) {
        trait = case$traits[i]
        loci = if (trait == 9L) c("f", "g") else c("a", "d", "l")[match(trait, c(1L, 4L, 16L))]
        dosage = setNames(as.numeric(geno[case$row, loci]), loci)
        h = 2 * (dosage / 2) * (1 - dosage / 2)
        # No observed dosage variation means an undefined slope, returned as zero.
        expected[[i]] = matrix(0, length(loci), 1)
      }
      expect_equal(case$gp$alpha, expected, tolerance = 1e-6)
    }
  })

  # ---- Average effects versus observed two-locus regressions ----

  # There is a distinction between AlphaSimR's marginal-reference alpha and
  # regression coefficients estimated from the observed joint genotype distribution.
  # See vignettes/QuanGen.Rmd: "Interpretation under linkage disequilibrium"
  # and "Relationship to classical theory" in vignette("QuanGen", package = "AlphaSimR").
  # AlphaSimR averages over the partner locus using its marginal frequencies.
  # A separate regression of total GV on one dosage also captures associated
  # effects at the partner locus; those slopes cannot simply be added to form BV.
  # Joint regression fits both dosages together and gives the observed additive
  # least-squares projection, which need not equal AlphaSimR's parameterization.

  # Traits 7 and 10 share the same partial-LD genotypes, without and with AA
  # Trait 9 has the same effects as trait 10, but its loci are independent
  # Expected joint slopes: AD/partial LD agrees, ADE/LE agrees, ADE/partial LD differs
  regressionTraits = c(7, 9, 10)
  myJointAlpha = mySeparateAlpha = myReferenceResidualCov = vector("list", length = length(regressionTraits))
  # These fixtures pair gametes independently, so dosage covariance is twice
  # gametic D and each dosage variance is 2*p*q
  # Under this construction dominance residuals are uncorrelated with both dosages
  # The AA residual gives Cov(x_i, N) = aa*2*D*(1-2*p_i)
  # Thus Cov(X, G) = Var(X)*alpha + Cov(X, N)
  # Separate slopes divide by each dosage variance; joint slopes remove the
  # associations between predictors by solving the two-locus normal equations
  for (case in seq_along(regressionTraits)) {
    trait = regressionTraits[case]
    loci = names(myAlpha[[trait]])
    gameticD = mean(haplo[, loci[1]]*haplo[, loci[2]]) - prod(p[loci])
    dosageCov = matrix(2*gameticD, nrow = length(loci), ncol = length(loci),
                       dimnames = list(loci, loci))
    diag(dosageCov) = 2*p[loci]*q[loci]
    myReferenceResidualCov[[case]] = aa[trait]*2*gameticD*(1-2*p[loci])
    myJointAlpha[[case]] = myAlpha[[trait]] +
      as.vector(solve(dosageCov, b = myReferenceResidualCov[[case]]))
    mySeparateAlpha[[case]] = (as.vector(dosageCov %*% myAlpha[[trait]]) +
      myReferenceResidualCov[[case]]) / diag(dosageCov)
  }

  jointAlpha = separateAlpha = referenceResidualCov = jointResidualCov = vector("list", length = length(regressionTraits))
  jointRank = integer(length(regressionTraits))
  for (case in seq_along(regressionTraits)) {
    trait = regressionTraits[case]
    loci = names(myAlpha[[trait]])
    dosages = sweep(geno[, loci, drop = FALSE], MARGIN = 2, STATS = 2*p[loci], FUN = "-")
    values = myGv[, trait] - mean(myGv[, trait])
    # Centered response and predictors require no intercept
    # Repeated individuals supply the observed genotype weights, so this equals
    # weighted least squares on genotype classes with their observed frequencies
    fit = lm.fit(dosages, y = values)
    jointAlpha[[case]] = fit$coefficients
    jointRank[case] = fit$rank
    separateAlpha[[case]] = colMeans(dosages*values)/colMeans(dosages^2)
    referenceResidual = values - as.vector(dosages %*% myAlpha[[trait]])
    referenceResidualCov[[case]] = colMeans(dosages*referenceResidual)
    jointResidualCov[[case]] = colMeans(dosages*fit$residuals)
  }

  # Caveats: agreement for trait 7 relies on its particular joint distribution.
  # HWE marginals alone do not guarantee agreement for arbitrary two-locus samples.
  # Perfect LD and fixed predictors cannot identify two separate joint slopes,
  # so traits 8, 11 and 17 are deliberately excluded from this comparison.
  # This is a dosage-only projection, not a full regression including dominance
  # and interaction predictors, nor a replacement for AlphaSimR's decomposition.
  # Nonzero residual covariances are intentional here; the later variance tests
  # must retain cross-component terms to reconstruct the same total genetic variance.

  test_that("genParam: Average effects versus observed two-locus regressions", {
    for (case in seq_along(regressionTraits)) {
      expect_equal(jointRank[case], 2)
      expect_equal(jointAlpha[[case]], myJointAlpha[[case]], tolerance = 1e-6)
      expect_equal(separateAlpha[[case]], mySeparateAlpha[[case]], tolerance = 1e-6)
      expect_equal(referenceResidualCov[[case]], myReferenceResidualCov[[case]], tolerance = 1e-6)
      expect_equal(unname(jointResidualCov[[case]]), c(0, 0), tolerance = 1e-6)
    }
    for (case in 1:2) {
      trait = regressionTraits[case]
      expect_equal(as.vector(gp$alpha[[trait]]), unname(jointAlpha[[case]]), tolerance = 1e-6)
    }
    # Trait 10 retains the marginal-reference alpha rather than the joint fit
    expect_equal(as.vector(gp$alpha[[10]]) - unname(jointAlpha[[3]]),
                 unname(myAlpha[[10]] - myJointAlpha[[3]]), tolerance = 1e-6)
  })

  # ---- Breeding value (idealized population) ----

  # Falconer (1961), p. 121: sum the average effects of the alleles carried
  # https://archive.org/details/introductiontoqu0000falc/page/121
  # For diploids this is (dosage-2p)*alpha_HW at each locus, summed over loci
  # Evaluate the HWE-reference effects on the actual individuals, including non-HWE traits
  # These values still have mean zero because both references use the same p
  x = sweep(geno, MARGIN = 2, STATS = 2*p, FUN = "-")
  myBv_HW = matrix(0, nrow=nInd, ncol=nTrt)

  # Single locus traits 1-5
  for (trait in 1:5) {
    locus = singleLoci[trait]
    myBv_HW[, trait] = x[, locus]*myAlpha_HW[[trait]][locus]
  }

  # Multiple locus AD traits: sum both loci using their own substitution effects
  myBv_HW[, 6] = x[, "f"]*myAlpha_HW[[6]]["f"] + x[, "g"]*myAlpha_HW[[6]]["g"]
  myBv_HW[, 7] = x[, "h"]*myAlpha_HW[[7]]["h"] + x[, "i"]*myAlpha_HW[[7]]["i"]
  myBv_HW[, 8] = x[, "j"]*myAlpha_HW[[8]]["j"] + x[, "k"]*myAlpha_HW[[8]]["k"]
  myBv_HW[, 12] = x[, "a"]*myAlpha_HW[[12]]["a"] + x[, "b"]*myAlpha_HW[[12]]["b"]
  myBv_HW[, 14] = x[, "b"]*myAlpha_HW[[14]]["b"] + x[, "c"]*myAlpha_HW[[14]]["c"]

  # ADE counterparts: the same dosage sums, with AA already included in alpha_HW
  myBv_HW[, 9] = x[, "f"]*myAlpha_HW[[9]]["f"] + x[, "g"]*myAlpha_HW[[9]]["g"]
  myBv_HW[, 10] = x[, "h"]*myAlpha_HW[[10]]["h"] + x[, "i"]*myAlpha_HW[[10]]["i"]
  myBv_HW[, 11] = x[, "j"]*myAlpha_HW[[11]]["j"] + x[, "k"]*myAlpha_HW[[11]]["k"]
  myBv_HW[, 13] = x[, "a"]*myAlpha_HW[[13]]["a"] + x[, "b"]*myAlpha_HW[[13]]["b"]
  myBv_HW[, 15] = x[, "b"]*myAlpha_HW[[15]]["b"] + x[, "c"]*myAlpha_HW[[15]]["c"]

  # Fixed locus l contributes zero centered dosage in traits 16&17
  # Trait 16 remains zero; in trait 17 AA shifts b's slope to 1.6
  myBv_HW[, 17] = x[, "b"]*myAlpha_HW[[17]]["b"]

  myBv_HW[, 18] = myBv_HW[, 9] + myBv_HW[, 10]

  myBv_HW[, 19] = x[, c("m", "n", "o", "p")] %*% myAlpha_HW[[19]]
  myBv_HW[, 20] = myBv_HW[, 17]

  # genParam does not return bv_HW, so construct it from its alpha_HW output
  # Locus names in the expectations specify the column order for each trait
  bv_HW = matrix(0, nrow=nInd, ncol=nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha_HW[[trait]])
    bv_HW[, trait] = x[, loci, drop=FALSE] %*% gp$alpha_HW[[trait]]
  }

  test_that("genParam: Breeding value (idealized population)", {
    for (trait in 1:nTrt) {
      expect_equal(bv_HW[, trait], myBv_HW[, trait], tolerance=1e-6)
      expect_equal(mean(bv_HW[, trait]), 0, tolerance=1e-6)
    }
    # The above tests all traits; below we just test against an published example
    # Falconer (1961), p. 121, Example 7.5: values at dosage 0, 1, 2
    expect_equal(bv_HW[match(0:2, table = geno[, "a"]), 1], c(-4.32, -1.92, 0.48), tolerance=1e-6)
    expect_equal(bv_HW[match(0:2, table = geno[, "b"]), 2], c(-4.32, -0.72, 2.88), tolerance=1e-6)
    # Same p and HWE slope for c, despite its heterozygote deficit
    expect_equal(bv_HW[match(0:2, table = geno[, "c"]), 3], c(-4.32, -0.72, 2.88), tolerance=1e-6)
    # Fully inbred e has only the two homozygotes; fixed l has no variation
    expect_equal(bv_HW[match(c(0, 2), table = geno[, "e"]), 5], c(-4.32, 2.88), tolerance=1e-6)
    expect_equal(bv_HW[, 16], rep(0, times = nInd), tolerance=1e-6)
  })

  # ---- Breeding value (observed population) ----

  # Use the observed-frequency alpha in the same centered dosage sum
  # Falconer (1985), p. 341: additive value from regression on dosage
  # https://www.cambridge.org/core/services/aop-cambridge-core/content/view/26DFA92B3BA3EA92CD76847BFA21C5C8/S0016672300022825a.pdf#page=5
  # Under non-random mating this need not have the usual progeny interpretation
  myBv = matrix(0, nrow = nInd, ncol = nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha[[trait]])
    myBv[, trait] = x[, loci, drop = FALSE] %*% myAlpha[[trait]]
  }

  test_that("genParam: Breeding value (observed population)", {
    summary = bv(pop, simParam = SP)
    expect_equal(unname(summary), myBv, tolerance = 1e-6)
    expect_identical(colnames(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$bv[, trait]), myBv[, trait], tolerance=1e-6)
      expect_equal(mean(gp$bv[, trait]), 0, tolerance=1e-6)
    }
  })

  test_that("genParam: Breeding value (observed population) - one individual", {
    for (case in singleCases) {
      # Centering on a single individual leaves no observed deviation.
      expected = matrix(0, 1, length(case$traits),
                        dimnames = list(NULL, case$SP$traitNames))
      expect_equal(case$gp$bv, expected, tolerance = 1e-6)
      expect_equal(bv(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Inbreeding depression load (observed population) ----

  # A locally-derived quantity, not a genParam output
  # IDL = sum((dosage-2p)*(alpha-alpha_HW)) = BV_observed - BV_HW
  # This centered individual quantity is distinct from population mean depression
  # Antonios et al. (2025), pp. 2-3: individual genetic inbreeding load
  # https://doi.org/10.1186/s12711-024-00945-z

  # Expected alpha difference: subtract a+d(q-p) from a+d(q-p)*(1-F)/(1+F)
  # This gives -2F*d(q-p)/(1+F)

  myIdl = matrix(0, nrow = nInd, ncol = nTrt)

  # Single locus traits 1-5: only non-HWE loci c and e have non-zero IDL
  # Traits 1&2 have F = 0; trait 4 also has d = 0
  myIdl[, 3] = x[, "c"]*(-2*F["c"]*d["c"]*(q["c"]-p["c"])/(1+F["c"]))
  myIdl[, 5] = x[, "e"]*(-2*F["e"]*d["e"]*(q["e"]-p["e"])/(1+F["e"]))

  # Traits 6-13: HWE marginals give zero IDL, with or without LD and AA
  # Traits 14&15: b contributes zero, so reuse the correction for c
  # AA cancels between the two alphas because mean dosage is unchanged
  myIdl[, 14] = myIdl[, 3]
  myIdl[, 15] = myIdl[, 3]

  # Traits 16&17: fixed l has zero centered dosage and b has HWE frequencies
  # Both remain zero without evaluating the undefined F at l

  # Trait 19: its own dominance effects determine the non-HWE correction
  myIdl[, 19] = x[, c("m", "n", "o", "p")] %*%
    (-2*F[c("m", "n", "o", "p")]*d[c("m", "n", "o", "p")]*(q[c("m", "n", "o", "p")]-p[c("m", "n", "o", "p")])/(1+F[c("m", "n", "o", "p")]))

  # Compute IDL from the returned alphas, using the same locus order as BV
  alpha_diff = Map(`-`, gp$alpha, gp$alpha_HW)
  idl = matrix(0, nrow = nInd, ncol = nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha[[trait]])
    idl[, trait] = x[, loci, drop = FALSE] %*% alpha_diff[[trait]]
  }

  test_that("genParam: Inbreeding depression load", {
    for (trait in 1:nTrt) {
      expect_equal(idl[, trait], myIdl[, trait], tolerance=1e-6)
      expect_equal(unname(gp$bv[, trait]) - bv_HW[, trait], myIdl[, trait], tolerance=1e-6)
      expect_equal(mean(idl[, trait]), 0, tolerance=1e-6)
    }
  })

  # ---- Dominance deviation (idealized population) ----

  # Falconer (1961), p. 123: deviations at dosages 0, 1, 2
  # https://archive.org/details/introductiontoqu0000falc/page/123
  # DD_HW = c(-2*p^2*d, 2*p*q*d, -2*q^2*d)
  # Sum the locus residuals using each trait's own named dominance effects
  myDd_HW = matrix(0, nrow = nInd, ncol = nTrt)
  myMeanDd_HW = numeric(nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha[[trait]])
    for (locus in loci) {
      deviations = c(-2*p[locus]^2, 2*p[locus]*q[locus], -2*q[locus]^2) * d[locus]
      myDd_HW[, trait] = myDd_HW[, trait] + deviations[geno[, locus]+1]
    }
    # Mean over actual individuals, not over hypothetical HWE frequencies
    myMeanDd_HW[trait] = sum(d[loci]*(H[loci]-2*p[loci]*q[loci]))
  }
  # AA changes the linear slope but leaves these dominance residuals unchanged
  # Traits 3, 14&15 have mean -0.72; fully inbred trait 5 has mean -0.96
  # Trait 19 uses its own dominance effects in the mean correction above
  # Remaining means are zero, including fixed-locus traits 16, 17&20

  # genParam does not return dd_HW
  # Undo observed centering and the alpha correction in its returned dd
  dd_HW = matrix(0, nrow = nInd, ncol = nTrt)
  for (trait in 1:nTrt) {
    dd_HW[, trait] = gp$dd[, trait] + idl[, trait] + myMeanDd_HW[trait]
  }

  test_that("genParam: Dominance deviation (idealized population)", {
    for (trait in 1:nTrt) {
      expect_equal(dd_HW[, trait], myDd_HW[, trait], tolerance = 1e-6)
      expect_equal(mean(dd_HW[, trait]), myMeanDd_HW[trait], tolerance = 1e-6)
    }
    # Falconer (1961), Example 7.6 (starts p. 122, numerical table p. 123)
    # https://archive.org/details/introductiontoqu0000falc/page/123
    # Pygmy-gene example: trait 1 uses q = 0.1, trait 2 uses q = 0.4
    # Reverse the book's ++, +pg, pgpg order to dosages 0, 1, 2 (pgpg, +pg, ++)
    expect_equal(dd_HW[match(0:2, table = geno[, "a"]), 1], c(-3.24, 0.36, -0.04), tolerance = 1e-6)
    expect_equal(dd_HW[match(0:2, table = geno[, "b"]), 2], c(-1.44, 0.96, -0.64), tolerance = 1e-6)
  })

  # ---- Dominance deviation (observed population) ----

  # Center the HWE-reference residuals on the observed population and remove
  # the extra linear contribution already assigned to observed BV through IDL
  # DD = DD_HW - mean_observed(DD_HW) - IDL
  myDd = matrix(0, nrow = nInd, ncol = nTrt)
  for (trait in 1:nTrt) {
    myDd[, trait] = myDd_HW[, trait] - myMeanDd_HW[trait] - myIdl[, trait]
  }

  test_that("genParam: Dominance deviation (observed population)", {
    summary = dd(pop, simParam = SP)
    expect_equal(unname(summary), myDd, tolerance = 1e-6)
    expect_identical(colnames(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$dd[, trait]), myDd[, trait], tolerance = 1e-6)
      expect_equal(mean(gp$dd[, trait]), 0, tolerance = 1e-6)
    }
    # Trait 4 has d = 0; trait 5 has only homozygotes; trait 16 is fixed
    # Each has zero observed DD
    for (trait in c(4, 5, 16)) {
      expect_equal(unname(gp$dd[, trait]), rep(0, times = nInd), tolerance = 1e-6)
    }
  })

  test_that("genParam: Dominance deviation (observed population) - one individual", {
    for (case in singleCases) {
      # Centering on a single individual leaves no observed deviation.
      expected = matrix(0, 1, length(case$traits),
                        dimnames = list(NULL, case$SP$traitNames))
      expect_equal(case$gp$dd, expected, tolerance = 1e-6)
      expect_equal(dd(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Additive-by-additive epistatic deviation (idealized population) ----

  # Falconer (1961), pp. 126-128: epistatic residual after marginal effects
  # https://archive.org/details/introductiontoqu0000falc/page/126
  # https://archive.org/details/introductiontoqu0000falc/page/128
  # For our AA coding, removing the mean and both linear terms gives aa*x1*x2
  # Here x is centered dosage; the idealized reference assumes HWE and LE
  myAa_HW = matrix(0, nrow = nInd, ncol = nTrt)
  myAa_HW[, 9] = aa[9]*x[, "f"]*x[, "g"]
  myAa_HW[, 10] = aa[10]*x[, "h"]*x[, "i"]
  myAa_HW[, 11] = aa[11]*x[, "j"]*x[, "k"]
  myAa_HW[, 13] = aa[13]*x[, "a"]*x[, "b"]
  myAa_HW[, 15] = aa[15]*x[, "b"]*x[, "c"]
  myAa_HW[, 17] = aa[17]*x[, "b"]*x[, "l"]
  myAa_HW[, 18] = myAa_HW[, 9] + myAa_HW[, 10]
  myAa_HW[, 19] = aa[19]*x[, "m"]*x[, "n"]
  myAa_HW[, 20] = aa[20]*x[, "l"]*x[, "b"]
  # Traits without AA remain zero; fixed l also makes trait 17's residual zero
  # Under LE the product has mean zero; observed LD can give a nonzero mean
  myMeanAa_HW = colMeans(myAa_HW)

  # genParam does not return aa_HW, so reconstruct the remaining deviation
  aa_HW = matrix(0, nrow = nInd, ncol = nTrt)
  for (trait in 1:nTrt) {
    aa_HW[, trait] = gp$gv[, trait] - gp$mu_HW[trait] - bv_HW[, trait] - dd_HW[, trait]
  }

  test_that("genParam: Additive-by-additive epistatic deviation (idealized population)", {
    for (trait in 1:nTrt) {
      expect_equal(aa_HW[, trait], myAa_HW[, trait], tolerance = 1e-6)
      expect_equal(mean(aa_HW[, trait]), myMeanAa_HW[trait], tolerance = 1e-6)
    }
  })

  # ---- Additive-by-additive epistatic deviation (observed population) ----

  # Subtract the observed mean of the centered-dosage product
  # Changes in the dominance correction cancel between BV and DD
  myAa = matrix(0, nrow = nInd, ncol = nTrt)
  for (trait in 1:nTrt) {
    myAa[, trait] = myAa_HW[, trait] - myMeanAa_HW[trait]
  }

  test_that("genParam: Additive-by-additive epistatic deviation (observed population)", {
    summary = aa(pop, simParam = SP)
    expect_equal(unname(summary), myAa, tolerance = 1e-6)
    expect_identical(colnames(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$aa[, trait]), myAa[, trait], tolerance = 1e-6)
      expect_equal(mean(gp$aa[, trait]), 0, tolerance = 1e-6)
    }
  })

  test_that("genParam: Additive-by-additive epistatic deviation (observed population) - one individual", {
    for (case in singleCases) {
      # Centering on a single individual leaves no observed deviation.
      expected = matrix(0, 1, length(case$traits),
                        dimnames = list(NULL, case$SP$traitNames))
      expect_equal(case$gp$aa, expected, tolerance = 1e-6)
      expect_equal(aa(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Imprinting deviation (idealized population) ----

  test_that("genParam: Imprinting deviation (idealized population)", {
    skip("Imprinting traits and a genParam imprinting output are not implemented")
  })

  # ---- Imprinting deviation (observed population) ----

  test_that("genParam: Imprinting deviation (observed population)", {
    skip("Imprinting traits and a genParam imprinting output are not implemented")
  })

  # ---- Non-additive deviation (idealized population) ----

  # Combine the dominance and AA deviations already tested above
  myNd_HW = myDd_HW + myAa_HW
  nd_HW = matrix(0, nrow = nInd, ncol = nTrt)
  for (trait in 1:nTrt) {
    nd_HW[, trait] = gp$gv[, trait] - gp$mu_HW[trait] - bv_HW[, trait]
  }

  test_that("genParam: Non-additive deviation (idealized population)", {
    for (trait in 1:nTrt) {
      expect_equal(nd_HW[, trait], myNd_HW[, trait], tolerance = 1e-6)
      expect_equal(mean(nd_HW[, trait]), myMeanDd_HW[trait] + myMeanAa_HW[trait], tolerance = 1e-6)
    }
  })

  # ---- Non-additive deviation (observed population) ----

  myNd = myDd + myAa

  test_that("genParam: Non-additive deviation (observed population)", {
    summary = nd(pop, simParam = SP)
    expect_equal(unname(summary), myNd, tolerance = 1e-6)
    expect_identical(colnames(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$nd[, trait]), myNd[, trait], tolerance = 1e-6)
      expect_equal(mean(gp$nd[, trait]), 0, tolerance = 1e-6)
      # Complete the population-centered decomposition G = mean + BV + DD + AA
      expect_equal(unname(gp$gv[, trait]),
                   gp$mu[trait] + myBv[, trait] + myDd[, trait] + myAa[, trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Non-additive deviation (observed population) - one individual", {
    for (case in singleCases) {
      # Centering on a single individual leaves no observed deviation.
      expected = matrix(0, 1, length(case$traits),
                        dimnames = list(NULL, case$SP$traitNames))
      expect_equal(case$gp$nd, expected, tolerance = 1e-6)
      expect_equal(nd(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Variance and reference population notes ----

  # Population references follow the QG vignette:
  # V_observed = varX = genicVarX + covX_HW + covX_L: observed population
  # V_LE = genicVarX + covX_HW: linkage-equilibrium (LE) population
  # V_HW = genicVarX: idealized population (HWE and LE)
  # covX_HW = V_LE - V_HW
  # covX_L = V_observed - V_LE
  # For the trait we also have covariances between components within each trait
  # covAD_L, covAAA_L, covDAA_L and covAN_L

  # ---- Total genetic variance (observed population) ----

  # Use the already tested genetic values, centered on each trait's observed mean
  # Equal individual weights retain observed genotype frequencies and LD
  # Cov(G_i, G_j) = mean((G_i-mean(G_i))*(G_j-mean(G_j)))
  # On the diagonal this is Var(G_i)
  # We use denominator nInd (via mean()), not nInd-1
  myGvCentered = matrix(0, nrow = nInd, ncol = nTrt)
  for (trait in 1:nTrt) {
    myGvCentered[, trait] = myGv[, trait] - mean(myGv[, trait])
  }
  myVarG = matrix(0, nrow = nTrt, ncol = nTrt)
  for (trait1 in 1:nTrt) {
    for (trait2 in 1:nTrt) {
      myVarG[trait1, trait2] = mean(myGvCentered[, trait1] * myGvCentered[, trait2])
    }
  }
  # Includes every pair of traits, whether they share loci, have associated loci, or
  # have zero variance and zero covariance with every other trait (Trait 16)

  test_that("genParam: Total genetic variance (observed population)", {
    summary = varG(pop)
    expect_equal(unname(summary), myVarG, tolerance = 1e-6)
    expect_identical(dimnames(summary), list(SP$traitNames, SP$traitNames))
    expect_equal(dim(gp$varG), c(nTrt, nTrt))
    expect_equal(rownames(gp$varG), colnames(gp$gv))
    expect_equal(colnames(gp$varG), colnames(gp$gv))
    for (trait1 in 1:nTrt) {
      for (trait2 in 1:nTrt) {
        expect_equal(unname(gp$varG[trait1, trait2]), myVarG[trait1, trait2], tolerance = 1e-6)
      }
    }
    expect_equal(unname(gp$varG[16, ]), rep(0, times = nTrt), tolerance = 1e-6)
    expect_equal(unname(gp$varG[, 16]), rep(0, times = nTrt), tolerance = 1e-6)
  })

  test_that("genParam: Total genetic variance (observed population) - one individual", {
    for (case in singleCases) {
      expected = matrix(0, length(case$traits), length(case$traits),
                        dimnames = list(case$SP$traitNames, case$SP$traitNames))
      expect_equal(case$gp$varG, expected, tolerance = 1e-6)
      expect_equal(varG(case$pop), expected, tolerance = 1e-6)
    }
  })

  # ---- Total genetic variance (linkage-equilibrium population) ----

  # Named AA pairs and their input effects for the reference genotype grids
  myPairLoci = myPairEffects = vector("list", length = nTrt)
  for (trait in c(9, 10, 11, 13, 15, 17)) {
    myPairLoci[[trait]] = matrix(names(myAlpha[[trait]]), nrow = 1)
    myPairEffects[[trait]] = aa[trait]
  }
  myPairLoci[[18]] = rbind(c("f", "g"), c("h", "i"))
  myPairEffects[[18]] = rep(aa[18], times = 2)
  myPairLoci[[19]] = rbind(c("m", "n"), c("o", "p"))
  myPairEffects[[19]] = c(aa[19], 0)
  myPairLoci[[20]] = matrix(c("l", "b"), nrow = 1)
  myPairEffects[[20]] = aa[20]

  # This is genic variance with observed genotype frequencies
  # Keep observed genotype frequencies at each locus, but make loci independent
  # Enumerate 3^nLoci combinations: 3, 9 or 81 rows for one, two or four loci
  # This includes combinations absent from the observed population under LD
  # Store only the genotype grid, weights and genetic values for the two references
  myRefGeno = myRefFreq_LE = myRefFreq_HW = myRefGv = vector("list", length = nTrt)
  myVarG_LE = numeric(nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha[[trait]])
    genotypes = as.matrix(expand.grid(setNames(rep(list(0:2), times = length(loci)), nm = loci)))
    freq_LE = freq_HW = rep(1, times = nrow(genotypes))
    for (locus in loci) {
      freq_LE = freq_LE * c(Q[locus], H[locus], P[locus])[genotypes[, locus]+1]
      freq_HW = freq_HW * c(Q_HW[locus], H_HW[locus], P_HW[locus])[genotypes[, locus]+1]
    }
    # The input genetic-value model is unchanged between reference populations
    values = as.vector(intercept + (genotypes-1) %*% a[loci] + (genotypes==1) %*% d[loci])
    for (pair in seq_along(myPairEffects[[trait]])) {
      partners = myPairLoci[[trait]][pair, ]
      values = values + myPairEffects[[trait]][pair] *
        (genotypes[, partners[1]]-1)*(genotypes[, partners[2]]-1)
    }
    myRefGeno[[trait]] = genotypes
    myRefFreq_LE[[trait]] = freq_LE
    myRefFreq_HW[[trait]] = freq_HW
    myRefGv[[trait]] = values
    # Center on the LE mean, which can differ from the observed mean with AA and LD
    mean_LE = sum(freq_LE * values)
    myVarG_LE[trait] = sum(freq_LE * (values - mean_LE)^2)
  }
  # Single-locus traits retain their observed variance, including fully inbred e
  # Fixed l receives zero weight on its unobserved genotypes in both references

  test_that("genParam: Total genetic variance (LE population)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarG[trait] + gp$covG_HW[trait]),
                   myVarG_LE[trait], tolerance = 1e-6)
    }
  })

  # ---- Total genetic variance (idealized population) ----

  # This is genic variance with HWE genotype frequencies
  # Unlike LE, this reference includes hypothetical heterozygotes for fully inbred e
  # Falconer (1961), p. 136: genetic variance from genotype values and frequencies
  # https://archive.org/details/introductiontoqu0000falc/page/136
  myVarG_HW = numeric(nTrt)
  for (trait in 1:nTrt) {
    values = myRefGv[[trait]]
    frequencies = myRefFreq_HW[[trait]]
    mean_HW = sum(frequencies * values)
    myVarG_HW[trait] = sum(frequencies * (values - mean_HW)^2)
  }

  test_that("genParam: Total genetic variance (idealized population)", {
    summary = genicVarG(pop, simParam = SP)
    expect_equal(unname(summary), myVarG_HW, tolerance = 1e-6)
    expect_identical(names(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarG[trait]), myVarG_HW[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Total genetic variance (idealized population) - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      for (i in seq_along(case$traits)) {
        trait = case$traits[i]
        loci = if (trait == 9L) c("f", "g") else c("a", "d", "l")[match(trait, c(1L, 4L, 16L))]
        dosage = setNames(as.numeric(geno[case$row, loci]), loci)
        h = 2 * (dosage / 2) * (1 - dosage / 2)
        # Reweight this section's genotype grid using the individual's allele frequencies.
        grid = myRefGeno[[trait]]
        weights = rep(1, nrow(grid))
        for (locus in loci) {
          weights = weights * dbinom(grid[, locus], size = 2, prob = dosage[locus] / 2)
        }
        values = myRefGv[[trait]]
        expected[i] = sum(weights * (values - sum(weights * values))^2)
      }
      expect_equal(case$gp$genicVarG, expected, tolerance = 1e-6)
      expect_equal(genicVarG(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Total genetic variance (reference-population adjustments) ----

  # Test each adjustment against independently calculated reference variances
  # covG_HW changes HWE + LE to observed-marginal LE
  # covG_L changes LE to the observed population, including cross-component terms
  myCovG_HW = myVarG_LE - myVarG_HW
  myCovG_L = diag(myVarG) - myVarG_LE

  test_that("genParam: Total genetic variance (reference-population adjustments)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covG_HW[trait]), myCovG_HW[trait], tolerance = 1e-6)
      expect_equal(unname(gp$covG_L[trait]), myCovG_L[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Total genetic variance (reference-population adjustments) - one individual", {
    for (case in singleCases) {
      # Each observed marginal is fixed at one dosage, so LE variance is also zero.
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covG_HW, -case$gp$genicVarG, tolerance = 1e-6)
      expect_equal(case$gp$covG_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Additive genetic variance (observed population) ----

  # Use the already tested breeding values, and they are centered on zero by definition
  # Mean products give population variances and covariances with denominator nInd
  myVarA = matrix(0, nrow = nTrt, ncol = nTrt)
  for (trait1 in 1:nTrt) {
    for (trait2 in 1:nTrt) {
      myVarA[trait1, trait2] = mean(myBv[, trait1] * myBv[, trait2])
    }
  }

  test_that("genParam: Additive genetic variance (observed population)", {
    summary = varA(pop, simParam = SP)
    expect_equal(unname(summary), myVarA, tolerance = 1e-6)
    expect_identical(dimnames(summary), list(SP$traitNames, SP$traitNames))
    expect_equal(dim(gp$varA), c(nTrt, nTrt))
    expect_equal(rownames(gp$varA), colnames(gp$bv))
    expect_equal(colnames(gp$varA), colnames(gp$bv))
    for (trait1 in 1:nTrt) {
      for (trait2 in 1:nTrt) {
        expect_equal(unname(gp$varA[trait1, trait2]), myVarA[trait1, trait2], tolerance = 1e-6)
      }
    }
    expect_equal(unname(gp$varA[16, ]), rep(0, times = nTrt), tolerance = 1e-6)
    expect_equal(unname(gp$varA[, 16]), rep(0, times = nTrt), tolerance = 1e-6)
  })

  test_that("genParam: Additive genetic variance (observed population) - one individual", {
    for (case in singleCases) {
      expected = matrix(0, length(case$traits), length(case$traits),
                        dimnames = list(case$SP$traitNames, case$SP$traitNames))
      expect_equal(case$gp$varA, expected, tolerance = 1e-6)
      expect_equal(varA(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Additive genetic variance (linkage-equilibrium population) ----

  # Keep observed marginal genotype frequencies and alpha, but remove LD
  # Reuse the genotype grids and weights from total genetic variance
  # BV is the sum of centered dosages multiplied by the locus-specific alphas
  myVarA_LE = numeric(nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha[[trait]])
    centeredDosages = sweep(myRefGeno[[trait]], MARGIN = 2, STATS = 2*p[loci], FUN = "-")
    values = as.vector(centeredDosages %*% myAlpha[[trait]])
    frequencies = myRefFreq_LE[[trait]]
    mean_LE = sum(frequencies*values)
    myVarA_LE[trait] = sum(frequencies*(values-mean_LE)^2)
  }
  # Unlike the observed variance, this excludes covariance between different loci

  test_that("genParam: Additive genetic variance (linkage-equilibrium population)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarA[trait] + gp$covA_HW[trait]), myVarA_LE[trait], tolerance = 1e-6)
    }
  })

  # ---- Additive genetic variance (idealized population) ----

  # Use HWE + LE weights and alpha_HW on the same centered genotype dosages
  # Both frequencies and substitution effects can change from the LE reference
  # Falconer (1961), p. 136: the single-locus result is 2*p*q*alpha_HW^2
  # https://archive.org/details/introductiontoqu0000falc/page/136
  myVarA_HW = numeric(nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha_HW[[trait]])
    centeredDosages = sweep(myRefGeno[[trait]], MARGIN = 2, STATS = 2*p[loci], FUN = "-")
    values = as.vector(centeredDosages %*% myAlpha_HW[[trait]])
    frequencies = myRefFreq_HW[[trait]]
    mean_HW = sum(frequencies*values)
    myVarA_HW[trait] = sum(frequencies*(values-mean_HW)^2)
  }
  # Fully inbred e now has hypothetical heterozygotes; fixed l remains fixed

  test_that("genParam: Additive genetic variance (idealized population)", {
    summary = genicVarA(pop, simParam = SP)
    expect_equal(unname(summary), myVarA_HW, tolerance = 1e-6)
    expect_identical(names(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarA[trait]), myVarA_HW[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Additive genetic variance (idealized population) - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      for (i in seq_along(case$traits)) {
        trait = case$traits[i]
        loci = if (trait == 9L) c("f", "g") else c("a", "d", "l")[match(trait, c(1L, 4L, 16L))]
        dosage = setNames(as.numeric(geno[case$row, loci]), loci)
        h = 2 * (dosage / 2) * (1 - dosage / 2)
        slope = a[loci] + d[loci] * (1 - dosage)
        if (trait == 9L) slope = slope + aa[trait] * rev(dosage - 1)
        expected[i] = sum(h * slope^2)
      }
      expect_equal(case$gp$genicVarA, expected, tolerance = 1e-6)
      expect_equal(genicVarA(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Additive genetic variance (reference-population adjustments) ----

  # Compare each returned adjustment with independent reference variances
  # The HW adjustment includes frequency and alpha changes; the LD adjustment
  # captures associations between loci with observed marginal alphas held fixed
  myCovA_HW = myVarA_LE - myVarA_HW
  myCovA_L = diag(myVarA) - myVarA_LE

  test_that("genParam: Additive genetic variance (reference-population adjustments)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covA_HW[trait]), myCovA_HW[trait], tolerance = 1e-6)
      expect_equal(unname(gp$covA_L[trait]), myCovA_L[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Additive genetic variance (reference-population adjustments) - one individual", {
    for (case in singleCases) {
      # Each observed marginal is fixed at one dosage, so LE variance is also zero.
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covA_HW, -case$gp$genicVarA, tolerance = 1e-6)
      expect_equal(case$gp$covA_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Dominance genetic variance (observed population) ----

  # The already tested dominance deviations have zero observed mean
  # Mean products give population variances and between-trait covariances
  myVarD = matrix(0, nrow = nTrt, ncol = nTrt)
  for (trait1 in 1:nTrt) {
    for (trait2 in 1:nTrt) {
      myVarD[trait1, trait2] = mean(myDd[, trait1] * myDd[, trait2])
    }
  }

  test_that("genParam: Dominance genetic variance (observed population)", {
    summary = varD(pop, simParam = SP)
    expect_equal(unname(summary), myVarD, tolerance = 1e-6)
    expect_identical(dimnames(summary), list(SP$traitNames, SP$traitNames))
    expect_equal(dim(gp$varD), c(nTrt, nTrt))
    expect_equal(rownames(gp$varD), colnames(gp$dd))
    expect_equal(colnames(gp$varD), colnames(gp$dd))
    for (trait1 in 1:nTrt) {
      for (trait2 in 1:nTrt) {
        expect_equal(unname(gp$varD[trait1, trait2]), myVarD[trait1, trait2], tolerance = 1e-6)
      }
    }
    # No dominance effect, full inbreeding and fixation respectively
    for (trait in c(4, 5, 16)) {
      expect_equal(unname(gp$varD[trait, ]), rep(0, times = nTrt), tolerance = 1e-6)
      expect_equal(unname(gp$varD[, trait]), rep(0, times = nTrt), tolerance = 1e-6)
    }
  })

  test_that("genParam: Dominance genetic variance (observed population) - one individual", {
    for (case in singleCases) {
      expected = matrix(0, length(case$traits), length(case$traits),
                        dimnames = list(case$SP$traitNames, case$SP$traitNames))
      expect_equal(case$gp$varD, expected, tolerance = 1e-6)
      expect_equal(varD(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Dominance genetic variance (linkage-equilibrium population) ----

  # Keep observed marginal genotype frequencies, but remove associations between loci
  # At each locus remove the mean and linear dosage effect from d*I(dosage==1)
  # The additive and marginal AA terms are linear, so leave no dominance residual
  # Cov(d*I(dosage==1), dosage) = d*H*(1-2*p)
  myRefDd_LE = vector("list", length = nTrt)
  myVarD_LE = numeric(nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha[[trait]])
    genotypes = myRefGeno[[trait]]
    values = numeric(nrow(genotypes))
    for (locus in loci) {
      dosageVariance = sum(c(Q[locus], H[locus], P[locus])*(0:2-2*p[locus])^2)
      dominanceSlope = 0
      if (dosageVariance > 0) {
        dominanceSlope = d[locus]*H[locus]*(1 - 2*p[locus])/dosageVariance
      }
      values = values + d[locus] * ((genotypes[, locus]==1) - H[locus]) -
        (genotypes[, locus] - 2*p[locus]) * dominanceSlope
    }
    myRefDd_LE[[trait]] = values
    frequencies = myRefFreq_LE[[trait]]
    mean_LE = sum(frequencies * values)
    myVarD_LE[trait] = sum(frequencies * (values - mean_LE)^2)
  }
  # Fixed loci have zero dosage variance and zero residual on their only genotype
  # Fully inbred e also has zero DD, because no observed genotype expresses d

  test_that("genParam: Dominance genetic variance (linkage-equilibrium population)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarD[trait] + gp$covD_HW[trait]), myVarD_LE[trait], tolerance = 1e-6)
    }
  })

  # ---- Dominance genetic variance (idealized population) ----

  # Reuse the HWE dominance deviations at dosages 0, 1, 2 from the values section
  # Weight their sums across loci by HWE & LE genotype frequencies
  # Falconer (1961), p. 136: the single-locus result is (2*p*q*d)^2
  # https://archive.org/details/introductiontoqu0000falc/page/136
  myRefDd_HW = vector("list", length = nTrt)
  myVarD_HW = numeric(nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha_HW[[trait]])
    genotypes = myRefGeno[[trait]]
    values = numeric(nrow(genotypes))
    for (locus in loci) {
      deviations = c(-2*p[locus]^2*d[locus],
                      2*p[locus]*q[locus]*d[locus],
                     -2*q[locus]^2*d[locus])
      values = values + deviations[genotypes[, locus] + 1]
    }
    myRefDd_HW[[trait]] = values
    frequencies = myRefFreq_HW[[trait]]
    mean_HW = sum(frequencies * values)
    myVarD_HW[trait] = sum(frequencies * (values - mean_HW)^2)
  }
  # Fully inbred e now has hypothetical heterozygotes and nonzero dominance variance
  # Fixed l remains fixed; trait 17 retains only b's dominance variance

  test_that("genParam: Dominance genetic variance (idealized population)", {
    summary = genicVarD(pop, simParam = SP)
    expect_equal(unname(summary), myVarD_HW, tolerance = 1e-6)
    expect_identical(names(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarD[trait]), myVarD_HW[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Dominance genetic variance (idealized population) - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      for (i in seq_along(case$traits)) {
        trait = case$traits[i]
        loci = if (trait == 9L) c("f", "g") else c("a", "d", "l")[match(trait, c(1L, 4L, 16L))]
        dosage = setNames(as.numeric(geno[case$row, loci]), loci)
        h = 2 * (dosage / 2) * (1 - dosage / 2)
        expected[i] = sum((h * d[loci])^2)
      }
      expect_equal(case$gp$genicVarD, expected, tolerance = 1e-6)
      expect_equal(genicVarD(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Dominance genetic variance (reference-population adjustments) ----

  # Compare each returned adjustment with independent reference variances
  myCovD_HW = myVarD_LE - myVarD_HW
  myCovD_L = diag(myVarD) - myVarD_LE

  test_that("genParam: Dominance genetic variance (reference-population adjustments)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covD_HW[trait]), myCovD_HW[trait], tolerance = 1e-6)
      expect_equal(unname(gp$covD_L[trait]), myCovD_L[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Dominance genetic variance (reference-population adjustments) - one individual", {
    for (case in singleCases) {
      # Each observed marginal is fixed at one dosage, so LE variance is also zero.
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covD_HW, -case$gp$genicVarD, tolerance = 1e-6)
      expect_equal(case$gp$covD_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Additive-by-additive epistatic genetic variance (observed population) ----

  # The already tested AA deviations have zero observed mean
  # Mean products give population variances and between-trait covariances
  myVarAA = matrix(0, nrow = nTrt, ncol = nTrt)
  for (trait1 in 1:nTrt) {
    for (trait2 in 1:nTrt) {
      myVarAA[trait1, trait2] = mean(myAa[, trait1] * myAa[, trait2])
    }
  }

  test_that("genParam: Additive-by-additive epistatic genetic variance (observed population)", {
    summary = varAA(pop, simParam = SP)
    expect_equal(unname(summary), myVarAA, tolerance = 1e-6)
    expect_identical(dimnames(summary), list(SP$traitNames, SP$traitNames))
    expect_equal(dim(gp$varAA), c(nTrt, nTrt))
    expect_equal(rownames(gp$varAA), colnames(gp$aa))
    expect_equal(colnames(gp$varAA), colnames(gp$aa))
    for (trait1 in 1:nTrt) {
      for (trait2 in 1:nTrt) {
        expect_equal(unname(gp$varAA[trait1, trait2]), myVarAA[trait1, trait2], tolerance = 1e-6)
      }
    }
    # No AA effect, or an AA pair with fixed partner l (trait 17)
    for (trait in c(1:8, 12, 14, 16, 17, 20)) {
      expect_equal(unname(gp$varAA[trait, ]), rep(0, times = nTrt), tolerance = 1e-6)
      expect_equal(unname(gp$varAA[, trait]), rep(0, times = nTrt), tolerance = 1e-6)
    }
  })

  test_that("genParam: Additive-by-additive epistatic genetic variance (observed population) - one individual", {
    for (case in singleCases) {
      expected = matrix(0, length(case$traits), length(case$traits),
                        dimnames = list(case$SP$traitNames, case$SP$traitNames))
      expect_equal(case$gp$varAA, expected, tolerance = 1e-6)
      expect_equal(varAA(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Additive-by-additive epistatic genetic variance (linkage-equilibrium population) ----

  # Reuse the centered-dosage product from the AA deviation section
  # https://archive.org/details/introductiontoqu0000falc/page/126
  # https://archive.org/details/introductiontoqu0000falc/page/128
  # Under LE, independent centered dosages give a zero mean product
  # Use observed marginal genotype frequencies, not the observed joint frequencies
  myRefAa_LE = vector("list", length = nTrt)
  myVarAA_LE = numeric(nTrt)
  for (trait in 1:nTrt) {
    loci = names(myAlpha[[trait]])
    genotypes = myRefGeno[[trait]]
    values = numeric(nrow(genotypes))
    centeredDosages = sweep(genotypes, MARGIN = 2, STATS = 2*p[loci], FUN = "-")
    for (pair in seq_along(myPairEffects[[trait]])) {
      partners = myPairLoci[[trait]][pair, ]
      values = values + myPairEffects[[trait]][pair] *
        centeredDosages[, partners[1]] * centeredDosages[, partners[2]]
    }
    frequencies = myRefFreq_LE[[trait]]
    mean_LE = sum(frequencies * values)
    myRefAa_LE[[trait]] = values - mean_LE
    myVarAA_LE[trait] = sum(frequencies * (values - mean_LE)^2)
  }
  # This removes the observed LD centering used by myAa
  # Fixed l contributes zero on every genotype with nonzero reference weight

  test_that("genParam: Additive-by-additive epistatic genetic variance (linkage-equilibrium population)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarAA[trait] + gp$covAA_HW[trait]), myVarAA_LE[trait], tolerance = 1e-6)
    }
  })

  # ---- Additive-by-additive epistatic genetic variance (idealized population) ----

  # The AA centered-dosage product is unchanged at the same allele frequencies
  # Unlike BV and DD, only its weights change between the two LE references
  # Independence gives Var(AA) = aa^2 * Var(dosage1) * Var(dosage2)
  # Under HWE each dosage variance is 2*p*q
  myRefAa_HW = vector("list", length = nTrt)
  myVarAA_HW = numeric(nTrt)
  for (trait in 1:nTrt) {
    values = myRefAa_LE[[trait]]
    frequencies = myRefFreq_HW[[trait]]
    mean_HW = sum(frequencies * values)
    myRefAa_HW[[trait]] = values - mean_HW
    myVarAA_HW[trait] = sum(frequencies * (values - mean_HW)^2)
  }
  # Mixed-HWE trait 15 changes variance; fixed-partner trait 17 still has zero AA variance

  test_that("genParam: Additive-by-additive epistatic genetic variance (idealized population)", {
    summary = genicVarAA(pop, simParam = SP)
    expect_equal(unname(summary), myVarAA_HW, tolerance = 1e-6)
    expect_identical(names(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarAA[trait]), myVarAA_HW[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Additive-by-additive epistatic genetic variance (idealized population) - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      for (i in seq_along(case$traits)) {
        trait = case$traits[i]
        loci = if (trait == 9L) c("f", "g") else c("a", "d", "l")[match(trait, c(1L, 4L, 16L))]
        dosage = setNames(as.numeric(geno[case$row, loci]), loci)
        h = 2 * (dosage / 2) * (1 - dosage / 2)
        if (trait == 9L) expected[i] = prod(h) * aa[trait]^2
      }
      expect_equal(case$gp$genicVarAA, expected, tolerance = 1e-6)
      expect_equal(genicVarAA(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Additive-by-additive epistatic genetic variance (reference-population adjustments) ----

  # Compare each returned adjustment with independent reference variances
  myCovAA_HW = myVarAA_LE - myVarAA_HW
  myCovAA_L = diag(myVarAA) - myVarAA_LE

  test_that("genParam: Additive-by-additive epistatic genetic variance (reference-population adjustments)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covAA_HW[trait]), myCovAA_HW[trait], tolerance = 1e-6)
      expect_equal(unname(gp$covAA_L[trait]), myCovAA_L[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Additive-by-additive epistatic genetic variance (reference-population adjustments) - one individual", {
    for (case in singleCases) {
      # Each observed marginal is fixed at one dosage, so LE variance is also zero.
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covAA_HW, -case$gp$genicVarAA, tolerance = 1e-6)
      expect_equal(case$gp$covAA_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Imprinting genetic variance (observed population) ----

  test_that("genParam: Imprinting genetic variance (observed population)", {
    skip("TODO: add tests for this section")
  })

  # ---- Imprinting genetic variance (linkage-equilibrium population) ----

  test_that("genParam: Imprinting genetic variance (linkage-equilibrium population)", {
    skip("TODO: add tests for this section")
  })

  # ---- Imprinting genetic variance (idealized population) ----

  test_that("genParam: Imprinting genetic variance (idealized population)", {
    skip("TODO: add tests for this section")
  })

  # ---- Non-additive genetic variance (observed population) ----

  # The already tested myNd = myDd + myAa has zero observed mean
  # Mean products include dominance-AA covariance and between-trait covariances
  myVarN = matrix(0, nrow = nTrt, ncol = nTrt)
  for (trait1 in 1:nTrt) {
    for (trait2 in 1:nTrt) {
      myVarN[trait1, trait2] = mean(myNd[, trait1] * myNd[, trait2])
    }
  }

  test_that("genParam: Non-additive genetic variance (observed population)", {
    summary = varN(pop, simParam = SP)
    expect_equal(unname(summary), myVarN, tolerance = 1e-6)
    expect_identical(dimnames(summary), list(SP$traitNames, SP$traitNames))
    expect_equal(dim(gp$varN), c(nTrt, nTrt))
    expect_equal(rownames(gp$varN), colnames(gp$nd))
    expect_equal(colnames(gp$varN), colnames(gp$nd))
    for (trait1 in 1:nTrt) {
      for (trait2 in 1:nTrt) {
        expect_equal(unname(gp$varN[trait1, trait2]), myVarN[trait1, trait2], tolerance = 1e-6)
      }
    }
    # No non-additive effect, full inbreeding and fixation respectively
    for (trait in c(4, 5, 16)) {
      expect_equal(unname(gp$varN[trait, ]), rep(0, times = nTrt), tolerance = 1e-6)
      expect_equal(unname(gp$varN[, trait]), rep(0, times = nTrt), tolerance = 1e-6)
    }
  })

  test_that("genParam: Non-additive genetic variance (observed population) - one individual", {
    for (case in singleCases) {
      expected = matrix(0, length(case$traits), length(case$traits),
                        dimnames = list(case$SP$traitNames, case$SP$traitNames))
      expect_equal(case$gp$varN, expected, tolerance = 1e-6)
      expect_equal(varN(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Non-additive genetic variance (linkage-equilibrium population) ----

  # Combine dominance and AA deviations on the existing reference genotype grids
  # Keep observed marginal genotype frequencies, but remove associations between loci
  myVarN_LE = numeric(nTrt)
  for (trait in 1:nTrt) {
    values = myRefDd_LE[[trait]] + myRefAa_LE[[trait]]
    frequencies = myRefFreq_LE[[trait]]
    mean_LE = sum(frequencies * values)
    myVarN_LE[trait] = sum(frequencies * (values - mean_LE)^2)
  }

  test_that("genParam: Non-additive genetic variance (linkage-equilibrium population)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarN[trait] + gp$covN_HW[trait]), myVarN_LE[trait], tolerance = 1e-6)
    }
  })

  # ---- Non-additive genetic variance (idealized population) ----

  # Combine HWE-reference dominance and AA deviations, weighted under HWE & LE
  myVarN_HW = numeric(nTrt)
  for (trait in 1:nTrt) {
    values = myRefDd_HW[[trait]] + myRefAa_HW[[trait]]
    frequencies = myRefFreq_HW[[trait]]
    mean_HW = sum(frequencies * values)
    myVarN_HW[trait] = sum(frequencies * (values - mean_HW)^2)
  }

  test_that("genParam: Non-additive genetic variance (idealized population)", {
    summary = genicVarN(pop, simParam = SP)
    expect_equal(unname(summary), myVarN_HW, tolerance = 1e-6)
    expect_identical(names(summary), SP$traitNames)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$genicVarN[trait]), myVarN_HW[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Non-additive genetic variance (idealized population) - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      for (i in seq_along(case$traits)) {
        trait = case$traits[i]
        loci = if (trait == 9L) c("f", "g") else c("a", "d", "l")[match(trait, c(1L, 4L, 16L))]
        dosage = setNames(as.numeric(geno[case$row, loci]), loci)
        h = 2 * (dosage / 2) * (1 - dosage / 2)
        expected[i] = sum((h * d[loci])^2)
        if (trait == 9L) expected[i] = expected[i] + prod(h) * aa[trait]^2
      }
      expect_equal(case$gp$genicVarN, expected, tolerance = 1e-6)
      expect_equal(genicVarN(case$pop, simParam = case$SP), expected, tolerance = 1e-6)
    }
  })

  # ---- Non-additive genetic variance (reference-population adjustments) ----

  # Compare each returned adjustment with independent reference variances
  myCovN_HW = myVarN_LE - myVarN_HW
  myCovN_L = diag(myVarN) - myVarN_LE

  test_that("genParam: Non-additive genetic variance (reference-population adjustments)", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covN_HW[trait]), myCovN_HW[trait], tolerance = 1e-6)
      expect_equal(unname(gp$covN_L[trait]), myCovN_L[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Non-additive genetic variance (reference-population adjustments) - one individual", {
    for (case in singleCases) {
      # Each observed marginal is fixed at one dosage, so LE variance is also zero.
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covN_HW, -case$gp$genicVarN, tolerance = 1e-6)
      expect_equal(case$gp$covN_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Covariance between breeding values and dominance deviations ----

  # The already tested observed deviations have zero mean
  # Mean products give within-trait population covariances (without a factor of two)
  # Single-locus regression values and dominance residuals are orthogonal;
  # however associations between loci can give a nonzero covariance
  myCovAD_L = numeric(nTrt)
  for (trait in 1:nTrt) {
    myCovAD_L[trait] = mean(myBv[, trait] * myDd[, trait])
  }

  test_that("genParam: Covariance between breeding values and dominance deviations", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covAD_L[trait]), myCovAD_L[trait], tolerance = 1e-6)
    }
    # Here is one example. Trait 5 has zero observed DD, but HWE-reference values
    # need not be orthogonal in the observed population.
    expect_equal(myCovAD_L[5], 0, tolerance = 1e-6)
    expect_true(mean(myBv_HW[, 5] * (myDd_HW[, 5] - mean(myDd_HW[, 5]))) > 0)
  })

  test_that("genParam: Covariance between breeding values and dominance deviations - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covAD_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Covariance between breeding values and additive-by-additive epistatic deviations ----

  # The already tested observed deviations have zero mean
  # Mean products give within-trait population covariances (without a factor of two)
  # Without AA effects, or with a fixed AA partner, the covariance is zero;
  # however associations between loci can make BV and AA covary
  myCovAAA_L = numeric(nTrt)
  for (trait in 1:nTrt) {
    myCovAAA_L[trait] = mean(myBv[, trait] * myAa[, trait])
  }

  test_that("genParam: Covariance between breeding values and additive-by-additive epistatic deviations", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covAAA_L[trait]), myCovAAA_L[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Covariance between breeding values and additive-by-additive epistatic deviations - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covAAA_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Covariance between dominance and additive-by-additive epistatic deviations ----

  # The already tested observed deviations have zero mean
  # Mean products give within-trait population covariances (without a factor of two)
  # Dominance and AA can covary in the observed population and
  # this covariance enters non-additive variance with a factor of two
  myCovDAA_L = numeric(nTrt)
  for (trait in 1:nTrt) {
    myCovDAA_L[trait] = mean(myDd[, trait] * myAa[, trait])
  }

  test_that("genParam: Covariance between dominance and additive-by-additive epistatic deviations", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covDAA_L[trait]), myCovDAA_L[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Covariance between dominance and additive-by-additive epistatic deviations - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covDAA_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Covariance between breeding values and non-additive deviations ----

  # The already tested observed deviations have zero mean
  # Mean products give within-trait population covariances (without a factor of two)
  # Reuse myNd = myDd + myAa to include both cross-component covariances
  myCovAN_L = numeric(nTrt)
  for (trait in 1:nTrt) {
    myCovAN_L[trait] = mean(myBv[, trait] * myNd[, trait])
  }

  test_that("genParam: Covariance between breeding values and non-additive deviations", {
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covAN_L[trait]), myCovAN_L[trait], tolerance = 1e-6)
    }
  })

  test_that("genParam: Covariance between breeding values and non-additive deviations - one individual", {
    for (case in singleCases) {
      expected = setNames(numeric(length(case$traits)), case$SP$traitNames)
      expect_equal(case$gp$covAN_L, expected, tolerance = 1e-6)
    }
  })

  # ---- Genetic variance decomposition reconstruction ----

  # Cross-component matrices include every pair of traits, not just the diagonal
  # Cov(A_i, N_j) need not equal Cov(A_j, N_i), so use C + t(C), not 2*C
  myCovAD = crossprod(myBv, y = myDd)/nInd
  myCovAAA = crossprod(myBv, y = myAa)/nInd
  myCovDAA = crossprod(myDd, y = myAa)/nInd
  myCovAN = crossprod(myBv, y = myNd)/nInd

  test_that("genParam: Genetic variance decomposition reconstruction", {
    expect_equal(unname(gp$varN), unname(gp$varD) + unname(gp$varAA) +
                   myCovDAA + t(myCovDAA), tolerance = 1e-6)
    expect_equal(unname(gp$varG), unname(gp$varA) + unname(gp$varN) +
                   myCovAN + t(myCovAN), tolerance = 1e-6)
    expect_equal(unname(gp$varG), unname(gp$varA) + unname(gp$varD) + unname(gp$varAA) +
                   myCovAD + t(myCovAD) + myCovAAA + t(myCovAAA) +
                   myCovDAA + t(myCovDAA), tolerance = 1e-6)
    for (trait in 1:nTrt) {
      expect_equal(unname(gp$covAN_L[trait]), unname(gp$covAD_L[trait] + gp$covAAA_L[trait]), tolerance = 1e-6)
      # Reference components are orthogonal under LE, including non-HWE marginals
      expect_equal(myVarG_LE[trait], myVarA_LE[trait] + myVarD_LE[trait] + myVarAA_LE[trait], tolerance = 1e-6)
      expect_equal(myVarG_HW[trait], myVarA_HW[trait] + myVarD_HW[trait] + myVarAA_HW[trait], tolerance = 1e-6)
      expect_equal(myVarN_LE[trait], myVarD_LE[trait] + myVarAA_LE[trait], tolerance = 1e-6)
      expect_equal(myVarN_HW[trait], myVarD_HW[trait] + myVarAA_HW[trait], tolerance = 1e-6)
      # Total LD adjustments already contain twice the cross-component covariances
      expect_equal(unname(gp$covG_L[trait]), unname(gp$covA_L[trait] + gp$covD_L[trait] + gp$covAA_L[trait] +
                   2*(gp$covAD_L[trait] + gp$covAAA_L[trait] + gp$covDAA_L[trait])), tolerance = 1e-6)
      expect_equal(unname(gp$covN_L[trait]), unname(gp$covD_L[trait] + gp$covAA_L[trait] +
                   2*gp$covDAA_L[trait]), tolerance = 1e-6)
    }
  })

})

