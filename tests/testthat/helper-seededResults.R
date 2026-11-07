# Scenarios for test-seededResults.R. Each one sets a seed, runs a short
# simulation with its own SimParam, and returns a named list of compact
# observations: exact checksums of genotypes and records of identifiers,
# and a few real-valued summaries. The expected values are stored in
# test-seededResults.R. When a change moves seeded results on purpose,
# regenerate them by running, from the package root with the package
# loaded,
#   testthat::source_test_helpers("tests/testthat", env = globalenv())
#   dput(lapply(seededScenarios, function(f) f()),
#        control = c("keepNA", "keepInteger", "niceNames", "digits17"))
# and copy the result over the stored values.
#
# MaCS founders are not used: the coalescent compares continuous values so
# often that a last-bit difference in a math library between platforms
# could change its output. quickHaplo draws from R's generator, and
# crossing draws from dqrng with Boost distributions, which are the same on
# every platform.

# An order-sensitive checksum of an integer matrix. Every term is an
# integer below 2^53, so the sum is exact whatever the platform.
seededChecksum = function(X) {
  X = as.vector(as.matrix(X))
  sum(as.numeric(X) * seq_along(X))
}

seededScenarios = list(

  # Founders, trait effects, sexes and phenotypes
  founders = function() {
    set.seed(9101)
    founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 50)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$setSexes("yes_rand")
    SP$addTraitA(10)
    SP$addTraitAD(10, meanDD = 0.2, varDD = 0.1)
    SP$addTraitADE(10, meanDD = 0.2, relAA = 0.5)
    SP$addTraitAG(10, varGxE = 0.5, varEnv = 0.2)
    SP$addSnpChip(10)
    SP$setVarE(h2 = c(0.5, 0.5, 0.5, 0.5))
    pop = newPop(founderPop, simParam = SP)
    list(haplo = seededChecksum(pullSegSiteHaplo(founderPop)),
         qtl = vapply(1:4, function(t)
           seededChecksum(pullQtlGeno(pop, trait = t, simParam = SP)), numeric(1)),
         snp = seededChecksum(pullSnpGeno(pop, simParam = SP)),
         sex = paste(pop@sex, collapse = ""),
         addEff = vapply(SP$traits, function(tr) sum(tr@addEff), numeric(1)),
         # The founders' genetic values sum to zero by construction, so a
         # few individuals' values are kept instead of the sum
         gv = as.vector(gv(pop)[1:3, ]),
         pheno = unname(colSums(pheno(pop))))
  },

  # Crossing, selection, doubled haploids, selfing and attrition
  crossing = function() {
    set.seed(9102)
    founderPop = quickHaplo(nInd = 20, nChr = 2, segSites = 50)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$setSexes("yes_sys")
    SP$addTraitAD(20, meanDD = 0.2)
    SP$setVarE(h2 = 0.4)
    pop = newPop(founderPop, simParam = SP)
    cross = randCross(pop, nCrosses = 30, nProgeny = 2, simParam = SP)
    selected = selectCross(cross, nFemale = 10, nMale = 5, nCrosses = 20,
                           simParam = SP)
    SP2 = SimParam$new(founderPop)
    SP2$nThreads = 1L
    SP2$addTraitA(20)
    SP2$setVarE(h2 = 0.4)
    pop2 = newPop(founderPop, simParam = SP2)
    dh = makeDH(pop2[1:5], nDH = 2, simParam = SP2)
    selfed = self(pop2[1:5], nProgeny = 2, simParam = SP2)
    reps = setPheno(selfed, reps = 3, simParam = SP2)
    kept = attrition(dh, p = 0.3)
    list(cross = seededChecksum(pullSegSiteHaplo(cross, simParam = SP)),
         crossParents = paste(cross@mother, cross@father, collapse = ","),
         selected = seededChecksum(pullSegSiteHaplo(selected, simParam = SP)),
         selectedParents = paste(selected@mother, selected@father, collapse = ","),
         dh = seededChecksum(pullSegSiteHaplo(dh, simParam = SP2)),
         selfed = seededChecksum(pullSegSiteHaplo(selfed, simParam = SP2)),
         kept = paste(kept@id, collapse = ","),
         pheno = sum(pheno(reps)))
  },

  # Pedigrees dropped through a map population and through a population
  pedigree = function() {
    set.seed(9103)
    founderPop = quickHaplo(nInd = 6, nChr = 2, segSites = 40)
    id = as.character(1:12)
    mother = c(NA, NA, NA, NA, "1", "3", "5", "5", "7", "7", "9", "9")
    father = c(NA, NA, NA, NA, "2", "4", "6", "6", "8", "8", "10", "10")
    mapDrop = pedigreeCross(founderPop, id, mother, father)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    pop = newPop(founderPop, simParam = SP)
    popDrop = pedigreeCross(pop, id, mother, father, simParam = SP)
    list(mapDrop = seededChecksum(pullSegSiteHaplo(mapDrop)),
         popDrop = seededChecksum(pullSegSiteHaplo(popDrop, simParam = SP)))
  },

  # Recombination tracking, mutation, and an alternative AD trait, whose
  # effects come from the C++ generator
  tracking = function() {
    set.seed(9104)
    founderPop = quickHaplo(nInd = 10, nChr = 2, segSites = 40)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$setTrackRec(TRUE)
    SP$altAddTraitAD(20, varD = 0.2, inbrDepr = 0.5, silent = TRUE)
    pop = randCross(newPop(founderPop, simParam = SP), nCrosses = 10,
                    simParam = SP)
    pop2 = randCross(pop, nCrosses = 10, simParam = SP)
    mutated = mutate(pop2, mutRate = 0.02, simParam = SP)
    list(ibd = seededChecksum(pullIbdHaplo(pop2, simParam = SP)),
         mutated = seededChecksum(pullSegSiteHaplo(mutated, simParam = SP)),
         addEff = sum(SP$traits[[1]]@addEff),
         domEff = sum(SP$traits[[1]]@domEff))
  },

  # Autotetraploid meiosis with quadrivalents, and changes of ploidy
  polyploid = function() {
    set.seed(9105)
    founderPop = quickHaplo(nInd = 10, nChr = 2, segSites = 40, ploidy = 4L)
    SP = SimParam$new(founderPop)
    SP$nThreads = 1L
    SP$quadProb = 0.5
    pop = newPop(founderPop, simParam = SP)
    cross = randCross(pop, nCrosses = 10, nProgeny = 2, simParam = SP)
    reduced = reduceGenome(cross[1:5], simParam = SP)
    doubled = doubleGenome(reduced, simParam = SP)
    list(cross = seededChecksum(pullSegSiteHaplo(cross, simParam = SP)),
         reduced = seededChecksum(pullSegSiteHaplo(reduced, simParam = SP)),
         doubled = seededChecksum(pullSegSiteHaplo(doubled, simParam = SP)))
  }
)
