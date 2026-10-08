context("hybrids")

founderPop = newMapPop(list(c(0)),
                       list(matrix(c(1,1,0,0),
                                   nrow=4,ncol=1)))

test_that("hybridCross",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  pop = newPop(founderPop,simParam=SP)
  #2x2
  hybrid = hybridCross(pop,pop,simParam=SP)
  expect_equal(hybrid@nInd,4L)
  hybrid = hybridCross(pop,pop,returnHybridPop=TRUE,
                       simParam=SP)
  expect_equal(hybrid@nInd,4L)
  #2x1
  hybrid = hybridCross(pop,pop[1],simParam=SP)
  expect_equal(hybrid@nInd,2L)
  hybrid = hybridCross(pop,pop[1],returnHybridPop=TRUE,
                       simParam=SP)
  expect_equal(hybrid@nInd,2L)
  #1x1
  hybrid = hybridCross(pop[1],pop[1],simParam=SP)
  expect_equal(hybrid@nInd,1L)
  hybrid = hybridCross(pop[1],pop[1],returnHybridPop=TRUE,
                       simParam=SP)
  expect_equal(hybrid@nInd,1L)
})

test_that("calcGCA",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=c(0,0),var=c(1,1))
  SP$setVarE(varE=c(1,1))
  pop = newPop(founderPop,simParam=SP)
  #2x2
  hybrid = hybridCross(pop,pop,returnHybridPop=TRUE,
                       simParam=SP)
  GCA = calcGCA(hybrid)
  expect_equal(nrow(GCA$GCAf),2L)
  expect_equal(nrow(GCA$GCAm),2L)
  expect_equal(nrow(GCA$SCA),4L)
  #2x1
  hybrid = hybridCross(pop,pop[1],returnHybridPop=TRUE,
                       simParam=SP)
  GCA = calcGCA(hybrid)
  expect_equal(nrow(GCA$GCAf),2L)
  expect_equal(nrow(GCA$GCAm),1L)
  expect_equal(nrow(GCA$SCA),2L)
  #1x2
  hybrid = hybridCross(pop[1],pop,returnHybridPop=TRUE,
                       simParam=SP)
  GCA = calcGCA(hybrid)
  expect_equal(nrow(GCA$GCAf),1L)
  expect_equal(nrow(GCA$GCAm),2L)
  expect_equal(nrow(GCA$SCA),2L)
  #1x1
  hybrid = hybridCross(pop[1],pop[1],returnHybridPop=TRUE,
                       simParam=SP)
  GCA = calcGCA(hybrid)
  expect_equal(nrow(GCA$GCAf),1L)
  expect_equal(nrow(GCA$GCAm),1L)
  expect_equal(nrow(GCA$SCA),1L)
})

# ---------------------------------------------------------------------------
# 1. hybridCross with map populations
# ---------------------------------------------------------------------------
#
# A pair of map populations means no simulation has been set up yet, and the
# hybrids are meant to become the founders of one. The check used below
# holds exactly: a hybrid of two fully inbred diploid parents carries one
# copy of each parent's only haplotype, so its dosage at every locus is the
# mean of its parents' dosages.

# Fully inbred haplotypes, split into female and male sets
hybridMaps = function(nFemale=3, nMale=2, nChr=2, segSites=20, seed=12001){
  set.seed(seed)
  inbreds = quickHaplo(nInd=nFemale+nMale, nChr=nChr, segSites=segSites,
                       inbred=TRUE)
  return(list(females=inbreds[seq_len(nFemale)],
              males=inbreds[nFemale+seq_len(nMale)]))
}

# A NamedMapPop, which is the only map class carrying ids. Nothing exported
# builds one from a MapPop, so it is assembled from the slots here.
asNamedMapPop = function(mapPop, id){
  return(new("NamedMapPop",
             id = as.character(id),
             mother = rep(NA_character_, mapPop@nInd),
             father = rep(NA_character_, mapPop@nInd),
             nInd = mapPop@nInd,
             nChr = mapPop@nChr,
             ploidy = mapPop@ploidy,
             nLoci = mapPop@nLoci,
             geno = mapPop@geno,
             genMap = mapPop@genMap,
             centromere = mapPop@centromere,
             inbred = mapPop@inbred))
}

test_that("MAP a pair of MapPops returns a MapPop of their hybrids", {
  d = hybridMaps()
  out = hybridCross(d$females, d$males, nThreads=1L)
  expect_true(isMapPop(out))
  expect_false(isNamedMapPop(out))
  expect_equal(out@nInd, 6L)
  expect_false(out@inbred)
  expect_equal(out@genMap, d$females@genMap)
  expect_equal(out@centromere, d$females@centromere)

  # A testcross runs every male within each female
  crossPlan = cbind(rep(1:3, each=2), rep(1:2, 3))
  fGeno = pullSegSiteGeno(d$females, nThreads=1L)
  mGeno = pullSegSiteGeno(d$males, nThreads=1L)
  expected = (fGeno[crossPlan[,1],] + mGeno[crossPlan[,2],])/2
  expect_equal(unname(pullSegSiteGeno(out, nThreads=1L)), unname(expected))
})

test_that("MAP a designed crossPlan is followed", {
  d = hybridMaps()
  crossPlan = cbind(c(3,1), c(1,2))
  out = hybridCross(d$females, d$males, crossPlan=crossPlan, nThreads=1L)
  expect_equal(out@nInd, 2L)
  fGeno = pullSegSiteGeno(d$females, nThreads=1L)
  mGeno = pullSegSiteGeno(d$males, nThreads=1L)
  expected = (fGeno[crossPlan[,1],] + mGeno[crossPlan[,2],])/2
  expect_equal(unname(pullSegSiteGeno(out, nThreads=1L)), unname(expected))
})

test_that("NAMED a pair of NamedMapPops returns named hybrids", {
  d = hybridMaps()
  females = asNamedMapPop(d$females, c("F1","F2","F3"))
  males = asNamedMapPop(d$males, c("M1","M2"))
  out = hybridCross(females, males, nThreads=1L)
  expect_true(isNamedMapPop(out))
  expect_equal(out@id, c("F1_M1","F1_M2","F2_M1","F2_M2","F3_M1","F3_M2"))
  expect_equal(out@mother, rep(c("F1","F2","F3"), each=2))
  expect_equal(out@father, rep(c("M1","M2"), 3))
  expect_false(out@inbred)
})

test_that("NAMED a NamedMapPop crossed to a MapPop returns a MapPop", {
  d = hybridMaps()
  females = asNamedMapPop(d$females, c("F1","F2","F3"))
  out = hybridCross(females, d$males, nThreads=1L)
  expect_true(isMapPop(out))
  expect_false(isNamedMapPop(out))
})

test_that("MAP the returned population starts a simulation", {
  d = hybridMaps()
  females = asNamedMapPop(d$females, c("F1","F2","F3"))
  males = asNamedMapPop(d$males, c("M1","M2"))
  out = hybridCross(females, males, nThreads=1L)
  SP = SimParam$new(out)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=5)
  pop = newPop(out, simParam=SP)
  expect_equal(pop@nInd, 6L)
  expect_equal(pop@id, out@id)
  expect_equal(pop@mother, out@mother)
})

# ---------------------------------------------------------------------------
# 2. Repeated hybrid ids
# ---------------------------------------------------------------------------

test_that("ID every copy of a repeated cross gets a letter code", {
  d = hybridMaps()
  females = asNamedMapPop(d$females, c("F1","F2","F3"))
  males = asNamedMapPop(d$males, c("M1","M2"))
  crossPlan = cbind(c(1,2,1,1), c(1,1,1,1))
  out = hybridCross(females, males, crossPlan=crossPlan, nThreads=1L)
  expect_equal(out@id, c("F1_M1_a","F2_M1","F1_M1_b","F1_M1_c"))
  expect_false(anyDuplicated(out@id)>0L)
  # Parentage is not altered by the renaming
  expect_equal(out@mother, c("F1","F2","F1","F1"))
  expect_equal(out@father, rep("M1", 4))
})

test_that("ID a letter code already in use is skipped", {
  id = c("A_B", "A_B", "A_B_a")
  expect_equal(AlphaSimR:::uniqueHybridId(id), c("A_B_b","A_B_c","A_B_a"))
  expect_equal(AlphaSimR:::uniqueHybridId(c("x","y")), c("x","y"))
})

test_that("ID letter codes follow z with aa", {
  letterCode = AlphaSimR:::letterCode
  expect_equal(letterCode(1L), "a")
  expect_equal(letterCode(26L), "z")
  expect_equal(letterCode(27L), "aa")
  expect_equal(letterCode(52L), "az")
  expect_equal(letterCode(53L), "ba")
  expect_equal(letterCode(702L), "zz")
  expect_equal(letterCode(703L), "aaa")
  id = AlphaSimR:::uniqueHybridId(rep("A_B", 28))
  expect_equal(id[c(1,26,27,28)], c("A_B_a","A_B_z","A_B_aa","A_B_ab"))
})

# ---------------------------------------------------------------------------
# 3. Arguments that cannot be honoured
# ---------------------------------------------------------------------------

test_that("ARGS returnHybridPop=TRUE is refused for map populations", {
  d = hybridMaps()
  expect_error(hybridCross(d$females, d$males, returnHybridPop=TRUE,
                           nThreads=1L),
               "returnHybridPop")
})

test_that("ARGS a Pop cannot be crossed to a map population", {
  d = hybridMaps()
  SP = SimParam$new(d$females)
  SP$nThreads = 1L
  pop = newPop(d$females, simParam=SP)
  expect_error(hybridCross(pop, d$males, simParam=SP),
               "both be map populations")
  expect_error(hybridCross(d$females, pop, simParam=SP),
               "both be map populations")
})

test_that("ARGS map populations must share a genetic map", {
  d = hybridMaps()
  set.seed(12002)
  other = quickHaplo(nInd=2, nChr=2, segSites=20, genLen=2, inbred=TRUE)
  expect_error(hybridCross(d$females, other, nThreads=1L), "genMap")
  fewer = quickHaplo(nInd=2, nChr=2, segSites=10, inbred=TRUE)
  expect_error(hybridCross(d$females, fewer, nThreads=1L), "nLoci")
})

test_that("ARGS simParam is ignored with a warning for map populations", {
  d = hybridMaps()
  SP = SimParam$new(d$females)
  SP$nThreads = 1L
  expect_warning(out <- hybridCross(d$females, d$males, simParam=SP,
                                    nThreads=1L),
                 "simParam is ignored")
  expect_true(isMapPop(out))
  # The ids handed out come from the temporary SimParam, not the user's
  expect_equal(SP$lastId, 0L)
})

test_that("ARGS recombination settings are checked and refused for a Pop", {
  d = hybridMaps()
  out = hybridCross(d$females, d$males, v=1, p=0.5, quadProb=0,
                    nThreads=1L)
  expect_equal(out@nInd, 6L)
  # A misspelling that is a prefix of a named argument, such as nThread,
  # would be partially matched to it rather than reaching ...
  expect_error(hybridCross(d$females, d$males, quadprob=0, nThreads=1L),
               "Unused arguments")
  expect_error(hybridCross(d$females, d$males, p=2, nThreads=1L),
               "p must be between zero and one")

  SP = SimParam$new(d$females)
  SP$nThreads = 1L
  pop = newPop(d$females, simParam=SP)
  expect_error(hybridCross(pop, pop, v=1, simParam=SP),
               "can only be set when females and males")
  # Refusing them leaves the user's SimParam as it was
  expect_equal(SP$v, SimParam$new(d$females)$v)
})

test_that("REFUSE a hybrid crossPlan outside the populations",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  pop = newPop(founderPop,simParam=SP)
  # returnHybridPop=TRUE does not go through makeCross2, and the C++ code
  # indexes the genotypes with crossPlan directly
  bad = list(cbind(c(1,0),c(1,1)),
             cbind(c(1,5),c(1,1)),
             cbind(c(1,1),c(1,5)),
             cbind(c(1,NA),c(1,1)),
             c(1,2))
  for(returnHybridPop in c(TRUE,FALSE)){
    for(crossPlan in bad){
      expect_error(hybridCross(pop,pop,crossPlan=crossPlan,
                               returnHybridPop=returnHybridPop,
                               simParam=SP),
                   "Invalid crossPlan")
    }
    expect_error(hybridCross(pop,pop,crossPlan=cbind("1","9"),
                             returnHybridPop=returnHybridPop,
                             simParam=SP),
                 "Failed to match supplied IDs")
  }
  # A crossPlan of IDs works for a HybridPop, as it does for a Pop
  hybrid = hybridCross(pop,pop,crossPlan=cbind(c("1","2"),c("2","1")),
                       returnHybridPop=TRUE,simParam=SP)
  expect_equal(hybrid@id,c("1_2","2_1"))
  # The C++ code checks the indexes as well
  expect_error(AlphaSimR:::getHybridGv(SP$traits[[1]],pop,c(1L,0L),
                                       pop,c(1L,1L),1L),
               "Invalid parent index")
  expect_error(AlphaSimR:::getHybridGv(SP$traits[[1]],pop,c(1L,5L),
                                       pop,c(1L,1L),1L),
               "Invalid parent index")
})
