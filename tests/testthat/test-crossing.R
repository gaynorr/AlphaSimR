context("crossing")

founderPop = newMapPop(list(c(0)),
                       list(matrix(c(1,1,0,0),
                                   nrow=4,ncol=1)))

test_that("makeCross",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop,simParam=SP)
  crossPlan = cbind(rep(1,10),rep(2,10))
  expect_error(makeCross(pop=pop,crossPlan=matrix(c("3", "4"), ncol = 2),simParam=SP),
               regexp = "Failed to match supplied IDs")
  expect_error(makeCross(pop=pop,crossPlan=matrix(c(3, 4), ncol = 2),simParam=SP),
               regexp = "Invalid crossPlan")
  #Match by number
  pop1 = makeCross(pop=pop,crossPlan=crossPlan,simParam=SP)
  expect_equal(SP$lastId,12L)
  expect_equal(unname(meanG(pop1)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),rep(1L,10L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),rep(2L,10L))
  #Match by id
  crossPlan = matrix(as.character(crossPlan),ncol=2)
  pop2 = makeCross(pop=pop,crossPlan=crossPlan,simParam=SP)
  expect_equal(SP$lastId,22L)
  expect_equal(unname(meanG(pop2)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),rep(1L,20L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),rep(2L,20L))
})

test_that("makeCross2",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop,simParam=SP)
  crossPlan = cbind(rep(1,10),rep(2,10))
  #Match by number
  pop1 = makeCross2(females=pop,males=pop,crossPlan=crossPlan,simParam=SP)
  expect_equal(SP$lastId,12L)
  expect_equal(unname(meanG(pop1)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),rep(1L,10L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),rep(2L,10L))
  #Match by id
  crossPlan = matrix(as.character(crossPlan),ncol=2)
  pop2 = makeCross2(females=pop,males=pop,crossPlan=crossPlan,simParam=SP)
  expect_equal(SP$lastId,22L)
  expect_equal(unname(meanG(pop2)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),rep(1L,20L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),rep(2L,20L))
})

test_that("randCross",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  SP$setTrackPed(TRUE)
  SP$setSexes("yes_sys")
  pop = newPop(founderPop,simParam=SP)
  crossPlan = cbind(rep(1,10),rep(2,10))
  pop1 = randCross(pop=pop,nCrosses=1,nProgeny=10,simParam=SP)
  expect_equal(SP$lastId,12L)
  expect_equal(unname(meanG(pop1)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),rep(2L,10L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),rep(1L,10L))
})

test_that("randCross2",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  SP$setTrackPed(TRUE)
  SP$setSexes("yes_sys")
  pop = newPop(founderPop,simParam=SP)
  crossPlan = cbind(rep(1,10),rep(2,10))
  pop1 = randCross2(females=pop,males=pop,nCrosses=1,nProgeny=10,simParam=SP)
  expect_equal(SP$lastId,12L)
  expect_equal(unname(meanG(pop1)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),rep(2L,10L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),rep(1L,10L))
})

test_that("self",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop,simParam=SP)
  crossPlan = cbind(rep(1,10),rep(2,10))
  pop1 = self(pop=pop,nProgeny=1,simParam=SP)
  expect_equal(SP$lastId,4L)
  expect_equal(unname(meanG(pop1)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),c(1L,2L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),c(1L,2L))
  expect_equal(unname(SP$pedigree[-(1:2),3L]),c(0L,0L))
})

test_that("makeDH",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop,simParam=SP)
  crossPlan = cbind(rep(1,10),rep(2,10))
  pop1 = makeDH(pop=pop,nDH=1,simParam=SP)
  expect_equal(SP$lastId,4L)
  expect_equal(unname(meanG(pop1)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),c(1L,2L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),c(1L,2L))
  expect_equal(unname(SP$pedigree[-(1:2),3L]),c(1L,1L))
})

test_that("selectCross",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  SP$setTrackPed(TRUE)
  SP$setSexes("yes_sys")
  pop = newPop(founderPop,simParam=SP)
  pop1 = selectCross(pop=pop,nFemale=1,nMale=1,use="rand",nCrosses=2,simParam=SP)
  expect_equal(SP$lastId,4L)
  expect_equal(unname(meanG(pop1)),0,tolerance=1e-6)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),c(2L,2L))
  expect_equal(unname(SP$pedigree[-(1:2),2L]),c(1L,1L))
})

# A few large families of full sibs, so truncation on genetic value
# concentrates on the best families and the restriction binds
ocsSetup = function(sexes="no", nProgeny=20, ploidy=2L){
  founder = quickHaplo(nInd=40, nChr=2, segSites=100, ploidy=ploidy)
  SP = SimParam$new(founder)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=20)
  SP$addSnpChip(nSnpPerChr=50)
  if(sexes!="no") SP$setSexes(sexes)
  pop = newPop(founder, simParam=SP)
  pop = randCross(pop, nCrosses=10, nProgeny=nProgeny, simParam=SP)
  list(SP=SP, pop=pop)
}

# Expected fixation of a set of parents contributing equally
ocsFix = function(Z, take){
  mean(colMeans(Z[take,,drop=FALSE])^2)
}

test_that("MATCH plain truncation when the restriction does not bind",{
  set.seed(1)
  s = ocsSetup()
  set.seed(2)
  pop1 = selectCross(s$pop, nInd=50, nCrosses=40, use="gv",
                     simParam=s$SP)
  set.seed(2)
  pop2 = suppressMessages(
    selectCross(s$pop, nInd=50, nCrosses=40, use="gv", restrInbr=TRUE,
                inbrTarget=1, inbrType="absolute", simParam=s$SP))
  expect_identical(pop1@mother, pop2@mother)
  expect_identical(pop1@father, pop2@father)
  # The IDs differ because the second call follows the first
  expect_identical(unname(pullSegSiteGeno(pop1, simParam=s$SP)),
                   unname(pullSegSiteGeno(pop2, simParam=s$SP)))
})

test_that("MEET a relative target without sexes",{
  set.seed(3)
  s = ocsSetup()
  Z = pullSnpGeno(s$pop, simParam=s$SP) - 1
  Ft = mean(colMeans(Z)^2)
  target = Ft + 0.01*(1-Ft)
  expect_message(
    take <- selectOCS(s$pop, nInd=100, use="gv", returnPop=FALSE,
                      simParam=s$SP),
    regexp="selectOCS: kappa")
  expect_length(take, 100L)
  expect_false(anyDuplicated(take) > 0)
  expect_lte(ocsFix(Z, take), target)
  # Truncation alone misses the target, so the penalty was needed
  trunc = selectInd(s$pop, nInd=100, use="gv", returnPop=FALSE,
                    simParam=s$SP)
  expect_gt(ocsFix(Z, trunc), target)

  # The same target given as an absolute value selects the same parents
  take2 = suppressMessages(
    selectOCS(s$pop, nInd=100, use="gv", inbrTarget=target,
              inbrType="absolute", returnPop=FALSE, simParam=s$SP))
  expect_identical(take, take2)

  # Returning a population gives the same individuals in the same order
  pop2 = suppressMessages(
    selectOCS(s$pop, nInd=100, use="gv", simParam=s$SP))
  expect_identical(pop2@id, s$pop@id[take])
})

test_that("MEET a relative target with sexes, pooling the sexes",{
  set.seed(4)
  s = ocsSetup("yes_sys", nProgeny=30)
  Z = pullSnpGeno(s$pop, simParam=s$SP) - 1
  Ft = mean(colMeans(Z)^2)
  target = Ft + 0.01*(1-Ft)
  take = suppressMessages(
    selectOCS(s$pop, nFemale=60, nMale=60, use="gv", returnPop=FALSE,
              simParam=s$SP))
  # Females come first, then males
  female = take[1:60]
  male = take[61:120]
  expect_true(all(s$pop@sex[female]=="F"))
  expect_true(all(s$pop@sex[male]=="M"))
  poolFix = function(female, male){
    mean(((colMeans(Z[female,,drop=FALSE]) +
           colMeans(Z[male,,drop=FALSE]))/2)^2)
  }
  expect_lte(poolFix(female, male), target)
  # Truncation alone misses the target, so the penalty was needed
  truncF = selectInd(s$pop, nInd=60, use="gv", sex="F",
                     returnPop=FALSE, simParam=s$SP)
  truncM = selectInd(s$pop, nInd=60, use="gv", sex="M",
                     returnPop=FALSE, simParam=s$SP)
  expect_gt(poolFix(truncF, truncM), target)

  # selectCross crosses exactly the individuals selectOCS selects, and
  # with balance every one of them is used
  pop2 = suppressMessages(
    selectCross(s$pop, nFemale=60, nMale=60, nCrosses=60, use="gv",
                restrInbr=TRUE, simParam=s$SP))
  expect_setequal(pop2@mother, s$pop@id[female])
  expect_setequal(pop2@father, s$pop@id[male])
})

test_that("USE QTL genotypes to measure fixation",{
  set.seed(5)
  s = ocsSetup()
  Z = pullQtlGeno(s$pop, simParam=s$SP) - 1
  Ft = mean(colMeans(Z)^2)
  take = suppressMessages(
    selectOCS(s$pop, nInd=100, use="gv", useQtl=TRUE, returnPop=FALSE,
              simParam=s$SP))
  expect_lte(ocsFix(Z, take), Ft + 0.01*(1-Ft))
})

test_that("AGREE between the precomputed and direct penalty kernels",{
  set.seed(9)
  for(ploidy in c(2L, 4L)){
    s = ocsSetup(ploidy=ploidy)
    Y =2*pullSnpGeno(s$pop, simParam=s$SP) - ploidy
    Z = Y/ploidy
    kerK = AlphaSimR:::.ocsKernel(Y, ploidy, useK=TRUE)
    kerD = AlphaSimR:::.ocsKernel(Y, ploidy, useK=FALSE)
    take = c(3, 17, 40, 41, 150)
    takeM = c(5, 60, 199)
    # Both work in whole numbers, so they agree exactly
    expect_identical(kerK$zu0, kerD$zu0)
    expect_identical(kerK$self, kerD$self)
    expect_equal(kerK$self, rowSums(Z^2))
    expect_identical(kerK$zu(take), kerD$zu(take))
    expect_identical(kerK$fix(take), kerD$fix(take))
    expect_identical(kerK$fixPool(take, takeM), kerD$fixPool(take, takeM))
    # and match the definitions in terms of Z
    u = colMeans(Z[take,])
    expect_equal(kerK$zu(take), drop(Z%*%u))
    expect_equal(kerK$fix(take), mean(u^2))
    uPool = (u + colMeans(Z[takeM,]))/2
    expect_equal(kerK$fixPool(take, takeM), mean(uPool^2))
    m = gv(s$pop)[,1]
    expect_identical(AlphaSimR:::.ocsTrunc(0.05, kerK, m, 50),
                     AlphaSimR:::.ocsTrunc(0.05, kerD, m, 50))
  }
})

test_that("MEET a relative target in a tetraploid",{
  set.seed(1)
  s = ocsSetup(ploidy=4L)
  # Dosages run from 0 to 4, so they are scaled to run from -1 to 1
  Z = pullSnpGeno(s$pop, simParam=s$SP)/2 - 1
  Ft = mean(colMeans(Z)^2)
  target = Ft + 0.01*(1-Ft)
  # A random set of forty tetraploid parents costs about 1/(4*40) of
  # the heterozygosity, leaving room within the target for selection
  expect_silent(take <- suppressMessages(
    selectOCS(s$pop, nInd=40, use="gv", returnPop=FALSE,
              simParam=s$SP)))
  expect_lte(ocsFix(Z, take), target)
  # Truncation alone misses the target, so the penalty was needed
  trunc = selectInd(s$pop, nInd=40, use="gv", returnPop=FALSE,
                    simParam=s$SP)
  expect_gt(ocsFix(Z, trunc), target)
})

test_that("WARN when the restriction is approximate or not met",{
  set.seed(6)
  s = ocsSetup("yes_sys")
  expect_warning(suppressMessages(
    selectCross(s$pop, nInd=100, nCrosses=50, use="gv", restrInbr=TRUE,
                simParam=s$SP)),
    regexp="only approximate")
  # The warning is about random crossing, so selectOCS does not give it
  expect_silent(suppressMessages(
    selectOCS(s$pop, nInd=100, use="gv", simParam=s$SP)))
  # Zero expected fixation needs the two parents to have an allele
  # frequency of exactly 0.5 at every locus, which is not found here
  expect_warning(suppressMessages(
    selectCross(s$pop, nFemale=1, nMale=1, nCrosses=2, use="gv",
                restrInbr=TRUE, inbrTarget=0, inbrType="absolute",
                simParam=s$SP)),
    regexp="target was not met")
})

test_that("REFUSE invalid restricted inbreeding settings",{
  set.seed(7)
  s = ocsSetup()
  expect_error(selectCross(s$pop, nInd=50, nCrosses=20, use="gv",
                           restrInbr=TRUE, inbrType="foo",
                           simParam=s$SP),
               regexp="inbrType")
  expect_error(selectCross(s$pop, nInd=50, nCrosses=20, use="gv",
                           restrInbr=TRUE, inbrTarget=-0.1,
                           simParam=s$SP),
               regexp="inbrTarget")
})

test_that("REFUSE a MultiPop or HybridPop in selectOCS",{
  set.seed(8)
  s = ocsSetup()
  expect_error(selectOCS(newMultiPop(s$pop, s$pop), nInd=10, use="gv",
                         simParam=s$SP),
               regexp="MultiPop")
  hybrid = hybridCross(s$pop[1:5], s$pop[6:10], returnHybridPop=TRUE,
                       simParam=s$SP)
  expect_error(selectOCS(hybrid, nInd=2, use="gv", simParam=s$SP),
               regexp="HybridPop")
})

test_that("selectOP",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop,simParam=SP)
  pop1 = selectOP(pop=pop,nInd=1,nSeeds=2,use="gv",simParam=SP)
  expect_equal(SP$lastId,4L)
  expect_equal(unname(meanG(pop1)),0,tolerance=1e-6)
  tmp = abs(SP$pedigree[-(1:2),1L]-SP$pedigree[-(1:2),2L])
  expect_equal(unname(tmp),c(1L,1L))
})

test_that("REFUSE self parents outside the population",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop,simParam=SP)
  # The C++ code indexes the genotypes with parents directly
  for(bad in list(0, -1, 3, NA, c(1, 3))){
    expect_error(self(pop,parents=bad,simParam=SP),
                 regexp = "Invalid parents")
  }
  # Nothing was made, so the IDs have not moved on
  expect_equal(SP$lastId,2L)
  pop1 = self(pop,parents=c(2,1),nProgeny=2,simParam=SP)
  expect_equal(pop1@nInd,4L)
  expect_equal(unname(SP$pedigree[-(1:2),1L]),c(2L,2L,1L,1L))
})

test_that("REFUSE invalid parents in the C++ crossing code",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop,simParam=SP)
  doCross = function(mother, father){
    AlphaSimR:::cross(pop@geno, mother, pop@geno, father,
                      SP$femaleMap, SP$maleMap, FALSE, 2L, 2L,
                      SP$v, SP$p, SP$femaleCentromere, SP$maleCentromere,
                      SP$quadProb, 1L)
  }
  expect_error(doCross(c(1L,0L), c(1L,2L)), regexp = "Invalid parent index")
  expect_error(doCross(c(1L,2L), c(1L,5L)), regexp = "Invalid parent index")
  expect_error(doCross(c(1L,-1L), c(1L,2L)), regexp = "Invalid parent index")
  expect_error(doCross(c(1L,2L), 1L), regexp = "same length")
})
