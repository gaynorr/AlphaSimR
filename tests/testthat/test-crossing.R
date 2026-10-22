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
ocsSetup = function(sexes="no", nProgeny=20, ploidy=2L, trackPed=FALSE){
  founder = quickHaplo(nInd=40, nChr=2, segSites=100, ploidy=ploidy)
  SP = SimParam$new(founder)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=20)
  SP$addSnpChip(nSnpPerChr=50)
  if(sexes!="no") SP$setSexes(sexes)
  if(trackPed) SP$setTrackPed(TRUE)
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

# Expected fixation measured by the pedigree: the average kinship of a set
# of parents contributing equally
ocsPedFix = function(K, take){
  mean(K[take,take,drop=FALSE])
}

test_that("USE the pedigree to measure fixation",{
  set.seed(10)
  for(ploidy in c(2L, 4L)){
    s = ocsSetup(ploidy=ploidy, trackPed=TRUE)
    K = pedigreeA(s$pop, simParam=s$SP)/ploidy
    Ft = ocsPedFix(K, seq_len(s$pop@nInd))
    target = Ft + 0.01*(1-Ft)
    take = suppressMessages(
      selectOCS(s$pop, nInd=40, use="gv", usePed=TRUE, returnPop=FALSE,
                simParam=s$SP))
    expect_length(take, 40L)
    expect_false(anyDuplicated(take) > 0)
    expect_lte(ocsPedFix(K, take), target)
    trunc = selectInd(s$pop, nInd=40, use="gv", returnPop=FALSE,
                      simParam=s$SP)
    # Whenever truncation misses the target, the penalty has to change the
    # selection
    if(ocsPedFix(K, trunc) > target){
      expect_false(setequal(take, trunc))
    }
    # For the diploid, truncation misses the target by a wide margin, so the
    # penalty is known to be exercised. For the tetraploid the margin is
    # thin, and the genotypes, and so which individuals truncation picks,
    # are not identical across platforms: macOS put truncation under the
    # target where Linux and Windows put it over.
    if(ploidy==2L){
      expect_gt(ocsPedFix(K, trunc), target)
    }
  }
})

test_that("USE the pedigree with sexes, pooling the sexes",{
  set.seed(11)
  s = ocsSetup("yes_sys", nProgeny=30, trackPed=TRUE)
  K = pedigreeA(s$pop, simParam=s$SP)/2
  Ft = ocsPedFix(K, seq_len(s$pop@nInd))
  target = Ft + 0.01*(1-Ft)
  # Each sex supplies half of the gametes
  poolFix = function(female, male){
    c = c(rep(1/(2*length(female)), length(female)),
          rep(1/(2*length(male)), length(male)))
    take = c(female, male)
    drop(crossprod(c, K[take,take]%*%c))
  }
  take = suppressMessages(
    selectOCS(s$pop, nFemale=30, nMale=30, use="gv", usePed=TRUE,
              returnPop=FALSE, simParam=s$SP))
  female = take[1:30]
  male = take[31:60]
  expect_true(all(s$pop@sex[female]=="F"))
  expect_true(all(s$pop@sex[male]=="M"))
  expect_lte(poolFix(female, male), target)
  truncF = selectInd(s$pop, nInd=30, use="gv", sex="F",
                     returnPop=FALSE, simParam=s$SP)
  truncM = selectInd(s$pop, nInd=30, use="gv", sex="M",
                     returnPop=FALSE, simParam=s$SP)
  expect_gt(poolFix(truncF, truncM), target)

  # selectCross passes usePed on and crosses the same individuals
  pop2 = suppressMessages(
    selectCross(s$pop, nFemale=30, nMale=30, nCrosses=30, use="gv",
                restrInbr=TRUE, usePed=TRUE, simParam=s$SP))
  expect_setequal(pop2@mother, s$pop@id[female])
  expect_setequal(pop2@father, s$pop@id[male])
})

test_that("MEASURE every set of unrelated founders alike by the pedigree",{
  set.seed(12)
  founder = quickHaplo(nInd=10, nChr=1, segSites=20)
  SP = SimParam$new(founder)
  SP$nThreads = 1L
  SP$addTraitA(10)
  SP$setTrackPed(TRUE)
  pop = newPop(founder, simParam=SP)
  # Five unrelated, outbred diploids always have an average kinship of
  # 1/(2*5), so that target is met by truncation and a lower one is not
  # met by any set
  take = suppressMessages(
    selectOCS(pop, nInd=5, use="gv", inbrTarget=0.1,
              inbrType="absolute", usePed=TRUE, returnPop=FALSE,
              simParam=SP))
  expect_setequal(take, selectInd(pop, nInd=5, use="gv",
                                  returnPop=FALSE, simParam=SP))
  expect_warning(suppressMessages(
    selectOCS(pop, nInd=5, use="gv", inbrTarget=0.09,
              inbrType="absolute", usePed=TRUE, simParam=SP)),
    regexp="target was not met")
})

test_that("SPREAD greedy pedigree selection across full-sib families",{
  # Three families of four outbred diploid full sibs from unrelated
  # parents, with the best merit in family 1 and the worst in family 3
  fam = rep(1:3, each=4)
  K = ifelse(outer(fam, fam, "=="), 1/4, 0)
  diag(K) = 1/2
  m = c(10:7, 6:3, 2:-1)
  ker = AlphaSimR:::.ocsKernelK(K, 1, 1)
  # Without a penalty the selection is truncation
  ans = AlphaSimR:::.ocsGreedy(0, ker, m, 3L)
  expect_identical(ans$take, 1:3)
  expect_equal(ans$F, mean(K[1:3,1:3]))
  # With a large penalty each family gives one, its best, which is the
  # least related set of three
  ans = AlphaSimR:::.ocsGreedy(1e6, ker, m, 3L)
  expect_identical(sort(ans$take), c(1L, 5L, 9L))
  expect_equal(ans$F, 1/6)

  # With sexes, females are odd and males even within each family, so
  # one of each sex is taken from different families
  female = which(seq_along(m)%%2L==1L)
  male = which(seq_along(m)%%2L==0L)
  ans = AlphaSimR:::.ocsGreedySex(1e6, ker, m, female, male, 1L, 1L)
  expect_identical(ans$take[1], 1L)
  expect_true(ans$take[2] %in% male[fam[male]!=1L])
  expect_equal(ans$F, ker$fixPool(ans$take[1], ans$take[2]))
  # A free candidate with a missing merit is still taken when no other
  # candidate is free
  ans = AlphaSimR:::.ocsGreedy(0, ker, c(m[-12], NA), 12L)
  expect_setequal(ans$take, 1:12)
})

test_that("REFUSE the pedigree without tracking or with useQtl",{
  set.seed(13)
  s = ocsSetup()
  expect_error(selectOCS(s$pop, nInd=10, use="gv", usePed=TRUE,
                         simParam=s$SP),
               regexp="trackPed")
  expect_error(selectCross(s$pop, nInd=10, nCrosses=5, use="gv",
                           restrInbr=TRUE, usePed=TRUE, simParam=s$SP),
               regexp="trackPed")
  s = ocsSetup(trackPed=TRUE)
  expect_error(selectOCS(s$pop, nInd=10, use="gv", usePed=TRUE,
                         useQtl=TRUE, simParam=s$SP),
               regexp="cannot both be TRUE")
})

test_that("AGREE between the genotype kernel and a kinship kernel",{
  set.seed(14)
  for(ploidy in c(2L, 4L)){
    s = ocsSetup(ploidy=ploidy)
    Y = 2*pullSnpGeno(s$pop, simParam=s$SP) - ploidy
    Z = Y/ploidy
    nLoci = ncol(Z)
    kerY = AlphaSimR:::.ocsKernel(Y, ploidy)
    # The kinship that the genotypes estimate, given as a pedigree's
    # would be
    kerA = AlphaSimR:::.ocsKernelK(tcrossprod(Z)/nLoci, 1, 1)
    take = c(3, 17, 40, 41, 150)
    takeM = c(5, 60, 199)
    # The penalties are per locus rather than summed over loci
    expect_equal(kerA$nLoci, 1)
    expect_equal(kerA$zu0*nLoci, kerY$zu0)
    expect_equal(kerA$self*nLoci, kerY$self)
    expect_equal(kerA$zu(take)*nLoci, kerY$zu(take))
    expect_equal(kerA$fix(take), kerY$fix(take))
    expect_equal(kerA$fixPool(take, takeM), kerY$fixPool(take, takeM))
    # so the same selection needs a penalty nLoci times larger
    m = gv(s$pop)[,1]
    ansY = AlphaSimR:::.ocsTrunc(0.05, kerY, m, 50)
    ansA = AlphaSimR:::.ocsTrunc(0.05*nLoci, kerA, m, 50)
    expect_identical(ansA$take, ansY$take)
    expect_equal(ansA$F, ansY$F)
  }
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
                      SP$femaleMap, SP$maleMap,
                      SP$femaleGenLen, SP$maleGenLen, FALSE, 2L, 2L,
                      SP$v, SP$p, SP$femaleCentromere, SP$maleCentromere,
                      SP$quadProb, 1L)
  }
  expect_error(doCross(c(1L,0L), c(1L,2L)), regexp = "Invalid parent index")
  expect_error(doCross(c(1L,2L), c(1L,5L)), regexp = "Invalid parent index")
  expect_error(doCross(c(1L,-1L), c(1L,2L)), regexp = "Invalid parent index")
  expect_error(doCross(c(1L,2L), 1L), regexp = "same length")
})

test_that("All noninterfering quadrivalent crossovers match the Poisson pathway", {
  set.seed(6016)
  founders = quickHaplo(nInd=4,nChr=2,segSites=40,ploidy=4)
  sp = SimParam$new(founders)
  sp$nThreads = 1L
  sp$quadProb = 1
  pop = newPop(founders,simParam=sp)
  sp$p = 1
  sp$v = 2.6
  set.seed(6017)
  allType2 = makeCross(pop,crossPlan=matrix(c(1,2,3,4),2,2),simParam=sp)
  # With no interference, the equivalent type-1 renewal process has v=1.
  sp$p = 0
  sp$v = 1
  set.seed(6017)
  poisson = makeCross(pop,crossPlan=matrix(c(1,2,3,4),2,2),simParam=sp)
  expect_equal(allType2@ploidy,4L)
  expect_equal(allType2@nInd,2L)
  expect_identical(allType2@geno,poisson@geno)
})

test_that("Mixed quadrivalent pathways preserve fixed alleles on long chromosomes", {
  # Long arms exhaust the initial batch of gamma draws in both pathways.
  h = matrix(rep(c(0,1),each=16),nrow=16,ncol=40)
  founders = newMapPop(list(seq(0.1,30,length.out=40)),list(h),ploidy=4,genLen=31)
  sp = SimParam$new(founders)
  sp$nThreads = 1L
  sp$quadProb = 1
  sp$p = 0.5
  set.seed(6020)
  pop = newPop(founders,simParam=sp)
  kids = makeCross(pop,cbind(1:4,4:1),nProgeny=4,simParam=sp)
  expect_equal(unname(pullSegSiteHaplo(kids,simParam=sp)),
               matrix(rep(c(0,1),each=64),nrow=64,ncol=40))
})

test_that("Meiosis kernels validate dimensions and normalize thread requests", {
  set.seed(6021)
  f = quickHaplo(2,1,10)
  sp = SimParam$new(f)
  crossArgs = list(motherGeno=f@geno,mother=1L,fatherGeno=f@geno,father=2L,
    femaleMap=sp$femaleMap,maleMap=sp$maleMap,
    femaleGenLen=sp$femaleGenLen,maleGenLen=sp$maleGenLen,
    trackRec=FALSE,motherPloidy=2L,fatherPloidy=2L,v=sp$v,p=sp$p,
    motherCentromere=sp$femaleCentromere,fatherCentromere=sp$maleCentromere,
    quadProb=sp$quadProb,nThreads=1L)
  dhArgs = list(geno=f@geno,nDH=c(1L,1L),genMap=f@genMap,genLen=f@genLen,
                v=sp$v,p=sp$p,trackRec=FALSE,nThreads=1L)
  rgArgs = list(geno=f@geno,nProgeny=c(1L,1L),genMap=f@genMap,genLen=f@genLen,
                v=sp$v,p=sp$p,trackRec=FALSE,ploidy=2L,
                centromere=f@centromere,quadProb=sp$quadProb,nThreads=1L)
  badCross = crossArgs
  badCross$fatherGeno = c(f@geno,f@geno)
  expect_error(do.call(AlphaSimR:::cross,badCross),
               "same chromosomes")
  expect_error(do.call(AlphaSimR:::cross,modifyList(crossArgs,list(femaleGenLen=numeric()))),
               "genLen must")
  expect_error(do.call(AlphaSimR:::createDH2,modifyList(dhArgs,list(genLen=numeric()))),
               "genLen must")
  expect_error(do.call(AlphaSimR:::createDH2,modifyList(dhArgs,list(nDH=1L))),
               "nDH must")
  expect_error(do.call(AlphaSimR:::createReducedGenome,modifyList(rgArgs,list(genLen=numeric()))),
               "genLen must")
  expect_error(do.call(AlphaSimR:::createReducedGenome,modifyList(rgArgs,list(nProgeny=1L))),
               "nProgeny must")
  for(case in list(list(fun=AlphaSimR:::cross,args=crossArgs),
                  list(fun=AlphaSimR:::createDH2,args=dhArgs),
                  list(fun=AlphaSimR:::createReducedGenome,args=rgArgs))){
    set.seed(6022)
    ref = do.call(case$fun,case$args)
    for(threads in c(0L,4L)){
      set.seed(6022)
      expect_identical(do.call(case$fun,modifyList(case$args,list(nThreads=threads))),ref)
    }
  }
})

test_that("Meiosis kernels check maps, centromeres and ploidy before the parallel region", {
  # Everything below is read on worker threads through raw pointers, where a
  # length that does not match would be read out of bounds rather than
  # reported
  set.seed(6023)
  f = quickHaplo(2,1,10)
  sp = SimParam$new(f)
  crossArgs = list(motherGeno=f@geno,mother=1L,fatherGeno=f@geno,father=2L,
    femaleMap=sp$femaleMap,maleMap=sp$maleMap,
    femaleGenLen=sp$femaleGenLen,maleGenLen=sp$maleGenLen,
    trackRec=FALSE,motherPloidy=2L,fatherPloidy=2L,v=sp$v,p=sp$p,
    motherCentromere=sp$femaleCentromere,fatherCentromere=sp$maleCentromere,
    quadProb=sp$quadProb,nThreads=1L)
  dhArgs = list(geno=f@geno,nDH=c(1L,1L),genMap=f@genMap,genLen=f@genLen,
                v=sp$v,p=sp$p,trackRec=FALSE,nThreads=1L)
  rgArgs = list(geno=f@geno,nProgeny=c(1L,1L),genMap=f@genMap,genLen=f@genLen,
                v=sp$v,p=sp$p,trackRec=FALSE,ploidy=2L,
                centromere=f@centromere,quadProb=sp$quadProb,nThreads=1L)
  # Replaces arguments whole. modifyList would merge a list argument, such
  # as a map or genotypes, into the one already there.
  setArgs = function(args, new){
    args[names(new)] = new
    return(args)
  }
  twoMaps = c(f@genMap,f@genMap)
  longMap = list(seq(0,1,length.out=100))

  expect_error(do.call(AlphaSimR:::cross,setArgs(crossArgs,list(maleMap=twoMaps))),
               "genMap must")
  expect_error(do.call(AlphaSimR:::cross,setArgs(crossArgs,list(femaleMap=longMap))),
               "more loci")
  # 20 loci take three bytes, 10 take two
  g = quickHaplo(2,1,20)
  expect_error(do.call(AlphaSimR:::cross,setArgs(crossArgs,list(fatherGeno=g@geno))),
               "same loci")
  expect_error(do.call(AlphaSimR:::cross,setArgs(crossArgs,list(motherPloidy=4L))),
               "ploidy does not match")
  expect_error(do.call(AlphaSimR:::createDH2,setArgs(dhArgs,list(genMap=twoMaps))),
               "genMap must")
  expect_error(do.call(AlphaSimR:::createDH2,setArgs(dhArgs,list(genMap=longMap))),
               "more loci")
  expect_error(do.call(AlphaSimR:::createReducedGenome,setArgs(rgArgs,list(genMap=twoMaps))),
               "genMap must")
  expect_error(do.call(AlphaSimR:::createReducedGenome,setArgs(rgArgs,list(ploidy=4L))),
               "ploidy does not match")
  # A diploid never reads its centromeres, so they need not be given
  expect_silent(do.call(AlphaSimR:::cross,
                        setArgs(crossArgs,list(motherCentromere=numeric(),
                                                  fatherCentromere=numeric()))))

  # A tetraploid can form quadrivalents, which read them
  f4 = quickHaplo(2,1,10,ploidy=4)
  sp4 = SimParam$new(f4)
  crossArgs4 = setArgs(crossArgs,list(motherGeno=f4@geno,fatherGeno=f4@geno,
    motherPloidy=4L,fatherPloidy=4L,quadProb=1))
  expect_error(do.call(AlphaSimR:::cross,setArgs(crossArgs4,list(fatherCentromere=numeric()))),
               "centromere must")
  rgArgs4 = setArgs(rgArgs,list(geno=f4@geno,ploidy=4L,quadProb=1))
  expect_error(do.call(AlphaSimR:::createReducedGenome,setArgs(rgArgs4,list(centromere=numeric()))),
               "centromere must")
  expect_silent(do.call(AlphaSimR:::cross,crossArgs4))
})

# Meiosis reads the parents' genotypes in place, from the same memory as the
# population passed in, and writes the progeny into arrays that are returned
# to R as they are. The blocks below check both ends of that.

test_that("KEEP the parents' genotypes unchanged by crossing",{
  set.seed(7101)
  # A deep copy, since an ordinary copy in R shares the parents' memory
  deepCopy = function(x) unserialize(serialize(x, NULL))
  founder2 = quickHaplo(nInd=20, nChr=2, segSites=37)
  SP = SimParam$new(founder2)
  SP$nThreads = 2L
  SP$setTrackRec(TRUE)
  pop = newPop(founder2, simParam=SP)
  before = deepCopy(pop@geno)
  invisible(makeCross(pop, cbind(1:10, 11:20), simParam=SP))
  invisible(randCross(pop, 30, simParam=SP))
  invisible(self(pop, nProgeny=2, simParam=SP))
  invisible(makeDH(pop, nDH=2, simParam=SP))
  invisible(makeCross2(pop[1:10], pop[11:20], cbind(1:10, 10:1),
                       simParam=SP))
  expect_identical(pop@geno, before)

  founder4 = quickHaplo(nInd=10, nChr=2, segSites=37, ploidy=4L)
  SP4 = SimParam$new(founder4)
  SP4$nThreads = 2L
  pop4 = newPop(founder4, simParam=SP4)
  before4 = deepCopy(pop4@geno)
  invisible(reduceGenome(pop4, nProgeny=3, simParam=SP4))
  invisible(randCross(pop4, 10, simParam=SP4))
  expect_identical(pop4@geno, before4)
})

test_that("PASS inbred parents' haplotypes to their progeny unchanged",{
  set.seed(7102)
  # An inbred parent's gametes are its haplotype whatever the crossovers,
  # so every haplotype of the progeny is known in advance. This checks that
  # each progeny's genotypes are written where they belong.
  founderIn = quickHaplo(nInd=70, nChr=3, segSites=29, inbred=TRUE)
  SP = SimParam$new(founderIn)
  SP$nThreads = 2L
  pop = newPop(founderIn, simParam=SP)
  H = pullSegSiteHaplo(pop, haplo=1, simParam=SP)
  plan = cbind(sample(70, 90, replace=TRUE), sample(70, 90, replace=TRUE))
  kids = makeCross(pop, plan, simParam=SP)
  expect_true(all(pullSegSiteHaplo(kids, haplo=1, simParam=SP) == H[plan[,1],]))
  expect_true(all(pullSegSiteHaplo(kids, haplo=2, simParam=SP) == H[plan[,2],]))
  dh = makeDH(pop, nDH=2, simParam=SP)
  expect_true(all(pullSegSiteHaplo(dh, haplo=1, simParam=SP) ==
                    H[rep(1:70, each=2),]))
  expect_true(all(pullSegSiteHaplo(dh, haplo=2, simParam=SP) ==
                    H[rep(1:70, each=2),]))
  selfed = self(pop[c(3, 9)], nProgeny=4, simParam=SP)
  expect_true(all(pullSegSiteGeno(selfed, simParam=SP) ==
                    2L*H[rep(c(3, 9), each=4),]))
})

test_that("REFUSE parent genotypes that are not packed raw arrays",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop,simParam=SP)
  doCross = function(geno){
    AlphaSimR:::cross(geno, 1L, pop@geno, 2L,
                      SP$femaleMap, SP$maleMap,
                      SP$femaleGenLen, SP$maleGenLen, FALSE, 2L, 2L,
                      SP$v, SP$p, SP$femaleCentromere, SP$maleCentromere,
                      SP$quadProb, 1L)
  }
  notRaw = lapply(pop@geno, function(x) array(as.integer(x), dim=dim(x)))
  expect_error(doCross(notRaw), regexp = "raw vectors")
  flat = lapply(pop@geno, as.vector)
  expect_error(doCross(flat), regexp = "three dimensional")
})
