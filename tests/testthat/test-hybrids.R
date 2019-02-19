context("hybrids")

founderPop = newMapPop(list(c(0)),
                       list(matrix(c(1,1,0,0),
                                   nrow=4,ncol=1)))

test_that("hybridCross",{
  SP = SimParam$new(founderPop=founderPop)
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
