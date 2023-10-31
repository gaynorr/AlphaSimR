context("addTrait")

#Population with 2 individuals, 1 chromosome and 1 QTL
#Population is fully inbred and p=q=0.5
founderPop = newMapPop(list(c(0)),
                       list(matrix(c(1,1,0,0),
                                   nrow=4,ncol=1)))

test_that("addTraitA",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),0,tolerance=1e-6)
  expect_equal(SP$varA,1,tolerance=1e-6)
  expect_equal(SP$varG,1,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(unname(c(ans$varA)),1,tolerance=1e-6)
  expect_equal(unname(c(ans$varD)),0,tolerance=1e-6)
  expect_equal(unname(c(ans$varG)),1,tolerance=1e-6)
  expect_equal(unname(ans$genicVarA),0.5,tolerance=1e-6)
  expect_equal(unname(ans$genicVarD),0,tolerance=1e-6)
  expect_equal(unname(ans$genicVarG),0.5,tolerance=1e-6)
})

test_that("addTraitAD",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitAD(nQtlPerChr=1,mean=0,var=1,meanDD=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(SP$traits[[1]]@domEff,1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),0,tolerance=1e-6)
  expect_equal(SP$varA,1,tolerance=1e-6)
  expect_equal(SP$varG,1,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(unname(c(ans$varA)),1,tolerance=1e-6)
  expect_equal(unname(c(ans$varD)),0,tolerance=1e-6)
  expect_equal(unname(c(ans$varG)),1,tolerance=1e-6)
  expect_equal(unname(ans$genicVarA),0.5,tolerance=1e-6)
  expect_equal(unname(ans$genicVarD),0.25,tolerance=1e-6)
  expect_equal(unname(ans$genicVarG),0.75,tolerance=1e-6)
})

test_that("addTraitAG",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitAG(nQtlPerChr=1,mean=0,var=1,varEnv=1,varGxE=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@gxeEff),1,tolerance=1e-6)
  expect_equal(SP$traits[[1]]@envVar,1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@gxeInt-1),0,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),0,tolerance=1e-6)
})

test_that("addTraitADG",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitADG(nQtlPerChr=1,mean=0,var=1,meanDD=1,varEnv=1,varGxE=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(SP$traits[[1]]@domEff,1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@gxeEff),1,tolerance=1e-6)
  expect_equal(SP$traits[[1]]@envVar,1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@gxeInt-1),0,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),0,tolerance=1e-6)
})

#Population with 2 individuals, 1 chromosome and 2 QTL
#Population is fully inbred and p=q=0.5
founderPop = newMapPop(list(c(0,0)),
                       list(cbind(c(1,1,0,0),c(1,1,0,0))))

test_that("addTraitAE",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitAE(nQtlPerChr=2,mean=0,var=1,relAA=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(SP$varA,1,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(unname(c(ans$varA)),1,tolerance=1e-6)
})

test_that("addTraitADE",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitADE(nQtlPerChr=2,mean=0,var=1,meanDD=1,relAA=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(SP$varA,1,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(unname(c(ans$varA)),1,tolerance=1e-6)
})

test_that("addTraitAEG",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitAEG(nQtlPerChr=2,mean=0,var=1,varEnv=1,varGxE=1,relAA=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(SP$varA,1,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(unname(c(ans$varA)),1,tolerance=1e-6)
})

test_that("addTraitADEG",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitADEG(nQtlPerChr=2,mean=0,var=1,meanDD=1,varEnv=1,varGxE=1,relAA=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(SP$varA,1,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(unname(c(ans$varA)),1,tolerance=1e-6)
})
