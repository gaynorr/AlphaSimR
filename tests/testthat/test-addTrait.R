context("addTrait")

#Population with 2 individuals, 1 chromosome and 1 QTL
#Population is fully inbred and p=q=0.5
founderPop = newMapPop(list(c(0)),
                       list(matrix(c(1,1,0,0),
                                   nrow=4,ncol=1)))


test_that("addTraitA",{
  SP = SimParam$new(founderPop=founderPop)
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),1,tolerance=1e-6)
  expect_equal(SP$varA,1,tolerance=1e-6)
  expect_equal(SP$varG,1,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(c(ans$varA),1,tolerance=1e-6)
  expect_equal(c(ans$varD),0,tolerance=1e-6)
  expect_equal(c(ans$varG),1,tolerance=1e-6)
  expect_equal(ans$genicVarA,1,tolerance=1e-6)
  expect_equal(ans$genicVarD,0,tolerance=1e-6)
  expect_equal(ans$genicVarG,1,tolerance=1e-6)
})

test_that("addTraitAD",{
  SP = SimParam$new(founderPop=founderPop)
  SP$addTraitAD(nQtlPerChr=1,mean=0,var=1,meanDD=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(SP$traits[[1]]@domEff,1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),1,tolerance=1e-6)
  expect_equal(SP$varA,1,tolerance=1e-6)
  expect_equal(SP$varG,1,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(c(ans$varA),1,tolerance=1e-6)
  expect_equal(c(ans$varD),0,tolerance=1e-6)
  expect_equal(c(ans$varG),1,tolerance=1e-6)
  expect_equal(ans$genicVarA,1,tolerance=1e-6)
  expect_equal(ans$genicVarD,0,tolerance=1e-6)
  expect_equal(ans$genicVarG,1,tolerance=1e-6)
})

test_that("addTraitAG",{
  SP = SimParam$new(founderPop=founderPop)
  SP$addTraitAG(nQtlPerChr=1,mean=0,var=1,varEnv=1,varGxE=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@gxeEff),1,tolerance=1e-6)
  expect_equal(SP$traits[[1]]@envVar,1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@gxeInt-1),1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),1,tolerance=1e-6)
})

test_that("addTraitADG",{
  SP = SimParam$new(founderPop=founderPop)
  SP$addTraitADG(nQtlPerChr=1,mean=0,var=1,meanDD=1,varEnv=1,varGxE=1)
  pop = newPop(founderPop,simParam=SP)
  expect_equal(abs(SP$traits[[1]]@addEff),1,tolerance=1e-6)
  expect_equal(SP$traits[[1]]@domEff,1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@gxeEff),1,tolerance=1e-6)
  expect_equal(SP$traits[[1]]@envVar,1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@gxeInt-1),1,tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),1,tolerance=1e-6)
})
