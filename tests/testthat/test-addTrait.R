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

#Population with 4 individuals, 1 chromosome and 1 QTL and p=q=0.5
founderPop = newMapPop(list(c(0)),
                       list(matrix(c(0,0,0,1,1,0,1,1),
                                   nrow=8,ncol=1)))

test_that("addTraitAI",{
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitAI(nQtlPerChr=1,mean=0,var=1,meanID=1)
  pop = newPop(founderPop,simParam=SP)
  # Just one locus and complete alias between a and i, will give varG=2, so
  # scaled effect will be sqrt(varG)
  expect_equal(abs(SP$traits[[1]]@addEff),sqrt(2),tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@impEff),sqrt(2),tolerance=1e-6)
  expect_equal(abs(SP$traits[[1]]@intercept),0,tolerance=1e-6)
  expect_equal(SP$varA,1,tolerance=1e-6)
  expect_equal(SP$varG,2,tolerance=1e-6)
  ans = genParam(pop,simParam=SP)
  expect_equal(unname(c(ans$varA)),1,tolerance=1e-6)
  expect_equal(unname(c(ans$varI)),1,tolerance=1e-6)
  expect_equal(unname(c(ans$varG)),2,tolerance=1e-6)
  expect_equal(unname(ans$genicVarA),1,tolerance=1e-6)
  expect_equal(unname(ans$genicVarI),1,tolerance=1e-6)
  expect_equal(unname(ans$genicVarG),2,tolerance=1e-6)
  gv_a = c(-SP$traits[[1]]@addEff,0,0,SP$traits[[1]]@addEff)
  gv_i = c(0,SP$traits[[1]]@impEff,-SP$traits[[1]]@impEff,0)
  idM = c(+SP$traits[[1]]@impEff, 0, 0, -SP$traits[[1]]@impEff)
  # +2pi, (p-q)i, (p-q)i, -2qi --> with p=q=0.5 --> +i, 0, 0, -i
  idP = c(-SP$traits[[1]]@impEff, 0, 0, +SP$traits[[1]]@impEff)
  # -2pi, (q-p)i, (q-p)i, +2qi --> with p=q=0.5 --> -i, 0, 0, +i
  expect_equal(unname(c(ans$gv_a)),gv_a,tolerance=1e-6)
  expect_equal(unname(c(ans$gv_i)),gv_i,tolerance=1e-6)
  expect_equal(unname(c(ans$gv)),gv_a+gv_i,tolerance=1e-6)
  expect_equal(unname(c(ans$bv)),gv_a,tolerance=1e-6)
  expect_equal(unname(c(ans$idM)),idM,tolerance=1e-6)
  expect_equal(unname(c(ans$bvM)),unname(c(ans$bv))+idM,tolerance=1e-6)
  expect_equal(unname(c(ans$bvP)),unname(c(ans$bv))+idP,tolerance=1e-6)
  # Test logic of how bv() and id() interact with sex
  expect_equal(bv(pop),ans$bv,tolerance=1e-6)
  expect_equal(bvM(pop),ans$bvM,tolerance=1e-6)
  expect_equal(bvP(pop),ans$bvP,tolerance=1e-6)
  pop@sex[] = "M"
  expect_equal(bv(pop),ans$bvP,tolerance=1e-6)
  expect_equal(unname(c(id(pop))),idP,tolerance=1e-6)
  pop@sex[] = "F"
  expect_equal(bv(pop),ans$bvM,tolerance=1e-6)
  expect_equal(unname(c(id(pop))),idM,tolerance=1e-6)
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
