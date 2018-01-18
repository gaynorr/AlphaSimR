# These tests will fail with a small probability and/or 
# they are slower tests, so they are skipped on CRAN.
# If any of these test fail, rerun the tests multiple times.
# Frequent failures indicate a problem.
context("statistics")

test_that("addError",{
  skip_on_cran()
  gv = matrix(0,nrow=10000,ncol=2)
  varE = c(1,1)
  pheno = AlphaSimR:::addError(gv=gv,varE=varE,reps=1)
  expect_equal(var(pheno),diag(varE),tol=0.1)
  varE = diag(2)
  pheno = AlphaSimR:::addError(gv=gv,varE=varE,reps=1)
  expect_equal(var(pheno),varE,tol=0.1)
  pheno = AlphaSimR:::addError(gv=gv,varE=varE,reps=4)
  expect_equal(var(pheno),varE/4,tol=0.1)
  varE = 0.5*diag(2)+0.5
  pheno = AlphaSimR:::addError(gv=gv,varE=varE,reps=1)
  expect_equal(var(pheno),varE,tol=0.1)
  pheno = AlphaSimR:::addError(gv=gv,varE=varE,reps=4)
  expect_equal(var(pheno),varE/4,tol=0.1)
  varE = 1.5*diag(2)-0.5
  pheno = AlphaSimR:::addError(gv=gv,varE=varE,reps=1)
  expect_equal(var(pheno),varE,tol=0.1)
  pheno = AlphaSimR:::addError(gv=gv,varE=varE,reps=4)
  expect_equal(var(pheno),varE/4,tol=0.1)
  varE = matrix(c(1,0.5,-0.5,1),ncol=2)
  expect_error(AlphaSimR:::addError(gv=gv,varE=varE,reps=1))
})

test_that("Univariate GS",{
  # All mixed model solvers should give nearly the
  # same prediction
  skip_on_cran()
  capture_output((founderPop = runMacs(nInd=10,nChr=10,segSites=100)))
  SP = SimParam$new(founderPop)
  SP$addTraitA(100,0,1)
  SP$setVarE(varE=10)
  pop = newPop(founderPop,simParam=SP)
  pop = randCross(pop,1000,simParam=SP)
  y = pop@pheno
  X = matrix(1,nrow=pop@nInd)
  Z = diag(pop@nInd)
  M = pullQtlGeno(pop,simParam=SP)
  G = calcG(M)
  M = scale(M,scale=FALSE)
  K = diag(ncol(M))
  # GBLUP models
  ansAni = solveAniModel(y,X,G)
  ansAni = c(ansAni$u)
  ansAniUVM = solveUVM(y,X,Z,G)
  ansAniUVM = c(ansAniUVM$u)
  ansAniMVM = solveMVM(y,X,Z,G)
  ansAniMVM = c(ansAniMVM$u)
  ansAniMKM = solveMKM(y,X,list(Z),list(G))
  ansAniMKM = c(ansAniMKM$u[[1]])
  # Compare models
  # tolerance based on simularity of models
  expect_equal(ansAni,ansAniUVM,tol=1e-4)
  expect_equal(ansAni,ansAniMVM,tol=1e-2)
  expect_equal(ansAni,ansAniMKM,tol=1e-3)
  # RR-BLUP models
  ansRR = solveRRBLUP(y,X,M)
  ansRR = c(M%*%ansRR$u)
  ansRRUVM = solveUVM(y,X,M,K)
  ansRRUVM = c(M%*%ansRRUVM$u)
  ansRRMVM = solveMVM(y,X,M,K)
  ansRRMVM = c(M%*%ansRRMVM$u)
  ansRRMKM = solveMKM(y,X,list(M),list(K))
  ansRRMKM = c(M%*%ansRRMKM$u[[1]])
  # Compare models
  expect_equal(ansRR,ansRRUVM,tol=1e-4)
  expect_equal(ansRR,ansRRMVM,tol=1e-2)
  expect_equal(ansRR,ansRRMKM,tol=1e-3)
  
  # Compare different types of models
  expect_equal(ansAni,ansRR,tol=1e-4)
  expect_equal(ansAniUVM,ansRRUVM,tol=1e-4)
  expect_equal(ansAniMVM,ansRRMVM,tol=1e-4)
  expect_equal(ansAniMKM,ansRRMKM,tol=1e-4)
})
