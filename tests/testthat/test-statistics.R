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
