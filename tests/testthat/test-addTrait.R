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

# Correlation arguments must be valid correlation matrices. A covariance
# matrix used to be accepted silently, and for gamma distributed effects it
# distorted the marginal distribution.
test_that("REJECT correlation arguments that are not correlation matrices",{
  checkCorMat = AlphaSimR:::checkCorMat
  R = matrix(c(1,0.5,0.5,1),nrow=2)
  expect_silent(checkCorMat(R,2))
  expect_silent(checkCorMat(diag(3),3))
  expect_error(checkCorMat(1,1),"numeric matrix")
  expect_error(checkCorMat(matrix("a"),1),"numeric matrix")
  expect_error(checkCorMat(R,3),"3 by 3")
  expect_error(checkCorMat(matrix(c(1,NA,NA,1),nrow=2),2),"missing")
  expect_error(checkCorMat(matrix(c(1,0.5,0.4,1),nrow=2),2),"symmetric")
  expect_error(checkCorMat(2*R,2),"diagonal")
  expect_error(checkCorMat(matrix(c(1,1.5,1.5,1),nrow=2),2),"outside")

  # Every correlation argument of every addTrait function is checked
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  bad = 2*R
  args = list(addTraitA=c("corA"),
              addTraitAD=c("corA","corDD"),
              addTraitAG=c("corA","corGxE"),
              addTraitADG=c("corA","corDD","corGxE"),
              addTraitAE=c("corA","corAA"),
              addTraitADE=c("corA","corDD","corAA"),
              addTraitAEG=c("corA","corAA","corGxE"),
              addTraitADEG=c("corA","corDD","corAA","corGxE"))
  for(f in names(args)){
    for(a in args[[f]]){
      call = list(nQtlPerChr=2,mean=c(0,0),var=c(1,1))
      if(grepl("G",f)) call$varGxE = c(1,1)
      call[[a]] = bad
      expect_error(do.call(SP[[f]],call),paste(a,"must have ones"),
                   info=paste(f,a))
    }
  }
})

test_that("LATENT correlation leaves normal traits and zero targets alone",{
  R = matrix(c(1,0.5,0,
               0.5,1,-0.3,
               0,-0.3,1),nrow=3)
  # All normal: nothing to adjust
  expect_identical(AlphaSimR:::gammaLatentCorr(R,rep(FALSE,3),rep(1,3)),R)
  # Gamma traits: zero stays zero and the adjustment moves away from zero
  Rz = AlphaSimR:::gammaLatentCorr(R,c(TRUE,TRUE,FALSE),c(0.5,0.5,1))
  expect_equal(Rz[1,3],0)
  expect_equal(diag(Rz),rep(1,3))
  expect_true(isSymmetric(Rz))
  expect_true(Rz[1,2]>0.5)
  expect_true(Rz[2,3]< -0.3)
  # Two traits with the same shape can be perfectly correlated
  R1 = matrix(1,nrow=2,ncol=2)
  expect_equal(AlphaSimR:::gammaLatentCorr(R1,c(TRUE,TRUE),c(0.2,0.2)),R1)
})

test_that("WARN when a correlation cannot be reached",{
  # A normal trait and a gamma trait with shape 0.2 correlate at most 0.78
  R = matrix(c(1,0.9,0.9,1),nrow=2)
  expect_warning(Rz <- AlphaSimR:::gammaLatentCorr(R,c(FALSE,TRUE),c(1,0.2)),
                 "cannot be reached")
  expect_equal(Rz[1,2],1)
})

# The gamma transform shrinks correlations towards zero, by 0.15 at shape 0.2
# and a target of 0.6. With 2e5 loci the sampling standard deviation of the
# correlation is about 0.003, so a tolerance of 0.02 separates the corrected
# and uncorrected methods while almost never failing by chance.
test_that("GAMMA effects reach the requested correlation",{
  skip_on_cran()
  set.seed(4417)
  nLoci = 200000L
  qtlLoci = new("LociMap",nLoci=nLoci,lociPerChr=nLoci,lociLoc=seq_len(nLoci))
  # The adjusted latent matrix must itself be a valid correlation matrix for
  # the target to be reached exactly. With a correlation of -0.4 between the
  # first and third traits it is not, and transMat's smoothing leaves the
  # result about 0.03 away.
  R = matrix(c(1,0.6,-0.3,
               0.6,1,0.2,
               -0.3,0.2,1),nrow=3)
  gamma = c(TRUE,TRUE,FALSE)
  shape = c(0.2,0.5,1)

  addEff = AlphaSimR:::sampAddEff(qtlLoci=qtlLoci,nTraits=3,corr=R,
                                  gamma=gamma,shape=shape)
  expect_lt(max(abs(cor(addEff)-R)),0.02)

  epiEff = AlphaSimR:::sampEpiEff(qtlLoci=qtlLoci,nTraits=3,corr=R,
                                  gamma=gamma,shape=shape,relVar=rep(1,3))
  expect_lt(max(abs(cor(epiEff)-R)),0.02)
})
