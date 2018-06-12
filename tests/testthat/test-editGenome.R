context("editGenome")

test_that("editGenome",{
  genMap = list(c(0))
  founderPop = trackHaploPop(genMap=genMap,nInd=2,inbred=TRUE)
  SP = SimParam$new(founderPop=founderPop)
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  pop = newPop(founderPop,simParam=SP)
  pop = editGenome(pop=pop,ind=1,chr=1,segSites=1,allele=1,
                   simParam=SP)
  expect_equal(c(varG(pop)),0,tolerance=1e-6)
  pop = randCross(pop=pop,nCrosses=10,simParam=SP)
  expect_equal(c(varG(pop)),0,tolerance=1e-6)
})
