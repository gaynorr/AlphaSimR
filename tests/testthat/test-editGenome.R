context("editGenome")

test_that("editGenome",{
  founderPop = newMapPop(list(c(0)),
                         list(matrix(c(0,0,1,1),
                                     nrow=4,ncol=1)))
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=1,mean=0,var=1)
  pop = newPop(founderPop,simParam=SP)
  pop = editGenome(pop=pop,ind=1,chr=1,segSites=1,allele=1,
                   simParam=SP)
  expect_equal(c(varG(pop)),0,tolerance=1e-6)
  pop = randCross(pop=pop,nCrosses=10,simParam=SP)
  expect_equal(c(varG(pop)),0,tolerance=1e-6)
})

test_that("GENOME EDITING forwards local simulation parameters", {
  f = miscFunctionFixture()
  withWrongGlobalMiscSP(f$sp, {
    editSP = SimParam$new(f$founder)
    editSP$nThreads = 1L
    editSP$addTraitA(10)
    pop = newPop(f$founder, simParam=editSP)
    calls = 0L
    editSP$subsetMiscPop = function(pop, i, simParam){
      expect_identical(simParam, editSP)
      expect_identical(i, 1L)
      calls <<- calls + 1L
      list()
    }
    edited = editGenomeTopQtl(pop, ind=1, nQtl=1, simParam=editSP)
    expect_identical(edited@nInd, pop@nInd)
    expect_identical(calls, 1L)
  })
})

test_that("editGenomeTopQtl accepts a trait name",{
  founderPop = quickHaplo(nInd=10, nChr=2, segSites=20)
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=5, mean=c(0,0), var=c(1,1), corA=diag(2),
               name=c("Yield","Height"))
  pop = newPop(founderPop, simParam=SP)

  byNumber = editGenomeTopQtl(pop, ind=1:3, nQtl=2, trait=2,
                              simParam=SP)
  byName = editGenomeTopQtl(pop, ind=1:3, nQtl=2, trait="Height",
                            simParam=SP)
  expect_identical(pullQtlGeno(byName, trait=2, simParam=SP),
                   pullQtlGeno(byNumber, trait=2, simParam=SP))
  expect_error(editGenomeTopQtl(pop, ind=1, nQtl=1, trait="Weight",
                                simParam=SP),
               "did not match any trait names")
})
