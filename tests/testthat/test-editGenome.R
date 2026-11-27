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

test_that("editGenomeTopQtl decreases trait values and skips fixed QTL",{
  # 3 QTL on chr 1 with known dosages for 2 individuals:
  # ind 1 (inbred): dosages 2, 0, 2
  # ind 2 (mixed):  dosages 0, 1, 2
  haplo = rbind(
    c(1, 0, 1),
    c(1, 0, 1),
    c(0, 1, 1),
    c(0, 0, 1)
  )
  founderPop = newMapPop(list(c(0, 0.5, 1)), list(haplo))
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$manAddTrait(new("TraitA", nLoci=3L, lociPerChr=3L, lociLoc=1:3,
                     addEff=c(3, -2, 1), intercept=0))
  pop = newPop(founderPop, simParam=SP)

  # Decreasing ind 1 by 2 QTL must edit locus 1 (2 -> 0, eff +3) and
  # locus 2 (0 -> 2, eff -2), leaving locus 3 unchanged.
  dec1 = editGenomeTopQtl(pop, ind=1, nQtl=2, increase=FALSE, simParam=SP)
  expect_equal(unname(pullQtlGeno(dec1, simParam=SP)[1, ]), c(0L, 2L, 2L))
  expect_lt(gv(dec1)[1, 1], gv(pop)[1, 1])

  # Decreasing ind 2 by 2 QTL must skip locus 1 (already 0 for eff +3) and
  # edit locus 2 (1 -> 2, eff -2) and locus 3 (2 -> 0, eff +1).
  dec2 = editGenomeTopQtl(pop, ind=2, nQtl=2, increase=FALSE, simParam=SP)
  expect_equal(unname(pullQtlGeno(dec2, simParam=SP)[2, ]), c(0L, 2L, 0L))
  expect_lt(gv(dec2)[2, 1], gv(pop)[2, 1])
})

test_that("editGenomeTopQtl handles polyploid dosages up to ploidy",{
  # Tetraploid with dosages 2, 3, 4 across 3 QTL (all positive effects)
  haplo = rbind(
    c(1, 1, 1),
    c(1, 1, 1),
    c(0, 1, 1),
    c(0, 0, 1)
  )
  founderPop = newMapPop(list(c(0, 0.5, 1)), list(haplo), ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$manAddTrait(new("TraitA", nLoci=3L, lociPerChr=3L, lociLoc=1:3,
                     addEff=c(3, 2, 1), intercept=0))
  pop = newPop(founderPop, simParam=SP)

  # Increasing by 2 QTL must edit loci 1 (dosage 2 -> 4) and 2 (dosage 3 -> 4),
  # not skip them because dosage >= 2.
  inc = editGenomeTopQtl(pop, ind=1, nQtl=2, increase=TRUE, simParam=SP)
  expect_equal(unname(pullQtlGeno(inc, simParam=SP)[1, ]), c(4L, 4L, 4L))
  expect_gt(gv(inc)[1, 1], gv(pop)[1, 1])

  # Decreasing by 1 QTL must edit locus 1 (dosage 2 -> 0).
  dec = editGenomeTopQtl(pop, ind=1, nQtl=1, increase=FALSE, simParam=SP)
  expect_equal(unname(pullQtlGeno(dec, simParam=SP)[1, ]), c(0L, 3L, 4L))
  expect_lt(gv(dec)[1, 1], gv(pop)[1, 1])
})
