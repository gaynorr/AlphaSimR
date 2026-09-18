test_that("SimParam nThreads validates values and NULL resets to default", {
  founder <- quickHaplo(nInd = 2, nChr = 2, segSites = 4)
  SP <- SimParam$new(founder)
  SP$nThreads <- 1L
  pop <- newPop(founder, simParam = SP)
  expect_equal(SP$nThreads, 1L)

  SP$nThreads <- NULL
  expect_equal(SP$nThreads, getNumThreads())
  SP$nThreads <- 1L
  expect_silent(pullSegSiteGeno(pop, simParam = SP))

  expect_error(
    SP$nThreads <- 0,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- 0L,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- 1.5,
    regexp = "single positive integer or NULL to reset"
  )
  expect_error(
    SP$nThreads <- NA_integer_,
    regexp = "single positive integer or NULL to reset"
  )
})

test_that("setRecombRatio sets centromeres for both sexes",{
  founderPop = quickHaplo(nInd=4,nChr=2,segSites=8,ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L

  #Centromeres are shared before a sex-specific map is set
  expect_equal(SP$centromere,rep(0.5,2))

  SP$setRecombRatio(2) #Twice as much recombination in females

  #Both sexes must have a centromere for every chromosome
  expect_equal(length(SP$femaleCentromere),2L)
  expect_equal(length(SP$maleCentromere),2L)
  expect_equal(length(SP$centromere),2L)

  #Centromeres are scaled with the map they sit on
  expect_equal(SP$femaleCentromere,rep(2/3,2))
  expect_equal(SP$maleCentromere,rep(1/3,2))

  #A centromere must lie on its own map
  expect_true(all(SP$femaleCentromere<=sapply(SP$femaleMap,max)))
  expect_true(all(SP$maleCentromere<=sapply(SP$maleMap,max)))

  #The sex-average is unchanged by setting a ratio
  expect_equal(SP$centromere,rep(0.5,2))
  expect_equal(unname(sapply(SP$genMap,max)),rep(1,2))
})

test_that("setRecombRatio does not drift when called repeatedly",{
  founderPop = quickHaplo(nInd=4,nChr=2,segSites=8,ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L

  SP$setRecombRatio(2)
  femaleCentromere = SP$femaleCentromere
  maleCentromere = SP$maleCentromere
  femaleLen = sapply(SP$femaleMap,max)
  maleLen = sapply(SP$maleMap,max)

  SP$setRecombRatio(2)

  #The maps are rebuilt from the sex-average, so the centromeres must be too
  expect_equal(SP$femaleCentromere,femaleCentromere)
  expect_equal(SP$maleCentromere,maleCentromere)
  expect_equal(sapply(SP$femaleMap,max),femaleLen)
  expect_equal(sapply(SP$maleMap,max),maleLen)

  #The centromere stays at the same relative position on its map
  expect_equal(unname(SP$femaleCentromere/sapply(SP$femaleMap,max)),rep(0.5,2))
  expect_equal(unname(SP$maleCentromere/sapply(SP$maleMap,max)),rep(0.5,2))
})

test_that("sex-specific maps keep chromosome names",{
  founderPop = quickHaplo(nInd=4,nChr=2,segSites=8,ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L

  chrNames = names(SP$genMap)
  expect_equal(length(chrNames),2L)

  SP$setRecombRatio(2)

  #The averaged map must keep the names of the sex-specific maps,
  #because getGenMap uses them to build its chr column
  expect_equal(names(SP$genMap),chrNames)
  expect_equal(names(SP$femaleMap),chrNames)
  expect_equal(names(SP$maleMap),chrNames)

  map = getGenMap(SP)
  expect_true("chr" %in% colnames(map))
  expect_equal(nrow(map),16L)
  expect_equal(unique(map$chr),chrNames)
})

test_that("sex-specific maps work with quadrivalent pairing",{
  founderPop = quickHaplo(nInd=4,nChr=2,segSites=8,ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$quadProb = 1 #Always form quadrivalents
  SP$setRecombRatio(2)
  pop = newPop(founderPop,simParam=SP)

  #Reaches findQuadrivalentCO, which indexes the centromere vectors
  progeny = randCross(pop,nCrosses=4,simParam=SP)
  expect_equal(progeny@nInd,4L)
  expect_equal(progeny@ploidy,4L)
})

test_that("reduceGenome uses the centromeres of the map it recombines on",{
  founderPop = quickHaplo(nInd=4,nChr=2,segSites=8,ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$quadProb = 1 #Always form quadrivalents
  SP$setRecombRatio(4) #Male map is shorter than the female centromere
  pop = newPop(founderPop,simParam=SP)

  #The male centromere must be taken from the male map, not the female one
  expect_true(all(SP$maleCentromere<=sapply(SP$maleMap,max)))
  expect_true(any(SP$femaleCentromere>sapply(SP$maleMap,max)))

  female = reduceGenome(pop,nProgeny=1,useFemale=TRUE,simParam=SP)
  male = reduceGenome(pop,nProgeny=1,useFemale=FALSE,simParam=SP)
  expect_equal(female@ploidy,2L)
  expect_equal(male@ploidy,2L)
  expect_equal(male@nInd,4L)
})
