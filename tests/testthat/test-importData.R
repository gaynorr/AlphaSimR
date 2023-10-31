context("importData")

test_that("importTrait",{
  # Create haplotype data
  haplo = rbind(c(1,1,0,1,0),
                c(1,1,0,1,0),
                c(0,1,1,0,0),
                c(0,1,1,0,0))
  colnames(haplo) = letters[1:5]
  
  # Create genetic map
  genMap = data.frame(markerName=letters[1:5],
                      chromosome=c(1,1,1,2,2),
                      position=c(0,0.5,1,0.15,0.4))
  
  # Create pedigree
  ped = data.frame(id=c("a","b"),
                   mother=c(0,0),
                   father=c(0,0))
  
  # Generate an external trait
  myTrait = data.frame(marker = c("a","c","d"),
                       a = c(1,-1,1))
  
  founderPop = importHaplo(haplo = haplo, 
                           genMap = genMap,
                           ploidy = 2L,
                           ped = ped)
  
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L
  
  # Import trait
  SP$importTrait(markerNames = myTrait$marker,
                 addEff = myTrait$a,
                 name = "myTrait")
  
  expect_equal(SP$traits[[1]]@addEff, myTrait$a, tolerance=1e-6)
  
  expect_equal(SP$traits[[1]]@intercept, 0, tolerance=1e-6)
  
  pop = newPop(founderPop,simParam=SP)
  
  expect_equal(unname(pop@gv[1,1]), 3, tolerance=1e-6)
  
  expect_equal(unname(pop@gv[2,1]), -3, tolerance=1e-6)
})

