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

# Founder population shared by the epistasis tests below. Chromosome 1
# carries markers a, b and c, and chromosome 2 carries d, e and f.
makeEpiFounder = function(){
  haplo = rbind(c(1,0,1,0,1,1),
                c(0,1,1,1,0,0),
                c(1,1,0,0,1,0),
                c(0,0,1,1,1,1),
                c(1,0,0,1,0,1),
                c(0,1,0,0,1,1),
                c(1,1,1,0,0,0),
                c(0,0,0,1,1,0))
  colnames(haplo) = letters[1:6]
  genMap = data.frame(markerName=letters[1:6],
                      chromosome=c(1,1,1,2,2,2),
                      position=c(0,0.3,0.7,0.1,0.5,0.9))
  ped = data.frame(id=c("w","x","y","z"),
                   mother=c(0,0,0,0),
                   father=c(0,0,0,0))
  importHaplo(haplo=haplo, genMap=genMap, ploidy=2L, ped=ped)
}

test_that("importTrait REORDERED markers and pairs spanning chromosomes give a TraitAE",{
  founderPop = makeEpiFounder()
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L

  # Markers are deliberately out of map order, and every pair spans
  # the two chromosomes
  markers = c("e","a","c","f","b","d")
  addEff = c(e=0.5, a=1, c=-1, f=2, b=-0.5, d=0.25)
  pairs = rbind(c("e","a"),
                c("c","f"),
                c("b","d"))
  aa = c(0.3, -0.7, 1.1)
  SP$importTrait(markerNames=markers,
                 addEff=unname(addEff),
                 epiPairs=pairs,
                 epiEff=aa,
                 intercept=2,
                 name="imported")

  # The same trait built by hand, with loci in map order a-f and pair
  # members given as indices into that order
  mapOrder = letters[1:6]
  manual = new("TraitAE",
               nLoci=6L,
               lociPerChr=c(3L,3L),
               lociLoc=c(1L,2L,3L,1L,2L,3L),
               addEff=unname(addEff[mapOrder]),
               epiEff=rbind(c(5,1,aa[1]),
                            c(3,6,aa[2]),
                            c(2,4,aa[3])),
               intercept=2,
               name="manual")
  SP$manAddTrait(manual)

  imported = SP$traits[[1]]
  expect_true(is(imported, "TraitAE"))
  expect_false(is(imported, "TraitADE"))
  expect_identical(imported@lociPerChr, manual@lociPerChr)
  expect_identical(imported@lociLoc, manual@lociLoc)
  expect_equal(imported@addEff, manual@addEff)
  expect_equal(imported@epiEff, manual@epiEff)

  # Each pair must still name the markers it was given with
  expect_identical(cbind(mapOrder[imported@epiEff[,1]],
                         mapOrder[imported@epiEff[,2]]),
                   pairs)

  pop = newPop(founderPop, simParam=SP)
  expect_equal(unname(pop@gv[,1]), unname(pop@gv[,2]))
})

test_that("importTrait MULTIPLE traits with dominance give TraitADE",{
  founderPop = makeEpiFounder()
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L

  markers = c("d","b","f","a","e","c")
  addEff = cbind(c(1,2,3,4,5,6), c(-1,0.5,0,2,-2,1))
  domEff = cbind(c(0.1,0.2,0.3,0.4,0.5,0.6), c(0,-0.3,0.2,0,0.1,0.4))
  pairs = data.frame(first=c("a","f","b"),
                     second=c("d","c","e"))
  aa = cbind(c(0.5,-1,2), c(1,0,-0.25))
  SP$importTrait(markerNames=markers,
                 addEff=addEff,
                 domEff=domEff,
                 epiPairs=pairs,
                 epiEff=aa,
                 name=c("T1","T2"))
  expect_identical(SP$nTraits, 2L)

  # Map order a-f corresponds to these rows of the inputs
  ord = match(letters[1:6], markers)
  pairIndex = cbind(c(1,6,2), c(4,3,5))
  for(i in 1:2){
    manual = new("TraitADE",
                 nLoci=6L,
                 lociPerChr=c(3L,3L),
                 lociLoc=c(1L,2L,3L,1L,2L,3L),
                 addEff=addEff[ord,i],
                 domEff=domEff[ord,i],
                 epiEff=cbind(pairIndex, aa[,i]),
                 intercept=0,
                 name=paste0("manual",i))
    SP$manAddTrait(manual)

    imported = SP$traits[[i]]
    expect_true(is(imported, "TraitADE"))
    expect_identical(imported@name, paste0("T",i))
    expect_equal(imported@addEff, manual@addEff)
    expect_equal(imported@domEff, manual@domEff)
    expect_equal(imported@epiEff, manual@epiEff)
  }

  pop = newPop(founderPop, simParam=SP)
  expect_equal(unname(pop@gv[,1]), unname(pop@gv[,3]))
  expect_equal(unname(pop@gv[,2]), unname(pop@gv[,4]))
})

test_that("importTrait FACTOR inputs match character inputs",{
  founderPop = makeEpiFounder()
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L

  markers = c("c","a","d","b")
  pairs = rbind(c("a","d"), c("c","b"))
  SP$importTrait(markerNames=markers, addEff=c(1,2,3,4),
                 epiPairs=pairs, epiEff=c(0.5,-0.5))
  SP$importTrait(markerNames=factor(markers), addEff=c(1,2,3,4),
                 epiPairs=data.frame(x=factor(pairs[,1]),
                                     y=factor(pairs[,2])),
                 epiEff=c(0.5,-0.5))
  expect_equal(SP$traits[[2]]@epiEff, SP$traits[[1]]@epiEff)
  expect_equal(SP$traits[[2]]@addEff, SP$traits[[1]]@addEff)
})

test_that("importTrait WITHOUT epiPairs keeps the A and AD trait classes",{
  founderPop = makeEpiFounder()
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L

  SP$importTrait(markerNames=c("b","a"), addEff=c(1,2))
  SP$importTrait(markerNames=c("b","a"), addEff=c(1,2), domEff=c(0.5,0))
  expect_identical(class(SP$traits[[1]])[1], "TraitA")
  expect_identical(class(SP$traits[[2]])[1], "TraitAD")
})

test_that("importTrait REJECTS invalid epistasis inputs",{
  founderPop = makeEpiFounder()
  SP = SimParam$new(founderPop=founderPop)
  SP$nThreads = 1L

  markers = c("a","b","c","d")
  addEff = c(1,2,3,4)
  pairs = rbind(c("a","b"), c("c","d"))

  expect_error(SP$importTrait(markers, addEff, epiPairs=pairs),
               "supplied together")
  expect_error(SP$importTrait(markers, addEff, epiEff=c(1,1)),
               "supplied together")
  expect_error(SP$importTrait(markers, addEff,
                              epiPairs=pairs[,1,drop=FALSE], epiEff=c(1,1)),
               "two-column")
  expect_error(SP$importTrait(markers, addEff,
                              epiPairs=rbind(c(1,2), c(3,4)), epiEff=c(1,1)),
               "two-column")
  expect_error(SP$importTrait(markers, addEff,
                              epiPairs=rbind(c("a","b"), c("c","z")),
                              epiEff=c(1,1)),
               "must be in markerNames")
  expect_error(SP$importTrait(markers, addEff,
                              epiPairs=rbind(c("a","b"), c("c",NA)),
                              epiEff=c(1,1)),
               "must be in markerNames")
  # Overlapping pairs
  expect_error(SP$importTrait(markers, addEff,
                              epiPairs=rbind(c("a","b"), c("b","c")),
                              epiEff=c(1,1)),
               "exactly one pair")
  # A locus paired with itself
  expect_error(SP$importTrait(markers, addEff,
                              epiPairs=rbind(c("a","a"), c("c","d")),
                              epiEff=c(1,1)),
               "exactly one pair")
  # A locus left without a partner
  expect_error(SP$importTrait(markers, addEff,
                              epiPairs=rbind(c("a","b")),
                              epiEff=1),
               "exactly one pair")
  expect_error(SP$importTrait(c("a","b","c"), c(1,2,3),
                              epiPairs=rbind(c("a","b")),
                              epiEff=1),
               "exactly one pair")
  expect_error(SP$importTrait(markers, addEff,
                              epiPairs=pairs, epiEff=c(1,1,1)),
               "one row for each pair")
  expect_error(SP$importTrait(markers, cbind(addEff, addEff),
                              epiPairs=pairs, epiEff=c(1,1)),
               "one column for each trait")
  expect_error(SP$importTrait(c("a","b","a","d"), addEff),
               "must be unique")

  # Nothing was added by the failed calls
  expect_identical(SP$nTraits, 0L)
})

