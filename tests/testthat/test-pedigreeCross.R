context("pedigreeCross")

# pedigreeCross takes a pedigree and builds the individuals it describes.
#
# Three rules shape everything below. id, mother and father are character
# vectors. An unknown parent is NA or one of the unknownParent codes, and a
# parent given any other value names an individual. A parent named but
# without a row of its own is added to the front of the pedigree as a
# founder before anything else happens, when it is used twice or more or
# when matchID is TRUE; otherwise it is treated as unknown.
#
# The structural check used throughout is Mendelian consistency, which holds
# exactly and needs no reference implementation: a diploid child takes one
# allele from each parent, so its dosage cannot be below the number of
# parents fixed for the 1 allele, nor above the number carrying one at all.

pedSP = function(nInd=4, nChr=2, segSites=20, seed=11001, trackRec=FALSE){
  set.seed(seed)
  founderPop = quickHaplo(nInd=nInd, nChr=nChr, segSites=segSites)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  if(trackRec){
    SP$setTrackRec(TRUE)
  }
  SP$addTraitA(nQtlPerChr=5)
  SP$setVarE(h2=0.5)
  return(list(map=founderPop, SP=SP))
}

# A NamedMapPop, which is the only map class carrying ids. Nothing exported
# builds one from a MapPop, so it is assembled from the slots here.
asNamedMapPop = function(mapPop, id){
  return(new("NamedMapPop",
             id = as.character(id),
             mother = rep(NA_character_, mapPop@nInd),
             father = rep(NA_character_, mapPop@nInd),
             nInd = mapPop@nInd,
             nChr = mapPop@nChr,
             ploidy = mapPop@ploidy,
             nLoci = mapPop@nLoci,
             geno = mapPop@geno,
             genMap = mapPop@genMap,
             centromere = mapPop@centromere,
             inbred = mapPop@inbred))
}

# child dosage must lie within what the two parents could have passed on
expect_mendelian = function(child, mother, father, label=""){
  low = (mother==2) + (father==2)
  high = (mother>0) + (father>0)
  expect_true(all(child>=low),
              info=paste(label, "child carries fewer 1 alleles than its parents forced"))
  expect_true(all(child<=high),
              info=paste(label, "child carries a 1 allele neither parent had"))
}

# The bound above only holds for a direct cross. Once a line is selfed or
# doubled, a heterozygous locus can go to either homozygote, so the only
# thing that still holds is that a locus where both ancestors were fixed for
# the same allele stays fixed.
expect_descendant = function(child, p1, p2, label=""){
  bothZero = (p1==0) & (p2==0)
  bothTwo = (p1==2) & (p2==2)
  expect_true(all(child[bothZero]==0),
              info=paste(label, "descendant gained a 1 allele no ancestor had"))
  expect_true(all(child[bothTwo]==2),
              info=paste(label, "descendant lost a 1 allele both ancestors were fixed for"))
}

# The pedigree used by the help page: a biparental cross then a chain of selfs
biparentalPed = function(){
  return(list(id = as.character(1:10),
              mother = c(NA, NA, "1", as.character(3:9)),
              father = c(NA, NA, "2", as.character(3:9))))
}

# ---------------------------------------------------------------------------
# PART 1  Building a pedigree, matchID = FALSE
# ---------------------------------------------------------------------------

test_that("a fully specified pedigree builds the population it describes", {
  d = pedSP()
  pop = newPop(d$map, simParam=d$SP)
  p = biparentalPed()

  set.seed(101)
  out = pedigreeCross(pop, p$id, p$mother, p$father, simParam=d$SP)

  expect_true(isPop(out))
  expect_equal(nInd(out), 10L)
  expect_equal(out@id, p$id)
  expect_equal(out@mother, p$mother)
  expect_equal(out@father, p$father)
  expect_true(isTRUE(validObject(out, test=TRUE)))

  geno = pullSegSiteGeno(out, simParam=d$SP)
  for(i in 3:10){
    m = match(p$mother[i], p$id)
    f = match(p$father[i], p$id)
    expect_mendelian(geno[i,], geno[m,], geno[f,], label=paste("individual", i))
  }
})

test_that("pedigreeCross is reproducible from a seed", {
  p = biparentalPed()

  a = pedSP(seed=11011)
  popA = newPop(a$map, simParam=a$SP)
  set.seed(202)
  outA = pedigreeCross(popA, p$id, p$mother, p$father, simParam=a$SP)

  b = pedSP(seed=11011)
  popB = newPop(b$map, simParam=b$SP)
  set.seed(202)
  outB = pedigreeCross(popB, p$id, p$mother, p$father, simParam=b$SP)

  expect_equal(unname(pullSegSiteHaplo(outA, simParam=a$SP)),
               unname(pullSegSiteHaplo(outB, simParam=b$SP)))
})

test_that("a pedigree needing more founders than supplied is refused", {
  d = pedSP(nInd=2)
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=as.character(1:4),
                             mother=c(NA,NA,NA,"1"),
                             father=c(NA,NA,NA,"2"),
                             simParam=d$SP),
               "founders")
})

test_that("half founders each get a founder genome of their own", {
  # An individual with one parent known and one unknown needs a genome for
  # the unknown parent, and no two half founders may share one
  d = pedSP(nInd=8, nChr=1, segSites=40, seed=11021, trackRec=TRUE)
  pop = newPop(d$map, simParam=d$SP)

  id     = c("1","2","A","B")
  mother = c(NA, NA, NA, NA)
  father = c(NA, NA, "1", "1")

  set.seed(303)
  out = pedigreeCross(pop, id, mother, father, simParam=d$SP)
  expect_equal(nInd(out), 4L)

  ibd = pullIbdHaplo(out, simParam=d$SP)
  founderOf = function(hap) return((hap+1L)%/%2L)
  rows = match(c("A_1","A_2","B_1","B_2"), rownames(ibd))
  expect_false(any(is.na(rows)))
  # One shared father plus two distinct mothers
  expect_equal(length(unique(founderOf(c(ibd[rows,])))), 3L)
})

test_that("DH and nSelf are applied to the individuals they name", {
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)

  id = as.character(1:4)
  mother = c(NA, NA, "1", "1")
  father = c(NA, NA, "2", "2")

  set.seed(404)
  out = pedigreeCross(pop, id, mother, father,
                      DH=c(FALSE,FALSE,TRUE,FALSE),
                      nSelf=c(0,0,0,2),
                      simParam=d$SP)

  expect_equal(nInd(out), 4L)
  geno = pullSegSiteGeno(out, simParam=d$SP)
  expect_true(all(geno[3,] %in% c(0,2)))
  expect_descendant(geno[3,], geno[1,], geno[2,], label="doubled haploid")
  expect_descendant(geno[4,], geno[1,], geno[2,], label="selfed individual")
})

test_that("the input vectors are checked against each other", {
  # Enumerated against checkPedigreeInput in the CHECKPED block below
  d = pedSP()
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("1","1"), mother=c(NA,NA),
                             father=c(NA,NA), simParam=d$SP),
               "duplicates")
  expect_error(pedigreeCross(pop, id=c("1","2"), mother=NA_character_,
                             father=c(NA,NA), simParam=d$SP),
               "length\\(mother\\)")
})

test_that("pedigreeCross refuses a simulation that uses sexes", {
  set.seed(11031)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setSexes("yes_sys")
  SP$addTraitA(nQtlPerChr=5)
  pop = newPop(founderPop, simParam=SP)

  expect_error(pedigreeCross(pop, id=c("1","2","3"), mother=c(NA,NA,"1"),
                             father=c(NA,NA,"2"), simParam=SP),
               "sex")
})

test_that("a pedigree given out of order still sorts", {
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)

  id     = c("5","4","3","2","1")
  mother = c("4","3","2","1",NA)
  father = c("4","3","2","1",NA)

  out = pedigreeCross(pop, id, mother, father, simParam=d$SP)
  expect_equal(nInd(out), 5L)
  expect_equal(out@id, id)

  geno = pullSegSiteGeno(out, simParam=d$SP)
  for(i in 1:4){
    m = match(mother[i], id)
    expect_mendelian(geno[i,], geno[m,], geno[m,], label=paste("row", i))
  }
})

test_that("an unsortable pedigree is reported", {
  d = pedSP()
  pop = newPop(d$map, simParam=d$SP)

  # A and B are each other's parents
  expect_error(pedigreeCross(pop, id=c("A","B"), mother=c("B","A"),
                             father=c("B","A"), simParam=d$SP),
               "cycle")
  expect_error(pedigreeCross(pop, id=c("A","B"), mother=c("B","A"),
                             father=c("B","A"), simParam=d$SP),
               "A")

  # Running out of passes is a different failure
  expect_error(pedigreeCross(pop, id=c("5","4","3","2","1"),
                             mother=c("4","3","2","1",NA),
                             father=c("4","3","2","1",NA),
                             maxCycle=2, simParam=d$SP),
               "maxCycle")
})

# ---------------------------------------------------------------------------
# PART 2  Extending the pedigree backwards
# ---------------------------------------------------------------------------

test_that("EXTEND a parent used once is dropped and its child is a founder", {
  # "3" names two parents that appear nowhere else. Neither carries a
  # relationship to anything, so both are treated as unknown and "3" is a
  # founder rather than the cross of two added rows.
  d = pedSP(nInd=4, seed=12001)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id="3", mother="2", father="1", simParam=d$SP)

  expect_equal(nInd(out), 1L)
  expect_equal(out@id, "3")
  expect_equal(out@mother, NA_character_)
  expect_equal(out@father, NA_character_)

  # The genome came whole from founderPop rather than being made
  geno = pullSegSiteGeno(out, simParam=d$SP)
  founders = pullSegSiteGeno(pop, simParam=d$SP)
  isCopy = apply(founders, 1L, function(x) all(x==geno[1,]))
  expect_true(any(isCopy))
})

test_that("EXTEND a parent used more than once is added", {
  # "mum" and "dad" are each the parent of two individuals, so dropping
  # them would leave two full sibs unrelated
  d = pedSP(nInd=8, seed=12011)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id=c("a","b"),
                      mother=c("mum","mum"),
                      father=c("dad","dad"),
                      simParam=d$SP)

  # The added rows are sorted, and come before the supplied pedigree
  expect_equal(out@id, c("dad","mum","a","b"))
  expect_equal(nInd(out), 4L)

  geno = pullSegSiteGeno(out, simParam=d$SP)
  expect_mendelian(geno[3,], geno[2,], geno[1,], label="a")
  expect_mendelian(geno[4,], geno[2,], geno[1,], label="b")
})

test_that("EXTEND a half sib group keeps the parent they share", {
  # The shared mother is added, the two fathers are used once each and are
  # dropped, which leaves the two half sibs as half founders
  d = pedSP(nInd=8, seed=12013)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id=c("a","b"),
                      mother=c("mum","mum"),
                      father=c("dad1","dad2"),
                      simParam=d$SP)

  expect_equal(out@id, c("mum","a","b"))
  expect_equal(out@mother, c(NA,"mum","mum"))
  expect_equal(out@father, c(NA_character_, NA_character_, NA_character_))
})

test_that("EXTEND a name used as both parents of one individual is added", {
  # Two uses, but only one individual. Dropping it would turn a self into
  # an outcross of two unrelated founders.
  d = pedSP(nInd=4, seed=12016)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id="kid", mother="p", father="p", simParam=d$SP)

  expect_equal(out@id, c("p","kid"))
  expect_equal(nInd(out), 2L)

  geno = pullSegSiteGeno(out, simParam=d$SP)
  expect_mendelian(geno[2,], geno[1,], geno[1,], label="self")
})

test_that("EXTEND a pedigree that needs no extension is left alone", {
  d = pedSP(nInd=4, seed=12021)
  pop = newPop(d$map, simParam=d$SP)
  p = biparentalPed()

  out = pedigreeCross(pop, p$id, p$mother, p$father, simParam=d$SP)
  expect_equal(out@id, p$id)
  expect_equal(nInd(out), 10L)
})

test_that("EXTEND zero is a name unless unknownParent says otherwise", {
  d = pedSP(nInd=4, seed=12031)
  pop = newPop(d$map, simParam=d$SP)

  # Used twice, so "0" names an individual and is added
  out = pedigreeCross(pop, id=c("3","4"), mother=c("0","0"),
                      father=c(NA,NA), simParam=d$SP)
  expect_true("0" %in% out@id)
  expect_equal(out@mother[out@id=="3"], "0")

  # Saying that "0" means unknown makes both rows founders instead
  out = pedigreeCross(pop, id=c("3","4"), mother=c("0","0"),
                      father=c(NA,NA), unknownParent="0", simParam=d$SP)
  expect_equal(out@id, c("3","4"))
  expect_true(all(is.na(out@mother)))
  expect_true(all(is.na(out@father)))
})

test_that("EXTEND several unknown codes can be given at once", {
  d = pedSP(nInd=4, seed=12036)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id=c("a","b"), mother=c("0","a"),
                      father=c("", "a"),
                      unknownParent=c("0",""), simParam=d$SP)

  expect_equal(out@id, c("a","b"))
  expect_equal(out@mother, c(NA,"a"))
  expect_equal(out@father, c(NA,"a"))
})

test_that("EXTEND added founders count towards the founder budget", {
  # "1" and "2" are each used twice, so both are added, and the two added
  # rows need two founder genomes where only one is on offer
  d = pedSP(nInd=1, seed=12041)
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("3","4"), mother=c("2","2"),
                             father=c("1","1"), simParam=d$SP),
               "founders")
})

# ---------------------------------------------------------------------------
# PART 3  matchID = TRUE
# ---------------------------------------------------------------------------

test_that("MATCHID a matched individual is copied from founderPop", {
  d = pedSP(nInd=4, seed=13001)
  pop = newPop(d$map, simParam=d$SP)
  expect_equal(pop@id, as.character(1:4))
  p = biparentalPed()

  out = pedigreeCross(pop, p$id, p$mother, p$father, matchID=TRUE,
                      simParam=d$SP)

  # "1" to "4" are matched, and 5 to 10 descend from them
  expect_equal(nInd(out), 10L)
  expect_equal(out@id, p$id)

  # The matched individuals are copies, so their genotypes are identical
  fromPed = pullSegSiteGeno(out, simParam=d$SP)[1:2,,drop=FALSE]
  fromPop = pullSegSiteGeno(pop[c("1","2")], simParam=d$SP)
  expect_equal(unname(fromPed), unname(fromPop))
})

test_that("MATCHID the ancestors of a match are skipped", {
  # "2" is in founderPop, so its own parents are not simulated and do not
  # appear in the result
  d = pedSP(nInd=4, seed=13011)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id=c("A","B","2","kid"),
                      mother=c(NA,NA,"A","2"),
                      father=c(NA,NA,"B","2"),
                      matchID=TRUE, simParam=d$SP)

  expect_equal(nInd(out), 2L)
  expect_equal(out@id, c("2","kid"))
  # The pedigree it reports is the one it was given
  expect_equal(out@mother, c("A","2"))

  geno = pullSegSiteGeno(out, simParam=d$SP)
  fromPop = pullSegSiteGeno(pop["2"], simParam=d$SP)
  expect_equal(unname(geno[1,]), unname(fromPop[1,]))
  expect_mendelian(geno[2,], geno[1,], geno[1,], label="kid")
})

test_that("MATCHID an unmatched pedigree is an error", {
  d = pedSP(nInd=4, seed=13021)
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("x","y","z"),
                             mother=c(NA,NA,"x"),
                             father=c(NA,NA,"y"),
                             matchID=TRUE, simParam=d$SP),
               "matches")
})

test_that("MATCHID an individual that cannot be reached is an error", {
  # "2" is matched, but X and Y hang off nothing that can be made
  d = pedSP(nInd=4, seed=13031)
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("2","X","Y"),
                             mother=c(NA,NA,"X"),
                             father=c(NA,NA,"X"),
                             matchID=TRUE, simParam=d$SP),
               "X")
})

test_that("MATCHID extension and matching work together", {
  # "1" and "2" are named as parents without rows of their own. They are
  # added by the extension and then matched against founderPop.
  d = pedSP(nInd=4, seed=13041)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id="kid", mother="1", father="2",
                      matchID=TRUE, simParam=d$SP)

  expect_equal(nInd(out), 3L)
  expect_equal(out@id, c("1","2","kid"))

  geno = pullSegSiteGeno(out, simParam=d$SP)
  par = pullSegSiteGeno(pop[c("1","2")], simParam=d$SP)
  expect_equal(unname(geno[1:2,]), unname(par))
  expect_mendelian(geno[3,], geno[1,], geno[2,], label="kid")
})

test_that("MATCHID reuses an id rather than refusing it", {
  # A row whose name is already in founderPop is matched, which is the
  # point of matchID, rather than being rejected as a clash
  d = pedSP(nInd=4, seed=13051)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id=c("1","2","3"),
                      mother=c(NA,NA,"1"),
                      father=c(NA,NA,"2"),
                      matchID=TRUE, simParam=d$SP)

  expect_equal(nInd(out), 3L)
  # "3" is in founderPop, so it is a copy and not the cross the pedigree
  # describes
  geno = pullSegSiteGeno(out, simParam=d$SP)
  fromPop = pullSegSiteGeno(pop["3"], simParam=d$SP)
  expect_equal(unname(geno[3,]), unname(fromPop[1,]))
})

test_that("MATCHID needs no founder budget", {
  # Nothing is sampled at random, so a founderPop with one individual is
  # enough as long as the pedigree reaches it
  d = pedSP(nInd=1, seed=13061)
  pop = newPop(d$map, simParam=d$SP)
  expect_equal(pop@id, "1")

  out = pedigreeCross(pop, id=c("1","kid"),
                      mother=c(NA,"1"),
                      father=c(NA,"1"),
                      matchID=TRUE, simParam=d$SP)
  expect_equal(nInd(out), 2L)
})

# ---------------------------------------------------------------------------
# PART 4  Map populations
# ---------------------------------------------------------------------------

test_that("MAP a MapPop means no simulation has been set up yet", {
  d = pedSP(nInd=4, seed=14001)
  p = biparentalPed()
  expect_false(exists("SP", envir=globalenv(), inherits=FALSE))

  set.seed(505)
  out = pedigreeCross(d$map, p$id, p$mother, p$father)

  expect_true(isNamedMapPop(out))
  expect_false(isPop(out))
  expect_equal(nInd(out), 10L)
  expect_equal(out@id, p$id)
  expect_equal(out@genMap, d$map@genMap)
  expect_true(isTRUE(validObject(out, test=TRUE)))
  expect_false(exists("SP", envir=globalenv(), inherits=FALSE))

  geno = pullSegSiteGeno(out)
  for(i in 3:10){
    m = match(p$mother[i], p$id)
    f = match(p$father[i], p$id)
    expect_mendelian(geno[i,], geno[m,], geno[f,], label=paste("individual", i))
  }
})

test_that("MAP the returned map population starts a simulation", {
  d = pedSP(nInd=4, seed=14011)
  p = biparentalPed()

  set.seed(515)
  out = pedigreeCross(d$map, p$id, p$mother, p$father)

  SP = SimParam$new(out)
  SP$nThreads = 1L
  SP$addTraitA(nQtlPerChr=5)
  SP$setVarE(h2=0.5)
  pop = newPop(out, simParam=SP)

  expect_true(isPop(pop))
  expect_equal(pop@id, p$id)
  expect_true(all(is.finite(c(gv(pop)))))
})

test_that("MAP a supplied simParam is ignored with a warning", {
  d = pedSP(nInd=4, seed=14021)
  p = biparentalPed()

  expect_warning(pedigreeCross(d$map, p$id, p$mother, p$father,
                               simParam=d$SP),
                 "ignored")
})

test_that("MAP a NamedMapPop can be matched on", {
  d = pedSP(nInd=4, seed=14031)
  named = asNamedMapPop(d$map, id=c("F1","F2","F3","F4"))

  out = pedigreeCross(named, id=c("F1","F2","kid"),
                      mother=c(NA,NA,"F1"),
                      father=c(NA,NA,"F2"),
                      matchID=TRUE)

  expect_true(isNamedMapPop(out))
  expect_equal(nInd(out), 3L)
  expect_equal(out@id, c("F1","F2","kid"))

  geno = pullSegSiteGeno(out)
  expect_mendelian(geno[3,], geno[1,], geno[2,], label="kid")
})

test_that("MAP matchID on a plain MapPop says why it cannot work", {
  d = pedSP(nInd=4, seed=14041)
  p = biparentalPed()

  expect_error(pedigreeCross(d$map, p$id, p$mother, p$father,
                             matchID=TRUE),
               "matchID")
})

# ---------------------------------------------------------------------------
# PART 5  Recombination settings passed through ...
# ---------------------------------------------------------------------------

test_that("RECOMB not passing v, p and quadProb uses the SimParam defaults", {
  d = pedSP(nInd=4, nChr=1, segSites=60, seed=15001)
  ped = biparentalPed()

  set.seed(2101)
  bare = pedigreeCross(d$map, ped$id, ped$mother, ped$father)
  set.seed(2101)
  spelled = pedigreeCross(d$map, ped$id, ped$mother, ped$father,
                          v=2.6, p=0, quadProb=0)

  expect_equal(unname(pullSegSiteHaplo(bare)),
               unname(pullSegSiteHaplo(spelled)))
})

test_that("RECOMB v and p reach the meiosis", {
  d = pedSP(nInd=4, nChr=1, segSites=60, seed=15011)
  ped = biparentalPed()

  set.seed(2102)
  kosambi = pedigreeCross(d$map, ped$id, ped$mother, ped$father)
  set.seed(2102)
  haldane = pedigreeCross(d$map, ped$id, ped$mother, ped$father, v=1)
  expect_false(isTRUE(all.equal(unname(pullSegSiteHaplo(kosambi)),
                                unname(pullSegSiteHaplo(haldane)))))

  set.seed(2103)
  noPath = pedigreeCross(d$map, ped$id, ped$mother, ped$father)
  set.seed(2103)
  withPath = pedigreeCross(d$map, ped$id, ped$mother, ped$father, p=0.5)
  expect_false(isTRUE(all.equal(unname(pullSegSiteHaplo(noPath)),
                                unname(pullSegSiteHaplo(withPath)))))
})

test_that("RECOMB quadProb reaches the meiosis of an autopolyploid", {
  skip_on_cran()
  set.seed(15021)
  mapPop = quickHaplo(nInd=6, nChr=1, segSites=60, ploidy=4L)

  id     = as.character(1:4)
  mother = c(NA, NA, "1", "3")
  father = c(NA, NA, "2", "2")

  set.seed(2104)
  bivalent = pedigreeCross(mapPop, id, mother, father, quadProb=0)
  set.seed(2104)
  quadrivalent = pedigreeCross(mapPop, id, mother, father, quadProb=1)

  expect_equal(bivalent@ploidy, 4L)
  expect_equal(nInd(quadrivalent), 4L)
  expect_false(isTRUE(all.equal(unname(pullSegSiteHaplo(bivalent)),
                                unname(pullSegSiteHaplo(quadrivalent)))))
})

test_that("RECOMB the arguments are checked and a Pop refuses them", {
  # Enumerated against checkRecombArgs in the CHECKARGS block below
  d = pedSP(nInd=4, seed=15031)
  ped = biparentalPed()
  cross = function(...) pedigreeCross(d$map, ped$id, ped$mother, ped$father, ...)

  # maxCycles is not a prefix of maxCycle, so R does not partial match it
  expect_error(cross(maxCycles=2), "Unused arguments")
  expect_error(cross(v=0), "greater than zero")
  expect_error(cross(p=1.5), "between zero and one")

  # The edges of the ranges run the whole cross, not just the check
  expect_equal(nInd(cross(p=0)), 10L)
  expect_equal(nInd(cross(p=1)), 10L)

  pop = newPop(d$map, simParam=d$SP)
  expect_error(pedigreeCross(pop, ped$id, ped$mother, ped$father,
                             simParam=d$SP, v=1),
               "map population")
  expect_equal(d$SP$v, 2.6)
})

# ---------------------------------------------------------------------------
# PART 6  checkRecombArgs on its own
#
# pedigreeCross cannot test an unnamed argument at all: R matches a bare
# value positionally to matchID long before it reaches the dots.
# ---------------------------------------------------------------------------

checkArgs = function(...) AlphaSimR:::checkRecombArgs(list(...))

test_that("CHECKARGS nothing to check is nothing to complain about", {
  expect_equal(AlphaSimR:::checkRecombArgs(list()), list())
  expect_equal(checkArgs(v=1.5), list(v=1.5))
  expect_equal(checkArgs(v=2.6, p=0.25, quadProb=0.5),
               list(v=2.6, p=0.25, quadProb=0.5))
})

test_that("CHECKARGS only the three recombination names are allowed", {
  expect_error(checkArgs(maxCycles=2), "Unused arguments")
  expect_error(checkArgs(V=1), "Unused arguments")
  expect_error(checkArgs(v=1, nonsense=2), "nonsense")
  expect_error(AlphaSimR:::checkRecombArgs(list(2.6)), "Unused arguments")
  expect_error(AlphaSimR:::checkRecombArgs(list(v=1, 2)), "Unused arguments")
})

test_that("CHECKARGS each value is a single finite number", {
  expect_error(checkArgs(v="a"), "single number")
  expect_error(checkArgs(v=TRUE), "single number")
  expect_error(checkArgs(v=c(1,2)), "single number")
  expect_error(checkArgs(v=numeric(0)), "single number")
  expect_error(checkArgs(p=NA_real_), "single number")
  expect_error(checkArgs(p=Inf), "single number")
  expect_error(checkArgs(quadProb="a"), "quadProb")
})

test_that("CHECKARGS the ranges are enforced at their edges", {
  expect_error(checkArgs(v=0), "greater than zero")
  expect_error(checkArgs(p=-0.1), "between zero and one")
  expect_error(checkArgs(quadProb=2), "between zero and one")

  expect_equal(checkArgs(v=1e-8)$v, 1e-8)
  expect_equal(checkArgs(p=0)$p, 0)
  expect_equal(checkArgs(p=1)$p, 1)
  expect_equal(checkArgs(quadProb=1)$quadProb, 1)
})

# ---------------------------------------------------------------------------
# PART 7  checkPedigreeInput on its own
# ---------------------------------------------------------------------------

checkPed = function(id, mother, father, DH=NULL, nSelf=NULL,
                    unknownParent=NA_character_, matchID=FALSE){
  return(AlphaSimR:::checkPedigreeInput(id=id, mother=mother, father=father,
                                        DH=DH, nSelf=nSelf,
                                        unknownParent=unknownParent,
                                        matchID=matchID))
}

test_that("CHECKPED the pedigree comes back coerced to character", {
  r = checkPed(id=1:3, mother=c(NA,NA,1), father=c(NA,NA,2))

  expect_equal(r$id, c("1","2","3"))
  expect_equal(r$mother, c(NA,NA,"1"))
  expect_equal(r$father, c(NA,NA,"2"))
  expect_true(is.character(r$id))
  expect_equal(r$DH, c(FALSE,FALSE,FALSE))
  expect_equal(r$nSelf, c(0L,0L,0L))
  expect_true(is.integer(r$nSelf))
  expect_equal(r$nAdded, 0L)
})

test_that("CHECKPED a missing name used once is treated as unknown", {
  r = checkPed(id="3", mother="2", father="1")

  expect_equal(r$id, "3")
  expect_equal(r$mother, NA_character_)
  expect_equal(r$father, NA_character_)
  expect_equal(r$nAdded, 0L)
  expect_equal(r$DH, FALSE)
  expect_equal(r$nSelf, 0L)
})

test_that("CHECKPED a missing name used twice is added", {
  r = checkPed(id=c("3","4"), mother=c("2","2"), father=c("1","1"))

  expect_equal(r$id, c("1","2","3","4"))
  expect_equal(r$mother, c(NA_character_, NA_character_, "2", "2"))
  expect_equal(r$father, c(NA_character_, NA_character_, "1", "1"))
  expect_equal(r$nAdded, 2L)
  # The added founders are neither selfed nor doubled
  expect_equal(r$DH, rep(FALSE, 4))
  expect_equal(r$nSelf, rep(0L, 4))
})

test_that("CHECKPED uses are counted across mother and father together", {
  # One individual, but two uses, so the name is kept and the individual
  # comes back as a self rather than as a founder
  r = checkPed(id="kid", mother="p", father="p")

  expect_equal(r$id, c("p","kid"))
  expect_equal(r$mother, c(NA_character_, "p"))
  expect_equal(r$father, c(NA_character_, "p"))
  expect_equal(r$nAdded, 1L)
})

test_that("CHECKPED a name that already has a row is never counted", {
  # "3" is used twice, but it is a row of the pedigree, so nothing is added
  # and the singleton "2" is still dropped
  r = checkPed(id=c("3","4"), mother=c("2","3"), father=c(NA,"3"))

  expect_equal(r$id, c("3","4"))
  expect_equal(r$mother, c(NA_character_, "3"))
  expect_equal(r$father, c(NA_character_, "3"))
  expect_equal(r$nAdded, 0L)
})

test_that("CHECKPED matchID adds every missing name", {
  # A name has to have a row before it can be matched, so the extension is
  # not limited when matchID is TRUE
  r = checkPed(id="3", mother="2", father="1", matchID=TRUE)

  expect_equal(r$id, c("1","2","3"))
  expect_equal(r$mother, c(NA_character_, NA_character_, "2"))
  expect_equal(r$father, c(NA_character_, NA_character_, "1"))
  expect_equal(r$nAdded, 2L)
})

test_that("CHECKPED extension keeps the supplied per individual vectors", {
  r = checkPed(id=c("3","4"), mother=c("2","2"), father=c("1","1"),
               DH=c(TRUE,FALSE), nSelf=c(0,2))

  expect_equal(r$id, c("1","2","3","4"))
  expect_equal(r$nAdded, 2L)
  # The added rows take defaults, the supplied rows keep their values
  expect_equal(r$DH, c(FALSE,FALSE,TRUE,FALSE))
  expect_equal(r$nSelf, c(0L,0L,0L,2L))
})

test_that("CHECKPED only NA is an unknown parent by default", {
  # "0" names an individual unless it is given as unknownParent
  r = checkPed(id=c("A","B"), mother=c("0","0"), father=c(NA,NA))

  expect_equal(r$id, c("0","A","B"))
  expect_equal(r$mother, c(NA_character_, "0", "0"))
  expect_equal(r$nAdded, 1L)
})

test_that("CHECKPED unknownParent recodes the parent vectors", {
  r = checkPed(id=c("A","B"), mother=c("0","0"), father=c(NA,NA),
               unknownParent="0")

  expect_equal(r$id, c("A","B"))
  expect_equal(r$mother, c(NA_character_, NA_character_))
  expect_equal(r$father, c(NA_character_, NA_character_))
  expect_equal(r$nAdded, 0L)
})

test_that("CHECKPED unknownParent takes more than one code", {
  r = checkPed(id=c("A","B","C"), mother=c("0","","-9"),
               father=c("0","","A"),
               unknownParent=c("0","","-9"))

  expect_equal(r$id, c("A","B","C"))
  expect_equal(r$mother, c(NA_character_, NA_character_, NA_character_))
  expect_equal(r$father, c(NA_character_, NA_character_, "A"))
})

test_that("CHECKPED unknownParent is coerced like the pedigree", {
  r = checkPed(id=c("A","B"), mother=c("0","0"), father=c(NA,NA),
               unknownParent=0)

  expect_equal(r$id, c("A","B"))
  expect_true(all(is.na(r$mother)))
})

test_that("CHECKPED NA stays unknown whatever unknownParent is", {
  r = checkPed(id=c("A","B"), mother=c("0",NA), father=c(NA,"A"),
               unknownParent="0")

  expect_equal(r$id, c("A","B"))
  expect_equal(r$mother, c(NA_character_, NA_character_))
  expect_equal(r$father, c(NA_character_, "A"))
})

test_that("CHECKPED an id that is an unknown code is refused", {
  expect_error(checkPed(id=c("0","A"), mother=c(NA,"0"), father=c(NA,NA),
                        unknownParent="0"),
               "unknownParent")
})

test_that("CHECKPED an id of NA is refused", {
  expect_error(checkPed(id=c("1",NA), mother=c(NA,NA), father=c(NA,NA)),
               "id can not contain NA")
})

test_that("CHECKPED the vectors are checked against each other", {
  expect_error(checkPed(id=c("1","1"), mother=c(NA,NA), father=c(NA,NA)),
               "duplicates")
  expect_error(checkPed(id=c("1","2"), mother=NA_character_,
                        father=c(NA,NA)),
               "length\\(mother\\)")
  expect_error(checkPed(id=c("1","2"), mother=c(NA,NA), father=NA_character_),
               "length\\(father\\)")
  expect_error(checkPed(id=c("1","2"), mother=c(NA,NA), father=c(NA,NA),
                        DH=TRUE),
               "length\\(DH\\)")
  expect_error(checkPed(id=c("1","2"), mother=c(NA,NA), father=c(NA,NA),
                        nSelf=0),
               "length\\(nSelf\\)")
  expect_error(checkPed(id=character(0), mother=character(0),
                        father=character(0)),
               "empty")
})

test_that("CHECKPED nSelf and DH are checked", {
  with3 = function(...) checkPed(id=c("1","2","3"),
                                 mother=c(NA,NA,"1"),
                                 father=c(NA,NA,"2"), ...)

  expect_error(with3(nSelf=c(0,0,-1)), "nSelf")
  expect_error(with3(nSelf=c(0,0,NA)), "nSelf")
  expect_error(with3(nSelf=c(0,0,"a")), "nSelf")
  expect_error(with3(DH=c(FALSE,FALSE,NA)), "DH")

  expect_equal(with3(nSelf=c(0,1,2))$nSelf, c(0L,1L,2L))
  expect_equal(with3(DH=c(TRUE,FALSE,TRUE))$DH, c(TRUE,FALSE,TRUE))
})

# ---------------------------------------------------------------------------
# PART 8  resolveFounders on its own
#
# Founder allocation is where every bug in this function has been. The
# blocks above reach it through a whole simulation; these reach it directly,
# with no population, no SimParam and no genotypes. The pedigree it is given
# has already been extended, so a parent is either a row or NA.
# ---------------------------------------------------------------------------

resolve = function(id, mother, father, founderIds=character(0),
                   nFounderInd=0L, matchID=FALSE){
  return(AlphaSimR:::resolveFounders(id=id, mother=mother, father=father,
                                     founderIds=founderIds,
                                     nFounderInd=nFounderInd,
                                     matchID=matchID))
}

test_that("RESOLVE a biparental pedigree names two founders", {
  r = resolve(id=c("1","2","3"),
              mother=c(NA,NA,"1"),
              father=c(NA,NA,"2"),
              nFounderInd=5L)

  expect_true(all(r$build))
  expect_equal(r$motherPed, c(NA_integer_,NA_integer_,1L))
  expect_equal(r$fatherPed, c(NA_integer_,NA_integer_,2L))
  expect_equal(sum(!is.na(r$founderRowFP)), 2L)
  expect_true(all(is.na(r$motherFP)))
})

test_that("RESOLVE every half founder gets a genome of its own", {
  set.seed(18001)
  r = resolve(id=c("1","2","A","B"),
              mother=c(NA,NA,NA,NA),
              father=c(NA,NA,"1","1"),
              nFounderInd=8L)

  expect_false(any(is.na(r$motherFP[3:4])))
  expect_false(r$motherFP[3]==r$motherFP[4])
  taken = c(r$founderRowFP, r$motherFP, r$fatherFP)
  taken = taken[!is.na(taken)]
  expect_equal(length(taken), 4L)
  expect_equal(anyDuplicated(taken), 0L)
})

test_that("RESOLVE half founders count towards the founder budget", {
  expect_error(resolve(id=c("1","2","A","B"),
                       mother=c(NA,NA,NA,NA),
                       father=c(NA,NA,"1","1"),
                       nFounderInd=3L),
               "founders")
})

test_that("RESOLVE matchID copies a match and skips its ancestors", {
  r = resolve(id=c("A","B","2","kid"),
              mother=c(NA,NA,"A","2"),
              father=c(NA,NA,"B","2"),
              founderIds=as.character(1:4), matchID=TRUE)

  expect_equal(r$build, c(FALSE,FALSE,TRUE,TRUE))
  # "2" is the second individual of founderPop
  expect_equal(r$founderRowFP[3], 2L)
  # kid is crossed from its parents, not copied
  expect_true(is.na(r$founderRowFP[4]))
  expect_equal(r$motherPed[4], 3L)
  # matchID never draws an anonymous founder
  expect_true(all(is.na(r$motherFP)))
  expect_true(all(is.na(r$fatherFP)))
})

test_that("RESOLVE matchID needs at least one match", {
  expect_error(resolve(id=c("x","y"), mother=c(NA,NA), father=c(NA,NA),
                       founderIds=c("1","2"), matchID=TRUE),
               "matches")
})

test_that("RESOLVE matchID refuses an individual it cannot reach", {
  expect_error(resolve(id=c("2","X","Y"),
                       mother=c(NA,NA,"X"),
                       father=c(NA,NA,"X"),
                       founderIds=as.character(1:4), matchID=TRUE),
               "X")
})

# ---------------------------------------------------------------------------
# PART 9  Generations, and building one generation at a time
#
# A generation number is one more than the larger of its parents', so it is
# the longest path back to a founder. It decides both which individuals are
# made together and the order they are made in, so a change to it would move
# every seeded result.
# ---------------------------------------------------------------------------

sortGen = function(id, mother, father, maxCycle=NULL){
  return(AlphaSimR:::sortPed(id=id, mother=mother, father=father,
                             maxCycle=maxCycle))
}

test_that("SORTPED an individual with no known parent is generation one", {
  g = sortGen(id=c("a","b"), mother=c(NA,NA), father=c(NA,NA))
  expect_equal(g, c(1L,1L))
  expect_true(is.integer(g))
})

test_that("SORTPED a generation is one past the later of its parents", {
  # "4" has a founder for a mother and a second generation father, so it is
  # third rather than second
  g = sortGen(id=c("1","2","3","4"),
              mother=c(NA,NA,"1","1"),
              father=c(NA,NA,"2","3"))
  expect_equal(g, c(1L,1L,2L,3L))
})

test_that("SORTPED an unknown parent counts as nothing", {
  # A half founder is one generation past its one known parent
  g = sortGen(id=c("1","2"), mother=c(NA,"1"), father=c(NA,NA))
  expect_equal(g, c(1L,2L))
})

test_that("SORTPED the order the pedigree arrives in does not matter", {
  id = c("1","2","3","4")
  mother = c(NA,NA,"1","3")
  father = c(NA,NA,"2","3")
  forward = sortGen(id, mother, father)

  back = sortGen(rev(id), rev(mother), rev(father))
  expect_equal(rev(back), forward)
})

test_that("SORTPED a deep pedigree sorts", {
  n = 60L
  id = as.character(1:n)
  parent = c(NA, as.character(1:(n-1)))
  expect_equal(sortGen(id, parent, parent), 1:n)
})

test_that("SORTPED the default bound follows the size of the pedigree", {
  # A chain is the deepest a pedigree of a given size can be, so the number
  # of individuals is always enough. This one used to need maxCycle raising
  # by hand, because it is more than 100 deep.
  n = 150L
  id = as.character(1:n)
  parent = c(NA, as.character(1:(n-1)))
  expect_equal(sortGen(id, parent, parent), 1:n)
})

test_that("SORTPED a pedigree deeper than a given maxCycle is reported", {
  id = as.character(1:5)
  parent = c(NA, as.character(1:4))
  expect_error(sortGen(id, parent, parent, maxCycle=3), "maxCycle")
})

test_that("SORTPED a cycle is reported and names the individuals in it", {
  expect_error(sortGen(id=c("A","B"), mother=c("B","A"), father=c("B","A")),
               "cycle")
  expect_error(sortGen(id=c("A","B"), mother=c("B","A"), father=c("B","A")),
               "A, B")
})

test_that("SORTPED a cycle is told apart from running out of passes", {
  # Both leave individuals unassigned, and the two are worth distinguishing
  expect_error(sortGen(id=c("A","B"), mother=c("B","A"), father=c("B","A"),
                       maxCycle=100),
               "cycle")
})

test_that("BATCH a generation mixing copies and crosses keeps its order", {
  # "3" is matched and copied, "kid" is crossed, and both are second
  # generation, so they are made in the same batch
  d = pedSP(nInd=4, seed=14001)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id=c("1","2","3","kid"),
                      mother=c(NA,NA,"1","1"),
                      father=c(NA,NA,"2","2"),
                      matchID=TRUE, simParam=d$SP)

  expect_equal(out@id, c("1","2","3","kid"))

  geno = pullSegSiteGeno(out, simParam=d$SP)
  fromPop = pullSegSiteGeno(pop, simParam=d$SP)
  expect_equal(unname(geno["3",]), unname(fromPop["3",]))
  expect_mendelian(geno["kid",], geno["1",], geno["2",], label="kid")
})

test_that("BATCH a generation selfed different amounts keeps its order", {
  d = pedSP(nInd=6, seed=14011)
  pop = newPop(d$map, simParam=d$SP)

  id = c("p1","p2","a","b","c")
  mother = c(NA,NA,"p1","p1","p1")
  father = c(NA,NA,"p2","p2","p2")

  out = pedigreeCross(pop, id, mother, father,
                      nSelf=c(0,0,0,2,1), simParam=d$SP)

  expect_equal(out@id, id)

  geno = pullSegSiteGeno(out, simParam=d$SP)
  # "a" is a direct cross, so the Mendelian bound holds. "b" and "c" have
  # been selfed, so only the weaker bound does.
  expect_mendelian(geno["a",], geno["p1",], geno["p2",], label="a")
  expect_descendant(geno["b",], geno["p1",], geno["p2",], label="b")
  expect_descendant(geno["c",], geno["p1",], geno["p2",], label="c")
})

test_that("BATCH a generation with only some doubled haploids keeps its order", {
  d = pedSP(nInd=6, seed=14021)
  pop = newPop(d$map, simParam=d$SP)

  id = c("p1","p2","a","b","c")
  mother = c(NA,NA,"p1","p1","p1")
  father = c(NA,NA,"p2","p2","p2")

  out = pedigreeCross(pop, id, mother, father,
                      DH=c(FALSE,FALSE,TRUE,FALSE,TRUE), simParam=d$SP)

  expect_equal(out@id, id)

  geno = pullSegSiteGeno(out, simParam=d$SP)
  # The two doubled haploids are homozygous everywhere and "b", which was
  # left out of the same call, is not
  expect_false(any(geno["a",]==1))
  expect_false(any(geno["c",]==1))
  expect_true(any(geno["b",]==1))
})

test_that("BATCH selfing and doubled haploids combine in one generation", {
  d = pedSP(nInd=6, seed=14031)
  pop = newPop(d$map, simParam=d$SP)

  id = c("p1","p2","a","b")
  mother = c(NA,NA,"p1","p1")
  father = c(NA,NA,"p2","p2")

  out = pedigreeCross(pop, id, mother, father,
                      nSelf=c(0,0,2,0),
                      DH=c(FALSE,FALSE,TRUE,TRUE),
                      simParam=d$SP)

  expect_equal(out@id, id)

  geno = pullSegSiteGeno(out, simParam=d$SP)
  expect_false(any(geno["a",]==1))
  expect_false(any(geno["b",]==1))
  expect_descendant(geno["a",], geno["p1",], geno["p2",], label="a")
})
