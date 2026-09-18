context("pedigreeCross")

# Issue #79. pedigreeCross had no test coverage at all, so this file is in
# two parts.
#
#   PART 1 pins the behavior that exists today and must survive the change.
#   PART 2 specifies the four things issue #79 asks for. Those tests fail
#          until the feature is written, which is the point of them.
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
             mother = rep("0", mapPop@nInd),
             father = rep("0", mapPop@nInd),
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

# The Mendelian bound above only holds for a direct cross. Once a line is
# selfed or doubled, a heterozygous locus can go to either homozygote, so
# the only thing that still holds over any number of generations is that a
# locus where both ancestors were fixed for the same allele stays fixed.
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
              mother = as.character(c(0,0,1,3:9)),
              father = as.character(c(0,0,2,3:9))))
}

# ---------------------------------------------------------------------------
# PART 1  Behavior that exists today and must not change
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

  # Every non founder is consistent with the parents the pedigree names
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

test_that("matchID uses the named founders from the population", {
  # founderPop holds exactly the two founders. A larger one would collide
  # with pedigree rows 3 and 4, which is now an error in its own right.
  d = pedSP(nInd=2)
  pop = newPop(d$map, simParam=d$SP)
  # newPop numbers individuals from one, so the founders are "1" and "2"
  expect_equal(pop@id[1:2], c("1","2"))

  p = biparentalPed()
  out = pedigreeCross(pop, p$id, p$mother, p$father, matchID=TRUE,
                      simParam=d$SP)

  expect_equal(nInd(out), 10L)
  expect_equal(out@id, p$id)

  # The two founders are the founderPop individuals of those names, copied
  # rather than generated, so their genotypes match exactly
  fromPed = pullSegSiteGeno(out, simParam=d$SP)[1:2,,drop=FALSE]
  fromPop = pullSegSiteGeno(pop[c("1","2")], simParam=d$SP)
  expect_equal(unname(fromPed), unname(fromPop))
})

test_that("matchID reports founders that are absent from the population", {
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("x","y","z"),
                             mother=c("0","0","x"),
                             father=c("0","0","y"),
                             matchID=TRUE, simParam=d$SP),
               "missing")
})

test_that("a pedigree needing more founders than supplied is refused", {
  d = pedSP(nInd=2)
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("1","2","3","4"),
                             mother=c("0","0","0","1"),
                             father=c("0","0","0","2"),
                             simParam=d$SP),
               "founders")
})

test_that("DH and nSelf are applied to the individuals they name", {
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)

  id = c("1","2","3","4")
  mother = c("0","0","1","1")
  father = c("0","0","2","2")

  set.seed(303)
  out = pedigreeCross(pop, id, mother, father,
                      DH=c(FALSE,FALSE,TRUE,FALSE),
                      nSelf=c(0,0,0,2),
                      simParam=d$SP)

  expect_equal(nInd(out), 4L)
  expect_equal(out@id, id)

  # A doubled haploid is homozygous everywhere
  geno = pullSegSiteGeno(out, simParam=d$SP)
  expect_true(all(geno[3,] %in% c(0,2)))
  expect_descendant(geno[3,], geno[1,], geno[2,], label="doubled haploid")

  # Individual 4 is the cross of 1 and 2 followed by two generations of
  # selfing, so only the weaker ancestry invariant applies to it
  expect_descendant(geno[4,], geno[1,], geno[2,], label="selfed individual")
})

test_that("the input vectors are checked against each other", {
  d = pedSP()
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("1","1"), mother=c("0","0"),
                             father=c("0","0"), simParam=d$SP),
               "duplicates")
  expect_error(pedigreeCross(pop, id=c("1","2"), mother=c("0"),
                             father=c("0","0"), simParam=d$SP),
               "length\\(mother\\)")
  expect_error(pedigreeCross(pop, id=c("1","2"), mother=c("0","0"),
                             father=c("0"), simParam=d$SP),
               "length\\(father\\)")
  expect_error(pedigreeCross(pop, id=c("1","2"), mother=c("0","0"),
                             father=c("0","0"), DH=TRUE, simParam=d$SP),
               "length\\(DH\\)")
  expect_error(pedigreeCross(pop, id=c("1","2"), mother=c("0","0"),
                             father=c("0","0"), nSelf=0, simParam=d$SP),
               "length\\(nSelf\\)")
})

test_that("pedigreeCross refuses a simulation that uses sexes", {
  set.seed(11021)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setSexes("yes_sys")
  SP$addTraitA(nQtlPerChr=5)
  pop = newPop(founderPop, simParam=SP)

  expect_error(pedigreeCross(pop, id=c("1","2","3"), mother=c("0","0","1"),
                             father=c("0","0","2"), simParam=SP),
               "sex")
})

test_that("a pedigree given out of order still sorts", {
  # sortPed makes repeated passes, so a pedigree listed youngest first needs
  # one pass per generation but must still come out right
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)

  id     = c("5","4","3","2","1")
  mother = c("4","3","2","1","0")
  father = c("4","3","2","1","0")

  out = pedigreeCross(pop, id, mother, father, simParam=d$SP)
  expect_equal(nInd(out), 5L)
  expect_equal(out@id, id)

  geno = pullSegSiteGeno(out, simParam=d$SP)
  for(i in 1:4){
    m = match(mother[i], id)
    expect_mendelian(geno[i,], geno[m,], geno[m,], label=paste("row", i))
  }
})

test_that("a single absent parent is already looked up in founderPop", {
  # This is the asymmetry issue #79 ask 4 is about. "2" has no pedigree row
  # but is found in founderPop, and one missing parent is handled today.
  # Two missing parents are not, which is what the ISSUE79 block below
  # specifies.
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)

  out = pedigreeCross(pop, id=c("1","kid"),
                      mother=c("0","1"),
                      father=c("0","2"),
                      matchID=TRUE, simParam=d$SP)

  expect_equal(nInd(out), 2L)
  expect_equal(out@id, c("1","kid"))

  geno = pullSegSiteGeno(out, simParam=d$SP)
  par2 = pullSegSiteGeno(pop["2"], simParam=d$SP)[1,]
  expect_mendelian(geno[2,], geno[1,], par2, label="kid")
})

test_that("an unsortable pedigree is an error", {
  d = pedSP()
  pop = newPop(d$map, simParam=d$SP)

  # A and B are each other's parents
  expect_error(pedigreeCross(pop, id=c("A","B"), mother=c("B","A"),
                             father=c("B","A"), simParam=d$SP))
})

# ---------------------------------------------------------------------------
# PART 2  Issue #79. These fail until the feature is written.
# ---------------------------------------------------------------------------

test_that("ISSUE79 pedigreeCross accepts a MapPop", {
  # Ask 1. A MapPop has no id slot, so only matchID=FALSE makes sense
  d = pedSP(nInd=4)
  p = biparentalPed()

  set.seed(404)
  out = pedigreeCross(d$map, p$id, p$mother, p$father, simParam=d$SP)

  expect_true(isPop(out))
  expect_equal(nInd(out), 10L)
  expect_equal(out@id, p$id)
  expect_true(isTRUE(validObject(out, test=TRUE)))

  geno = pullSegSiteGeno(out, simParam=d$SP)
  for(i in 3:10){
    m = match(p$mother[i], p$id)
    f = match(p$father[i], p$id)
    expect_mendelian(geno[i,], geno[m,], geno[f,], label=paste("individual", i))
  }
})

test_that("ISSUE79 pedigreeCross accepts a NamedMapPop with matchID", {
  # Ask 1. A NamedMapPop carries ids, so it can be matched against
  d = pedSP(nInd=4)
  named = asNamedMapPop(d$map, id=c("F1","F2","F3","F4"))

  out = pedigreeCross(named, id=c("F1","F2","kid"),
                      mother=c("0","0","F1"),
                      father=c("0","0","F2"),
                      matchID=TRUE, simParam=d$SP)

  expect_true(isPop(out))
  expect_equal(nInd(out), 3L)
  expect_equal(out@id, c("F1","F2","kid"))

  geno = pullSegSiteGeno(out, simParam=d$SP)
  expect_mendelian(geno[3,], geno[1,], geno[2,], label="kid")
})

test_that("ISSUE79 matchID on a MapPop says why it cannot work", {
  # Ask 1 and ask 2. A plain MapPop has no ids to match, and the error
  # should say so rather than failing on a missing slot
  d = pedSP(nInd=4)
  p = biparentalPed()

  expect_error(pedigreeCross(d$map, p$id, p$mother, p$father,
                             matchID=TRUE, simParam=d$SP),
               "matchID")
})

test_that("ISSUE79 population members can be founders without pedigree rows", {
  # Ask 4. "1" and "2" are in founderPop but have no rows of their own.
  # Today both parents being absent makes the child itself a founder and
  # the two names are discarded; they should be looked up instead.
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)
  expect_equal(pop@id[1:2], c("1","2"))

  out = pedigreeCross(pop, id="kid", mother="1", father="2",
                      matchID=TRUE, simParam=d$SP)

  expect_equal(nInd(out), 1L)
  expect_equal(out@id, "kid")

  kid = pullSegSiteGeno(out, simParam=d$SP)[1,]
  par = pullSegSiteGeno(pop[c("1","2")], simParam=d$SP)
  expect_mendelian(kid, par[1,], par[2,], label="implicit founder cross")
})

test_that("ISSUE79 replacing an existing id is an error", {
  # Ask 3. "3" names an individual in founderPop, but the pedigree gives it
  # parents, so the pedigree would generate a different individual under a
  # name already in use
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)
  expect_equal(pop@id[1:3], c("1","2","3"))

  expect_error(pedigreeCross(pop, id=c("1","2","3"),
                             mother=c("0","0","1"),
                             father=c("0","0","2"),
                             matchID=TRUE, simParam=d$SP),
               "3")

  # The same pedigree under a name that is free is fine
  out = pedigreeCross(pop, id=c("1","2","kid"),
                      mother=c("0","0","1"),
                      father=c("0","0","2"),
                      matchID=TRUE, simParam=d$SP)
  expect_equal(nInd(out), 3L)

  # And matchID=FALSE is unaffected, because ids are not matched at all
  set.seed(505)
  free = pedigreeCross(pop, id=c("1","2","3"),
                       mother=c("0","0","1"),
                       father=c("0","0","2"),
                       simParam=d$SP)
  expect_equal(nInd(free), 3L)
})

test_that("ISSUE79 a cycle is reported as a cycle and names the individuals", {
  # Ask 2
  d = pedSP()
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("A","B"), mother=c("B","A"),
                             father=c("B","A"), simParam=d$SP),
               "A")
  expect_error(pedigreeCross(pop, id=c("A","B"), mother=c("B","A"),
                             father=c("B","A"), simParam=d$SP),
               "cycle")

  # An individual that is its own parent is the same failure
  expect_error(pedigreeCross(pop, id=c("1","2","3"),
                             mother=c("0","0","3"),
                             father=c("0","0","2"),
                             simParam=d$SP),
               "3")
})

test_that("ISSUE79 running out of cycles is reported separately", {
  # Ask 2. The current message blames "maxGen", which is not an argument
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)

  id     = c("5","4","3","2","1")
  mother = c("4","3","2","1","0")
  father = c("4","3","2","1","0")

  expect_error(pedigreeCross(pop, id, mother, father, maxCycle=2,
                             simParam=d$SP),
               "maxCycle")
  # The same pedigree sorts when given enough passes
  expect_equal(nInd(pedigreeCross(pop, id, mother, father, maxCycle=100,
                                  simParam=d$SP)), 5L)
})

test_that("ISSUE79 an unknown parent is warned about", {
  # Ask 2, the incomplete pedigree case. "ghost" is neither a pedigree row
  # nor a founderPop individual, so it becomes a silent extra founder today
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)

  expect_warning(pedigreeCross(pop, id=c("1","2","3"),
                               mother=c("0","0","ghost"),
                               father=c("0","0","2"),
                               simParam=d$SP),
                 "ghost")

  # The conventional unknown marker is not worth a warning
  expect_warning(pedigreeCross(pop, id=c("1","2","3"),
                               mother=c("0","0","1"),
                               father=c("0","0","2"),
                               simParam=d$SP),
                 NA)
})

test_that("ISSUE79 DH and nSelf are validated", {
  # Ask 2
  d = pedSP(nInd=4)
  pop = newPop(d$map, simParam=d$SP)
  id = c("1","2","3")
  mother = c("0","0","1")
  father = c("0","0","2")

  expect_error(pedigreeCross(pop, id, mother, father,
                             nSelf=c(0,0,-1), simParam=d$SP),
               "nSelf")
  expect_error(pedigreeCross(pop, id, mother, father,
                             nSelf=c(0,0,NA), simParam=d$SP),
               "nSelf")
  expect_error(pedigreeCross(pop, id, mother, father,
                             DH=c(FALSE,FALSE,NA), simParam=d$SP),
               "DH")
})

test_that("ISSUE79 an empty pedigree is rejected cleanly", {
  # Ask 2. Today this reaches max() on an empty vector
  d = pedSP()
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=character(0), mother=character(0),
                             father=character(0), simParam=d$SP),
               "empty")
})


# ---------------------------------------------------------------------------
# PART 3  Issue #131, half founders. An individual with one parent known and
# one unknown needs a founder genome for the unknown parent, and no two of
# them may share one.
# ---------------------------------------------------------------------------

test_that("ISSUE131 the pedigree from the bug report works", {
  # Individual 3 has no mother and individual 4 has no father
  d = pedSP(nInd=4, nChr=1, segSites=10, seed=13101)
  pop = newPop(d$map, simParam=d$SP)

  id = as.character(1:12)
  mother = as.character(c(0,0,0,1,1,3:9))
  father = as.character(c(0,0,2,0,2,3:9))

  set.seed(606)
  out = pedigreeCross(pop, id, mother, father, simParam=d$SP)

  expect_equal(nInd(out), 12L)
  expect_equal(out@id, id)
  expect_equal(out@mother, mother)
  expect_equal(out@father, father)
  expect_true(isTRUE(validObject(out, test=TRUE)))
})

test_that("ISSUE131 several half founders of the same kind each get a genome", {
  # Every unknown parent is coded "0", so counting founders by name gives
  # one where three are needed, and the founder population is undersized
  d = pedSP(nInd=8, nChr=1, segSites=10, seed=13111)
  pop = newPop(d$map, simParam=d$SP)

  id     = c("1","2","A","B","C")
  mother = c("0","0","0","0","0")
  father = c("0","0","1","1","2")

  set.seed(707)
  out = pedigreeCross(pop, id, mother, father, simParam=d$SP)

  expect_equal(nInd(out), 5L)
  expect_equal(out@id, id)
  expect_true(isTRUE(validObject(out, test=TRUE)))
})

test_that("ISSUE131 two half founders never share a founder genome", {
  # A and B have the same father and each has an unknown mother. If the two
  # unknown mothers collapse onto one founder, the number of distinct
  # founders behind A and B drops from three to two. IBD makes that visible.
  d = pedSP(nInd=8, nChr=1, segSites=40, seed=13121, trackRec=TRUE)
  pop = newPop(d$map, simParam=d$SP)

  id     = c("1","2","A","B")
  mother = c("0","0","0","0")
  father = c("0","0","1","1")

  set.seed(808)
  out = pedigreeCross(pop, id, mother, father, simParam=d$SP)
  expect_equal(nInd(out), 4L)

  ibd = pullIbdHaplo(out, simParam=d$SP)
  # Founder haplotypes are numbered two per individual, in order
  founderOf = function(hap) return((hap+1L)%/%2L)
  rows = match(c("A_1","A_2","B_1","B_2"), rownames(ibd))
  expect_false(any(is.na(rows)))
  used = unique(founderOf(c(ibd[rows,])))

  # One shared father plus two distinct mothers
  expect_equal(length(used), 3L)
})

test_that("ISSUE131 NA and 0 both mean unknown", {
  d = pedSP(nInd=4, nChr=1, segSites=10, seed=13131)
  pop = newPop(d$map, simParam=d$SP)

  id     = c("1","2","A","B")
  mother = c("0", NA, "0", NA)
  father = c(NA, "0", "1", "2")

  set.seed(909)
  out = pedigreeCross(pop, id, mother, father, simParam=d$SP)

  # Two founder rows and two unknown mothers, so four founders in total
  expect_equal(nInd(out), 4L)
  expect_true(isTRUE(validObject(out, test=TRUE)))
})

test_that("ISSUE131 half founders are counted in the founder budget", {
  # Two founder rows and two half founders need four, not two
  d = pedSP(nInd=3, nChr=1, segSites=10, seed=13141)
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=as.character(1:5),
                             mother=as.character(c(0,0,0,1,1)),
                             father=as.character(c(0,0,2,0,2)),
                             simParam=d$SP),
               "founders")
})

test_that("ISSUE131 Gregor's partial pedigree works", {
  # The worked example from the issue thread, which needs exactly eight
  # founders: four founder rows and four half founders
  d = pedSP(nInd=8, nChr=1, segSites=10, seed=13151)
  pop = newPop(d$map, simParam=d$SP)

  id     = as.character(1:10)
  mother = as.character(c(0,0,1,0,1,3,3,0,0,NA))
  father = as.character(c(0,0,0,2,2,4,0,4,0,NA))

  set.seed(1010)
  out = pedigreeCross(pop, id, mother, father, simParam=d$SP)

  expect_equal(nInd(out), 10L)
  expect_equal(out@id, id)
  expect_true(isTRUE(validObject(out, test=TRUE)))
})

test_that("ISSUE131 matchID draws unknown parents from unnamed individuals", {
  d = pedSP(nInd=6, nChr=1, segSites=40, seed=13161, trackRec=TRUE)
  pop = newPop(d$map, simParam=d$SP)
  expect_equal(pop@id, as.character(1:6))

  id     = c("1","2","A")
  mother = c("0","0","0")
  father = c("0","0","1")

  set.seed(1111)
  out = pedigreeCross(pop, id, mother, father, matchID=TRUE, simParam=d$SP)
  expect_equal(nInd(out), 3L)

  ibd = pullIbdHaplo(out, simParam=d$SP)
  founderOf = function(hap) return((hap+1L)%/%2L)
  rows = match(c("A_1","A_2"), rownames(ibd))
  used = sort(unique(founderOf(c(ibd[rows,]))))

  # The father is individual 1, and the unknown mother has to be one of the
  # individuals the pedigree does not name, so not individual 2
  expect_equal(length(used), 2L)
  expect_true(1L %in% used)
  expect_false(2L %in% used)
  expect_true(all(used %in% c(1L,3L,4L,5L,6L)))
})

test_that("ISSUE131 matchID needs unnamed individuals to draw from", {
  # Every individual in founderPop is named by the pedigree, so there is
  # nobody left to be the unknown mother
  d = pedSP(nInd=2, nChr=1, segSites=10, seed=13171)
  pop = newPop(d$map, simParam=d$SP)

  expect_error(pedigreeCross(pop, id=c("1","2","A"),
                             mother=c("0","0","0"),
                             father=c("0","0","1"),
                             matchID=TRUE, simParam=d$SP),
               "founders")
})
