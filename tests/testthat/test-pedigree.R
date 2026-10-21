context("pedigree")

# The relationship functions are checked against Henderson's tabular method
# written out in R without any rule for doubled haploids. A doubled haploid
# is instead written as a long chain of selfing: each generation of selfing
# keeps an individual's relationships to everyone else and halves its
# distance from full inbreeding, so after 40 generations it differs from a
# doubled haploid by about 1e-12. That makes the reference independent of
# the doubled haploid rule it checks.

# 1. Helpers ----

refPedA = function(mother, father){
  n = length(mother)
  A = matrix(0, n, n)
  for(i in seq_len(n)){
    s = mother[i]
    d = father[i]
    for(j in seq_len(i-1L)){
      a = 0
      if(s>0) a = a + A[j,s]/2
      if(d>0) a = a + A[j,d]/2
      A[i,j] = A[j,i] = a
    }
    A[i,i] = if(s>0 && d>0) 1 + A[s,d]/2 else 1
  }
  A
}

# Replaces every doubled haploid with nSelf generations of selfing from its
# parent, and every inbred founder with a founder followed by nSelf
# generations of selfing. Returns the new pedigree and where each original
# individual is in it.
expandDH = function(mother, father, isDH, nSelf=40L){
  n = length(mother)
  newMother = newFather = integer(0)
  pos = integer(n)
  for(i in seq_len(n)){
    s = if(mother[i]>0) pos[mother[i]] else 0L
    d = if(father[i]>0) pos[father[i]] else 0L
    if(isDH[i]==1L){
      if(s==0L){
        newMother = c(newMother, 0L)
        newFather = c(newFather, 0L)
        s = length(newMother)
      }
      for(k in seq_len(nSelf)){
        newMother = c(newMother, s)
        newFather = c(newFather, s)
        s = length(newMother)
      }
    }else{
      newMother = c(newMother, s)
      newFather = c(newFather, d)
    }
    pos[i] = length(newMother)
  }
  list(mother=newMother, father=newFather, pos=pos)
}

refPedADH = function(mother, father, isDH){
  ex = expandDH(mother, father, isDH)
  refPedA(ex$mother, ex$father)[ex$pos, ex$pos]
}

tripToMat = function(trip){
  M = matrix(0, trip$nrow, trip$ncol)
  M[cbind(trip$i, trip$j)] = trip$x
  M
}

# The tabular method for ploidy v, written from the pairs of alleles: of
# the choose(v,2) pairs in an individual, choose(v/2,2) lie within each
# parent's gamete and v^2/4 pair one allele from each gamete
refPedAPloidy = function(mother, father, inbred, v){
  n = length(mother)
  A = matrix(0, n, n)
  Fv = numeric(n)
  for(i in seq_len(n)){
    s = mother[i]
    d = father[i]
    for(j in seq_len(i-1L)){
      a = 0
      if(s>0) a = a + A[j,s]/2
      if(d>0) a = a + A[j,d]/2
      A[i,j] = A[j,i] = a
    }
    if(s==0 && d==0){
      Fv[i] = if(inbred[i]==1L) 1 else 0
    }else{
      Fs = if(s>0) Fv[s] else 0
      Fd = if(d>0) Fv[d] else 0
      asd = if(s>0 && d>0) A[s,d] else 0
      Fv[i] = (choose(v/2,2)*(Fs+Fd) + (v/2)^2*asd/v)/choose(v,2)
    }
    A[i,i] = 1 + (v-1)*Fv[i]
  }
  A
}

# The relationship matrix implied by the inverse and its basis, W*A_b*W'
ainvToA = function(ans){
  W = tripToMat(ans$W)
  W %*% solve(tripToMat(ans$Ainv)) %*% t(W)
}

# A random pedigree with crosses, selfs, doubled haploids and, optionally,
# inbred founders
randPed = function(nFounder, nInd, pDH=0, pSelf, pInbred=0, seed){
  set.seed(seed)
  n = nFounder + nInd
  mother = father = isDH = integer(n)
  isDH[seq_len(nFounder)] = as.integer(runif(nFounder)<pInbred)
  for(i in (nFounder+1L):n){
    u = runif(1)
    if(u<pDH){
      mother[i] = father[i] = sample.int(i-1L, 1L)
      isDH[i] = 1L
    }else if(u<(pDH+pSelf)){
      mother[i] = father[i] = sample.int(i-1L, 1L)
    }else{
      p = sample.int(i-1L, 2L)
      mother[i] = p[1]
      father[i] = p[2]
    }
  }
  list(mother=mother, father=father, isDH=isDH)
}

# A simulation with crosses, doubled haploids and selfs, tracking the
# pedigree
pedSim = function(seed, inbred=FALSE){
  set.seed(seed)
  founderPop = quickHaplo(nInd=6, nChr=1, segSites=10, inbred=inbred)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  pop1 = newPop(founderPop, simParam=SP)
  pop2 = randCross(pop1, nCrosses=4, nProgeny=2, simParam=SP)
  pop3 = makeDH(pop2, nDH=2, simParam=SP)
  pop4 = self(pop3, simParam=SP)
  pop5 = randCross(c(pop3, pop4), nCrosses=6, simParam=SP)
  pop6 = makeDH(pop5, nDH=1, simParam=SP)
  list(SP=SP, pop1=pop1, pop2=pop2, pop3=pop3, pop4=pop4, pop5=pop5,
       pop6=pop6)
}

# 2. Values worked by hand ----

test_that("VALUE crosses, selfs and doubled haploids match hand calculation", {
  # 1 and 2 are founders, 3 is their cross, 4 a self of 3, 5 a doubled
  # haploid of 3, 6 a doubled haploid of 5, and 7 a cross of 5 and 4
  mother = c(0L, 0L, 1L, 3L, 3L, 5L, 5L)
  father = c(0L, 0L, 2L, 3L, 3L, 5L, 4L)
  isDH =   c(0L, 0L, 0L, 0L, 1L, 1L, 0L)
  A = AlphaSimR:::calcPedA(mother, father, isDH)

  expect_equal(diag(A)[1:6], c(1, 1, 1, 1.5, 2, 2))
  expect_equal(A[1,3], 0.5)
  # A self and a doubled haploid relate to others as their parent does
  expect_equal(A[4,1:3], A[3,1:3])
  expect_equal(A[5,1:4], A[3,1:4])
  # A doubled haploid of a doubled haploid is a copy of it
  expect_equal(A[6,], A[5,])
  # The cross of 5 and 4 has F of half their relationship, which is that
  # of 3 with a self of itself, 1
  expect_equal(A[5,4], 1)
  expect_equal(A[7,7], 1.5)

  F = AlphaSimR:::calcPedInbreeding(mother, father, isDH)
  expect_equal(c(F), diag(A)-1)
})

test_that("VALUE inbred founders are fully inbred and unrelated", {
  mother = c(0L, 0L, 1L)
  father = c(0L, 0L, 2L)
  isDH =   c(1L, 1L, 0L)
  A = AlphaSimR:::calcPedA(mother, father, isDH)
  expect_equal(A, matrix(c(2, 0, 1,
                           0, 2, 1,
                           1, 1, 1), 3, 3))
  F = AlphaSimR:::calcPedInbreeding(mother, father, isDH)
  expect_equal(c(F), c(1, 1, 0))
})

test_that("VALUE polyploid selfs and crosses match hand calculation", {
  # 1 and 2 are founders, 3 their cross, 4 a self of 3, 5 a cross of 3 and 4
  mother = c(0L, 0L, 1L, 3L, 3L)
  father = c(0L, 0L, 2L, 3L, 4L)
  isDH = integer(5)
  # Tetraploid: a self of a non-inbred parent has F of 1/6, and the cross
  # of 3 and 4 has F of (1/6 + 1)/6. AGHmatrix gives 1.583333 for it.
  A4 = AlphaSimR:::calcPedA(mother, father, isDH, ploidy=4L)
  expect_equal(diag(A4), c(1, 1, 1, 1.5, 1+3*7/36))
  expect_equal(c(AlphaSimR:::calcPedInbreeding(mother, father, isDH,
                                               ploidy=4L)),
               c(0, 0, 0, 1/6, 7/36))
  # Hexaploid: F of 0.1 for the self and 0.12 for the cross
  A6 = AlphaSimR:::calcPedA(mother, father, isDH, ploidy=6L)
  expect_equal(diag(A6), c(1, 1, 1, 1.5, 1.6))
  # Off-diagonals are the same at every ploidy
  A2 = AlphaSimR:::calcPedA(mother, father, isDH)
  expect_equal(A4[upper.tri(A4)], A2[upper.tri(A2)])
})

# 3. Against the reference ----

test_that("REFERENCE relationships match the tabular method with DH as selfing", {
  for(seed in 1:5){
    ped = randPed(nFounder=5, nInd=40, pDH=0.25, pSelf=0.2, pInbred=0.4,
                  seed=8100+seed)
    A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH)
    ref = refPedADH(ped$mother, ped$father, ped$isDH)
    expect_equal(A, ref, tolerance=1e-9)
  }
})

test_that("REFERENCE inbreeding coefficients are the diagonal of A less one", {
  for(seed in 1:5){
    ped = randPed(nFounder=5, nInd=60, pDH=0.2, pSelf=0.3, pInbred=0.4,
                  seed=8200+seed)
    A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH)
    F = AlphaSimR:::calcPedInbreeding(ped$mother, ped$father, ped$isDH)
    expect_equal(c(F), diag(A)-1, tolerance=1e-12)
  }
})

test_that("REFERENCE Colleau's method gives the block of the tabular matrix", {
  for(seed in 1:5){
    ped = randPed(nFounder=5, nInd=60, pDH=0.25, pSelf=0.2, pInbred=0.4,
                  seed=8700+seed)
    A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH)
    take = sample.int(length(ped$mother), 20L)
    B = AlphaSimR:::calcPedASub(ped$mother, ped$father, ped$isDH, take)
    expect_equal(B, A[take, take], tolerance=1e-12)
    expect_true(isSymmetric(B, tol=0))
  }
})

test_that("REFERENCE polyploid relationships, inbreeding and blocks", {
  for(v in c(4L, 6L)){
    for(seed in 1:3){
      ped = randPed(nFounder=5, nInd=60, pDH=0, pSelf=0.3, pInbred=0.4,
                    seed=8900+10*v+seed)
      A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH, ploidy=v)
      expect_equal(A, refPedAPloidy(ped$mother, ped$father, ped$isDH, v),
                   tolerance=1e-12)
      F = AlphaSimR:::calcPedInbreeding(ped$mother, ped$father, ped$isDH,
                                        ploidy=v)
      expect_equal(c(F), (diag(A)-1)/(v-1), tolerance=1e-12)
      take = sample.int(length(ped$mother), 20L)
      expect_equal(AlphaSimR:::calcPedASub(ped$mother, ped$father, ped$isDH,
                                           take, ploidy=v),
                   A[take, take], tolerance=1e-12)
      ans = AlphaSimR:::calcPedAinvTriplet(ped$mother, ped$father, ped$isDH,
                                           ploidy=v)
      expect_equal(ainvToA(ans), A, tolerance=1e-9)
    }
  }
})

# 4. Inverse ----

test_that("INVERSE is the inverse of A when there are no clones", {
  # Without doubled haploids or inbred founders no Mendelian sampling
  # variance is zero, so the basis is the whole pedigree
  ped = randPed(nFounder=5, nInd=60, pDH=0, pSelf=0.3, seed=8301)
  ans = AlphaSimR:::calcPedAinvTriplet(ped$mother, ped$father, ped$isDH)
  expect_equal(ans$basis, seq_along(ped$mother))
  A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH)
  expect_equal(tripToMat(ans$Ainv) %*% A, diag(nrow(A)), tolerance=1e-10)
})

test_that("INVERSE reproduces A through the basis with doubled haploids", {
  for(seed in 1:5){
    ped = randPed(nFounder=5, nInd=60, pDH=0.25, pSelf=0.2, pInbred=0.4,
                  seed=8400+seed)
    ans = AlphaSimR:::calcPedAinvTriplet(ped$mother, ped$father, ped$isDH)
    A = AlphaSimR:::calcPedA(ped$mother, ped$father, ped$isDH)
    expect_equal(ainvToA(ans), A, tolerance=1e-9)
    # The inverse is symmetric and the basis has no zero variances
    Ainv = tripToMat(ans$Ainv)
    expect_equal(Ainv, t(Ainv))
    expect_true(all(ans$D[ans$basis]>1e-10))
  }
})

test_that("INVERSE leaves clones and F1s of inbred parents out of the basis", {
  # 1 and 2 are inbred founders, 3 their F1, 4 a doubled haploid of 1, 5 a
  # self of 4, and 6 a doubled haploid of 3
  mother = c(0L, 0L, 1L, 1L, 4L, 3L)
  father = c(0L, 0L, 2L, 1L, 4L, 3L)
  isDH =   c(1L, 1L, 0L, 1L, 0L, 1L)
  ans = AlphaSimR:::calcPedAinvTriplet(mother, father, isDH)
  expect_equal(ans$basis, c(1L, 2L, 6L))
  W = tripToMat(ans$W)
  expect_equal(W[3,], c(0.5, 0.5, 0))
  expect_equal(W[4,], c(1, 0, 0))
  expect_equal(W[5,], c(1, 0, 0))
  expect_equal(ainvToA(ans), AlphaSimR:::calcPedA(mother, father, isDH))
  expect_equal(AlphaSimR:::calcPedASub(mother, father, isDH, 6:1),
               AlphaSimR:::calcPedA(mother, father, isDH)[6:1, 6:1])
})

# 5. From a simulation ----

test_that("SIMULATION pedigreeA matches the reference on the recorded pedigree", {
  for(inbred in c(FALSE, TRUE)){
    sim = pedSim(seed=8501, inbred=inbred)
    ped = sim$SP$pedigree
    ref = refPedADH(ped[,1], ped[,2], ped[,3])
    for(pop in sim[-1]){
      A = pedigreeA(pop, simParam=sim$SP)
      expect_equal(unname(A), ref[pop@iid, pop@iid], tolerance=1e-9)
      expect_identical(dimnames(A), list(pop@id, pop@id))
      F = pedigreeInbreeding(pop, simParam=sim$SP)
      expect_equal(unname(F), diag(ref)[pop@iid]-1, tolerance=1e-9)
      expect_identical(names(F), pop@id)
    }
  }
})

test_that("SIMULATION doubled haploids are fully inbred and copy their parent", {
  sim = pedSim(seed=8502)
  F = pedigreeInbreeding(sim$pop3, simParam=sim$SP)
  expect_true(all(F==1))
  both = c(sim$pop2, sim$pop3)
  A = pedigreeA(both, simParam=sim$SP)
  # The doubled haploids of pop2's first individual are its first two
  dh = sim$pop2@nInd + 1:2
  expect_equal(unname(A[dh[1], 1:sim$pop2@nInd]), unname(A[1, 1:sim$pop2@nInd]))
  expect_equal(unname(A[dh[2], 1:sim$pop2@nInd]), unname(A[1, 1:sim$pop2@nInd]))
})

test_that("PRUNE subsets and repeats give the same relationships", {
  sim = pedSim(seed=8503)
  all = c(sim$pop4, sim$pop5)
  Aall = pedigreeA(all, simParam=sim$SP)
  take = c(3L, 1L, 12L, 3L)
  Asub = pedigreeA(all[take], simParam=sim$SP)
  expect_equal(unname(Asub), unname(Aall[take, take]))
})

test_that("PRUNE resetPed keeps the pedigree usable", {
  sim = pedSim(seed=8504)
  sim$SP$resetPed(sim$pop2@iid[sim$pop2@nInd])
  pop3 = makeDH(sim$pop2, nDH=1, simParam=sim$SP)
  F = pedigreeInbreeding(pop3, simParam=sim$SP)
  expect_true(all(F==1))
  A = pedigreeA(c(sim$pop2, pop3), simParam=sim$SP)
  n = sim$pop2@nInd
  expect_equal(unname(A[n+1:n, 1:n]), unname(A[1:n, 1:n]))
})

test_that("SIMULATION tetraploid pedigreeA matches the reference", {
  set.seed(8510)
  founderPop = quickHaplo(nInd=6, nChr=1, segSites=10, ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  pop1 = newPop(founderPop, simParam=SP)
  pop2 = randCross(pop1, nCrosses=4, nProgeny=2, simParam=SP)
  pop3 = self(pop2, simParam=SP)
  pop4 = randCross(c(pop2, pop3), nCrosses=6, simParam=SP)
  ped = SP$pedigree
  ref = refPedAPloidy(ped[,1], ped[,2], ped[,3], 4L)
  for(pop in list(pop2, pop3, pop4)){
    expect_equal(unname(pedigreeA(pop, simParam=SP)), ref[pop@iid, pop@iid],
                 tolerance=1e-12)
    expect_equal(unname(pedigreeInbreeding(pop, simParam=SP)),
                 (diag(ref)[pop@iid]-1)/3, tolerance=1e-12)
  }
})

# 6. Cache ----

test_that("CACHE known inbreeding coefficients are used and the rest computed", {
  ped = randPed(nFounder=5, nInd=80, pDH=0.2, pSelf=0.2, pInbred=0.3,
                seed=8801)
  F = c(AlphaSimR:::calcPedInbreeding(ped$mother, ped$father, ped$isDH))
  known = F
  known[sample.int(length(F), 40L)] = NA
  expect_equal(c(AlphaSimR:::calcPedInbreeding(ped$mother, ped$father,
                                               ped$isDH, known)), F)
  expect_equal(c(AlphaSimR:::calcPedInbreeding(ped$mother, ped$father,
                                               ped$isDH, rep(NA_real_, length(F)))), F)
  # A value given as known is taken as it is
  known = rep(NA_real_, length(F))
  known[1] = 0.25
  expect_equal(AlphaSimR:::calcPedInbreeding(ped$mother, ped$father,
                                             ped$isDH, known)[1], 0.25)
})

test_that("CACHE pedigreeInbreeding fills the SimParam's cache and reuses it", {
  sim = pedSim(seed=8802)
  cache = function() sim$SP$.__enclos_env__$private$.pedF
  F3 = pedigreeInbreeding(sim$pop3, simParam=sim$SP)
  # pop3 and its ancestors are cached, and later individuals are not
  expect_equal(cache()[sim$pop3@iid], unname(F3))
  expect_true(all(is.na(cache()[sim$pop6@iid])))
  # A cached value is reused rather than recomputed
  sim$SP$.__enclos_env__$private$.pedF[sim$pop3@iid[1]] = 0.123
  expect_equal(unname(pedigreeInbreeding(sim$pop3, simParam=sim$SP)[1]), 0.123)
})

test_that("CACHE resetPed cuts the cache back with the pedigree", {
  sim = pedSim(seed=8803)
  pedigreeInbreeding(sim$pop6, simParam=sim$SP)
  lastId = sim$pop2@iid[sim$pop2@nInd]
  sim$SP$resetPed(lastId)
  expect_true(length(sim$SP$.__enclos_env__$private$.pedF)<=lastId)
})

# 7. Errors ----

test_that("ERROR pedigree relationships need the pedigree tracked", {
  set.seed(8601)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop, simParam=SP)
  expect_error(pedigreeA(pop, simParam=SP), "trackPed")
  expect_error(pedigreeInbreeding(pop, simParam=SP), "trackPed")
})

test_that("ERROR tracking turned on late is refused", {
  set.seed(8602)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  pop = newPop(founderPop, simParam=SP)
  SP$setTrackPed(TRUE, force=TRUE)
  pop2 = randCross(pop, nCrosses=2, simParam=SP)
  expect_error(pedigreeA(pop2, simParam=SP), "does not cover")
})

test_that("ERROR odd ploidy levels are refused", {
  set.seed(8603)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10, ploidy=3L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$setTrackPed(TRUE)
  pop = newPop(founderPop, simParam=SP)
  expect_error(pedigreeA(pop, simParam=SP), "even ploidy")
  expect_error(pedigreeInbreeding(pop, simParam=SP), "even ploidy")
})

test_that("ERROR malformed pedigrees are refused", {
  expect_error(AlphaSimR:::calcPedA(c(0L, 2L), c(0L, 0L), c(0L, 0L)),
               "before its progeny")
  expect_error(AlphaSimR:::calcPedA(c(0L, 0L, 1L), c(0L, 0L, 2L),
                                    c(0L, 0L, 1L)),
               "same mother and father")
  expect_error(AlphaSimR:::calcPedInbreeding(c(0L, NA), c(0L, 1L),
                                             c(0L, 0L)),
               "before its progeny")
  expect_error(AlphaSimR:::calcPedAinvTriplet(c(0L, 1L), c(0L), c(0L, 0L)),
               "same length")
  expect_error(AlphaSimR:::calcPedASub(c(0L, 0L), c(0L, 0L), c(0L, 0L), 3L),
               "row of the pedigree")
  expect_error(AlphaSimR:::calcPedInbreeding(c(0L, 0L), c(0L, 0L),
                                             c(0L, 0L), c(0, 0, 0)),
               "one value per individual")
  expect_error(AlphaSimR:::calcPedA(c(0L, 1L), c(0L, 1L), c(0L, 1L),
                                    ploidy=4L),
               "only supported for diploids")
  expect_error(AlphaSimR:::calcPedA(0L, 0L, 0L, ploidy=3L), "even ploidy")
})
