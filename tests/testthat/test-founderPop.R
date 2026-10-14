context("Founder population edge cases")

test_that("REFUSE inconsistent haplotype dimensions", {
  expect_error(newMapPop(list(c(0.2,0.8),c(0.2,0.8)),
                         list(matrix(0,4,2),matrix(0,6,2))),
               "Number of rows must be equal")
  expect_error(newMapPop(list(c(0.2,0.8)),list(matrix(0,3,2))),
               "divisible by ploidy")
  expect_error(newMapPop(list(c(0.2,0.8)),list(matrix(0,4,3))),
               "segregating sites.*don't match")
})

test_that("Sampling inbred founders preserves source haplotypes and maps", {
  # Distinct haplotypes let us check sampling without replacement directly.
  h = diag(4)
  pop = newMapPop(list(chrA=c(0.1,0.3,0.5,0.8)),list(h),genLen=1)
  set.seed(6012)
  out = sampleHaplo(pop,nInd=4,inbred=TRUE,ploidy=4,replace=FALSE)
  sampled = pullSegSiteHaplo(out,nThreads=1)
  expect_true(out@inbred)
  expect_equal(out@ploidy,4L)
  for(i in seq_len(4)){
    rows = (i-1)*4+seq_len(4)
    expect_equal(unname(sampled[rows,,drop=FALSE]),
                 matrix(sampled[rows[1],],4,4,byrow=TRUE))
  }
  expect_equal(unname(sort(apply(sampled[seq(1,16,4),],1,paste,collapse=""))),
               sort(apply(h,1,paste,collapse="")))
  expect_equal(out@genMap,pop@genMap)
  expect_equal(out@genLen,pop@genLen)
  expect_error(sampleHaplo(pop,nInd=5,inbred=TRUE,replace=FALSE))
  expect_error(sampleHaplo(pop,nInd=3,inbred=FALSE,replace=FALSE))
})

test_that("Adding named sites updates inbred status only for heterozygotes", {
  pop = newMapPop(list(chrA=c(left=0.2,right=0.8)),
                  list(matrix(c(0,1,1,0),2,2)),inbred=TRUE,genLen=1)
  homo = addSegSite(pop,"homo","chrA",0.5,c(0,0,1,1))
  expect_true(homo@inbred)
  expect_equal(as.vector(pullMarkerGeno(homo,"homo",nThreads=1)),c(0,2))
  hetero = addSegSite(pop,"hetero","chrA",0.5,c(0,1,1,1),nThreads=1)
  expect_false(hetero@inbred)
  expect_equal(as.vector(pullMarkerGeno(hetero,"hetero",nThreads=1)),c(1,2))
  expect_equal(homo@genLen,pop@genLen)
  expect_equal(unname(homo@genMap[[1]]),c(0.2,0.5,0.8))
  expect_equal(pullMarkerGeno(homo,c("left","right"),nThreads=1),
               pullMarkerGeno(pop,c("left","right"),nThreads=1))
  expect_error(addSegSite(pop,"left",1,0.5,c(0,0,1,1),nThreads=1),
               "already in use")
})

test_that("MaCS rejects unknown presets and missing manual map lengths", {
  expect_error(runMacs(2,species="unknown",nThreads=1),"No rules for species")
  expect_error(runMacs(2,manualCommand="1000 -t 0.01",nThreads=1),
               "must define manualGenLen")
})

test_that("Small samples exercise all presets with valid full chromosome lengths", {
  lengths = c(GENERIC=1,CATTLE=round(2.8e9/30)*9.26e-9,WHEAT=1.43,MAIZE=2)
  for(species in names(lengths)){
    # Two haplotypes and five retained sites keep this suitable for routine CI.
    set.seed(6013)
    pop = runMacs(1,segSites=5,species=species,nThreads=1)
    expect_equal(pop@nLoci,5L,info=species)
    expect_equal(pop@genLen,unname(lengths[species]),info=species)
    expect_true(all(pop@genMap[[1]]>=0 & pop@genMap[[1]]<=pop@genLen))
    expect_equal(as.vector(colSums(pullSegSiteHaplo(pop,nThreads=1))),rep(1,5))
  }
})

test_that("Split commands scale with ploidy and sample inbred founders correctly", {
  for(inbred in c(FALSE,TRUE)){
    cmd = runMacs2(nInd=4,ploidy=4,inbred=inbred,Ne=100,bp=1000,
                   genLen=0.1,mutRate=1e-8,histNe=200,histGen=80,
                   split=40,returnCommand=TRUE)
    tokens = strsplit(cmd," +")[[1]]
    flag_args = function(flag,n) tokens[match(flag,tokens)+seq_len(n)]
    expect_equal(as.numeric(flag_args("-t",1)),8e-6)
    expect_equal(as.numeric(flag_args("-r",1)),0.08)
    expect_equal(as.numeric(flag_args("-eN",2)),c(0.1,2))
    expect_equal(as.numeric(flag_args("-I",3)),c(2,rep(if(inbred) 2 else 8,2)))
    expect_equal(as.numeric(flag_args("-ej",3)),c(0.050001,2,1))
  }
  set.seed(6014)
  pop = runMacs(2,segSites=5,split=100,nThreads=1)
  expect_equal(pop@nInd,2L)
  expect_equal(pop@nLoci,5L)
  expect_error(runMacs(3,inbred=TRUE,split=100,nThreads=1))
})
