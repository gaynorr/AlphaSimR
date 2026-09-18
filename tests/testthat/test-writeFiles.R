context("writeFiles")

# writePlink and writeRecords had no coverage. File writers fail quietly:
# a transposed matrix or an off by one in the header still produces a file,
# and nothing notices until something else tries to read it.
#
# These write into a temporary directory and read the result back, checking
# the shape of what was written and that the genotypes in the file are the
# genotypes in the population.

filePop = function(nInd=8, nChr=2, segSites=20, nQtl=5, nSnp=5, seed=9971,
                   sexes="no"){
  set.seed(seed)
  founderPop = quickHaplo(nInd=nInd, nChr=nChr, segSites=segSites)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  if(sexes != "no"){
    SP$setSexes(sexes)
  }
  SP$restrSegSites(minQtlPerChr=nQtl, minSnpPerChr=nSnp, overlap=FALSE)
  SP$addTraitA(nQtlPerChr=nQtl)
  SP$addSnpChip(nSnpPerChr=nSnp)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)
  return(list(pop=pop, SP=SP, nSnp=nChr*nSnp, nQtl=nChr*nQtl))
}

# A directory that exists only for one test. tempfile() is used rather than a
# name built from a random number, so the name cannot collide and the tests do
# not depend on the state of the random number generator.
withTempDir = function(code){
  dir = tempfile("asr_")
  dir.create(dir, showWarnings=FALSE, recursive=TRUE)
  on.exit(unlink(dir, recursive=TRUE), add=TRUE)
  force(code(dir))
}

test_that("writePlink writes a ped and a map of the right shape", {
  d = filePop(nInd=8, sexes="yes_sys")
  withTempDir(function(dir){
    base = file.path(dir, "plinkTest")
    writePlink(d$pop, baseName=base, simParam=d$SP)

    expect_true(file.exists(paste0(base, ".ped")))
    expect_true(file.exists(paste0(base, ".map")))

    ped = read.table(paste0(base, ".ped"), header=FALSE,
                     stringsAsFactors=FALSE, colClasses="character")
    map = read.table(paste0(base, ".map"), header=FALSE,
                     stringsAsFactors=FALSE)

    # One row per individual, six leading columns then two per marker
    expect_equal(nrow(ped), nInd(d$pop))
    expect_equal(ncol(ped), 6L + 2L*d$nSnp)

    # One row per marker, four columns
    expect_equal(nrow(map), d$nSnp)
    expect_equal(ncol(map), 4L)

    # The leading columns are the pedigree and the phenotype
    expect_equal(ped[,2], d$pop@id)
    expect_true(all(ped[,1] == "1"))
    expect_true(all(ped[,5] %in% c("0","1","2")))
    expect_equal(as.numeric(ped[,6]), unname(c(pheno(d$pop))),
                 tolerance=1e-6)

    # PLINK codes alleles as 1 and 2, so nothing else may appear
    alleles = unlist(ped[, 7:ncol(ped)])
    expect_true(all(alleles %in% c("1","2")))

    # The genotype in the file is the genotype in the population. PLINK
    # writes the two alleles of a marker side by side, and AlphaSimR's 1
    # allele is written as 2, so the dosage is the count of 2s.
    geno = pullSnpGeno(d$pop, simParam=d$SP)
    fromFile = matrix(0L, nrow=nInd(d$pop), ncol=d$nSnp)
    for(j in seq_len(d$nSnp)){
      cols = c(6L + 2L*j - 1L, 6L + 2L*j)
      fromFile[,j] = rowSums(ped[, cols] == "2")
    }
    expect_equal(unname(fromFile), unname(geno))
  })
})

test_that("writePlink's map matches the SNP map", {
  d = filePop()
  withTempDir(function(dir){
    base = file.path(dir, "mapTest")
    writePlink(d$pop, baseName=base, simParam=d$SP)
    map = read.table(paste0(base, ".map"), header=FALSE,
                     stringsAsFactors=FALSE)

    snpMap = getSnpMap(snpChip=1, simParam=d$SP)
    expect_equal(as.character(map[,1]), as.character(snpMap$chr))
    expect_equal(as.character(map[,2]), as.character(snpMap$id))
    expect_equal(as.numeric(map[,3]), snpMap$pos*100, tolerance=1e-8)
    expect_equal(as.numeric(map[,4]), as.numeric(snpMap$site))
  })
})

test_that("writePlink can write QTL instead of a SNP chip", {
  d = filePop()
  withTempDir(function(dir){
    base = file.path(dir, "qtlTest")
    writePlink(d$pop, baseName=base, useQtl=TRUE, simParam=d$SP)
    ped = read.table(paste0(base, ".ped"), header=FALSE,
                     stringsAsFactors=FALSE, colClasses="character")
    expect_equal(ncol(ped), 6L + 2L*d$nQtl)
  })
})

test_that("writePlink refuses a ploidy it cannot represent", {
  set.seed(9981)
  founderPop = quickHaplo(nInd=4, nChr=1, segSites=10, ploidy=4L)
  SP = SimParam$new(founderPop)
  SP$nThreads = 1L
  SP$addTraitA(5)
  SP$addSnpChip(5)
  SP$setVarE(h2=0.5)
  pop = setPheno(newPop(founderPop, simParam=SP), simParam=SP)

  withTempDir(function(dir){
    expect_error(writePlink(pop, baseName=file.path(dir, "x"), simParam=SP),
                 "ploidy")
  })
})

test_that("writeRecords writes the files it says it will", {
  d = filePop(nInd=8)
  withTempDir(function(dir){
    writeRecords(d$pop, dir=dir, simParam=d$SP)

    for(f in c("nMarkers.txt","markerType.txt","info.txt","gv.txt",
               "pheno.txt","genotype.txt")){
      expect_true(file.exists(file.path(dir, f)), info=f)
    }

    expect_equal(scan(file.path(dir,"nMarkers.txt"), integer(), quiet=TRUE),
                 d$nSnp)
    expect_equal(scan(file.path(dir,"markerType.txt"), character(),
                      quiet=TRUE), "SNP_1")

    info = read.table(file.path(dir,"info.txt"), header=TRUE,
                      stringsAsFactors=FALSE, colClasses="character")
    expect_equal(nrow(info), nInd(d$pop))
    expect_equal(info$id, d$pop@id)

    gvFile = read.table(file.path(dir,"gv.txt"), header=FALSE)
    expect_equal(nrow(gvFile), nInd(d$pop))
    expect_equal(gvFile[[1]], unname(c(gv(d$pop))), tolerance=1e-6)

    phenoFile = read.table(file.path(dir,"pheno.txt"), header=FALSE)
    expect_equal(phenoFile[[1]], unname(c(pheno(d$pop))), tolerance=1e-6)

    # Written by Armadillo in raw ascii, so one line per individual and one
    # space separated dosage per marker
    genoFile = read.table(file.path(dir,"genotype.txt"), header=FALSE)
    expect_equal(nrow(genoFile), nInd(d$pop))
    expect_equal(ncol(genoFile), d$nSnp)
    expect_equal(unname(as.matrix(genoFile)),
                 unname(pullSnpGeno(d$pop, simParam=d$SP)))
  })
})

test_that("writeRecords appends a second population to the first", {
  d = filePop(nInd=6)
  withTempDir(function(dir){
    writeRecords(d$pop, dir=dir, simParam=d$SP)
    writeRecords(d$pop, dir=dir, append=TRUE, simParam=d$SP)

    info = read.table(file.path(dir,"info.txt"), header=TRUE,
                      stringsAsFactors=FALSE, colClasses="character")
    expect_equal(nrow(info), 2L*nInd(d$pop))

    genoFile = read.table(file.path(dir,"genotype.txt"), header=FALSE)
    expect_equal(nrow(genoFile), 2L*nInd(d$pop))

    # The second block is the same population written again
    first = as.matrix(genoFile[1:nInd(d$pop),])
    second = as.matrix(genoFile[(nInd(d$pop)+1L):(2L*nInd(d$pop)),])
    expect_equal(unname(first), unname(second))
  })
})

test_that("writeRecords replaces the directory when told not to append", {
  d = filePop(nInd=6)
  withTempDir(function(dir){
    writeRecords(d$pop, dir=dir, simParam=d$SP)
    writeRecords(d$pop, dir=dir, append=FALSE, simParam=d$SP)

    info = read.table(file.path(dir,"info.txt"), header=TRUE,
                      stringsAsFactors=FALSE, colClasses="character")
    expect_equal(nrow(info), nInd(d$pop))
  })
})

test_that("writeRecords refuses to mix marker sets in one directory", {
  d = filePop(nInd=6, nSnp=5)
  withTempDir(function(dir){
    writeRecords(d$pop, dir=dir, simParam=d$SP)

    # The QTL set is a different number of markers and a different type, so
    # adding it to the same directory has to fail rather than corrupt it
    expect_error(writeRecords(d$pop, dir=dir, useQtl=TRUE, append=TRUE,
                              simParam=d$SP))
  })
})

test_that("writeRecords can skip genotypes and can add haplotypes", {
  d = filePop(nInd=6)

  withTempDir(function(dir){
    writeRecords(d$pop, dir=dir, snpChip=0, simParam=d$SP)
    expect_false(file.exists(file.path(dir,"genotype.txt")))
    expect_equal(scan(file.path(dir,"markerType.txt"), character(),
                      quiet=TRUE), "NULL")
  })

  withTempDir(function(dir){
    writeRecords(d$pop, dir=dir, includeHaplo=TRUE, simParam=d$SP)
    expect_true(file.exists(file.path(dir,"haplotype1.txt")))
    expect_true(file.exists(file.path(dir,"haplotype2.txt")))

    h1 = as.matrix(read.table(file.path(dir,"haplotype1.txt"),
                              header=FALSE))
    h2 = as.matrix(read.table(file.path(dir,"haplotype2.txt"),
                              header=FALSE))
    expect_equal(nrow(h1), nInd(d$pop))
    expect_equal(ncol(h1), d$nSnp)

    # The two haplotypes add up to the genotype
    geno = pullSnpGeno(d$pop, simParam=d$SP)
    expect_equal(unname(h1 + h2), unname(geno))

    # And each haplotype only ever holds one allele
    expect_true(all(h1 %in% c(0,1)))
    expect_true(all(h2 %in% c(0,1)))
  })
})
