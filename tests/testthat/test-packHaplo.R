context("Haplotype packing")

test_that("REFUSE empty haplotype input for inbred and outbred packing", {
  for(inbred in c(FALSE,TRUE)){
    expect_error(AlphaSimR:::packHaplo(matrix(as.raw(0),0,8),2L,inbred),
                 "packHaplo was given no genotype data",fixed=TRUE)
  }
})

test_that("REFUSE incomplete individuals when packing outbred haplotypes", {
  # Exercise the C++ guard directly: newMapPop rejects these inputs earlier.
  for(ploidy in c(2L,4L)){
    haplo = matrix(as.raw(0),ploidy+1L,8)
    expect_error(AlphaSimR:::packHaplo(haplo,ploidy,FALSE),
                 "Number of rows not a factor of ploidy",fixed=TRUE)
    # For inbred founders each row instead represents a whole individual.
    packed = AlphaSimR:::packHaplo(haplo,ploidy,TRUE)
    expect_equal(dim(packed),c(1L,ploidy,ploidy+1L))
    expect_true(all(packed==as.raw(0)))
  }
})
