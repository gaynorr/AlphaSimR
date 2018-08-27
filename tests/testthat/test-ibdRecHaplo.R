context("ibdRecHaplo")

# ---- Data ----

# A simple pedigree
pedigree = matrix(data = 0L, nrow = 10, ncol = 2)
pedigree[ 5, 1:2] = c(1L, 2L)
pedigree[ 6, 1:2] = c(3L, 4L)
pedigree[ 7, 1:2] = c(5L, 6L)
pedigree[ 8, 1:2] = c(5L, 7L)
pedigree[ 9, 1:2] = c(8L, 8L)
pedigree[10, 1:2] = c(8L, 9L)

nLociPerChr = c(300, 300)

# Recombinations (as they are stored in simParam$recHist)
recHist = vector(mode = "list", length = "10")
# 2 chromosomes
recHist[[5]] = recHist[[6]] = recHist[[7]] = recHist[[8]] =
  recHist[[9]] = recHist[[9]] = vector(mode = "list", length = "2")
# 2 gametes of a chromosome
recHist[[5]][[1]] = recHist[[5]][[2]] = 
  recHist[[6]][[1]] = recHist[[6]][[2]] =
  recHist[[7]][[1]] = recHist[[7]][[2]] = 
  recHist[[8]][[1]] = recHist[[8]][[2]] = 
  recHist[[9]][[1]] = recHist[[9]][[2]] = 
  recHist[[10]][[1]] = recHist[[10]][[2]] = vector(mode = "list", length = "2")

ind = 5
chr = 1; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 5, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(1L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(2L, 10L)
recHist[[ind]][[chr]][[par]][3, 1:2] = c(1L, 50L)
recHist[[ind]][[chr]][[par]][4, 1:2] = c(2L, 100L)
recHist[[ind]][[chr]][[par]][5, 1:2] = c(1L, 150L)
chr = 1; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 5, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(1L, 5L)
recHist[[ind]][[chr]][[par]][3, 1:2] = c(2L, 10L)
recHist[[ind]][[chr]][[par]][4, 1:2] = c(1L, 15L)
recHist[[ind]][[chr]][[par]][5, 1:2] = c(2L, 20L)
chr = 2; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 10, ncol = 2)
recHist[[ind]][[chr]][[par]][ 1, 1:2] = c(1L, 1L)
recHist[[ind]][[chr]][[par]][ 2, 1:2] = c(2L, 2L)
recHist[[ind]][[chr]][[par]][ 3, 1:2] = c(1L, 3L)
recHist[[ind]][[chr]][[par]][ 4, 1:2] = c(2L, 4L)
recHist[[ind]][[chr]][[par]][ 5, 1:2] = c(1L, 5L)
recHist[[ind]][[chr]][[par]][ 6, 1:2] = c(2L, 6L)
recHist[[ind]][[chr]][[par]][ 7, 1:2] = c(1L, 7L)
recHist[[ind]][[chr]][[par]][ 8, 1:2] = c(2L, 8L)
recHist[[ind]][[chr]][[par]][ 9, 1:2] = c(1L, 9L)
recHist[[ind]][[chr]][[par]][10, 1:2] = c(2L, 10L)
chr = 2; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 10, ncol = 2)
recHist[[ind]][[chr]][[par]][ 1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][ 2, 1:2] = c(1L, 3L)
recHist[[ind]][[chr]][[par]][ 3, 1:2] = c(2L, 5L)
recHist[[ind]][[chr]][[par]][ 4, 1:2] = c(1L, 7L)
recHist[[ind]][[chr]][[par]][ 5, 1:2] = c(2L, 9L)
recHist[[ind]][[chr]][[par]][ 6, 1:2] = c(1L, 11L)
recHist[[ind]][[chr]][[par]][ 7, 1:2] = c(2L, 13L)
recHist[[ind]][[chr]][[par]][ 8, 1:2] = c(1L, 15L)
recHist[[ind]][[chr]][[par]][ 9, 1:2] = c(2L, 17L)
recHist[[ind]][[chr]][[par]][10, 1:2] = c(1L, 19L)

ind = 6
chr = 1; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 4, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(1L, 150L)
recHist[[ind]][[chr]][[par]][3, 1:2] = c(2L, 200L)
recHist[[ind]][[chr]][[par]][4, 1:2] = c(1L, 250L)
chr = 1; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 5, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(1L, 9L)
recHist[[ind]][[chr]][[par]][3, 1:2] = c(2L, 11L)
recHist[[ind]][[chr]][[par]][4, 1:2] = c(1L, 50L)
recHist[[ind]][[chr]][[par]][5, 1:2] = c(2L, 100L)
chr = 2; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 11, ncol = 2)
recHist[[ind]][[chr]][[par]][ 1, 1:2] = c(1L, 1L)
recHist[[ind]][[chr]][[par]][ 2, 1:2] = c(2L, 10L)
recHist[[ind]][[chr]][[par]][ 3, 1:2] = c(1L, 20L)
recHist[[ind]][[chr]][[par]][ 4, 1:2] = c(2L, 30L)
recHist[[ind]][[chr]][[par]][ 5, 1:2] = c(1L, 40L)
recHist[[ind]][[chr]][[par]][ 6, 1:2] = c(2L, 50L)
recHist[[ind]][[chr]][[par]][ 7, 1:2] = c(1L, 60L)
recHist[[ind]][[chr]][[par]][ 8, 1:2] = c(2L, 70L)
recHist[[ind]][[chr]][[par]][ 9, 1:2] = c(1L, 80L)
recHist[[ind]][[chr]][[par]][10, 1:2] = c(2L, 90L)
recHist[[ind]][[chr]][[par]][11, 1:2] = c(1L, 100L)
chr = 2; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 6, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(1L, 50L)
recHist[[ind]][[chr]][[par]][3, 1:2] = c(2L, 100L)
recHist[[ind]][[chr]][[par]][4, 1:2] = c(1L, 150L)
recHist[[ind]][[chr]][[par]][5, 1:2] = c(2L, 200L)
recHist[[ind]][[chr]][[par]][6, 1:2] = c(1L, 250L)

ind = 7
chr = 1; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 6, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(1L, 6L)
recHist[[ind]][[chr]][[par]][3, 1:2] = c(2L, 10L)
recHist[[ind]][[chr]][[par]][4, 1:2] = c(1L, 45L)
recHist[[ind]][[chr]][[par]][5, 1:2] = c(2L, 110L)
recHist[[ind]][[chr]][[par]][6, 1:2] = c(1L, 140L)
chr = 1; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
chr = 2; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 6, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(1L, 6L)
recHist[[ind]][[chr]][[par]][3, 1:2] = c(2L, 10L)
recHist[[ind]][[chr]][[par]][4, 1:2] = c(1L, 45L)
recHist[[ind]][[chr]][[par]][5, 1:2] = c(2L, 110L)
recHist[[ind]][[chr]][[par]][6, 1:2] = c(1L, 140L)
chr = 2; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 7, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(1L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(2L, 50L)
recHist[[ind]][[chr]][[par]][3, 1:2] = c(1L, 75L)
recHist[[ind]][[chr]][[par]][4, 1:2] = c(2L, 99L)
recHist[[ind]][[chr]][[par]][5, 1:2] = c(1L, 101L)
recHist[[ind]][[chr]][[par]][6, 1:2] = c(2L, 149L)
recHist[[ind]][[chr]][[par]][7, 1:2] = c(1L, 151L)

ind = 8
chr = 1; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
chr = 1; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
chr = 2; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
chr = 2; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)

ind = 9
chr = 1; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(1L, 1L)
chr = 1; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
chr = 2; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(1L, 1L)
chr = 2; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)

ind = 10
chr = 1; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 2, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(1L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(2L, 50L)
chr = 1; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 2, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(1L, 100L)
chr = 2; par = 1
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 2, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(1L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(2L, 50L)
chr = 2; par = 2
recHist[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 2, ncol = 2)
recHist[[ind]][[chr]][[par]][1, 1:2] = c(2L, 1L)
recHist[[ind]][[chr]][[par]][2, 1:2] = c(1L, 100L)

# ---- getIbdRecHist ----

test_that("getIbdRecHist converts recHist & pedigree to ibdRecHist correctly", {
  # IBD recombinations (since the base generation)
  expect = recHist
  
  ind = 1
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = c(1L, 1L), nrow = 1, ncol = 2)
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = c(2L, 1L), nrow = 1, ncol = 2)
  expect[[ind]][[2]] = expect[[ind]][[1]]
  
  ind = 2
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = c(3L, 1L), nrow = 1, ncol = 2)
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = c(4L, 1L), nrow = 1, ncol = 2)
  expect[[ind]][[2]] = expect[[ind]][[1]]
  
  ind = 3
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = c(5L, 1L), nrow = 1, ncol = 2)
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = c(6L, 1L), nrow = 1, ncol = 2)
  expect[[ind]][[2]] = expect[[ind]][[1]]
  
  ind = 4
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = c(7L, 1L), nrow = 1, ncol = 2)
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = c(8L, 1L), nrow = 1, ncol = 2)
  expect[[ind]][[2]] = expect[[ind]][[1]]
  
  ind = 5
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[ind]][[chr]][[par]][1, 1:2] = c(1L, 1L)
  expect[[ind]][[chr]][[par]][2, 1:2] = c(2L, 10L)
  expect[[ind]][[chr]][[par]][3, 1:2] = c(1L, 50L)
  expect[[ind]][[chr]][[par]][4, 1:2] = c(2L, 100L)
  expect[[ind]][[chr]][[par]][5, 1:2] = c(1L, 150L)
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[ind]][[chr]][[par]][1, 1:2] = c(4L, 1L)
  expect[[ind]][[chr]][[par]][2, 1:2] = c(3L, 5L)
  expect[[ind]][[chr]][[par]][3, 1:2] = c(4L, 10L)
  expect[[ind]][[chr]][[par]][4, 1:2] = c(3L, 15L)
  expect[[ind]][[chr]][[par]][5, 1:2] = c(4L, 20L)
  chr = 2; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 10, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(1L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(2L, 2L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(1L, 3L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(2L, 4L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(1L, 5L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(2L, 6L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(1L, 7L)
  expect[[ind]][[chr]][[par]][ 8, 1:2] = c(2L, 8L)
  expect[[ind]][[chr]][[par]][ 9, 1:2] = c(1L, 9L)
  expect[[ind]][[chr]][[par]][10, 1:2] = c(2L, 10L)
  chr = 2; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 10, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(4L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(3L, 3L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(4L, 5L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(3L, 7L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(4L, 9L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(3L, 11L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(4L, 13L)
  expect[[ind]][[chr]][[par]][ 8, 1:2] = c(3L, 15L)
  expect[[ind]][[chr]][[par]][ 9, 1:2] = c(4L, 17L)
  expect[[ind]][[chr]][[par]][10, 1:2] = c(3L, 19L)
  
  ind = 6
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 4, ncol = 2)
  expect[[ind]][[chr]][[par]][1, 1:2] = c(6L, 1L)
  expect[[ind]][[chr]][[par]][2, 1:2] = c(5L, 150L)
  expect[[ind]][[chr]][[par]][3, 1:2] = c(6L, 200L)
  expect[[ind]][[chr]][[par]][4, 1:2] = c(5L, 250L)
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[ind]][[chr]][[par]][1, 1:2] = c(8L, 1L)
  expect[[ind]][[chr]][[par]][2, 1:2] = c(7L, 9L)
  expect[[ind]][[chr]][[par]][3, 1:2] = c(8L, 11L)
  expect[[ind]][[chr]][[par]][4, 1:2] = c(7L, 50L)
  expect[[ind]][[chr]][[par]][5, 1:2] = c(8L, 100L)
  chr = 2; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 11, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(5L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(6L, 10L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(5L, 20L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(6L, 30L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(5L, 40L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(6L, 50L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(5L, 60L)
  expect[[ind]][[chr]][[par]][ 8, 1:2] = c(6L, 70L)
  expect[[ind]][[chr]][[par]][ 9, 1:2] = c(5L, 80L)
  expect[[ind]][[chr]][[par]][10, 1:2] = c(6L, 90L)
  expect[[ind]][[chr]][[par]][11, 1:2] = c(5L, 100L)
  chr = 2; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 6, ncol = 2)
  expect[[ind]][[chr]][[par]][1, 1:2] = c(8L, 1L)
  expect[[ind]][[chr]][[par]][2, 1:2] = c(7L, 50L)
  expect[[ind]][[chr]][[par]][3, 1:2] = c(8L, 100L)
  expect[[ind]][[chr]][[par]][4, 1:2] = c(7L, 150L)
  expect[[ind]][[chr]][[par]][5, 1:2] = c(8L, 200L)
  expect[[ind]][[chr]][[par]][6, 1:2] = c(7L, 250L)
  
  ind = 7
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 12, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(4L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(3L, 5L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(1L, 6L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(4L, 10L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(3L, 15L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(4L, 20L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(2L, 45L)
  expect[[ind]][[chr]][[par]][ 8, 1:2] = c(1L, 50L)
  expect[[ind]][[chr]][[par]][ 9, 1:2] = c(2L, 100L)
  expect[[ind]][[chr]][[par]][10, 1:2] = c(4L, 110L)
  expect[[ind]][[chr]][[par]][11, 1:2] = c(2L, 140L)
  expect[[ind]][[chr]][[par]][12, 1:2] = c(1L, 150L)
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[ind]][[chr]][[par]][1, 1:2] = c(8L, 1L)
  expect[[ind]][[chr]][[par]][2, 1:2] = c(7L, 9L)
  expect[[ind]][[chr]][[par]][3, 1:2] = c(8L, 11L)
  expect[[ind]][[chr]][[par]][4, 1:2] = c(7L, 50L)
  expect[[ind]][[chr]][[par]][5, 1:2] = c(8L, 100L)
  chr = 2; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 16, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(4L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(3L, 3L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(4L, 5L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(2L, 6L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(1L, 7L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(2L, 8L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(1L, 9L)
  expect[[ind]][[chr]][[par]][ 8, 1:2] = c(4L, 10L)
  expect[[ind]][[chr]][[par]][ 9, 1:2] = c(3L, 11L)
  expect[[ind]][[chr]][[par]][10, 1:2] = c(4L, 13L)
  expect[[ind]][[chr]][[par]][11, 1:2] = c(3L, 15L)
  expect[[ind]][[chr]][[par]][12, 1:2] = c(4L, 17L)  
  expect[[ind]][[chr]][[par]][13, 1:2] = c(3L, 19L)  
  expect[[ind]][[chr]][[par]][14, 1:2] = c(2L, 45L)  
  expect[[ind]][[chr]][[par]][15, 1:2] = c(3L, 110L)  
  expect[[ind]][[chr]][[par]][16, 1:2] = c(2L, 140L)  
  chr = 2; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 15, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(5L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(6L, 10L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(5L, 20L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(6L, 30L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(5L, 40L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(7L, 50L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(6L, 75L)
  expect[[ind]][[chr]][[par]][ 8, 1:2] = c(5L, 80L)
  expect[[ind]][[chr]][[par]][ 9, 1:2] = c(6L, 90L)
  expect[[ind]][[chr]][[par]][10, 1:2] = c(7L, 99L)
  expect[[ind]][[chr]][[par]][11, 1:2] = c(8L, 100L)
  expect[[ind]][[chr]][[par]][12, 1:2] = c(5L, 101L)  
  expect[[ind]][[chr]][[par]][13, 1:2] = c(8L, 149L)  
  expect[[ind]][[chr]][[par]][14, 1:2] = c(7L, 150L)  
  expect[[ind]][[chr]][[par]][15, 1:2] = c(5L, 151L)
  
  ind = 8
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = expect[[5]][[chr]][[2]]
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = expect[[7]][[chr]][[2]]
  chr = 2; par = 1
  expect[[ind]][[chr]][[par]] = expect[[5]][[chr]][[2]]
  chr = 2; par = 2
  expect[[ind]][[chr]][[par]] = expect[[7]][[chr]][[2]]
  
  ind = 9
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = expect[[8]][[chr]][[1]]
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = expect[[8]][[chr]][[2]]
  chr = 2; par = 1
  expect[[ind]][[chr]][[par]] = expect[[8]][[chr]][[1]]
  chr = 2; par = 2
  expect[[ind]][[chr]][[par]] = expect[[8]][[chr]][[2]]
  
  ind = 10
  chr = 1; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 7, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(4L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(3L, 5L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(4L, 10L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(3L, 15L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(4L, 20L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(7L, 50L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(8L, 100L)
  chr = 1; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[ind]][[chr]][[par]][1, 1:2] = c(8L, 1L)
  expect[[ind]][[chr]][[par]][2, 1:2] = c(7L, 9L)
  expect[[ind]][[chr]][[par]][3, 1:2] = c(8L, 11L)
  expect[[ind]][[chr]][[par]][4, 1:2] = c(7L, 50L)
  expect[[ind]][[chr]][[par]][5, 1:2] = c(4L, 100L)
  chr = 2; par = 1
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 20, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(4L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(3L, 3L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(4L, 5L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(3L, 7L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(4L, 9L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(3L, 11L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(4L, 13L)
  expect[[ind]][[chr]][[par]][ 8, 1:2] = c(3L, 15L)
  expect[[ind]][[chr]][[par]][ 9, 1:2] = c(4L, 17L)
  expect[[ind]][[chr]][[par]][10, 1:2] = c(3L, 19L)
  expect[[ind]][[chr]][[par]][11, 1:2] = c(7L, 50L)
  expect[[ind]][[chr]][[par]][12, 1:2] = c(6L, 75L)  
  expect[[ind]][[chr]][[par]][13, 1:2] = c(5L, 80L)  
  expect[[ind]][[chr]][[par]][14, 1:2] = c(6L, 90L)  
  expect[[ind]][[chr]][[par]][15, 1:2] = c(7L, 99L)  
  expect[[ind]][[chr]][[par]][16, 1:2] = c(8L, 100L)  
  expect[[ind]][[chr]][[par]][17, 1:2] = c(5L, 101L)  
  expect[[ind]][[chr]][[par]][18, 1:2] = c(8L, 149L)  
  expect[[ind]][[chr]][[par]][19, 1:2] = c(7L, 150L)  
  expect[[ind]][[chr]][[par]][20, 1:2] = c(5L, 151L)  
  chr = 2; par = 2
  expect[[ind]][[chr]][[par]] = matrix(data = 0L, nrow = 11, ncol = 2)
  expect[[ind]][[chr]][[par]][ 1, 1:2] = c(5L, 1L)
  expect[[ind]][[chr]][[par]][ 2, 1:2] = c(6L, 10L)
  expect[[ind]][[chr]][[par]][ 3, 1:2] = c(5L, 20L)
  expect[[ind]][[chr]][[par]][ 4, 1:2] = c(6L, 30L)
  expect[[ind]][[chr]][[par]][ 5, 1:2] = c(5L, 40L)
  expect[[ind]][[chr]][[par]][ 6, 1:2] = c(7L, 50L)
  expect[[ind]][[chr]][[par]][ 7, 1:2] = c(6L, 75L)
  expect[[ind]][[chr]][[par]][ 8, 1:2] = c(5L, 80L)
  expect[[ind]][[chr]][[par]][ 9, 1:2] = c(6L, 90L)
  expect[[ind]][[chr]][[par]][10, 1:2] = c(7L, 99L)
  expect[[ind]][[chr]][[par]][11, 1:2] = c(3L, 100L)
  
  result = AlphaSimR:::getIbdRecHist(recHist     = recHist,
                                     pedigree    = pedigree,
                                     nLociPerChr = nLociPerChr)$ibdRecHist
  
  for (ind in 1:10) {
    for (chr in 1:2) {
      for (par in 1:2) {
        expect_identical(object = result[[ind]][[chr]][[par]],
                         expected = expect[[ind]][[chr]][[par]],
                         info = paste0("ind ", ind, " chr ", chr, " par ", par))
      }
    }
  }
})

# ---- pullIbdHaplo2 ----

test_that("pullIbdHaplo2 gets correct IBD haplotypes", {
  ibdRecHist = AlphaSimR:::getIbdRecHist(recHist     = recHist,
                                         pedigree    = pedigree,
                                         nLociPerChr = nLociPerChr)$ibdRecHist
  output = AlphaSimR:::getIbdHaplo(ibdRecHist  = ibdRecHist,
                                   individuals = 1L:nrow(pedigree),
                                   nLociPerChr = nLociPerChr)
  
  expect = matrix(data = 0L, nrow = 18, ncol = 600)
  
  # ind 1
  gam = 1
  expect[gam, ] = 1L
  expect_identical(object = output[gam, ],
                   expected = expect[gam, ],
                   info = paste0("ind 1, chr 1 & 2, gamete 1"))
  gam = 2
  expect[gam, ] = 2L
  expect_identical(object = output[gam, ],
                   expected = expect[gam, ],
                   info = paste0("ind 1, chr 1 & 2, gamete 2"))
  
  # ind 2
  gam = 3
  expect[gam, ] = 3L
  expect_identical(object = output[gam, ],
                   expected = expect[gam, ],
                   info = paste0("ind 2, chr 1 & 2, gamete 1"))
  gam = 4
  expect[gam, ] = 4L
  expect_identical(object = output[gam, ],
                   expected = expect[gam, ],
                   info = paste0("ind 2, chr 1 & 2, gamete 2"))
  
  # ind 3
  gam = 5
  expect[gam, ] = 5L
  expect_identical(object = output[gam, ],
                   expected = expect[gam, ],
                   info = paste0("ind 3, chr 1 & 2, gamete 1"))
  gam = 6
  expect[gam, ] = 6L
  expect_identical(object = output[gam, ],
                   expected = expect[gam, ],
                   info = paste0("ind 3, chr 1 & 2, gamete 2"))
  
  # ind 4
  gam = 7
  expect[gam, ] = 7L
  expect_identical(object = output[gam, ],
                   expected = expect[gam, ],
                   info = paste0("ind 4, chr 1 & 2, gamete 1"))
  gam = 8
  expect[gam, ] = 8L
  expect_identical(object = output[gam, ],
                   expected = expect[gam, ],
                   info = paste0("ind 4, chr 1 & 2, gamete 2"))
  
  # ind 5 - chr 1, maternal
  gam = 9
  expect[gam,   1:009] = 1L
  expect[gam,  10:049] = 2L
  expect[gam,  50:099] = 1L
  expect[gam, 100:149] = 2L
  expect[gam, 150:300] = 1L
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 5, chr 1, gamete 1"))
  # ind 5 - chr 1, paternal
  gam = 10
  expect[gam,  1:004] = 4L
  expect[gam,  5:009] = 3L
  expect[gam, 10:014] = 4L
  expect[gam, 15:019] = 3L
  expect[gam, 20:300] = 4L
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 5, chr 1, gamete 2"))
  # ind 5 - chr 2, maternal
  gam = 9
  expect[gam, 300 + 001:001] = 1L
  expect[gam, 300 + 002:002] = 2L
  expect[gam, 300 + 003:003] = 1L
  expect[gam, 300 + 004:004] = 2L
  expect[gam, 300 + 005:005] = 1L
  expect[gam, 300 + 006:006] = 2L
  expect[gam, 300 + 007:007] = 1L
  expect[gam, 300 + 008:008] = 2L
  expect[gam, 300 + 009:009] = 1L
  expect[gam, 300 + 010:300] = 2L
  expect_identical(object = output[gam, 301:600],
                   expected = expect[gam, 301:600],
                   info = paste0("ind 5, chr 2, gamete 1"))
  # ind 5 - chr 2, paternal
  gam = 10
  expect[gam, 300 + 001:002] = 4L
  expect[gam, 300 + 003:004] = 3L
  expect[gam, 300 + 005:006] = 4L
  expect[gam, 300 + 007:008] = 3L
  expect[gam, 300 + 009:010] = 4L
  expect[gam, 300 + 011:012] = 3L
  expect[gam, 300 + 013:014] = 4L
  expect[gam, 300 + 015:016] = 3L
  expect[gam, 300 + 017:018] = 4L
  expect[gam, 300 + 019:300] = 3L
  expect_identical(object = output[gam, 301:600],
                   expected = expect[gam, 301:600],
                   info = paste0("ind 5, chr 2, gamete 2"))
  
  # ind 6 - chr 1, maternal
  gam = 11
  expect[gam,   1:149] = 6L
  expect[gam, 150:199] = 5L
  expect[gam, 200:249] = 6L
  expect[gam, 250:300] = 5L
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 6, chr 1, gamete 1"))
  # ind 6 - chr 1, paternal
  gam = 12
  expect[gam,   1:008] = 8L
  expect[gam,   9:010] = 7L
  expect[gam,  11:049] = 8L
  expect[gam,  50:099] = 7L
  expect[gam, 100:300] = 8L
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 6, chr 1, gamete 2"))
  # ind 6 - chr 2, maternal
  gam = 11
  expect[gam, 300 + 001:009] = 5L
  expect[gam, 300 + 010:019] = 6L
  expect[gam, 300 + 020:029] = 5L
  expect[gam, 300 + 030:039] = 6L
  expect[gam, 300 + 040:049] = 5L
  expect[gam, 300 + 050:059] = 6L
  expect[gam, 300 + 060:069] = 5L
  expect[gam, 300 + 070:079] = 6L
  expect[gam, 300 + 080:089] = 5L
  expect[gam, 300 + 090:099] = 6L
  expect[gam, 300 + 100:300] = 5L
  expect_identical(object = output[gam, 301:600],
                   expected = expect[gam, 301:600],
                   info = paste0("ind 6, chr 2, gamete 1"))
  # ind 6 - chr 2, paternal
  gam = 12
  expect[gam, 300 + 001:049] = 8L
  expect[gam, 300 + 050:099] = 7L
  expect[gam, 300 + 100:149] = 8L
  expect[gam, 300 + 150:199] = 7L
  expect[gam, 300 + 200:249] = 8L
  expect[gam, 300 + 250:300] = 7L
  expect_identical(object = output[gam, 301:600],
                   expected = expect[gam, 301:600],
                   info = paste0("ind 6, chr 2, gamete 2"))
  
  # ind 7 - chr 1, maternal
  gam = 13
  expect[gam,   1:004] = 4L
  expect[gam,   5:005] = 3L
  expect[gam,   6:009] = 1L
  expect[gam,  10:014] = 4L
  expect[gam,  15:019] = 3L
  expect[gam,  20:044] = 4L
  expect[gam,  45:049] = 2L
  expect[gam,  50:099] = 1L
  expect[gam, 100:109] = 2L
  expect[gam, 110:139] = 4L
  expect[gam, 140:149] = 2L
  expect[gam, 150:300] = 1L
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 7, chr 1, gamete 1"))
  # ind 7 - chr 1, paternal
  gam = 14
  expect[gam,   1:008] = 8L
  expect[gam,   9:010] = 7L
  expect[gam,  11:049] = 8L
  expect[gam,  50:099] = 7L
  expect[gam, 100:300] = 8L
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 7, chr 1, gamete 2"))
  # ind 7 - chr 2, maternal
  gam = 13
  expect[gam, 300 + 001:002] = 4L
  expect[gam, 300 + 003:004] = 3L
  expect[gam, 300 + 005:005] = 4L
  expect[gam, 300 + 006:006] = 2L
  expect[gam, 300 + 007:007] = 1L
  expect[gam, 300 + 008:008] = 2L
  expect[gam, 300 + 009:009] = 1L
  expect[gam, 300 + 010:010] = 4L
  expect[gam, 300 + 011:012] = 3L
  expect[gam, 300 + 013:014] = 4L
  expect[gam, 300 + 015:016] = 3L
  expect[gam, 300 + 017:018] = 4L
  expect[gam, 300 + 019:044] = 3L
  expect[gam, 300 + 045:109] = 2L
  expect[gam, 300 + 110:139] = 3L
  expect[gam, 300 + 140:300] = 2L
  expect_identical(object = output[gam, 301:600],
                   expected = expect[gam, 301:600],
                   info = paste0("ind 7, chr 2, gamete 1"))
  # ind 7 - chr 2, paternal
  gam = 14
  expect[gam, 300 + 001:009] = 5L
  expect[gam, 300 + 010:019] = 6L
  expect[gam, 300 + 020:029] = 5L
  expect[gam, 300 + 030:039] = 6L
  expect[gam, 300 + 040:049] = 5L
  expect[gam, 300 + 050:074] = 7L
  expect[gam, 300 + 075:079] = 6L
  expect[gam, 300 + 080:089] = 5L
  expect[gam, 300 + 090:098] = 6L
  expect[gam, 300 + 099:099] = 7L
  expect[gam, 300 + 100:100] = 8L
  expect[gam, 300 + 101:148] = 5L
  expect[gam, 300 + 149:149] = 8L
  expect[gam, 300 + 150:150] = 7L
  expect[gam, 300 + 151:300] = 5L
  expect_identical(object = output[gam, 301:600],
                   expected = expect[gam, 301:600],
                   info = paste0("ind 7, chr 2, gamete 2"))
  
  # ind 8 - chr 1, maternal
  gam = 15
  expect[gam, 1:300] = expect[10, 1:300]
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 8, chr 1, gamete 1"))
  # ind 8 - chr 1, paternal
  gam = 16
  expect[gam, 1:300] = expect[14, 1:300]
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 8, chr 1, gamete 2"))
  # ind 8 - chr 2, maternal
  gam = 15
  expect[gam, 300 + 1:300] = expect[10, 300 + 1:300]
  expect_identical(object = output[gam, 300 + 1:300],
                   expected = expect[gam, 300 + 1:300],
                   info = paste0("ind 8, chr 2, gamete 1"))
  # ind 8 - chr 2, paternal
  gam = 16
  expect[gam, 300 + 1:300] = expect[14, 300 + 1:300]
  expect_identical(object = output[gam, 300 + 1:300],
                   expected = expect[gam, 300 + 1:300],
                   info = paste0("ind 8, chr 2, gamete 2"))
  
  # ind 9 - chr 1, maternal
  gam = 17
  expect[gam, 1:300] = expect[15, 1:300]
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 9, chr 1, gamete 1"))
  # ind 9 - chr 1, paternal
  gam = 18
  expect[gam, 1:300] = expect[16, 1:300]
  expect_identical(object = output[gam, 1:300],
                   expected = expect[gam, 1:300],
                   info = paste0("ind 9, chr 1, gamete 2"))
  # ind 8 - chr 2, maternal
  gam = 17
  expect[gam, 300 + 1:300] = expect[15, 300 + 1:300]
  expect_identical(object = output[gam, 300 + 1:300],
                   expected = expect[gam, 300 + 1:300],
                   info = paste0("ind 9, chr 2, gamete 1"))
  # ind 8 - chr 2, paternal
  gam = 18
  expect[gam, 300 + 1:300] = expect[16, 300 + 1:300]
  expect_identical(object = output[gam, 300 + 1:300],
                   expected = expect[gam, 300 + 1:300],
                   info = paste0("ind 9, chr 2, gamete 2"))
  
  # Test subsetting
  output = AlphaSimR:::getIbdHaplo(ibdRecHist  = ibdRecHist,
                                   individuals = 7L,
                                   nLociPerChr = nLociPerChr)
  # ind 7 - chr 1, maternal
  expect_identical(object = output[1, 1:300],
                   expected = expect[13, 1:300],
                   info = paste0("ind 7, chr 1, gamete 1, from subsetted call"))
  # ind 7 - chr 1, paternal
  expect_identical(object = output[2, 1:300],
                   expected = expect[14, 1:300],
                   info = paste0("ind 7, chr 1, gamete 2, from subsetted call"))
  # ind 7 - chr 2, maternal
  expect_identical(object = output[1, 301:600],
                   expected = expect[13, 301:600],
                   info = paste0("ind 7, chr 2, gamete 1, from subsetted call"))
  # ind 7 - chr 2, paternal
  expect_identical(object = output[2, 301:600],
                   expected = expect[14, 301:600],
                   info = paste0("ind 7, chr 2, gamete 2, from subsetted call"))
})
