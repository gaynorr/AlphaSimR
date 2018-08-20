context("ibdRecHaplo")

# ---- Data ----

# A simple pedigree
pedigree = matrix(data = 0L, nrow = 7, ncol = 2)
pedigree[5, 1:2] = c(1L, 2L)
pedigree[6, 1:2] = c(3L, 4L)
pedigree[7, 1:2] = c(5L, 6L)

nLociPerChr = c(300, 300)

# Recombinations (as they are stored in simParam$recHist)
recHist = vector(mode = "list", length = "7")
# 2 chromosomes
recHist[[5]] = recHist[[6]] = recHist[[7]] = vector(mode = "list", length = "2")
# 2 gametes of a chromosome
recHist[[5]][[1]] = recHist[[5]][[2]] = 
  recHist[[6]][[1]] = recHist[[6]][[2]] =
  recHist[[7]][[1]] = recHist[[7]][[2]] = vector(mode = "list", length = "2")
# ind 5
recHist[[5]][[1]][[1]] = matrix(data = 0L, nrow = 5, ncol = 2)
recHist[[5]][[1]][[1]][1, 1:2] = c(1L, 1L)
recHist[[5]][[1]][[1]][2, 1:2] = c(2L, 10L)
recHist[[5]][[1]][[1]][3, 1:2] = c(1L, 50L)
recHist[[5]][[1]][[1]][4, 1:2] = c(2L, 100L)
recHist[[5]][[1]][[1]][5, 1:2] = c(1L, 150L)
recHist[[5]][[1]][[2]] = matrix(data = 0L, nrow = 5, ncol = 2)
recHist[[5]][[1]][[2]][1, 1:2] = c(2L, 1L)
recHist[[5]][[1]][[2]][2, 1:2] = c(1L, 5L)
recHist[[5]][[1]][[2]][3, 1:2] = c(2L, 10L)
recHist[[5]][[1]][[2]][4, 1:2] = c(1L, 15L)
recHist[[5]][[1]][[2]][5, 1:2] = c(2L, 20L)
recHist[[5]][[2]][[1]] = matrix(data = 0L, nrow = 10, ncol = 2)
recHist[[5]][[2]][[1]][ 1, 1:2] = c(1L, 1L)
recHist[[5]][[2]][[1]][ 2, 1:2] = c(2L, 2L)
recHist[[5]][[2]][[1]][ 3, 1:2] = c(1L, 3L)
recHist[[5]][[2]][[1]][ 4, 1:2] = c(2L, 4L)
recHist[[5]][[2]][[1]][ 5, 1:2] = c(1L, 5L)
recHist[[5]][[2]][[1]][ 6, 1:2] = c(2L, 6L)
recHist[[5]][[2]][[1]][ 7, 1:2] = c(1L, 7L)
recHist[[5]][[2]][[1]][ 8, 1:2] = c(2L, 8L)
recHist[[5]][[2]][[1]][ 9, 1:2] = c(1L, 9L)
recHist[[5]][[2]][[1]][10, 1:2] = c(2L, 10L)
recHist[[5]][[2]][[2]] = matrix(data = 0L, nrow = 10, ncol = 2)
recHist[[5]][[2]][[2]][ 1, 1:2] = c(2L, 1L)
recHist[[5]][[2]][[2]][ 2, 1:2] = c(1L, 3L)
recHist[[5]][[2]][[2]][ 3, 1:2] = c(2L, 5L)
recHist[[5]][[2]][[2]][ 4, 1:2] = c(1L, 7L)
recHist[[5]][[2]][[2]][ 5, 1:2] = c(2L, 9L)
recHist[[5]][[2]][[2]][ 6, 1:2] = c(1L, 11L)
recHist[[5]][[2]][[2]][ 7, 1:2] = c(2L, 13L)
recHist[[5]][[2]][[2]][ 8, 1:2] = c(1L, 15L)
recHist[[5]][[2]][[2]][ 9, 1:2] = c(2L, 17L)
recHist[[5]][[2]][[2]][10, 1:2] = c(1L, 19L)
# ind 6
recHist[[6]][[1]][[1]] = matrix(data = 0L, nrow = 4, ncol = 2)
recHist[[6]][[1]][[1]][1, 1:2] = c(2L, 1L)
recHist[[6]][[1]][[1]][2, 1:2] = c(1L, 150L)
recHist[[6]][[1]][[1]][3, 1:2] = c(2L, 200L)
recHist[[6]][[1]][[1]][4, 1:2] = c(1L, 250L)
recHist[[6]][[1]][[2]] = matrix(data = 0L, nrow = 5, ncol = 2)
recHist[[6]][[1]][[2]][1, 1:2] = c(2L, 1L)
recHist[[6]][[1]][[2]][2, 1:2] = c(1L, 9L)
recHist[[6]][[1]][[2]][3, 1:2] = c(2L, 11L)
recHist[[6]][[1]][[2]][4, 1:2] = c(1L, 50L)
recHist[[6]][[1]][[2]][5, 1:2] = c(2L, 100L)
recHist[[6]][[2]][[1]] = matrix(data = 0L, nrow = 11, ncol = 2)
recHist[[6]][[2]][[1]][ 1, 1:2] = c(1L, 1L)
recHist[[6]][[2]][[1]][ 2, 1:2] = c(2L, 10L)
recHist[[6]][[2]][[1]][ 3, 1:2] = c(1L, 20L)
recHist[[6]][[2]][[1]][ 4, 1:2] = c(2L, 30L)
recHist[[6]][[2]][[1]][ 5, 1:2] = c(1L, 40L)
recHist[[6]][[2]][[1]][ 6, 1:2] = c(2L, 50L)
recHist[[6]][[2]][[1]][ 7, 1:2] = c(1L, 60L)
recHist[[6]][[2]][[1]][ 8, 1:2] = c(2L, 70L)
recHist[[6]][[2]][[1]][ 9, 1:2] = c(1L, 80L)
recHist[[6]][[2]][[1]][10, 1:2] = c(2L, 90L)
recHist[[6]][[2]][[1]][11, 1:2] = c(1L, 100L)
recHist[[6]][[2]][[2]] = matrix(data = 0L, nrow = 6, ncol = 2)
recHist[[6]][[2]][[2]][1, 1:2] = c(2L, 1L)
recHist[[6]][[2]][[2]][2, 1:2] = c(1L, 50L)
recHist[[6]][[2]][[2]][3, 1:2] = c(2L, 100L)
recHist[[6]][[2]][[2]][4, 1:2] = c(1L, 150L)
recHist[[6]][[2]][[2]][5, 1:2] = c(2L, 200L)
recHist[[6]][[2]][[2]][6, 1:2] = c(1L, 250L)
# ind 7
recHist[[7]][[1]][[1]] = matrix(data = 0L, nrow = 6, ncol = 2)
recHist[[7]][[1]][[1]][1, 1:2] = c(2L, 1L)
recHist[[7]][[1]][[1]][2, 1:2] = c(1L, 6L)
recHist[[7]][[1]][[1]][3, 1:2] = c(2L, 10L)
recHist[[7]][[1]][[1]][4, 1:2] = c(1L, 45L)
recHist[[7]][[1]][[1]][5, 1:2] = c(2L, 110L)
recHist[[7]][[1]][[1]][6, 1:2] = c(1L, 140L)
recHist[[7]][[1]][[2]] = matrix(data = 0L, nrow = 1, ncol = 2)
recHist[[7]][[1]][[2]][1, 1:2] = c(2L, 1L)
recHist[[7]][[2]][[1]] = matrix(data = 0L, nrow = 6, ncol = 2)
recHist[[7]][[2]][[1]][1, 1:2] = c(2L, 1L)
recHist[[7]][[2]][[1]][2, 1:2] = c(1L, 6L)
recHist[[7]][[2]][[1]][3, 1:2] = c(2L, 10L)
recHist[[7]][[2]][[1]][4, 1:2] = c(1L, 45L)
recHist[[7]][[2]][[1]][5, 1:2] = c(2L, 110L)
recHist[[7]][[2]][[1]][6, 1:2] = c(1L, 140L)
recHist[[7]][[2]][[2]] = matrix(data = 0L, nrow = 7, ncol = 2)
recHist[[7]][[2]][[2]][1, 1:2] = c(1L, 1L)
recHist[[7]][[2]][[2]][2, 1:2] = c(2L, 50L)
recHist[[7]][[2]][[2]][3, 1:2] = c(1L, 75L)
recHist[[7]][[2]][[2]][4, 1:2] = c(2L, 99L)
recHist[[7]][[2]][[2]][5, 1:2] = c(1L, 101L)
recHist[[7]][[2]][[2]][6, 1:2] = c(2L, 149L)
recHist[[7]][[2]][[2]][7, 1:2] = c(1L, 151L)

# ---- getIbdRecHist ----

test_that("getIbdRecHist converts recHist & pedigree to ibdRecHist correctly", {
  # IBD recombinations (since the base generation)
  expect = recHist
  # ind 1
  expect[[1]][[1]][[1]] = matrix(data = c(1L, 1L), nrow = 1, ncol = 2)
  expect[[1]][[1]][[2]] = matrix(data = c(2L, 1L), nrow = 1, ncol = 2)
  expect[[1]][[2]] = expect[[1]][[1]]
  # ind 2
  expect[[2]][[1]][[1]] = matrix(data = c(3L, 1L), nrow = 1, ncol = 2)
  expect[[2]][[1]][[2]] = matrix(data = c(4L, 1L), nrow = 1, ncol = 2)
  expect[[2]][[2]] = expect[[2]][[1]]
  # ind 3
  expect[[3]][[1]][[1]] = matrix(data = c(5L, 1L), nrow = 1, ncol = 2)
  expect[[3]][[1]][[2]] = matrix(data = c(6L, 1L), nrow = 1, ncol = 2)
  expect[[3]][[2]] = expect[[3]][[1]]
  # ind 4
  expect[[4]][[1]][[1]] = matrix(data = c(7L, 1L), nrow = 1, ncol = 2)
  expect[[4]][[1]][[2]] = matrix(data = c(8L, 1L), nrow = 1, ncol = 2)
  expect[[4]][[2]] = expect[[4]][[1]]
  # ind 5
  expect[[5]][[1]][[1]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[5]][[1]][[1]][1, 1:2] = c(1L, 1L)
  expect[[5]][[1]][[1]][2, 1:2] = c(2L, 10L)
  expect[[5]][[1]][[1]][3, 1:2] = c(1L, 50L)
  expect[[5]][[1]][[1]][4, 1:2] = c(2L, 100L)
  expect[[5]][[1]][[1]][5, 1:2] = c(1L, 150L)
  expect[[5]][[1]][[2]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[5]][[1]][[2]][1, 1:2] = c(4L, 1L)
  expect[[5]][[1]][[2]][2, 1:2] = c(3L, 5L)
  expect[[5]][[1]][[2]][3, 1:2] = c(4L, 10L)
  expect[[5]][[1]][[2]][4, 1:2] = c(3L, 15L)
  expect[[5]][[1]][[2]][5, 1:2] = c(4L, 20L)
  expect[[5]][[2]][[1]] = matrix(data = 0L, nrow = 10, ncol = 2)
  expect[[5]][[2]][[1]][ 1, 1:2] = c(1L, 1L)
  expect[[5]][[2]][[1]][ 2, 1:2] = c(2L, 2L)
  expect[[5]][[2]][[1]][ 3, 1:2] = c(1L, 3L)
  expect[[5]][[2]][[1]][ 4, 1:2] = c(2L, 4L)
  expect[[5]][[2]][[1]][ 5, 1:2] = c(1L, 5L)
  expect[[5]][[2]][[1]][ 6, 1:2] = c(2L, 6L)
  expect[[5]][[2]][[1]][ 7, 1:2] = c(1L, 7L)
  expect[[5]][[2]][[1]][ 8, 1:2] = c(2L, 8L)
  expect[[5]][[2]][[1]][ 9, 1:2] = c(1L, 9L)
  expect[[5]][[2]][[1]][10, 1:2] = c(2L, 10L)
  expect[[5]][[2]][[2]] = matrix(data = 0L, nrow = 10, ncol = 2)
  expect[[5]][[2]][[2]][ 1, 1:2] = c(4L, 1L)
  expect[[5]][[2]][[2]][ 2, 1:2] = c(3L, 3L)
  expect[[5]][[2]][[2]][ 3, 1:2] = c(4L, 5L)
  expect[[5]][[2]][[2]][ 4, 1:2] = c(3L, 7L)
  expect[[5]][[2]][[2]][ 5, 1:2] = c(4L, 9L)
  expect[[5]][[2]][[2]][ 6, 1:2] = c(3L, 11L)
  expect[[5]][[2]][[2]][ 7, 1:2] = c(4L, 13L)
  expect[[5]][[2]][[2]][ 8, 1:2] = c(3L, 15L)
  expect[[5]][[2]][[2]][ 9, 1:2] = c(4L, 17L)
  expect[[5]][[2]][[2]][10, 1:2] = c(3L, 19L)
  # ind 6
  expect[[6]][[1]][[1]] = matrix(data = 0L, nrow = 4, ncol = 2)
  expect[[6]][[1]][[1]][1, 1:2] = c(6L, 1L)
  expect[[6]][[1]][[1]][2, 1:2] = c(5L, 150L)
  expect[[6]][[1]][[1]][3, 1:2] = c(6L, 200L)
  expect[[6]][[1]][[1]][4, 1:2] = c(5L, 250L)
  expect[[6]][[1]][[2]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[6]][[1]][[2]][1, 1:2] = c(8L, 1L)
  expect[[6]][[1]][[2]][2, 1:2] = c(7L, 9L)
  expect[[6]][[1]][[2]][3, 1:2] = c(8L, 11L)
  expect[[6]][[1]][[2]][4, 1:2] = c(7L, 50L)
  expect[[6]][[1]][[2]][5, 1:2] = c(8L, 100L)
  expect[[6]][[2]][[1]] = matrix(data = 0L, nrow = 11, ncol = 2)
  expect[[6]][[2]][[1]][ 1, 1:2] = c(5L, 1L)
  expect[[6]][[2]][[1]][ 2, 1:2] = c(6L, 10L)
  expect[[6]][[2]][[1]][ 3, 1:2] = c(5L, 20L)
  expect[[6]][[2]][[1]][ 4, 1:2] = c(6L, 30L)
  expect[[6]][[2]][[1]][ 5, 1:2] = c(5L, 40L)
  expect[[6]][[2]][[1]][ 6, 1:2] = c(6L, 50L)
  expect[[6]][[2]][[1]][ 7, 1:2] = c(5L, 60L)
  expect[[6]][[2]][[1]][ 8, 1:2] = c(6L, 70L)
  expect[[6]][[2]][[1]][ 9, 1:2] = c(5L, 80L)
  expect[[6]][[2]][[1]][10, 1:2] = c(6L, 90L)
  expect[[6]][[2]][[1]][11, 1:2] = c(5L, 100L)
  expect[[6]][[2]][[2]] = matrix(data = 0L, nrow = 6, ncol = 2)
  expect[[6]][[2]][[2]][1, 1:2] = c(8L, 1L)
  expect[[6]][[2]][[2]][2, 1:2] = c(7L, 50L)
  expect[[6]][[2]][[2]][3, 1:2] = c(8L, 100L)
  expect[[6]][[2]][[2]][4, 1:2] = c(7L, 150L)
  expect[[6]][[2]][[2]][5, 1:2] = c(8L, 200L)
  expect[[6]][[2]][[2]][6, 1:2] = c(7L, 250L)
  # ind 7
  expect[[7]][[1]][[1]] = matrix(data = 0L, nrow = 12, ncol = 2)
  expect[[7]][[1]][[1]][ 1, 1:2] = c(4L, 1L)
  expect[[7]][[1]][[1]][ 2, 1:2] = c(3L, 5L)
  expect[[7]][[1]][[1]][ 3, 1:2] = c(1L, 6L)
  expect[[7]][[1]][[1]][ 4, 1:2] = c(4L, 10L)
  expect[[7]][[1]][[1]][ 5, 1:2] = c(3L, 15L)
  expect[[7]][[1]][[1]][ 6, 1:2] = c(4L, 20L)
  expect[[7]][[1]][[1]][ 7, 1:2] = c(2L, 45L)
  expect[[7]][[1]][[1]][ 8, 1:2] = c(1L, 50L)
  expect[[7]][[1]][[1]][ 9, 1:2] = c(2L, 100L)
  expect[[7]][[1]][[1]][10, 1:2] = c(4L, 110L)
  expect[[7]][[1]][[1]][11, 1:2] = c(2L, 140L)
  expect[[7]][[1]][[1]][12, 1:2] = c(1L, 150L)
  expect[[7]][[1]][[2]] = matrix(data = 0L, nrow = 5, ncol = 2)
  expect[[7]][[1]][[2]][1, 1:2] = c(8L, 1L)
  expect[[7]][[1]][[2]][2, 1:2] = c(7L, 9L)
  expect[[7]][[1]][[2]][3, 1:2] = c(8L, 11L)
  expect[[7]][[1]][[2]][4, 1:2] = c(7L, 50L)
  expect[[7]][[1]][[2]][5, 1:2] = c(8L, 100L)
  expect[[7]][[2]][[1]] = matrix(data = 0L, nrow = 16, ncol = 2)
  expect[[7]][[2]][[1]][ 1, 1:2] = c(4L, 1L)
  expect[[7]][[2]][[1]][ 2, 1:2] = c(3L, 3L)
  expect[[7]][[2]][[1]][ 3, 1:2] = c(4L, 5L)
  expect[[7]][[2]][[1]][ 4, 1:2] = c(2L, 6L)
  expect[[7]][[2]][[1]][ 5, 1:2] = c(1L, 7L)
  expect[[7]][[2]][[1]][ 6, 1:2] = c(2L, 8L)
  expect[[7]][[2]][[1]][ 7, 1:2] = c(1L, 9L)
  expect[[7]][[2]][[1]][ 8, 1:2] = c(4L, 10L)
  expect[[7]][[2]][[1]][ 9, 1:2] = c(3L, 11L)
  expect[[7]][[2]][[1]][10, 1:2] = c(4L, 13L)
  expect[[7]][[2]][[1]][11, 1:2] = c(3L, 15L)
  expect[[7]][[2]][[1]][12, 1:2] = c(4L, 17L)  
  expect[[7]][[2]][[1]][13, 1:2] = c(3L, 19L)  
  expect[[7]][[2]][[1]][14, 1:2] = c(2L, 45L)  
  expect[[7]][[2]][[1]][15, 1:2] = c(3L, 110L)  
  expect[[7]][[2]][[1]][16, 1:2] = c(2L, 140L)  
  expect[[7]][[2]][[2]] = matrix(data = 0L, nrow = 15, ncol = 2)
  expect[[7]][[2]][[2]][ 1, 1:2] = c(5L, 1L)
  expect[[7]][[2]][[2]][ 2, 1:2] = c(6L, 10L)
  expect[[7]][[2]][[2]][ 3, 1:2] = c(5L, 20L)
  expect[[7]][[2]][[2]][ 4, 1:2] = c(6L, 30L)
  expect[[7]][[2]][[2]][ 5, 1:2] = c(5L, 40L)
  expect[[7]][[2]][[2]][ 6, 1:2] = c(7L, 50L)
  expect[[7]][[2]][[2]][ 7, 1:2] = c(6L, 75L)
  expect[[7]][[2]][[2]][ 8, 1:2] = c(5L, 80L)
  expect[[7]][[2]][[2]][ 9, 1:2] = c(6L, 90L)
  expect[[7]][[2]][[2]][10, 1:2] = c(7L, 99L)
  expect[[7]][[2]][[2]][11, 1:2] = c(8L, 100L)
  expect[[7]][[2]][[2]][12, 1:2] = c(5L, 101L)  
  expect[[7]][[2]][[2]][13, 1:2] = c(8L, 149L)  
  expect[[7]][[2]][[2]][14, 1:2] = c(7L, 150L)  
  expect[[7]][[2]][[2]][15, 1:2] = c(5L, 151L)  
  result = AlphaSimR:::getIbdRecHist(recHist     = recHist,
                                     pedigree    = pedigree,
                                     nLociPerChr = nLociPerChr)$ibdRecHist
  for (ind in 1:7) {
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
  output = AlphaSimR:::getIbdHaplo2(ibdRecHist  = ibdRecHist,
                                    individuals = 1L:nrow(pedigree),
                                    nLociPerChr = nLociPerChr)
  expect = matrix(data = 0L, nrow = 14, ncol = 600)
  # ind 1
  expect[1, ] = 1L
  expect_identical(object = output[1, ],
                   expected = expect[1, ],
                   info = paste0("ind 1, chr 1 & 2, gamete 1"))
  expect[2, ] = 2L
  expect_identical(object = output[2, ],
                   expected = expect[2, ],
                   info = paste0("ind 1, chr 1 & 2, gamete 2"))
  # ind 2
  expect[3, ] = 3L
  expect_identical(object = output[3, ],
                   expected = expect[3, ],
                   info = paste0("ind 2, chr 1 & 2, gamete 1"))
  expect[4, ] = 4L
  expect_identical(object = output[4, ],
                   expected = expect[4, ],
                   info = paste0("ind 2, chr 1 & 2, gamete 2"))
  # ind 3
  expect[5, ] = 5L
  expect_identical(object = output[5, ],
                   expected = expect[5, ],
                   info = paste0("ind 3, chr 1 & 2, gamete 1"))
  expect[6, ] = 6L
  expect_identical(object = output[6, ],
                   expected = expect[6, ],
                   info = paste0("ind 3, chr 1 & 2, gamete 2"))
  # ind 4
  expect[7, ] = 7L
  expect_identical(object = output[7, ],
                   expected = expect[7, ],
                   info = paste0("ind 4, chr 1 & 2, gamete 1"))
  expect[8, ] = 8L
  expect_identical(object = output[8, ],
                   expected = expect[8, ],
                   info = paste0("ind 4, chr 1 & 2, gamete 2"))
  # ind 5 - chr 1, paternal
  expect[9,   1:009] = 1L
  expect[9,  10:049] = 2L
  expect[9,  50:099] = 1L
  expect[9, 100:149] = 2L
  expect[9, 150:300] = 1L
  expect_identical(object = output[9, 1:300],
                   expected = expect[9, 1:300],
                   info = paste0("ind 5, chr 1, gamete 1"))
  # ind 5 - chr 1, maternal
  expect[10,  1:004] = 4L
  expect[10,  5:009] = 3L
  expect[10, 10:014] = 4L
  expect[10, 15:019] = 3L
  expect[10, 20:300] = 4L
  expect_identical(object = output[10, 1:300],
                   expected = expect[10, 1:300],
                   info = paste0("ind 5, chr 1, gamete 2"))
  # ind 5 - chr 2, paternal
  expect[9, 300 + 001:001] = 1L
  expect[9, 300 + 002:002] = 2L
  expect[9, 300 + 003:003] = 1L
  expect[9, 300 + 004:004] = 2L
  expect[9, 300 + 005:005] = 1L
  expect[9, 300 + 006:006] = 2L
  expect[9, 300 + 007:007] = 1L
  expect[9, 300 + 008:008] = 2L
  expect[9, 300 + 009:009] = 1L
  expect[9, 300 + 010:300] = 2L
  expect_identical(object = output[9, 301:600],
                   expected = expect[9, 301:600],
                   info = paste0("ind 5, chr 2, gamete 1"))
  # ind 5 - chr 2, maternal
  expect[10, 300 + 001:002] = 4L
  expect[10, 300 + 003:004] = 3L
  expect[10, 300 + 005:006] = 4L
  expect[10, 300 + 007:008] = 3L
  expect[10, 300 + 009:010] = 4L
  expect[10, 300 + 011:012] = 3L
  expect[10, 300 + 013:014] = 4L
  expect[10, 300 + 015:016] = 3L
  expect[10, 300 + 017:018] = 4L
  expect[10, 300 + 019:300] = 3L
  expect_identical(object = output[10, 301:600],
                   expected = expect[10, 301:600],
                   info = paste0("ind 5, chr 2, gamete 2"))
  # ind 6 - chr 1, paternal
  expect[11,   1:149] = 6L
  expect[11, 150:199] = 5L
  expect[11, 200:249] = 6L
  expect[11, 250:300] = 5L
  expect_identical(object = output[11, 1:300],
                   expected = expect[11, 1:300],
                   info = paste0("ind 6, chr 1, gamete 1"))
  # ind 6 - chr 1, maternal
  expect[12,   1:008] = 8L
  expect[12,   9:010] = 7L
  expect[12,  11:049] = 8L
  expect[12,  50:099] = 7L
  expect[12, 100:300] = 8L
  expect_identical(object = output[12, 1:300],
                   expected = expect[12, 1:300],
                   info = paste0("ind 6, chr 1, gamete 2"))
  # ind 6 - chr 2, paternal
  expect[11, 300 + 001:009] = 5L
  expect[11, 300 + 010:019] = 6L
  expect[11, 300 + 020:029] = 5L
  expect[11, 300 + 030:039] = 6L
  expect[11, 300 + 040:049] = 5L
  expect[11, 300 + 050:059] = 6L
  expect[11, 300 + 060:069] = 5L
  expect[11, 300 + 070:079] = 6L
  expect[11, 300 + 080:089] = 5L
  expect[11, 300 + 090:099] = 6L
  expect[11, 300 + 100:300] = 5L
  expect_identical(object = output[11, 301:600],
                   expected = expect[11, 301:600],
                   info = paste0("ind 6, chr 2, gamete 1"))
  # ind 6 - chr 2, maternal
  expect[12, 300 + 001:049] = 8L
  expect[12, 300 + 050:099] = 7L
  expect[12, 300 + 100:149] = 8L
  expect[12, 300 + 150:199] = 7L
  expect[12, 300 + 200:249] = 8L
  expect[12, 300 + 250:300] = 7L
  expect_identical(object = output[12, 301:600],
                   expected = expect[12, 301:600],
                   info = paste0("ind 6, chr 2, gamete 2"))
  # ind 7 - chr 1, paternal
  expect[13,   1:004] = 4L
  expect[13,   5:005] = 3L
  expect[13,   6:009] = 1L
  expect[13,  10:014] = 4L
  expect[13,  15:019] = 3L
  expect[13,  20:044] = 4L
  expect[13,  45:049] = 2L
  expect[13,  50:099] = 1L
  expect[13, 100:109] = 2L
  expect[13, 110:139] = 4L
  expect[13, 140:149] = 2L
  expect[13, 150:300] = 1L
  expect_identical(object = output[13, 1:300],
                   expected = expect[13, 1:300],
                   info = paste0("ind 7, chr 1, gamete 1"))
  # ind 7 - chr 1, maternal
  expect[14,   1:008] = 8L
  expect[14,   9:010] = 7L
  expect[14,  11:049] = 8L
  expect[14,  50:099] = 7L
  expect[14, 100:300] = 8L
  expect_identical(object = output[14, 1:300],
                   expected = expect[14, 1:300],
                   info = paste0("ind 7, chr 1, gamete 2"))
  # ind 7 - chr 2, paternal
  expect[13, 300 + 001:002] = 4L
  expect[13, 300 + 003:004] = 3L
  expect[13, 300 + 005:005] = 4L
  expect[13, 300 + 006:006] = 2L
  expect[13, 300 + 007:007] = 1L
  expect[13, 300 + 008:008] = 2L
  expect[13, 300 + 009:009] = 1L
  expect[13, 300 + 010:010] = 4L
  expect[13, 300 + 011:012] = 3L
  expect[13, 300 + 013:014] = 4L
  expect[13, 300 + 015:016] = 3L
  expect[13, 300 + 017:018] = 4L
  expect[13, 300 + 019:044] = 3L
  expect[13, 300 + 045:109] = 2L
  expect[13, 300 + 110:139] = 3L
  expect[13, 300 + 140:300] = 2L
  expect_identical(object = output[13, 301:600],
                   expected = expect[13, 301:600],
                   info = paste0("ind 7, chr 2, gamete 1"))
  # ind 7 - chr 2, maternal
  expect[14, 300 + 001:009] = 5L
  expect[14, 300 + 010:019] = 6L
  expect[14, 300 + 020:029] = 5L
  expect[14, 300 + 030:039] = 6L
  expect[14, 300 + 040:049] = 5L
  expect[14, 300 + 050:074] = 7L
  expect[14, 300 + 075:079] = 6L
  expect[14, 300 + 080:089] = 5L
  expect[14, 300 + 090:098] = 6L
  expect[14, 300 + 099:099] = 7L
  expect[14, 300 + 100:100] = 8L
  expect[14, 300 + 101:148] = 5L
  expect[14, 300 + 149:149] = 8L
  expect[14, 300 + 150:150] = 7L
  expect[14, 300 + 151:300] = 5L
  expect_identical(object = output[14, 301:600],
                   expected = expect[14, 301:600],
                   info = paste0("ind 7, chr 2, gamete 2"))
  
  # Test subsetting
  output = AlphaSimR:::getIbdHaplo2(ibdRecHist  = ibdRecHist,
                                    individuals = 7L,
                                    nLociPerChr = nLociPerChr)
  # ind 7 - chr 1, paternal
  expect_identical(object = output[1, 1:300],
                   expected = expect[13, 1:300],
                   info = paste0("ind 7, chr 1, gamete 1, from subsetted call"))
  # ind 7 - chr 1, maternal
  expect_identical(object = output[2, 1:300],
                   expected = expect[14, 1:300],
                   info = paste0("ind 7, chr 1, gamete 2, from subsetted call"))
  # ind 7 - chr 2, paternal

  expect_identical(object = output[1, 301:600],
                   expected = expect[13, 301:600],
                   info = paste0("ind 7, chr 2, gamete 1, from subsetted call"))
  # ind 7 - chr 2, maternal
  expect_identical(object = output[2, 301:600],
                   expected = expect[14, 301:600],
                   info = paste0("ind 7, chr 2, gamete 2, from subsetted call"))
})