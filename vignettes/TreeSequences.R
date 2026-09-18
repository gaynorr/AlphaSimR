## ----eval=FALSE---------------------------------------------------------------
# library(AlphaSimR)
# 
# founderPop = quickHaplo(nInd=20, nChr=3, segSites=100)
# SP = SimParam$new(founderPop)
# SP$setTrackRec(TRUE)
# 
# founders = newPop(founderPop, simParam=SP)
# generation1 = randCross(founders, nCrosses=20, simParam=SP)
# generation2 = randCross(generation1, nCrosses=20, simParam=SP)
# 
# trees = asTreeSequence(generation2, simParam=SP)
# writeTreeSequence(trees, "breeding-program.trees")

## ----eval=FALSE---------------------------------------------------------------
# trees[[1]]$num_samples()
# trees[[1]]$num_trees()
# trees[[1]]$num_sites()
# trees[[1]]$num_mutations()

