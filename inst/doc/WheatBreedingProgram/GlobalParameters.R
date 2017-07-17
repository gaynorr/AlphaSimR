# Crossing->DH ------------------------------------------------------------

#Parents 
nParents = 50 #simulates an equal number of landraces
nParentsPYT = 20 #Relates to crossing block selection
nParentsAYT = 10 #Relate to crossing block selection

#Initial parents mean and variance
initMeanG = 50
initVarG = 10
initVarGE = 20

#Number of crosses per year
nCrosses = 100

#DH lines produced per cross
nDH = 100

#The maximum number of DH lines per cross to entry the PYT
#nCrosses*famMax must be greater than or equal to nPYT
famMax = 10

#Entries per yield trial
nPYT = 500
nAYT = 50
nEYT = 10

#Yield trial error variance, bushels per acre
#Relates to error variance for an entry mean
varE = 40

#Effective replication of yield trials
repHDRW = 4/9 #h2=0.1
repPYT = 1
repAYT = 4
repEYT = 8

#Number of QTL per chromosome
nQtl = 1000

#Number of SNP per chromosome
nSnp = 0
