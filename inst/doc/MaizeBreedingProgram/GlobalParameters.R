# Crossing->DH ------------------------------------------------------------

#Parents per heterotic pool
nParents = 50

#Number of crosses
nCrosses = 80

#DH lines produced per cross
nDH = 50

# Selection on GCA --------------------------------------------------------

#Number of inbreds per heterotic pool per stage
nInbred1 = nCrosses*nDH #Do not change
nInbred2 = 400
nInbred3 = 40

#Number of testers per heterotic pool per stage
#Values must be smaller than nElite
nTester1 = 1
nTester2 = 3

#Yield trial entries
nYT1 = nInbred1*nTester1 #Do not change
nYT2 = nInbred2*nTester2 #Do not change

# Selection on SCA --------------------------------------------------------

#Elite parents per heterotic pool
nElite = 5

#Elite YT size
nYT3 = nInbred3*nElite #Do not change
nYT4 = 20
nYT5 = 4


# Genetic and error terms -------------------------------------------------

#Number of QTL per chromosome
nQtl = 400

#Number of SNP per chromosome
nSnp = 0

#Heterotic pool split
nGenSplit = 100

#Initial inbred mean, bushels per acre
initMeanG = 70

#Initial inbred variance, bushels per acre
initVarG = 20
initVarGE = 40

#Degree of dominance
ddMean = 0.8
dd = rnorm(nQtl*10,mean=ddMean,sd=0.05)

#Yield trial error variance, bushels per acre
#Relates to error variance for an entry mean
varE = 60

#Yield trial effective replications for calculating error
#Roughly related to the number of locations relative to YT_1
repYT1 = 1
repYT2 = 2
repYT3 = 4
repYT4 = 8
repYT5 = 100
