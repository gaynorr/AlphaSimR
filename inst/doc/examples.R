# AlphaSimR Examples
# Works with version 0.3.X

# *Requires MaCS
# Which can be taken from AlphaSim on the AlphaGenes website
# http://www.alphagenes.roslin.ed.ac.uk/alphasuite-softwares/alphasim/


library(AlphaSimR)

# Set-up a simulation -----------------------------------------------------

# Generate initial haplotypes, replace macsPath
macsPath = "/Users/rgaynor/Documents/testPBS/macs"
FOUNDERPOP = runMacs(macs=macsPath, #"FOUNDERPOP" is the default for other functions
                     nInd=100, #Number of individuals
                     nChr=10, #Number of chromosomes
                     segSites=1100, #Sites for SNPs and QTLs
                     inbred=TRUE, #Are initial individuals inbred
                     species="TEST", #Species history, TEST is used for speed
                     split=100) #An optional split in the population, for create heterotic pools
rm(macsPath) #No longer needed

# Fill in parameters
SIMPARAM = createSimulation(FOUNDERPOP,
                            maxQtl=1000, #Number of sites for QTLs
                            maxSnp=100, #Number of sites for SNPs
                            snpQtlOverlap=FALSE, #Can QTLs also be SNPs
                            minSnpFreq = 0.1, #Sets minimum minor allow frequency
                            gender="no") #Determines if gender restricts matings

# Add a SNP chip, can be repeated for multiple SNP chips
SIMPARAM = addSnpChip(nSnpPerChr=100,simParam=SIMPARAM) #Must be less than or equal to maxSnp

# Add a trait, can be repeated for multiple traits
# This trait has both additive and dominance effect
SIMPARAM = addTraitAD(FOUNDERPOP,
                      nQtlPerChr=1000,
                      meanG=50, #Mean value of trait
                      varG=15, #Total genetic variance of the trait
                      domDegree=0.6,
                      simParam=SIMPARAM) #Degree of dominance for each QTL

# Set-up heterotic pools --------------------------------------------------

# Extract parents
femaleParents = newPop(FOUNDERPOP[1:50],simParam=SIMPARAM)
maleParents = newPop(FOUNDERPOP[51:100],simParam=SIMPARAM)
rm(FOUNDERPOP) #No longer needed

# View summary data -------------------------------------------------------

# View genotype data
View(pullSnpGeno(femaleParents,simParam=SIMPARAM))

# View haplotype data
View(pullSnpHaplo(femaleParents,simParam=SIMPARAM))

# View population mean and variance
meanG(femaleParents) #Mean genetic values for each trait
varG(femaleParents) #Total genetic variance for each trait

# Make crosses ------------------------------------------------------------

# Examine test crosses
hybrids = hybridCross(femaleParents,maleParents,
                      returnHybridPop=T,simParam=SIMPARAM) #Uses a test cross scheme by default
gca = calcGCA(hybrids,use="gv")
plot(femaleParents@gv,gca$females[,2],
     xlab="GV (per se)", ylab="GCA",
     main="Female Parents")
plot(maleParents@gv,gca$males[,2],
     xlab="GV (per se)", ylab="GCA",
     main="Male Parents")
hist(gca$SCA[,3],main="Hybrids",xlab="GV")
rm(hybrids,gca) #No longer needed

# Tricks for faster computation of hybrids
# Produce hybrids by formal crossing, i.e. produces a true population
system.time({hybrids1 = hybridCross(femaleParents,maleParents,returnHybridPop=F,simParam=SIMPARAM)})
# Produce hybrids using shortcuts, assumes parents are inbred and doesn't return genotype data
# This method won't allow for recalculation of phenotypes on traits with GxE
system.time({hybrids2 = hybridCross(femaleParents,maleParents,returnHybridPop=T,simParam=SIMPARAM)})
object.size(hybrids1)
object.size(hybrids2)
cor(hybrids1@gv,hybrids2@gv)
rm(hybrids1,hybrids2) #No longer need

# Create DH lines----

# Make bi-parental crosses
femaleF1 = randCross(pop=femaleParents,
                     nCrosses=50,
                     simParam=SIMPARAM)
maleF1 = randCross(pop=maleParents,
                   nCrosses=50,
                   simParam=SIMPARAM)

# Create F2s by selfing
femaleF2 = self(pop=femaleF1,nProgeny=10,simParam=SIMPARAM) #Creates 10 F2s per F1
maleF2 = self(pop=maleF1,nProgeny=10,simParam=SIMPARAM)
rm(femaleF2,maleF2) #No longer needed

# Make DH lines from F1s
femaleDH = makeDH(pop=femaleF1,
                  nDH=50, #Create 50 DH lines per F1
                  simParam=SIMPARAM) 
maleDH = makeDH(pop=maleF1,
                nDH=50,simParam=SIMPARAM)
rm(femaleF1,maleF1) #No longer needed

# Evaluate inbred lines----

# Add phenotype with error variance 10
femaleDH = setPheno(pop=femaleDH,varE=10,simParam=SIMPARAM)
maleDH = setPheno(pop=maleDH,varE=10,simParam=SIMPARAM)

# Select the best 500 on per se performance
femaleDH = selectInd(pop=femaleDH,nInd=500,simParam=SIMPARAM)
maleDH = selectInd(pop=maleDH,nInd=500,simParam=SIMPARAM)

# Genetic values of selected lines
hist(femaleDH@gv,main="Female DH",xlab="GV (per se)")
hist(maleDH@gv,main="Male DH",xlab="GV (per se)")

# Evaluate hybrids----

# Create testers
femaleTesters = selectInd(maleParents,3,use="gv",simParam=SIMPARAM) #Picks three testers from other heterotic pool
maleTesters = selectInd(femaleParents,3,use="gv",simParam=SIMPARAM)

# Create test crosses
femaleTC = hybridCross(femaleDH,femaleTesters,returnHybridPop=TRUE,
                       simParam=SIMPARAM)
maleTC = hybridCross(maleDH,maleTesters,returnHybridPop=TRUE,
                     simParam=SIMPARAM)

# Calculate true combining ability
femaleGCA = calcGCA(femaleTC,use="gv")
plot(femaleDH@gv,femaleGCA$females[,2],
     xlab="GV (per se)", ylab="GCA",
     main="Female DH")
maleGCA = calcGCA(maleTC,use="gv")
plot(maleDH@gv,maleGCA$females[,2], #Select females for GCA because testers were used as males
     xlab="GV (per se)", ylab="GCA",
     main="Male DH")
rm(femaleTC,maleTC,femaleGCA,maleGCA) #No longer needed

# Set DH phenotypes using test cross GCA
# For females with a phenotype with error variance of 10
femaleDH = setPhenoGCA(pop=femaleDH,testers=femaleTesters,
                       varE=10,inbred=TRUE,simParam=SIMPARAM)
# For males using genetic values, i.e. true values
maleDH = setPhenoGCA(pop=maleDH,testers=maleTesters,use="gv",
                     inbred=TRUE,simParam=SIMPARAM)
plot(maleDH@gv,maleDH@pheno,
     xlab="GV (per se)", ylab="GCA",
     main="Male DH")

