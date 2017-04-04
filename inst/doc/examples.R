# AlphaSimR Examples
# Works with version 0.1

# *Requires MaCS
# Which can be taken from AlphaSim on the AlphaGenes website
# http://www.alphagenes.roslin.ed.ac.uk/alphasuite-softwares/alphasim/

# Required R libraries:
# 1. devtools (to download from Bitbucket)
# 2. MASS (should come with R)
# 3. Rcpp (requires C++ compiler)
# 4. RcppArmadillo

# Install package ---------------------------------------------------------

library(devtools)
install_bitbucket("hickeyjohnteam/AlphaSimR")
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
SIMPARAM = createSimulation(maxQtl=1000, #Number of sites for QTLs
                            maxSnp=100, #Number of sites for SNPs
                            snpQtlOverlap=FALSE, #Can QTLs also be SNPs
                            minSnpFreq = 0.1, #Sets minimum minor allow frequency
                            gender="no") #Determines if gender restricts matings

# Add a SNP chip, can be repeated for multiple SNP chips
SIMPARAM = addSnpChip(nSnpPerChr=100) #Must be less than or equal to maxSnp

# Add a trait, can be repeated for multiple traits
# This trait has both additive and dominance effect
SIMPARAM = addTraitAD(nQtlPerChr=1000,
                      meanG=50, #Mean value of trait
                      varG=15, #Total genetic variance of the trait
                      domDegree=0.6) #Degree of dominance for each QTL

# Set-up heterotic pools --------------------------------------------------

# Extract parents
femaleParents = newPop(FOUNDERPOP[1:50],
                       id=paste0("G0_F",1:50)) #Giving lines a custom id
maleParents = newPop(FOUNDERPOP[51:100],
                     id=paste0("G0_M",1:50)) #Giving lines a custom id
rm(FOUNDERPOP) #No longer needed

# View summary data -------------------------------------------------------

# View genotype data
View(pullSnpGeno(femaleParents))

# View haplotype data
View(pullSnpHaplo(femaleParents))

# View population mean and variance
meanG(femaleParents) #Mean genetic values for each trait
varG(femaleParents) #Total genetic variance for each trait

# Make crosses ------------------------------------------------------------

# Examine test crosses
hybrids = hybridCross(femaleParents,maleParents,
                      returnHybridPop=T) #Uses a test cross scheme by default
gca = calcGCA(hybrids,useGv=TRUE)
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
system.time({hybrids1 = hybridCross(femaleParents,maleParents,returnHybridPop=F)})
# Produce hybrids using shortcuts, assumes parents are inbred and doesn't return genotype data
# This method won't allow for recalculation of phenotypes on traits with GxE
system.time({hybrids2 = hybridCross(femaleParents,maleParents,returnHybridPop=T)})
object.size(hybrids1)
object.size(hybrids2)
cor(hybrids1@gv,hybrids2@gv)
rm(hybrids1,hybrids2) #No longer need

# Create DH lines----

# Make bi-parental crosses
femaleF1 = randCross(pop=femaleParents,
                     nCrosses=50) #simParam=SIMPARAM by default
maleF1 = randCross(pop=maleParents,
                   nCrosses=50) #simParam=SIMPARAM by default

# Create F2s by selfing
femaleF2 = self(pop=femaleF1,nProgeny=10) #Creates 10 F2s per F1
maleF2 = self(pop=maleF1,nProgeny=10)
rm(femaleF2,maleF2) #No longer needed

# Make DH lines from F1s
femaleDH = makeDH(pop=femaleF1,
                  nDH=50, #Create 50 DH lines per F1
                  id=paste0("G1_F",1:(50*50))) #Optional for pedigree
maleDH = makeDH(pop=maleF1,
                nDH=50,
                id=paste0("G1_M",1:(50*50)))
rm(femaleF1,maleF1) #No longer needed

# Evaluate inbred lines----

# Add phenotype with error variance 10
femaleDH = setPheno(pop=femaleDH,varE=10)
maleDH = setPheno(pop=maleDH,varE=10)

# Select the best 500 on per se performance
femaleDH = selectInd(pop=femaleDH,nInd=500)
maleDH = selectInd(pop=maleDH,nInd=500)

# Genetic values of selected lines
hist(femaleDH@gv,main="Female DH",xlab="GV (per se)")
hist(maleDH@gv,main="Male DH",xlab="GV (per se)")

# Evaluate hybrids----

# Create testers
femaleTesters = selectInd(maleParents,3,useGv=TRUE) #Picks three testers from other heterotic pool
maleTesters = selectInd(femaleParents,3,useGv=TRUE)

# Create test crosses
femaleTC = hybridCross(femaleDH,femaleTesters,returnHybridPop=TRUE)
maleTC = hybridCross(maleDH,maleTesters,returnHybridPop=TRUE)

# Calculate true combining ability
femaleGCA = calcGCA(femaleTC,useGv=TRUE)
plot(femaleDH@gv,femaleGCA$females[,2],
     xlab="GV (per se)", ylab="GCA",
     main="Female DH")
maleGCA = calcGCA(maleTC,useGv=TRUE)
plot(maleDH@gv,maleGCA$females[,2], #Select females for GCA because testers were used as males
     xlab="GV (per se)", ylab="GCA",
     main="Male DH")
rm(femaleTC,maleTC,femaleGCA,maleGCA) #No longer needed

# Set DH phenotypes using test cross GCA
# For females with a phenotype with error variance of 10
femaleDH = setPhenoGCA(pop=femaleDH,testers=femaleTesters,varE=10,inbred=TRUE)
# For males using genetic values, i.e. true values
maleDH = setPhenoGCA(pop=maleDH,testers=maleTesters,useGv=TRUE,inbred=TRUE)
plot(maleDH@gv,maleDH@pheno,
     xlab="GV (per se)", ylab="GCA",
     main="Male DH")

