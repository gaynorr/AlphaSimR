# Generate initial haplotypes
FOUNDERPOP = runMacs(macs=macsPath,
                     nInd=nParents*2, 
                     nChr=10, 
                     segSites=nQtl+nSnp,
                     inbred=TRUE, 
                     species="MAIZE",
                     split=nGenSplit)
SIMPARAM = createSimulation(nQtl,nSnp)
if(nSnp>0){
  SIMPARAM = addSnpChip(nSnp)
}
SIMPARAM = addTraitADG(nQtl,initMeanG,initVarG,dd,initVarGE)

# Split heterotic pools to form initial parents
FemaleParents = newPop(FOUNDERPOP[1:nParents])
MaleParents = newPop(FOUNDERPOP[(nParents+1):(nParents*2)])

#Set hybrid parents for later yield trials
MaleElite = selectInd(MaleParents,nElite,use="gv")
FemaleElite = selectInd(FemaleParents,nElite,use="gv")

#Reverse order to keep best parent in longer
MaleElite = MaleElite[nElite:1]
FemaleElite = FemaleElite[nElite:1]

#Set initial testers for YT1 and YT2
#Requires nTesters to be smaller than nElite
MaleTester1 = MaleElite[1:nTester1]
FemaleTester1 = FemaleElite[1:nTester1]
MaleTester2 = MaleElite[1:nTester2]
FemaleTester2 = FemaleElite[1:nTester2]


