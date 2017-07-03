# Generate initial haplotypes
FOUNDERPOP = runMacs(macs=macsPath,
                     nInd=nParents, 
                     nChr=21, 
                     segSites=nQtl+nSnp,
                     inbred=TRUE, 
                     species="WHEAT")
SIMPARAM = createSimulation(nQtl,nSnp)
if(nSnp>0){
  SIMPARAM = addSnpChip(nSnp)
}
SIMPARAM = addTraitAG(nQtl,initMeanG,initVarG,initVarGE)

Parents = newPop(FOUNDERPOP)
#Add phenotype reflecting 2 years of evaluation in EYT
#ignoring GxY
Parents = setPheno(Parents,varE=varE,reps=repEYT*2)
rm(FOUNDERPOP)


