library(AlphaSimR)
macsPath = "/Users/rgaynor/Documents/testPBS/macs"
FOUNDERPOP = runMacs(macs=macsPath,nInd=100,nChr=10,segSites=1200,
                     inbred=FALSE,species="TEST")
rm(macsPath)
SIMPARAM = createSimulation(FOUNDERPOP,maxQtl=200,maxSnp=1000,gender="yes_sys")
SIMPARAM = addTraitA(FOUNDERPOP,nQtlPerChr=100,meanG=0,varG=1,simParam=SIMPARAM)
SIMPARAM = addSnpChip(1000,simParam=SIMPARAM)

pop = newPop(FOUNDERPOP,simParam=SIMPARAM)
unlink("GS",recursive=TRUE)
dir.create("GS")
writeRecords(pop,"GS",1,simParam=SIMPARAM)
for(gen in 1:9){
  pop = randCross(pop,100,simParam=SIMPARAM)
  writeRecords(pop,"GS",1,simParam=SIMPARAM)
}
#Predict
pop = randCross(pop,100,simParam=SIMPARAM)
gsModel = RRBLUP("GS",1,use="GV",simParam=SIMPARAM)
pop = setEBV(pop,gsModel)
cor(pop@gv[,1],pop@ebv[,1])
#Select on EBV
bestInd = selectInd(pop,10,use="ebv",simParam=SIMPARAM)

