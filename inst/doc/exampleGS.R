library(AlphaSimR)
macsPath = "/Users/rgaynor/Documents/testPBS/macs"
FOUNDERPOP = runMacs(macs=macsPath,nInd=100,nChr=10,segSites=1200,
                     inbred=FALSE,species="TEST")
rm(macsPath)
SIMPARAM = createSimulation(maxQtl=200,maxSnp=1000,gender="yes_sys")
SIMPARAM = addTraitA(nQtlPerChr=100,meanG=0,varG=1)
SIMPARAM = addSnpChip(1000)

pop = newPop(FOUNDERPOP)
unlink("GS",recursive=TRUE)
dir.create("GS")
writeRecords(pop,"GS",1)
for(gen in 1:9){
  pop = randCross(pop,100)
  writeRecords(pop,"GS",1)
}
#Predict
pop = randCross(pop,100)
gsModel = RRBLUP("GS",1,use="GV")
pop = setEBV(pop,gsModel)
cor(pop@gv[,1],pop@ebv[,1])
#Select on EBV
bestInd = selectInd(pop,10,use="ebv")

