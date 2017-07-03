#Set initial yield trials with unique individuals
#Ignoring GxY at this time
for(year in 1:7){
  cat("FillPipeline year:",year,"of 7\n")
  if(year<8){
    #Year 1
    F1 = randCross(Parents,nCrosses)
  }
  if(year<7){
    #Year 2
    DH = makeDH(F1,nDH)
  }
  if(year<6){
    #Year 3
    HDRW = setPheno(DH,varE=varE,reps=repHDRW)
  }
  if(year<5){
    #Year 4
    PYT = selectWithinFam(HDRW,famMax)
    PYT = selectInd(PYT,nPYT)
    PYT = setPheno(PYT,varE=varE,reps=repPYT)
  }
  if(year<4){
    #Year 5
    AYT = selectInd(PYT,nAYT)
    AYT = setPheno(AYT,varE=varE,reps=repAYT)
  }
  if(year<3){
    #Year 6
    EYT1 = selectInd(AYT,nEYT)
    EYT1 = setPheno(EYT1,varE=varE,reps=repEYT)
  }
  if(year<2){
    #Year 7
    EYT2 = setPheno(EYT1,varE=varE,reps=repEYT)
    EYT2@pheno = (EYT1@pheno+EYT2@pheno)/2
  }
  if(year<1){
    #Year 8
  }
}
