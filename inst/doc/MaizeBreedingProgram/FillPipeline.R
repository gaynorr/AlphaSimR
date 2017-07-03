#Set initial yield trials with unique individuals
for(year in 1:6){
  cat("FillPipeline year:",year,"of 6\n")
  if(year<7){
    #Year 1
    MaleF1 = randCross(MaleParents,nCrosses)
    FemaleF1 = randCross(FemaleParents,nCrosses)
    MaleDH = makeDH(MaleF1,nDH)
    FemaleDH = makeDH(FemaleF1,nDH)
  }
  if(year<6){
    #Year 2
    w=runif(1)
    MaleYT1 = setPhenoGCA(MaleDH,FemaleTester1,varE=varE,
                          reps=repYT1,inbred=T,w=w)
    FemaleYT1 = setPhenoGCA(FemaleDH,MaleTester1,varE=varE,
                            reps=repYT1,inbred=T,w=w)
  }
  if(year<5){
    #Year 3
    w = runif(1)
    MaleYT2 = selectInd(MaleYT1,nInbred2)
    FemaleYT2 = selectInd(FemaleYT1,nInbred2)
    MaleYT2 = setPhenoGCA(MaleYT2,FemaleTester2,varE=varE,
                          reps=repYT2,inbred=T,w=w)
    FemaleYT2 = setPhenoGCA(FemaleYT2,MaleTester2,varE=varE,
                            reps=repYT2,inbred=T,w=w)
  }
  if(year<4){
    #Year 4
    w = runif(1)
    MaleInbredYT3 = selectInd(MaleYT2,nInbred3)
    FemaleInbredYT3 = selectInd(FemaleYT2,nInbred3)
    MaleHybridYT3 = hybridCross(MaleInbredYT3,FemaleElite,
                                varE=varE,reps=repYT3,
                                returnHybridPop=F,w=w)
    FemaleHybridYT3 = hybridCross(FemaleInbredYT3,MaleElite,
                                  varE=varE,reps=repYT3,
                                  returnHybridPop=F,w=w)
  }
  if(year<3){
    #Year 5
    w = runif(1)
    MaleHybridYT4 = selectInd(MaleHybridYT3,nYT4)
    FemaleHybridYT4 = selectInd(FemaleHybridYT3,nYT4)
    MaleHybridYT4 = setPheno(MaleHybridYT4,varE=varE,
                             reps=repYT4,w=w)
    FemaleHybridYT4 = setPheno(FemaleHybridYT4,varE=varE,
                               reps=repYT4,w=w)
    MaleInbredYT4 = MaleInbredYT3[
      MaleInbredYT3@id%in%MaleHybridYT4@mother
      ]
    FemaleInbredYT4 = FemaleInbredYT3[
      FemaleInbredYT3@id%in%FemaleHybridYT4@mother
      ]
  }
  if(year<2){
    #Year 6
    w = runif(1)
    MaleHybridYT5 = selectInd(MaleHybridYT4,nYT5)
    FemaleHybridYT5 = selectInd(FemaleHybridYT4,nYT5)
    MaleHybridYT5 = setPheno(MaleHybridYT5,varE=varE,
                             reps=repYT5,w=w)
    FemaleHybridYT5 = setPheno(FemaleHybridYT5,varE=varE,
                               reps=repYT5,w=w)
    MaleInbredYT5 = MaleInbredYT4[
      MaleInbredYT4@id%in%MaleHybridYT5@mother
      ]
    FemaleInbredYT5 = FemaleInbredYT4[
      FemaleInbredYT4@id%in%FemaleHybridYT5@mother
      ]
  }
  if(year<1){
    #Year 7, release
  }
}
