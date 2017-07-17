#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data
w = runif(1)

#Year 7
#Release hybrid

#Year 6
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

#Year 5
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

#Year 4
MaleInbredYT3 = selectInd(MaleYT2,nInbred3)
FemaleInbredYT3 = selectInd(FemaleYT2,nInbred3)
MaleHybridYT3 = hybridCross(MaleInbredYT3,FemaleElite,
                            varE=varE,reps=repYT3,
                            returnHybridPop=F,w=w)
FemaleHybridYT3 = hybridCross(FemaleInbredYT3,MaleElite,
                              varE=varE,reps=repYT3,
                              returnHybridPop=F,w=w)

#Year 3
MaleYT2 = selectInd(MaleYT1,nInbred2)
FemaleYT2 = selectInd(FemaleYT1,nInbred2)
MaleYT2 = setPhenoGCA(MaleYT2,FemaleTester2,varE=varE,
                      reps=repYT2,inbred=T,w=w)
FemaleYT2 = setPhenoGCA(FemaleYT2,MaleTester2,varE=varE,
                        reps=repYT2,inbred=T,w=w)

#Year 2
MaleYT1 = setPhenoGCA(MaleDH,FemaleTester1,varE=varE,
                      reps=repYT1,inbred=T,w=w)
FemaleYT1 = setPhenoGCA(FemaleDH,MaleTester1,varE=varE,
                        reps=repYT1,inbred=T,w=w)

#Year 1
MaleF1 = randCross(MaleParents,nCrosses)
FemaleF1 = randCross(FemaleParents,nCrosses)
MaleDH = makeDH(MaleF1,nDH)
FemaleDH = makeDH(FemaleF1,nDH)

