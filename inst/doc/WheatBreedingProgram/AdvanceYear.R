#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data
#w is created outside the script

#Year 8
#Release variety

#Year 7
EYT2 = setPheno(EYT1,varE=varE,reps=repEYT,w=w)
EYT2@pheno = (EYT1@pheno+EYT2@pheno)/2

#Year 6
EYT1 = selectInd(AYT,nEYT)
EYT1 = setPheno(EYT1,varE=varE,reps=repEYT,w=w)

#Year 5
AYT = selectInd(PYT,nAYT)
AYT = setPheno(AYT,varE=varE,reps=repAYT,w=w)

#Year 4
PYT = selectWithinFam(HDRW,famMax)
PYT = selectInd(PYT,nPYT)
PYT = setPheno(PYT,varE=varE,reps=repPYT,w=w)

#Year 3
HDRW = setPheno(DH,varE=varE,reps=repHDRW,w=w)

#Year 2
DH = makeDH(F1,nDH)

#Year 1
F1 = randCross(Parents,nCrosses)

