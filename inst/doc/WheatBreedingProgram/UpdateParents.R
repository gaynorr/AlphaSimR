#Drop PYT parents
Parents = Parents[1:(nParents-nParentsPYT)]
#Add AYT parents
Parents = c(Parents,selectInd(AYT,nParentsAYT))
#Update parents in EYT with new phenotype
TmpPop = c(EYT2,EYT1)
take = !(Parents@id%in%TmpPop@id) #Parents not in EYT
Parents = c(Parents[take],TmpPop)
rm(TmpPop)
#Select best parents
Parents = selectInd(Parents,nParents-nParentsPYT)
#Add PYT parents
Parents = c(Parents,selectInd(PYT,nParentsPYT))
