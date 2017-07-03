#Sets new parents for inbreds
MaleParents = c(MaleElite,MaleInbredYT5,MaleInbredYT4)
nSelYT2 = nParents - MaleParents@nInd
if(nSelYT2>0){
  MaleParents = c(MaleParents,selectInd(MaleYT2,nSelYT2))
}else if(nSelYT2<0){
  MaleParents = MaleParents[1:nParents]
}

FemaleParents = c(FemaleElite,FemaleInbredYT5,FemaleInbredYT4)
nSelYT2 = nParents - FemaleParents@nInd
if(nSelYT2>0){
  FemaleParents = c(FemaleParents,selectInd(FemaleYT2,nSelYT2))
}else if(nSelYT2<0){
  FemaleParents = FemaleParents[1:nParents]
}
