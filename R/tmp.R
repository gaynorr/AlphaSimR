#' @title Add genetic values
#' 
#' @description Promotes class 'Pop' to class 'TraitPop'
#' 
#' @param pop an object of class 'Pop'
#' @param simParam an object of class 'SimParam'
#' 
#' @export
addGv = function(pop, simParam=SIMPARAM){
  stopifnot(class(pop)=="Pop")
  gv = lapply(simParam@traits,getGv,pop=pop,w=0)
  gv = do.call("cbind",gv)
  pheno = matrix(NA_real_,nrow=nrow(gv),ncol=ncol(gv))
  pop = new("TraitPop",pop,gv=gv,pheno=pheno,
            nTraits=simParam@nTraits)
  return(pop)
}

#' @title Add pedigree
#' 
#' @description Promotes class 'Pop' or 'TraitPop' to 'PedPop'
#' 
#' @param pop an object of class 'Pop' or 'TraitPop'
#' @param id a unique name for each individual
#' @param par1 the first/female parent for each individual
#' @param par2 the second/male parent for each individual
#' @param simParam an object of class 'SimParam'
#' 
#' @export
addPed = function(pop, id, par1, par2, simParam=SIMPARAM){
  if(class(pop)=="Pop"){
    pop = addGv(pop,simParam=simParam)
  }
  stopifnot(class(pop)=="TraitPop")
  pop = new("PedPop",pop,id=as.character(id),
            par1=as.character(par1),par2=as.character(par2))
  return(pop)
}