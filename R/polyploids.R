#' @title Create individuals with reduced ploidy
#' 
#' @description Creates new individuals from gametes. This function 
#' was created to model the creation of diploid potatoes from 
#' tetraploid potatoes. It can be used on any population with an 
#' even ploidy level. The newly created individuals will have half 
#' the ploidy level of the originals. The reduction can occur with 
#' or without genetic recombination.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nProgeny total number of progeny per individual
#' @param useFemale should female recombination rates be used. 
#' @param keepParents should previous parents be used for mother and 
#' father. 
#' @param simRecomb should genetic recombination be modeled.
#' @param simParam an object of 'SimParam' class
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Create individuals with reduced ploidy
#' pop2 = reduceGenome(pop, simParam=SP)
#' 
#' @export
reduceGenome = function(pop,nProgeny=1,useFemale=TRUE,keepParents=TRUE,
                        simRecomb=TRUE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(pop@ploidy%%2L){
    stop("You cannot reduce aneuploids")
  }
  
  if(simRecomb){
    if(useFemale){
      map = simParam$femaleMap
    }else{
      map = simParam$maleMap
    }
  }else{
    # Create dummy map with zero genetic distance
    map = vector("list",pop@nChr)
    for(i in 1:pop@nChr){
      map[[i]] = rep(0,pop@nLoci[i])
    }
    map = as.matrix(map)
  }
  
  tmp = createReducedGenome(pop@geno, nProgeny,
                            map,
                            simParam$v,
                            simParam$p,
                            simParam$isTrackRec,
                            pop@ploidy,
                            simParam$femaleCentromere,
                            simParam$quadProb,
                            simParam$nThreads)
  dim(tmp$geno) = NULL 
  
  rPop = new("RawPop",
             nInd=as.integer(pop@nInd*nProgeny),
             nChr=pop@nChr,
             ploidy=as.integer(pop@ploidy/2),
             nLoci=pop@nLoci,
             geno=tmp$geno)
  
  if(simParam$isTrackRec){
    hist = tmp$recHist
  }else{
    hist = NULL
  }
  
  if(keepParents){
    return(newPop(rawPop=rPop,
                  mother=rep(pop@mother,each=nProgeny),
                  father=rep(pop@father,each=nProgeny),
                  simParam=simParam,
                  iMother=rep(pop@iid,each=nProgeny),
                  iFather=rep(pop@iid,each=nProgeny),
                  femaleParentPop=pop,
                  maleParentPop=pop,
                  hist=hist
    ))
  }else{
    return(newPop(rawPop=rPop,
                  mother=rep(pop@id,each=nProgeny),
                  father=rep(pop@id,each=nProgeny),
                  simParam=simParam,
                  iMother=rep(pop@iid,each=nProgeny),
                  iFather=rep(pop@iid,each=nProgeny),
                  femaleParentPop=pop,
                  maleParentPop=pop,
                  hist=hist
    ))
  }
}

#' @title Double the ploidy of individuals
#' 
#' @description Creates new individuals with twice the ploidy. 
#' This function was created to model the formation of tetraploid 
#' potatoes from diploid potatoes. This function will work on any 
#' population.
#' 
#' @param pop an object of 'Pop' superclass
#' @param keepParents should previous parents be used for mother and 
#' father. 
#' @param simParam an object of 'SimParam' class
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Create individuals with doubled ploidy
#' pop2 = doubleGenome(pop, simParam=SP)
#' 
#' @export
doubleGenome = function(pop, keepParents=TRUE,
                        simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  geno = pop@geno
  for(i in 1:pop@nChr){
    geno[[i]] = geno[[i]][,rep(1:pop@ploidy,each=2),,drop=FALSE]
  }
  
  rPop = new("RawPop",
             nInd=as.integer(pop@nInd),
             nChr=pop@nChr,
             ploidy=2L*pop@ploidy,
             nLoci=pop@nLoci,
             geno=geno)
  
  if(simParam$isTrackRec){
    # Match haplotypes according to original ploidy
    hist = vector("list", pop@ploidy)
    for(i in 1:pop@ploidy){
      hist[[i]] = cbind(i, 1L)
    }
    
    # Double ploidy
    hist = rep(hist, each=2)
    
    # Rep for chromosomes and individuals
    hist = rep(list(hist), pop@nChr)
    hist = rep(list(hist), pop@nInd)
  }else{
    hist = NULL
  }
  
  if(keepParents){
    return(newPop(rawPop=rPop,
                  mother=pop@mother,
                  father=pop@father,
                  simParam=simParam,
                  iMother=pop@iid,
                  iFather=pop@iid,
                  femaleParentPop=pop,
                  maleParentPop=pop,
                  hist=hist
    ))
  }else{
    return(newPop(rawPop=rPop,
                  mother=pop@id,
                  father=pop@id,
                  simParam=simParam,
                  iMother=pop@iid,
                  iFather=pop@iid,
                  femaleParentPop=pop,
                  maleParentPop=pop,
                  hist=hist
    ))
  }
}

#' @title Combine genomes of individuals
#'
#' @description
#' This function is designed to model the pairing of gametes. The male
#' and female individuals are treated as gametes, so the ploidy of newly 
#' created individuals will be the sum of it parents.
#'
#' @param females an object of \code{\link{Pop-class}} for female parents.
#' @param males an object of \code{\link{Pop-class}} for male parents.
#' @param crossPlan a matrix with two column representing
#' female and male parents. Either integers for the position in
#' population or character strings for the IDs.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Cross individual 1 with individual 10
#' crossPlan = matrix(c(1,10), nrow=1, ncol=2)
#' pop2 = mergeGenome(pop, pop, crossPlan, simParam=SP)
#'
#' @export
mergeGenome = function(females,males,crossPlan,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],females@id),
                      match(crossPlan[,2],males@id))
    if(any(is.na(crossPlan))){
      stop("Failed to match supplied IDs")
    }
  }
  
  if((max(crossPlan[,1])>nInd(females)) | 
     (max(crossPlan[,2])>nInd(males)) |
     (min(crossPlan)<1L)){
    stop("Invalid crossPlan")
  }
  
  # Merge genotype data
  geno = vector("list", females@nChr)
  for(i in 1:females@nChr){
    geno[[i]] = array(as.raw(0), 
                      dim = c(dim(females@geno[[i]])[1],
                              females@ploidy+males@ploidy,
                              nrow(crossPlan)))
    for(j in 1:nrow(crossPlan)){
      # Add female gametes
      geno[[i]][,1:females@ploidy,j] = 
        females@geno[[i]][,,crossPlan[j,1]]
      
      # Add male gametes
      geno[[i]][,(females@ploidy+1):(females@ploidy+males@ploidy),j] = 
        males@geno[[i]][,,crossPlan[j,2]]
    }
  }
  
  rPop = new("RawPop",
             nInd=as.integer(nrow(crossPlan)),
             nChr=females@nChr,
             ploidy=females@ploidy+males@ploidy,
             nLoci=females@nLoci,
             geno=geno)
  
  
  if(simParam$isTrackRec){
    # Create history for haplotypes
    hist = vector("list", rPop@ploidy)
    
    # Add female contribution
    for(i in 1:females@ploidy){
      hist[[i]] = cbind(i, 1L)
    }
    
    # Add male contribution
    for(i in 1:males@ploidy){
      hist[[i+females@ploidy]] = cbind(i, 1L)
    }
    
    # Rep for chromosomes and individuals
    hist = rep(list(hist), rPop@nChr)
    hist = rep(list(hist), rPop@nInd)
  }else{
    hist = NULL
  }
  
  return(newPop(rawPop=rPop,
                mother=females@id[crossPlan[,1]],
                father=males@id[crossPlan[,2]],
                simParam=simParam,
                iMother=females@iid[crossPlan[,1]],
                iFather=males@iid[crossPlan[,2]],
                femaleParentPop=females,
                maleParentPop=males,
                hist=hist
  ))
}

