#' @title Create individuals with reduced ploidy
#' 
#' @description Creates new individuals from gametes. This function 
#' was created to model the creation of diploid potatoes from 
#' tetraploid potatoes. It can be used on any population with an 
#' even ploidy level. The newly created individuals will have half 
#' the ploidy level of the originals and they will first undergo
#' a single round of meiosis.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nProgeny total number of progeny per individual
#' @param useFemale should female recombination rates be used. 
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
#' #Create individuals with reduced ploidy
#' pop2 = reduceGenome(pop, simParam=SP)
#' 
#' @export
reduceGenome = function(pop,nProgeny=1,useFemale=TRUE,keepParents=TRUE,
                        simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(pop@ploidy%%2L){
    stop("You cannot reduce aneuploids")
  }
  if(useFemale){
    tmp = createReducedGenome(pop@geno, nProgeny,
                              simParam$femaleMap,
                              simParam$v,
                              simParam$isTrackRec,
                              pop@ploidy,
                              simParam$femaleCentromere,
                              simParam$quadProb,
                              simParam$nThreads)
  }else{
    tmp = createReducedGenome(pop@geno, nProgeny,
                              simParam$maleMap,
                              simParam$v,
                              simParam$isTrackRec,
                              pop@ploidy,
                              simParam$maleCentromere,
                              simParam$quadProb,
                              simParam$nThreads)
  }
  rPop = new("RawPop",
             nInd=as.integer(pop@nInd*nProgeny),
             nChr=pop@nChr,
             ploidy=as.integer(pop@ploidy/2),
             nLoci=pop@nLoci,
             geno=tmp$geno)
  if(simParam$isTrackRec){
    simParam$addToRec(tmp$recHist)
  }
  if(keepParents){
    return(newPop(rawPop=rPop,
                  mother=rep(pop@id,each=nProgeny),
                  father=rep(pop@id,each=nProgeny),
                  origM=rep(pop@mother,each=nProgeny),
                  origF=rep(pop@father,each=nProgeny),
                  isDH=FALSE,
                  simParam=simParam))
  }else{
    return(newPop(rawPop=rPop,
                  mother=rep(pop@id,each=nProgeny),
                  father=rep(pop@id,each=nProgeny),
                  isDH=FALSE,
                  simParam=simParam))
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
    geno[[i]] = geno[[i]][,rep(1:pop@ploidy,each=2),]
  }
  rPop = new("RawPop",
             nInd=as.integer(pop@nInd),
             nChr=pop@nChr,
             ploidy=2L*pop@ploidy,
             nLoci=pop@nLoci,
             geno=geno)
  if(keepParents){
    origM=pop@mother
    origF=pop@father
  }else{
    origM=pop@id
    origF=pop@id
  }
  if(simParam$isTrackPed){
    # Extract actual parents
    ped = simParam$ped
    id = as.numeric(pop@id)
    mother = ped[id,1]
    father = ped[id,2]
  }else{
    # Provide arbitrary parents (not actually used)
    mother = origM
    father = origF
  }
  if(simParam$isTrackRec){
    # Duplicate recombination histories
    oldHist = simParam$recHist
    newHist = vector("list", 2*pop@ploidy)
    newHist = rep(list(newHist), pop@nChr)
    newHist = rep(list(newHist), pop@nInd)
    for(i in 1:pop@nInd){
      for(j in 1:pop@nChr){
        k = 0
        for(l in 1:pop@ploidy){
          for(m in 1:2){
            k = k+1
            newHist[[i]][[j]][[k]] = 
              oldHist[[as.numeric(id[i])]][[j]][[l]]
          }
        }
      }
    }
    simParam$addToRec(newHist)
  }
  return(newPop(rawPop=rPop,
                mother=mother,
                father=father,
                origM=origM,
                origF=origF,
                isDH=TRUE,
                simParam=simParam))
}

