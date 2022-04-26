#' @title Import genetic map
#' 
#' @description
#' Formats a genetic map stored in a data.frame to
#' AlphaSimR's internal format. Map positions must be 
#' in Morgans. 
#' 
#' @param genMap genetic map as a data.frame. The first  
#' three columns must be: marker name, chromosome, and 
#' map position (Morgans). Marker name and chromosome are 
#' coerced using as.character.
#'
#' @return a list of named vectors
#' 
#' @examples 
#' genMap = data.frame(markerName=letters[1:5],
#'                     chromosome=c(1,1,1,2,2),
#'                     position=c(0,0.5,1,0.15,0.4))
#' 
#' asrMap = importGenMap(genMap=genMap)
#' 
#' str(asrMap)
#' 
#' @export
importGenMap = function(genMap){
  # Convert data type
  markerName = as.character(genMap[,1])
  chromosome = as.character(genMap[,2])
  position = as.numeric(genMap[,3])

  # Create list for map
  uniqueChr = unique(chromosome)
  genMap = vector("list", length=length(uniqueChr))
  names(genMap) = uniqueChr

  # Iterate through chromosomes
  for(i in 1:length(uniqueChr)){

    take = (chromosome==uniqueChr[i])
    tmpPos = position[take]
    tmpName = markerName[take]

    # Order and name
    take = order(tmpPos, decreasing=FALSE)
    tmpPos = tmpPos[take]
    names(tmpPos) = tmpName[take]

    genMap[[uniqueChr[i] ]] = tmpPos - tmpPos[1]
  }

  return(genMap)
}
 
#' @title Import diploid inbred genotypes
#' 
#' @description
#' Formats a genetic map stored in a data.frame to
#' AlphaSimR's internal format. Map positions must be 
#' in Morgans. 
#' 
#' @param geno a matrix of genotypes
#' @param genMap genetic map as a data.frame. The first  
#' three columns must be: marker name, chromosome, and 
#' map position (Morgans). Marker name and chromosome are 
#' coerced using as.character. See \link{importGenMap}
#' @param ped an optional pedigree for the supplied 
#' genotypes. The first three columns must be: id, 
#' mother, and father. All values are coerced using 
#' as.character.
#'
#' @return a \code{\link{MapPop-class}} if ped is NULL,
#' otherwise a \code{\link{NamedMapPop-class}}
#' 
#' @examples 
#' geno = rbind(c(2,2,0,2,0),
#'              c(0,2,2,0,0))
#' colnames(geno) = letters[1:5]
#' 
#' genMap = data.frame(markerName=letters[1:5],
#'                     chromosome=c(1,1,1,2,2),
#'                     position=c(0,0.5,1,0.15,0.4))
#' 
#' ped = data.frame(id=c("a","b"),
#'                  mother=c(0,0),
#'                  father=c(0,0))
#' 
#' founderPop = importInbredGeno(geno=geno,
#'                               genMap=genMap,
#'                               ped=ped)
#' 
#' @export
importInbredGeno = function(geno, genMap, ped=NULL){
  # Extract pedigree, if supplied
  if(!is.null(ped)){
    id = as.character(ped[,1])
    stopifnot(length(id)==nrow(geno),
              !any(duplicated(id)))
    mother = as.character(ped[,2])
    father = as.character(ped[,3])
  }
  
  genMap = importGenMap(genMap)
  
  # Get marker names
  if(is.data.frame(geno)){
    geno = as.matrix(geno)
  }
  markerName = colnames(geno)

  # Check marker coding and convert to haplotypes
  minGeno = min(geno)
  maxGeno = max(geno)
  stopifnot(minGeno >= (-1-1e-8) )

  if(minGeno < (0-1e-8) ){
    # Suspect -1,0,1 coding
    stopifnot(maxGeno <= (1+1e-8) )
    # Converting to 0,1 haplotypes with hard thresholds
    geno = matrix(as.integer( round( (geno+1)/2 ) ),
                  ncol=ncol(geno))
  }else{
    # Suspect 0,1,2 coding
    stopifnot(maxGeno <= (2+1e-8) )
    # Converting to 0,1 haplotypes with hard thresholds
    geno = matrix(as.integer( round( geno/2 ) ),
                  ncol=ncol(geno))
  }

  # Create haplotype list
  haplotypes = vector("list", length=length(genMap))

  # Order haplotypes by chromosome
  for(i in 1:length(genMap)){
    mapMarkers = names(genMap[[i]])
    take = match(mapMarkers, markerName)
    if(any(is.na(take))){
      genMap[[i]] = genMap[[i]][is.na(take)]
      stopifnot(length(genMap[[i]]) >= 1L)
      genMap[[i]] = genMap[[i]] - genMap[[i]]-genMap[[i]][1]
      take = na.omit(take)
    }
    haplotypes[[i]] = geno[,take]
  }

  founderPop = newMapPop(genMap=genMap,
                         haplotypes=haplotypes,
                         inbred=TRUE)

  if(!is.null(ped)){
    founderPop = new("NamedMapPop",
                     id=id,
                     mother=mother,
                     father=father,
                     founderPop)
  }

  return(founderPop)
}

# # importHaplo
# 
# 
# 
# # Convert to SimParam functions
# importSnpChip = function(markerName,
#                          chipName=NULL){
#   
# }
# 
# importTraitA = function(markerName, 
#                         addEff, 
#                         intercept=0, 
#                         traitName=NULL, 
#                         simParam=NULL){
#   if(is.null(simParam)){
#     simParam = get("SP",envir=.GlobalEnv)
#   }
#   
#   stopifnot(length(markerName)==length(addEff))
#   
#   # Extract genetic map and check if names are in map
#   genMap = simParam$genMap
#   genMapMarkerNames = unlist(lapply(genMap, names))
#   stopifnot(all(markerName%in%genMapMarkerNames))
#   
#   # Create trait variables
#   lociPerChr = integer(length(genMap))
#   addEffList = lociLoc = vector("list", length(genMap))
#   
#   # Loop through chromosomes
#   for(i in 1:length(genMap)){
#     
#     # Initialize variables
#     addEffList[[i]] = numeric()
#     lociLoc[[i]] = integer()
#     
#     # Find matches if they exist
#     take = match(names(genMap[[i]]), markerName)
#     lociPerChr[i] = length(na.omit(take))
#     if(lociPerChr[i]>0L){
#       lociLoc[[i]] = which(!is.na(take))
#       addEffList[[i]] = addEff[na.omit(take)]
#     }
#   }
#   addEff = unlist(addEffList)
#   lociLoc = unlist(lociLoc)
#   nLoci = sum(lociPerChr)
#   
#   # Create Trait
#   trait = new("TraitA", 
#               addEff=addEff,
#               intercept=as.numeric(intercept),
#               nLoci=nLoci,
#               lociPerChr=lociPerChr,
#               lociLoc=lociLoc)
#   
#   # Add trait to simParam
#   simParam$manAddTrait(trait)
#   
#   # Return nothing
#   invisible(NULL)
# }
# 
# importTraitAD = function(markerName, 
#                         addEff, 
#                         domEff,
#                         intercept=0, 
#                         traitName=NULL, 
#                         simParam=NULL){
#   if(is.null(simParam)){
#     simParam = get("SP",envir=.GlobalEnv)
#   }
#   
#   stopifnot(length(markerName)==length(addEff))
#   
#   # Extract genetic map and check if names are in map
#   genMap = simParam$genMap
#   genMapMarkerNames = unlist(lapply(genMap, names))
#   stopifnot(all(markerName%in%genMapMarkerNames))
#   
#   # Create trait variables
#   lociPerChr = integer(length(genMap))
#   addEffList = domEffList = lociLoc = 
#     vector("list", length(genMap))
#   
#   # Loop through chromosomes
#   for(i in 1:length(genMap)){
#     
#     # Initialize variables
#     addEffList[[i]] = domEffList[[i]] = numeric()
#     lociLoc[[i]] = integer()
#     
#     # Find matches if they exist
#     take = match(names(genMap[[i]]), markerName)
#     lociPerChr[i] = length(na.omit(take))
#     if(lociPerChr[i]>0L){
#       lociLoc[[i]] = which(!is.na(take))
#       addEffList[[i]] = addEff[na.omit(take)]
#       domEffList[[i]] = domEff[na.omit(take)]
#     }
#   }
#   addEff = unlist(addEffList)
#   domEff = unlist(domEffList)
#   lociLoc = unlist(lociLoc)
#   nLoci = sum(lociPerChr)
#   
#   # Create Trait
#   trait = new("TraitAD", 
#               addEff=addEff,
#               domEff=domEff,
#               intercept=as.numeric(intercept),
#               nLoci=nLoci,
#               lociPerChr=lociPerChr,
#               lociLoc=lociLoc)
#   
#   # Add trait to simParam
#   simParam$manAddTrait(trait)
#   
#   # Return nothing
#   invisible(NULL)
# }
