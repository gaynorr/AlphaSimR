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
 
#' @title Import inbred, diploid genotypes
#' 
#' @description
#' Formats the genotypes from inbred, diploid lines 
#' to an AlphaSimR population that can be used to 
#' initialize a simulation. An attempt is made to 
#' automatically detect 0,1,2 or -1,0,1 genotype coding. 
#' Heterozygotes or probabilistic genotypes are allowed, 
#' but will be coerced to the nearest homozygote. Pedigree 
#' information is optional and when provided will be 
#' passed to the population for easier identification 
#' in the simulation.
#' 
#' @param geno a matrix of genotypes
#' @param genMap genetic map as a data.frame. The first  
#' three columns must be: marker name, chromosome, and 
#' map position (Morgans). Marker name and chromosome are 
#' coerced using as.character. See \link{importGenMap}
#' @param ped an optional pedigree for the supplied 
#' genotypes. See details. 
#'
#' @details 
#' The optional pedigree can be a data.frame, matrix or a vector. 
#' If the object is a data.frame or matrix, the first three 
#' columns must include information in the following order: id, 
#' mother, and father. All values are coerced using 
#' as.character. If the object is a vector, it is assumed to only 
#' include the id. In this case, the mother and father will be set 
#' to "0" for all individuals.
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
    if(is.vector(ped)){
      id = as.character(ped)
      stopifnot(length(id)==nrow(geno),
                !any(duplicated(id)))
      mother = father = rep("0", length(id))
    }else{
      id = as.character(ped[,1])
      stopifnot(length(id)==nrow(geno),
                !any(duplicated(id)))
      mother = as.character(ped[,2])
      father = as.character(ped[,3])
    }
  }
  
  genMap = importGenMap(genMap)
  
  # Get marker names
  if(is.data.frame(geno)){
    geno = as.matrix(geno)
  }
  markerName = colnames(geno)

  # Check marker coding and convert to haplotypes
  if(is.raw(geno)){
    geno[geno==as.raw(1)] = as.raw(0) # For consistency with round
    geno[geno==as.raw(2)] = as.raw(1)
  }else{
    minGeno = min(geno)
    maxGeno = max(geno)
    stopifnot(minGeno >= (-1-1e-8) )
    
    if(minGeno < (0-1e-8) ){
      # Suspect -1,0,1 coding
      stopifnot(maxGeno <= (1+1e-8) )
      # Converting to 0,1 haplotypes with hard thresholds
      geno = matrix(as.raw( round( (geno+1)/2 ) ),
                    ncol=ncol(geno))
    }else{
      # Suspect 0,1,2 coding
      stopifnot(maxGeno <= (2+1e-8) )
      # Converting to 0,1 haplotypes with hard thresholds
      geno = matrix(as.raw( round( geno/2 ) ),
                    ncol=ncol(geno))
    }
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


#' @title Import haplotypes
#' 
#' @description
#' Formats haplotype in a matrix format to an 
#' AlphaSimR population that can be used to 
#' initialize a simulation. This function serves 
#' as wrapper for \code{\link{newMapPop}} that 
#' utilizes a more user friendly input format.
#' 
#' @param haplo a matrix of haplotypes
#' @param genMap genetic map as a data.frame. The first  
#' three columns must be: marker name, chromosome, and 
#' map position (Morgans). Marker name and chromosome are 
#' coerced using as.character. See \code{\link{importGenMap}}
#' @param ploidy ploidy level of the organism
#' @param ped an optional pedigree for the supplied 
#' genotypes. See details. 
#' 
#' @details 
#' The optional pedigree can be a data.frame, matrix or a vector. 
#' If the object is a data.frame or matrix, the first three 
#' columns must include information in the following order: id, 
#' mother, and father. All values are coerced using 
#' as.character. If the object is a vector, it is assumed to only 
#' include the id. In this case, the mother and father will be set 
#' to "0" for all individuals.
#'
#' @return a \code{\link{MapPop-class}} if ped is NULL,
#' otherwise a \code{\link{NamedMapPop-class}}
#' 
#' @examples 
#' haplo = rbind(c(1,1,0,1,0),
#'               c(1,1,0,1,0),
#'               c(0,1,1,0,0),
#'               c(0,1,1,0,0))
#' colnames(haplo) = letters[1:5]
#' 
#' genMap = data.frame(markerName=letters[1:5],
#'                     chromosome=c(1,1,1,2,2),
#'                     position=c(0,0.5,1,0.15,0.4))
#' 
#' ped = data.frame(id=c("a","b"),
#'                  mother=c(0,0),
#'                  father=c(0,0))
#' 
#' founderPop = importHaplo(haplo=haplo, 
#'                          genMap=genMap,
#'                          ploidy=2L,
#'                          ped=ped)
#' 
#' @export
importHaplo = function(haplo, genMap, ploidy=2L, ped=NULL){
  # Extract pedigree, if supplied
  if(!is.null(ped)){
    if(is.vector(ped)){
      id = as.character(ped)
      stopifnot(length(id)==(nrow(haplo)/ploidy),
                !any(duplicated(id)))
      mother = father = rep("0", length(id))
    }else{
      id = as.character(ped[,1])
      stopifnot(length(id)==(nrow(haplo)/ploidy),
                !any(duplicated(id)))
      mother = as.character(ped[,2])
      father = as.character(ped[,3])
    }
  }
  
  genMap = importGenMap(genMap)
  
  # Get marker names
  if(is.data.frame(haplo)){
    haplo = as.matrix(haplo)
  }
  markerName = colnames(haplo)
  
  # Convert haplotypes to raw
  haplo = matrix(as.raw(haplo), ncol=ncol(haplo))
  stopifnot(haplo==as.raw(0) | haplo==as.raw(1))
  
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
    haplotypes[[i]] = haplo[,take,drop=FALSE]
  }
  
  founderPop = newMapPop(genMap=genMap,
                         haplotypes=haplotypes,
                         ploidy=ploidy)
  
  if(!is.null(ped)){
    founderPop = new("NamedMapPop",
                     id=id,
                     mother=mother,
                     father=father,
                     founderPop)
  }
  
  return(founderPop)
}

