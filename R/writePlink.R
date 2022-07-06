#' @title Writes a Pop-class as PLINK files
#' 
#' @description
#' Writes a Pop-class to PLINK PED and MAP files. The arguments 
#' for this function were chosen for consistency with 
#' \code{\link{RRBLUP2}}. The base pair coordinate will the locus
#' position as stored in AlphaSimR and not an actual base pair 
#' position. This is because AlphaSimR doesn't track base pair 
#' positions, only relative positions for the loci used in the 
#' simulation. 
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param baseName basename for PED and MAP files.
#' @param traits an integer indicating the trait to write, a trait name, or a
#' function of the traits returning a single value.
#' @param use what to use for PLINK's phenotype field. Either phenotypes "pheno", 
#' genetic values "gv", estimated breeding values "ebv", breeding values "bv", 
#' or random values "rand".
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @examples 
#' \dontrun{
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$setSexes(sex="yes_rand")
#' SP$addTraitA(nQtlPerChr=10)
#' SP$addSnpChip(nSnpPerChr=5)
#' SP$setVarE(h2=0.5)
#' 
#' #Create population
#' pop = newPop(rawPop = founderPop)
#' 
#' # Write out PLINK files
#' writePlink(pop, baseName="test")
#' }
#' @export
writePlink = function(pop, baseName, traits=1, use="pheno", 
                      snpChip=1, useQtl=FALSE, simParam=NULL, 
                      ...){
  if(pop@ploidy!=2L){
    stop("writePlink() only supports ploidy=2")
  } 
  
  if(is.null(simParam)){ 
    simParam = get(x="SP", envir=.GlobalEnv)
  }
  
  # Pull "phenotype" data indicated by traits
  y = getResponse(pop=pop, trait=traits, use=use,
                  simParam=simParam, ...)
  
  # Pull QTL/SNP data indicated by snpChip and useQtl
  if(useQtl){
    H1 = pullQtlHaplo(pop=pop, trait=snpChip, 
                      haplo=1, asRaw=TRUE, 
                      simParam=simParam)
    
    H2 = pullQtlHaplo(pop=pop, trait=snpChip, 
                      haplo=2, asRaw=TRUE, 
                      simParam=simParam)
    
    map = getQtlMap(trait=snpChip, simParam=simParam)
  }else{
    H1 = pullSnpHaplo(pop=pop, snpChip=snpChip, 
                      haplo=1, asRaw=TRUE, 
                      simParam=simParam)
    
    H2 = pullSnpHaplo(pop=pop, snpChip=snpChip, 
                      haplo=2, asRaw=TRUE, 
                      simParam=simParam)
    
    map = getSnpMap(snpChip=snpChip, simParam=simParam)
  }
  
  ## Make .ped file
  
  # Format pop data for a .fam (first columns of .ped)
  # Format sex for PLINK
  sex = pop@sex
  sex[which(sex=="H")] = "0"
  sex[which(sex=="M")] = "1"
  sex[which(sex=="F")] = "2"
  
  # Determine within-family ID of father, "0" if not present
  father = pop@id[match(pop@father, pop@id)]
  father[is.na(father)] = "0"
  
  # Determine within-family ID of mother, "0" if not present
  mother = pop@id[match(pop@mother, pop@id)]
  mother[is.na(mother)] = "0"
  
  fam = rbind(rep("1", pop@nInd), # Family ID
              pop@id, # Within-family ID
              father, # Within-family ID of father
              mother, # Within-family ID of mother
              sex, # Sex
              as.character(c(y))) # Phenotype
  
  # Weave together haplotype data for writing to a file with 
  # the write function (requires a transposed matrix)
  H = unname(rbind(t(H1), t(H2)))
  H = H[c(matrix(1:nrow(H),nrow=2,byrow=T)),]
  
  # Free up some memory
  rm(H1, H2) 
  
  # Convert to characters for writing to file
  H = ifelse(H, "2", "1")
  
  # Append .fam
  H = rbind(fam,H)
  
  # Write .ped file
  write(H, file=paste0(baseName,".ped"), ncolumns=nrow(H))
  
  ## Make .map file
  map = rbind(map$chr, # Chromosome
              map$id, # Variant id
              as.character(map$pos*10), # Genetic map position (cM)
              as.character(map$site) # Physical map position
              )
  
  # Write .map file
  write(map, file=paste0(baseName,".map"), ncolumns=nrow(map))
  
  # Don't return anything
  return(invisible())
}
