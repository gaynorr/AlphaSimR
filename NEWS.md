# AlphaSimR 1.2.2

## Minor improvements and fixes

*added `getPed` to quick extract a population's pedigree

*added `getGenMap` to pull a genetic map in data.frame format

# AlphaSimR 1.2.1

## Minor improvements and fixes

*fixed bugs relating to `importData` functions

*fixed `writePlin`k errors and no longer requires equal length chromosomes

# AlphaSimR 1.2.0
  
## New features 

*added `importGenMap` to format genetic maps for AlphaSimR

*added `importInbredGeno` and `importHaplo` to make it easier to create a simulation from external data

*added `importSnpChip`, `importTrait` to `SimParam` to make it easier to manually define traits

*added `pullMarkerGeno` and `pullMarkerHaplo` to make it easier to extract genotypes and haplotypes of specific loci without defining a trait or SNP chip
    
## Minor improvements and fixes

*`reduceGenome`, `mergeGenome` and `doubleGenome` should really now work with pedigree and recombination tracking
    
## Known issues

*`pedigreeCross` fails without an appropriate warning for some incomplete pedigrees

# AlphaSimR 1.1.2
    
## Minor improvements and fixes

*added missing #ifdef _OPENMP to OCS.cpp

# AlphaSimR 1.1.1
    
## Minor improvements and fixes

*removed use of PI variable in C++ code which is compiler specific

# AlphaSimR 1.1.0
  
## New features

*added snpChip argument to `pullIbdHaplo` for backwards compatibility

*exposed internal mixed model solvers

*all selection functions now return a warning when there are not enough individuals
    
## Minor improvements and fixes

*fixed error in `pullIbdHaplo` when chr isn't NULL

*fixed an error with assigning 1 QTL and/or SNP

*changed geno slot from matrix to list to support future RcppArmadillo changes

*`doubleGenome` and `reduceGenome` now work with IBD tracking

# AlphaSimR 1.0.4

## Minor improvements and fixes

*fixed errors in implementation of Gamma Sprinkling model

# AlphaSimR 1.0.3

## Minor improvements and fixes

*fixed formatting error in genetic maps created by runMacs that broke genotype extraction functions

# AlphaSimR 1.0.2

## New features

*added h2 and H2 to `setPhenoGCA`

*`pullGeno` and `pullHaplo` functions now report marker names from the genetic map

# AlphaSimR 1.0.1

## Minor improvements and fixes

*removed lazyData field in DESCRIPTION

# AlphaSimR 1.0.0
  
*AlphaSimR manuscript has been published in G3 (citation added)
  
## New features

*changed to a Gamma Sprinkling model for crossovers, default is still a Gamma model

*change default interfernce parameter (v) to 2.6 to be consistent with the Kosambi mapping function (was 1, consistent with the Haldane mapping function)

*new internal id (iid) that allows user to freely change id slot in populations

*`runMacs2` now adjusts Ne for autopolyploids

*parent populations are now passed to `finalizePop`

*check added that throws an error when use of discontinued "gender" argument is detected

*added experimental `MegaPop-class`

# AlphaSimR 0.13.0

## New features

*references to gender have been changed to the more appropriate terms sex or sexes

*added misc slot to populations

*added `finalizePop` to `SimParam`

*added physical positions to `getSnpMap` and `getQtlMap`

*you can now use h2 and H2 to specify error variance in `setPheno`

*`SimParam$setVarE` now accepts a matrix for varE
    
## Minor improvements and fixes

*fixed a bug in `editGenome` when making multiple edits

*adding merging of centromere vector in `cChr`

# AlphaSimR 0.12.2
  
## New features

*GxE traits now default to random sampling of p-values
  
## Minor improvements and fixes

*fixed a bug in `restrSegSites`

# AlphaSimR 0.12.1

## Minor improvements and fixes

*fixed a bug in selection of segSites

# AlphaSimR 0.12.0
  
## New features

*changed output of `genParam` to match Bulmer, 1976

*nProgeny added to `makeCross` and `makeCross2`

*all `SimParam` documentation is now in `?SimParam`

*non-overlapping QTL and SNP is now the default

*new interface for `restrSegSites` in `SimParam`
  
## Minor improvements and fixes

*fixed subset by id for populations

*fixed major bug in `newMapPop`

# AlphaSimR 0.11.1
  
## New features

*switched to a circular design for the balance option in `randCross` and `randCross2`

*added `reduceGenome` and `doubleGenome` for changing ploidy levels

*added minSnpFreq to SimParam_addSnpChip for any reference population

*the `c` function now merges individuals for MapPop objects (was chromosomes before)

*the `cChr` function new merges chromosomes for MapPop objects 
    
## Minor improvements and fixes

*fixed broken SimParam_addStructuredSnpChip

*removed broken `pullMultipleSnpGeno` and `pullMultipleSnpHaplo`

*fixed broken `writePlink`

# AlphaSimR 0.11.0

## Breaking changes

*rework of `setEBV` (breaks some scripts)

## New features

*genotype data now stored as bits (was bytes)

*implemented a gamma model for crossover interference

*added the mutate function to model random mutations

*added a vignette explaining the biological model for traits

*GS models now handle polyploids

*heterogenous error variance is now optional in GS models (default is homogeneous error)

*improved gene drop functionality of pedigreeCross

*added keepParents option to makeDH and self (indirectly extends `selectFam` and `selectWithinFam`)

*added RRBLUP_SCA2

*set methods for the "show" function when applied to populations
    
## Minor improvements and fixes

*fixed a bug returning the first individual when selecting 0

*fixed error in recombination track when using `makeDH`

*fixed error causing epistatic effects to mask GxE effects

*fixed an error with `pullSegSiteGeno` and `pullSegSiteHaplo` with variable number of sites per chromosome

# AlphaSimR 0.10.0
  
## New features

*added traits with epistasis

*Max number of threads automatically detected

*added RRBLUP_D2

*added version tracking to `SimParam`

*removed `trackHaploPop` (super-ceded by `pullIbdHaplo`)

*added `fastRRBLUP`
  
## Minor improvements and fixes

*fixed faulty double crossover logic

*fixed broken `writePlink`

*fixed broken `pullIbdHaplo`

*`mergePops` no longer assumes diploidy

# AlphaSimR 0.9.0

## New features

*added support for autopolyploids

*added `RRBLUP_GCA2`

*`randCross2` can now "balance" crossing when not using gender
  
## Minor improvements and fixes

*fixed recombination tracking bug in `createDH2`

*removed bug in `setEBV` with append=TRUE

# AlphaSimR 0.8.2

## Minor improvements and fixes

*fixed ambiguous overloading in optimize.cpp

# AlphaSimR 0.8.1

## Minor improvements and fixes

*`setPheno` (not `setPhenoGCA`) passes the number of reps to populations

*fixed bug in `editGenomeTopQtl`

*fixed bug in `RRBLUP_D`

*fixed bug in `resetPop`

*fixed bug in SimParam_rescaleTraits

*removed unimplemented SimParam_restrSnpSites and SimParam_restrQtlSites

*add error message for no traits in `calcGCA`

# AlphaSimR 0.8.0

## New features

*added GxE traits with zero environmental variance

*faster trait scaling

*faster calculation of genetic values

*dsyevr now called via arma_fortran

*added OpenMP support

*parallelized `cross2`

*parallelized `runMacs`

*parallelized calculation of genetic values

*variance calculations now account for inbreeding
    
## Minor improvements and fixes

*fixes for male selection in `selectOP`

# AlphaSimR 0.7.1
  
## Minor improvements and fixes

*add fixEff to `setPhenoGCA`

# AlphaSimR 0.7.0

## New features

*added default `runMacs` option to return all segSites

*added ability to specify separate male and female genetic maps

*`pullGeno` and `pullHaplo` functions can now specify chromosomes

*added `RRBLUP2` for special GS cases

*improved speed by replacing Rcpp random number generators

*changed available MaCS species

*GS functions now use populations directly

*added `pullIbdHaplo`

*added `writePlink`
    
## Minor improvements and fixes

*fixed population sub-setting checks to prevent invalid selections

*fixed slow `calcGCA`

*fixed error in `addTraitAG` preventing multiple traits

*fixed bug with `mergePops` when merging ebv

*fixed bug in `setVarE` when using H2 and multiple traits

# AlphaSimR 0.6.1

## New features

*`selectFam` now handles half-sib families

*`selectWithinFa`m now handles half-sib families
    
## Minor improvements and fixes

*Removed restriction on varE=NULL in `setPhenoGCA`

# AlphaSimR 0.6.0

## New features

*Added NEWS file

*Added `selectOP` to model selection in open pollinating plants

*Added `runMacs2` as a wrapper for `runMacs`
    
## Minor improvements and fixes

*Fixed error when using H2 in SimParam_setVarE
    
