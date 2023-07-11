# AlphaSimR TBD

* Changed 'MegaPop' to 'MultiPop'

# AlphaSimR 1.4.2

* updated MaCS citation to https site

# AlphaSimR 1.4.1

* Changed citation to use `bibentry` instead of `citEntry`

# AlphaSimR 1.4.0

*fixed a bug in IBD tracking

*add `setFounderHap` to SimParam for applying custom haplotypes to founders

*added `addSnpChipByName` to SimParam for defining SNP chips by marker names

# AlphaSimR 1.3.4

*changed C++ using `sprintf` to use `snprintf`

# AlphaSimR 1.3.3

*fixed bug in calculation of genic variance

*fixed `importHaplo` not passing ploidy to `newMapPop`

*fixed bug with correlated error variances

# AlphaSimR 1.3.2

*fixed column name bug with multiple traits in `setEBV`

*fixed CTD caused by `runMacs` when too many segSites are requested

*fixed missing names in GV when using `resetPop`

*fixed bug in `importTrait`

*`popVar` now deals with matrices having 1 row

# AlphaSimR 1.3.1

*updated link to Gaynor, 2017

# AlphaSimR 1.3.0

*added ability to exclude loci by name in `SimParam$restrSegSites`

*`pullMarkerGeno` and `pullMarkerHaplo` now work with a MapPop class

*added `setMarkerHaplo` to manually change genotypes in a Pop or MapPop

*added `addSegSite` for manually adding segregating sites to a MapPop class

*`simParam$setCorE` has been deprecated in favor of a corE argument in `simParam$setVarE`

*`setPheno` now takes corE as an argument

*`setPheno` now allows the user to set phenotypes for a subset of traits

*add `newEmptyPop` to create populations with zero individuals

*removed reps slot from populations and heterogeneous residual variance GS models

*added h2, H2, and corE to `setPhenoGCA` and `setPhenoProgTest`

*the "EUROPEAN" species history was removed from `runMacs` due to lengthy runtime

# AlphaSimR 1.2.2

*added `getPed` to quick extract a population's pedigree

*added `getGenMap` to pull a genetic map in data.frame format

# AlphaSimR 1.2.1

*fixed bugs relating to `importData` functions

*fixed `writePlink` errors and no longer requires equal length chromosomes

# AlphaSimR 1.2.0

*added `importGenMap` to format genetic maps for AlphaSimR

*added `importInbredGeno` and `importHaplo` to make it easier to create a simulation from external data

*added `importSnpChip`, `importTrait` to `SimParam` to make it easier to manually define traits

*added `pullMarkerGeno` and `pullMarkerHaplo` to make it easier to extract genotypes and haplotypes of specific loci without defining a trait or SNP chip

*`reduceGenome`, `mergeGenome` and `doubleGenome` should really now work with pedigree and recombination tracking

# AlphaSimR 1.1.2

*added missing #ifdef _OPENMP to OCS.cpp

# AlphaSimR 1.1.1

*removed use of PI variable in C++ code due to it being compiler specific

# AlphaSimR 1.1.0

*added snpChip argument to `pullIbdHaplo` for backwards compatibility

*exposed internal mixed model solvers

*all selection functions now return a warning when there are not enough individuals

*fixed error in `pullIbdHaplo` when chr isn't NULL

*fixed an error with assigning 1 QTL and/or SNP

*changed geno slot from matrix to list to support future RcppArmadillo changes

*`doubleGenome` and `reduceGenome` now work with IBD tracking

# AlphaSimR 1.0.4

*fixed errors in implementation of Gamma Sprinkling model

# AlphaSimR 1.0.3

*fixed formatting error in genetic maps created by runMacs that broke genotype extraction functions

# AlphaSimR 1.0.2

*added h2 and H2 to `setPhenoGCA`

*`pullGeno` and `pullHaplo` functions now report marker names from the genetic map

# AlphaSimR 1.0.1

*removed lazyData field in DESCRIPTION

# AlphaSimR 1.0.0
  
*AlphaSimR manuscript has been published in G3 (citation added)

*changed to a Gamma Sprinkling model for crossovers, default is still a Gamma model

*change default interference parameter (v) to 2.6 to be consistent with the Kosambi mapping function (was 1, consistent with the Haldane mapping function)

*new internal id (iid) that allows user to freely change id slot in populations

*`runMacs2` now adjusts Ne for autopolyploids

*parent populations are now passed to `finalizePop`

*check added that throws an error when use of discontinued "gender" argument is detected

*added experimental `MegaPop-class`

# AlphaSimR 0.13.0

*references to gender have been changed to the more appropriate terms sex or sexes

*added misc slot to populations

*added `finalizePop` to `SimParam`

*added physical positions to `getSnpMap` and `getQtlMap`

*you can now use h2 and H2 to specify error variance in `setPheno`

*`SimParam$setVarE` now accepts a matrix for varE

*fixed a bug in `editGenome` when making multiple edits

*adding merging of centromere vector in `cChr`

# AlphaSimR 0.12.2

*GxE traits now default to random sampling of p-values

*fixed a bug in `restrSegSites`

# AlphaSimR 0.12.1

*fixed a bug in selection of segSites

# AlphaSimR 0.12.0

*changed output of `genParam` to match Bulmer, 1976

*nProgeny added to `makeCross` and `makeCross2`

*all `SimParam` documentation is now in `?SimParam`

*non-overlapping QTL and SNP is now the default

*new interface for `restrSegSites` in `SimParam`

*fixed subset by id for populations

*fixed major bug in `newMapPop`

# AlphaSimR 0.11.1

*switched to a circular design for the balance option in `randCross` and `randCross2`

*added `reduceGenome` and `doubleGenome` for changing ploidy levels

*added minSnpFreq to SimParam_addSnpChip for any reference population

*the `c` function now merges individuals for MapPop objects (was chromosomes before)

*the `cChr` function new merges chromosomes for MapPop objects 

*fixed broken SimParam_addStructuredSnpChip

*removed broken `pullMultipleSnpGeno` and `pullMultipleSnpHaplo`

*fixed broken `writePlink`

# AlphaSimR 0.11.0

*rework of `setEBV` (breaks some scripts)

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

*fixed a bug returning the first individual when selecting 0

*fixed error in recombination track when using `makeDH`

*fixed error causing epistatic effects to mask GxE effects

*fixed an error with `pullSegSiteGeno` and `pullSegSiteHaplo` with variable number of sites per chromosome

# AlphaSimR 0.10.0

*added traits with epistasis

*Max number of threads automatically detected

*added RRBLUP_D2

*added version tracking to `SimParam`

*removed `trackHaploPop` (super-ceded by `pullIbdHaplo`)

*added `fastRRBLUP`

*fixed faulty double crossover logic

*fixed broken `writePlink`

*fixed broken `pullIbdHaplo`

*`mergePops` no longer assumes diploidy

# AlphaSimR 0.9.0

*added support for autopolyploids

*added `RRBLUP_GCA2`

*`randCross2` can now "balance" crossing when not using gender

*fixed recombination tracking bug in `createDH2`

*removed bug in `setEBV` with append=TRUE

# AlphaSimR 0.8.2

*fixed ambiguous overloading in optimize.cpp

# AlphaSimR 0.8.1

*`setPheno` (not `setPhenoGCA`) passes the number of reps to populations

*fixed bug in `editGenomeTopQtl`

*fixed bug in `RRBLUP_D`

*fixed bug in `resetPop`

*fixed bug in SimParam_rescaleTraits

*removed unimplemented SimParam_restrSnpSites and SimParam_restrQtlSites

*add error message for no traits in `calcGCA`

# AlphaSimR 0.8.0

*added GxE traits with zero environmental variance

*faster trait scaling

*faster calculation of genetic values

*dsyevr now called via arma_fortran

*added OpenMP support

*parallelized `cross2`

*parallelized `runMacs`

*parallelized calculation of genetic values

*variance calculations now account for inbreeding

*fixes for male selection in `selectOP`

# AlphaSimR 0.7.1

*add fixEff to `setPhenoGCA`

# AlphaSimR 0.7.0

*added default `runMacs` option to return all segSites

*added ability to specify separate male and female genetic maps

*`pullGeno` and `pullHaplo` functions can now specify chromosomes

*added `RRBLUP2` for special GS cases

*improved speed by replacing Rcpp random number generators

*changed available MaCS species

*GS functions now use populations directly

*added `pullIbdHaplo`

*added `writePlink`

*fixed population sub-setting checks to prevent invalid selections

*fixed slow `calcGCA`

*fixed error in `addTraitAG` preventing multiple traits

*fixed bug with `mergePops` when merging ebv

*fixed bug in `setVarE` when using H2 and multiple traits

# AlphaSimR 0.6.1

*`selectFam` now handles half-sib families

*`selectWithinFa`m now handles half-sib families

*Removed restriction on varE=NULL in `setPhenoGCA`

# AlphaSimR 0.6.0

*Added NEWS file

*Added `selectOP` to model selection in open pollinating plants

*Added `runMacs2` as a wrapper for `runMacs`

*Fixed error when using H2 in SimParam_setVarE
    
