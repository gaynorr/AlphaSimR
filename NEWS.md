# AlphaSimR 2.1.0.9016

* Added `restrInbr` to `selectCross`, which selects parents with a restriction on the expected increase in fixation, an approximation to optimal contribution selection with equal contributions. Fixation is measured at a SNP chip or at a trait's QTL, chosen with `snpChip` and `useQtl` as in `RRBLUP`. The target is set with `inbrTarget` and `inbrType`, either relative to the current population or as an absolute maximum. The default is `FALSE`, which leaves `selectCross` unchanged. The new function `selectOCS` makes the same selection without crossing, so the selected individuals can be crossed with a custom crossing plan, for example with `makeCross`.

* Fixed several places where invalid input could read or write memory outside a population's genotypes, or outside internal MaCS structures, instead of giving an error. These now stop with an error:
  * `haplo` in `pullSnpHaplo`, `pullQtlHaplo`, `pullSegSiteHaplo` and `pullMarkerHaplo` must be `"all"` or a single whole number from 1 to the ploidy level.
  * A `LociMap`, such as a trait passed to `SimParam$manAddTrait`, whose loci fall outside their chromosome.
  * `parents` in `self` and `crossPlan` in `hybridCross` with `returnHybridPop=TRUE` must index individuals in the population. `hybridCross` now also accepts a `crossPlan` of IDs for `returnHybridPop=TRUE`, as `makeCross2` does.
  * `segSites` in `runMacs` must have length 1 or `nChr`.
  * Population IDs in a `runMacs` `manualCommand` (the `-n`, `-g` and `-m` options and the `-en`, `-eg`, `-es`, `-ej` and `-em` events) must name a population that exists, and the number of populations given to `-I` must be a positive whole number.

* Changed MaCS to pass its messages to R after the chromosomes have been simulated, instead of writing them from worker threads, which is not safe. A warning about the command is now shown once rather than once per chromosome, and when MaCS fails the error includes what MaCS reported.

* Standardized all vignettes and articles to use the shared `REFERENCES.bib` bibliography, with forward slashes on Windows to avoid Pandoc resource-path errors during package builds.

* Changed `pedigreeCross` to use `NA` for an unknown parent. A value of `0` now names an individual, as any other value does.

* Added an `unknownParent` argument to `pedigreeCross`, giving the values in `mother` and `father` that mean the parent is unknown, such as `"0"` or `""`. More than one may be given. `NA` is always treated as unknown.

* Changed `pedigreeCross` to extend a pedigree backwards, adding a founder row for any individual named as a parent but lacking a row of its own. With `matchID=FALSE`, only names used as a parent twice or more are added, because a name used once carries no relationship and its child can be a founder. With `matchID=TRUE` every missing name is added, since a name must have a row before it can be matched.

* Changed `pedigreeCross` with `matchID=TRUE` to return the matched individuals and their descendants. The ancestry above a match is no longer simulated, and an individual that is neither matched nor descended from a match is now an error.

* Changed `pedigreeCross` to accept a `MapPop` or `NamedMapPop` as meaning that no simulation has been set up yet. It then builds a temporary `SimParam` of its own and returns a `NamedMapPop` that can be passed to `SimParam$new`. The recombination settings `v`, `p` and `quadProb` may be passed through `...` in this case.

* Changed `hybridCross` to accept a pair of `MapPop` or `NamedMapPop` objects as meaning that no simulation has been set up yet, so that hybrids can be made to serve as the founder population of a hybrid breeding program. It then builds a temporary `SimParam` of its own and returns a `NamedMapPop` when both inputs are a `NamedMapPop`, and a `MapPop` otherwise. The recombination settings `v`, `p` and `quadProb` may be passed through `...` in this case, and `returnHybridPop` must be `FALSE`. Hybrid ids are [mother_id]_[father_id], and when a `crossPlan` repeats a cross every copy of the repeated id is given a further underscore and letter code, as in `A_B_a` and `A_B_b`.

* Added ID based subsetting to `NamedMapPop`, so `pop["a"]` works as it does for a `Pop`.

* Changed `pedigreeCross` to build a pedigree one generation at a time, making all of a generation's crosses in a single `makeCross2` call and batching its selfing and doubled haploid steps the same way. The number of calls into the crossing code now follows the depth of the pedigree rather than its size. Generation numbers are also assigned with a vectorised pass, so sorting a deep pedigree no longer loops over every individual in every pass. Results from a given seed differ from earlier versions, because the order in which random numbers are drawn has changed.

* Changed the default for `maxCycle` in `pedigreeCross` to `NULL`, which uses the number of individuals in the pedigree. That is the deepest a pedigree of that size can be, so the bound is never what stops a pedigree being sorted. A pedigree more than 100 generations deep previously had to have `maxCycle` raised by hand.

* Added a "Gene Drop Simulations" article covering `pedigreeCross`, with an emphasis on using an external pedigree and external genotypes.

* Changed the `addTrait` functions to give gamma distributed effects (`gamma=TRUE`) the requested correlation. Effects were sampled as correlated normal deviates and then transformed to a gamma distribution, which pulled the correlations towards zero, most strongly for small values of `shape`. For example, a requested correlation of 0.6 between two traits with `shape=0.2` gave about 0.45. The normal deviates are now sampled at an adjusted correlation that allows for the transform, so the effects match the requested correlation within sampling error. This applies to both additive and additive-by-additive effects. Traits whose effects follow different distributions, such as a normal trait and a gamma trait, cannot reach every correlation. A requested correlation beyond the attainable limit now gives a warning and uses the limit. Results from a given seed differ from earlier versions when `gamma=TRUE`.

* Added checks that `corA`, `corDD`, `corAA` and `corGxE` in the `addTrait` functions, and `corE` in `setPheno`, `SimParam$setVarE` and `SimParam$setCorE`, are valid correlation matrices. Each must be a numeric matrix with one row and column per trait, with no missing values, symmetric, with ones on the diagonal and all entries between -1 and 1. A covariance matrix was previously accepted without error and is now refused. A matrix that is not positive semi-definite is still smoothed with a warning, as before.

* Epistatic trait validity checks now reject malformed locus pairs, including repeated, missing, non-integer, non-finite, and out-of-range indices. This also applies to epistatic GxE traits.

* Added unit tests for `genParam`.

* Fixed the HWE substitution-effect guard for a fixed locus in an epistatic pair. It now sets that locus's effect to zero while preserving the other locus's effect.

* Fixed epistatic genic variance calculations to use each locus's own HWE substitution effect. Unequal effects or allele frequencies could previously produce incorrect genic variances and HWE adjustments that depended on epistatic pair order.

* Added vignette about recombination.

* Added vignette about quantitative genetic quantities.

* Added vignette about general and specific combining ability (MUST BE REVISED).

* Polished quantitative genetics terminology across the vignettes and documentation.

* Added vignette about genomic prediction and selection functionality.

* Added total non-additive deviations `nd` and associated variance `varN` for the breeding value parameterization and non-additive genotypic effect contributions `gv_n` to genetic values for the genotypic parameterization. These capture both dominance and epistasis.

* Parent average and Mendelian sampling had a bug when `use="bv"` and this option is currently not available.

* Fixed a bug in `SimParam$setRecombRatio` that left the male centromere positions empty. The male centromeres were scaled from a `NULL` starting value, which R silently returns as a zero length vector, so `SP$maleCentromere` and `SP$centromere` both became empty. This caused an out of bounds read in autopolyploid crosses using quadrivalent pairing.

* Fixed `SimParam$setRecombRatio` rescaling the centromeres from the sex-specific positions while rescaling the genetic maps from the sex-average. Repeated calls moved the centromeres relative to their maps. Both are now taken from the sex-average.

* Fixed a bug in `reduceGenome` that always used the female centromere positions, even when `useFemale=FALSE` selected the male genetic map.

* Fixed `SimParam$genMap` dropping chromosome names when sex-specific maps are in use. This removed the `chr` column from `getGenMap`, `getSnpMap` and `getQtlMap` output after any call to `setRecombRatio`, `switchFemaleMap` or `switchMaleMap`.

* Added vector support for nProgeny in `self`, `randCross`, `randCross2`, `makeCross`, and `makeCross2`

* Added MultiPop accessors (`[[`, `[`, `$`, and `names`), and their replacement methods (`[[<-`, `[<-`, `$<-`, and `names<-`).

* Added a `show()` method for MultiPop objects that displays nested structure, indices, names (if available) and number of items/individuals.

* Added `mergeMultiPops()` function to combine multiple Pop and/or MultiPop objects up to a defined `level`.

* Added `flattenMultiPop()` function to reduce the depth of a nested MultiPop object up to a defined `level`.

* Added `selectPop()` function to select a subset of populations from a MultiPop. This function uses the new internal `calcPopValue()` function to calculate a summary value for each population, which is then used for selection.

* Added `unnameMultiPop()` function to remove names from a MultiPop at specific `level`(s).

* Added `newEmptyMultiPop()` function to create an empty MultiPop. It has the same outcome as `newMultiPop()` with no arguments. To make this behavior consistent with `newEmptyPop()`, `newPop()` now also supports creating an empty population with defined ploidy.

* Update `mergePops()` to handle nested MultiPop objects.

* Added asPoisson() function.

* Added asLogNormal() function.

* Improved examples for asCategorical() function.

* Added `SimParam$finalizePheno` field, to provide a user-defined function applied to newly generated phenotypic values (similar to `SimParam$finalizePop`).

* Improved documentation for `SimParam$finalizePop` field.

* AlphaSimR now reports OpenMP support and the default thread count when attached in interactive sessions. Use `options(AlphaSimR.quiet = TRUE)` to silence this message.

* Added a short vignette explaining OpenMP support for parallelization.

* `SimParam$nThreads` now validates assignments. Setting it to `NULL` resets to `getNumThreads()`, and invalid values now fail with a clear error.

* Added optional `nThreads` arguments across OpenMP-enabled R functions and `SimParam` methods so thread counts can be controlled explicitly per call instead of only through `SimParam$nThreads` and is propagated across the package consistently.

* Clarified in function documentation that `simParam = NULL` uses the global `SP` object where applicable.

* Consolidated the use of RNG across the package to enable reproducibility. This is an internal change not visible to users.

* Made meiosis-related C++ RNG reproducible across serial and OpenMP execution by using `dqrng`. This is an internal change not visible to users, but will enable visible reproducibility.

* Fixed a reproducibility bug in `runMacs()` and `runMacs2()`: `set.seed()` can now reproduce MaCS founder simulations, including when chromosomes are simulated in parallel with OpenMP. This is an internal change not visible to users, but will enable visible reproducibility.

* Fixed a bug in `getGvE` with multiple traits.

* Added additional structure to help documents so that related functions will be shown in the "see also" section.

* Changed `mutate` into an S4 generic dispatching on its first argument. Anything that is not a population is passed to the next `mutate` on the search path, so attaching AlphaSimR after a package such as dplyr no longer takes `mutate` away from it. Attaching that package after AlphaSimR still masks this one, where `AlphaSimR::mutate` is needed.

* Added a regression test for `mutate()` assigning mutation sites across chromosomes.

* Performance optimization of functions in meiosis.cpp using Claude

* Performance optimization of function in MME.cpp using Claude

* Performance optimization of the MaCS code in algorithm.cpp, datastructures.cpp and simulator.cpp using Claude. When a limited number of segregating sites is requested, MaCS now samples them while simulating a chromosome instead of generating every site and discarding most of them afterwards. The sites retained are drawn the same way as before, but `runMacs` and `runMacs2` now draw random numbers in a different order, so a given seed will not reproduce founder populations made by earlier versions.

* Performance optimization of the genotype and haplotype extraction functions in getGeno.cpp using Claude. This affects the speed of `pullSnpGeno`, `pullQtlGeno`, `pullSegSiteGeno` and the matching haplotype functions, but not their output.

# AlphaSimR 2.1.0

* Changed R6 and methods from Depends to Imports to match current best practices for R packages

* Change order of call to `finalizePop` in `.newPop` to allow access to recombination tracking data

* Added `parentAverage` and `mendelianSampling` functions

* Fixed bug in `c` for RawPop, MapPop, and NamedMapPop

* Changed `popVar` to an R wrapper to automate casting of vectors to matrices

* Warn when misc lists don't match in length or names in `c(pop, pop2)`

* Corrected bibliography entry month from sept to sep

# AlphaSimR 2.0.0

* Added names to `SP$recHist`

* Added `asCategorical` to convert a normal (Gaussian) trait to an ordered categorical (threshold) trait

* Improved computational performance of simulations with multiple traits

* Added support for data.frames in SimParam genetic map switching functions

* Changed finalizePop function call in `.newPop` to pass simParam as an argument

* Updated version numbering to follow tidyverse format with a major version indicating backward compatibility has been broken

# AlphaSimR 1.6.1

* Fixed bug in `mergePops` and `[` (subset) methods - they were failing for populations that had a misc slot with a matrix - we now check if a misc slot element is a matrix and rbind them for `mergePops` and subset rows for `[` (assuming the first dimension represents individuals)

# AlphaSimR 1.6.0

* Exported `meanEBV` and added `varEBV` to complement `meanP`/`varP` and `meanG`/`varG`

* Changed all parameters of the CATTLE demographic model to exactly match Macleod et al. (2013) - specifically reducing the mutation rate from 2.5e-8 (from human literature) to 1.2e-8 (used in Macleod et al., 2013) and recombination rate from 1e-8 (generic) to 9.26e-9 (used in Macleod et al., 2013). These changes will reduce number of segregating sites to ~240K per chromosome for 100 samples and will run faster.

* changed misc slot in Pop class from a list organized as ind x nodes to a list organized as nodes x ind (this simplified code and increased speed)

* Removed `setMisc` and `getMisc` because the new misc slot structure makes it easy to set and get misc components with base R code

* Added `length` method for Pop class that returns number of individuals (like `nInd`)

* Added `length` method for MultiPop class that returns number of populations

* Fixed bug in quadrivalent pairing resulting in distribution of double reductions not respecting the centromere

# AlphaSimR 1.5.3

* Fixed bug in `SimParam$restrSegSites` with excluding sites at end of chromosome

# AlphaSimR 1.5.2

* Fix SimParam examples for CRAN

# AlphaSimR 1.5.1

* Deleted bad example code for `setMisc`

* Changed examples to use a single thread for CRAN testing this change is not shown in the documentation

# AlphaSimR 1.5.0

* Renamed `MegaPop` to `MultiPop`

* Fixed bug in `writePlink` to correctly export map positions in cM

* Fixed bug in `writeRecords` due to removed reps slot in pops

* Added `altAddTraitAD` for specifying traits with dominance effects using dominance variance and inbreeding depression

* Add miscPop slot to class `Pop`

# AlphaSimR 1.4.2

* Updated MaCS citation to https site

# AlphaSimR 1.4.1

* Changed citation to use `bibentry` instead of `citEntry`

# AlphaSimR 1.4.0

* Fixed a bug in IBD tracking

* Add `setFounderHap` to SimParam for applying custom haplotypes to founders

* Added `addSnpChipByName` to SimParam for defining SNP chips by marker names

# AlphaSimR 1.3.4

* Changed C++ using `sprintf` to use `snprintf`

# AlphaSimR 1.3.3

* Fixed bug in calculation of genic variance

* Fixed `importHaplo` not passing ploidy to `newMapPop`

* Fixed bug with correlated error variances

# AlphaSimR 1.3.2

* Fixed column name bug with multiple traits in `setEBV`

* Fixed CTD caused by `runMacs` when too many segSites are requested

* Fixed missing names in GV when using `resetPop`

* Fixed bug in `importTrait`

* `popVar` now deals with matrices having 1 row

# AlphaSimR 1.3.1

* Updated link to Gaynor, 2017

# AlphaSimR 1.3.0

* Added ability to exclude loci by name in `SimParam$restrSegSites`

* `pullMarkerGeno` and `pullMarkerHaplo` now work with a MapPop class

* Added `setMarkerHaplo` to manually change genotypes in a Pop or MapPop

* Added `addSegSite` for manually adding segregating sites to a MapPop class

* `simParam$setCorE` has been deprecated in favor of a corE argument in `simParam$setVarE`

* `setPheno` now takes corE as an argument

* `setPheno` now allows the user to set phenotypes for a subset of traits

* Add `newEmptyPop` to create populations with zero individuals

* Removed reps slot from populations and heterogeneous residual variance GS models

* Added h2, H2, and corE to `setPhenoGCA` and `setPhenoProgTest`

* The "EUROPEAN" species history was removed from `runMacs` due to lengthy runtime

# AlphaSimR 1.2.2

* Added `getPed` to quick extract a population's pedigree

* Added `getGenMap` to pull a genetic map in data.frame format

# AlphaSimR 1.2.1

* Fixed bugs relating to `importData` functions

* Fixed `writePlink` errors and no longer requires equal length chromosomes

# AlphaSimR 1.2.0

* Added `importGenMap` to format genetic maps for AlphaSimR

* Added `importInbredGeno` and `importHaplo` to make it easier to create a simulation from external data

* Added `importSnpChip`, `importTrait` to `SimParam` to make it easier to manually define traits

* Added `pullMarkerGeno` and `pullMarkerHaplo` to make it easier to extract genotypes and haplotypes of specific loci without defining a trait or SNP chip

* `reduceGenome`, `mergeGenome` and `doubleGenome` should really now work with pedigree and recombination tracking

# AlphaSimR 1.1.2

* Added missing #ifdef _OPENMP to OCS.cpp

# AlphaSimR 1.1.1

* Removed use of PI variable in C++ code due to it being compiler specific

# AlphaSimR 1.1.0

* Added snpChip argument to `pullIbdHaplo` for backward compatibility

* Exposed internal mixed model solvers

* All selection functions now return a warning when there are not enough individuals

* Fixed error in `pullIbdHaplo` when chr isn't NULL

* Fixed an error with assigning 1 QTL and/or SNP

* Changed geno slot from matrix to list to support future RcppArmadillo changes

* `doubleGenome` and `reduceGenome` now work with IBD tracking

# AlphaSimR 1.0.4

* Fixed errors in implementation of Gamma Sprinkling model

# AlphaSimR 1.0.3

* Fixed formatting error in genetic maps created by runMacs that broke genotype extraction functions

# AlphaSimR 1.0.2

* Added h2 and H2 to `setPhenoGCA`

* `pullGeno` and `pullHaplo` functions now report marker names from the genetic map

# AlphaSimR 1.0.1

* Removed lazyData field in DESCRIPTION

# AlphaSimR 1.0.0

* AlphaSimR manuscript has been published in G3 (citation added)

* Changed to a Gamma Sprinkling model for crossovers, default is still a Gamma model

* Change default interference parameter (v) to 2.6 to be consistent with the Kosambi mapping function (was 1, consistent with the Haldane mapping function)

* New internal id (iid) that allows user to freely change id slot in populations

* `runMacs2` now adjusts Ne for autopolyploids

* Parent populations are now passed to `finalizePop`

* Check added that throws an error when use of discontinued "gender" argument is detected

* Added experimental `MegaPop-class`

# AlphaSimR 0.13.0

* References to gender have been changed to the more appropriate terms sex or sexes

* Added misc slot to populations

* Added `finalizePop` to `SimParam`

* Added physical positions to `getSnpMap` and `getQtlMap`

* You can now use h2 and H2 to specify error variance in `setPheno`

* `SimParam$setVarE` now accepts a matrix for varE

* Fixed a bug in `editGenome` when making multiple edits

* Adding merging of centromere vector in `cChr`

# AlphaSimR 0.12.2

* GxE traits now default to random sampling of p-values

* Fixed a bug in `restrSegSites`

# AlphaSimR 0.12.1

* Fixed a bug in selection of segSites

# AlphaSimR 0.12.0

* Changed output of `genParam` to match Bulmer, 1976

* `nProgeny` added to `makeCross` and `makeCross2`

* All `SimParam` documentation is now in `?SimParam`

* Non-overlapping QTL and SNP is now the default

* New interface for `restrSegSites` in `SimParam`

* Fixed subset by id for populations

* Fixed major bug in `newMapPop`

# AlphaSimR 0.11.1

* Switched to a circular design for the balance option in `randCross` and `randCross2`

* Added `reduceGenome` and `doubleGenome` for changing ploidy levels

* Added minSnpFreq to SimParam_addSnpChip for any reference population

* The `c` function now merges individuals for MapPop objects (was chromosomes before)

* The `cChr` function new merges chromosomes for MapPop objects

* Fixed broken SimParam_addStructuredSnpChip

* Removed broken `pullMultipleSnpGeno` and `pullMultipleSnpHaplo`

* Fixed broken `writePlink`

# AlphaSimR 0.11.0

* Rework of `setEBV` (breaks some scripts)

* Genotype data now stored as bits (was bytes)

* Implemented a gamma model for crossover interference

* Added the mutate function to model random mutations

* Added a vignette explaining the biological model for traits

* GS models now handle polyploids

* Heterogenous error variance is now optional in GS models (default is homogeneous error)

* Improved gene drop functionality of pedigreeCross

* Added keepParents option to makeDH and self (indirectly extends `selectFam` and `selectWithinFam`)

* Added RRBLUP_SCA2

* Set methods for the "show" function when applied to populations

* Fixed a bug returning the first individual when selecting 0

* Fixed error in recombination track when using `makeDH`

* Fixed error causing epistatic effects to mask GxE effects

* Fixed an error with `pullSegSiteGeno` and `pullSegSiteHaplo` with variable number of sites per chromosome

# AlphaSimR 0.10.0

* Added traits with epistasis

* Max number of threads automatically detected

* Added RRBLUP_D2

* Added version tracking to `SimParam`

* Removed `trackHaploPop` (super-ceded by `pullIbdHaplo`)

* Added `fastRRBLUP`

* Fixed faulty double crossover logic

* Fixed broken `writePlink`

* Fixed broken `pullIbdHaplo`

* `mergePops` no longer assumes diploidy

# AlphaSimR 0.9.0

* Added support for autopolyploids

* Added `RRBLUP_GCA2`

* `randCross2` can now "balance" crossing when not using gender

* Fixed recombination tracking bug in `createDH2`

* Removed bug in `setEBV` with append=TRUE

# AlphaSimR 0.8.2

* Fixed ambiguous overloading in optimize.cpp

# AlphaSimR 0.8.1

* `setPheno` (not `setPhenoGCA`) passes the number of reps to populations

* Fixed bug in `editGenomeTopQtl`

* Fixed bug in `RRBLUP_D`

* Fixed bug in `resetPop`

* Fixed bug in SimParam_rescaleTraits

* Removed unimplemented SimParam_restrSnpSites and SimParam_restrQtlSites

* Add error message for no traits in `calcGCA`

# AlphaSimR 0.8.0

* Added GxE traits with zero environmental variance

* Faster trait scaling

* Faster calculation of genetic values

* `dsyevr` now called via arma_fortran

* Added OpenMP support

* Parallelized `cross2`

* Parallelized `runMacs`

* Parallelized calculation of genetic values

* Variance calculations now account for inbreeding

* Fixes for male selection in `selectOP`

# AlphaSimR 0.7.1

* Add fixEff to `setPhenoGCA`

# AlphaSimR 0.7.0

* Added default `runMacs` option to return all segSites

* Added ability to specify separate male and female genetic maps

* `pullGeno` and `pullHaplo` functions can now specify chromosomes

* Added `RRBLUP2` for special GS cases

* Improved speed by replacing Rcpp random number generators

* Changed available MaCS species

* GS functions now use populations directly

* Added `pullIbdHaplo`

* Added `writePlink`

* Fixed population sub-setting checks to prevent invalid selections

* Fixed slow `calcGCA`

* Fixed error in `addTraitAG` preventing multiple traits

* Fixed bug with `mergePops` when merging ebv

* Fixed bug in `setVarE` when using H2 and multiple traits

# AlphaSimR 0.6.1

* `selectFam` now handles half-sib families

* `selectWithinFam` now handles half-sib families

* Removed restriction on varE=NULL in `setPhenoGCA`

# AlphaSimR 0.6.0

* Added NEWS file

* Added `selectOP` to model selection in open pollinating plants

* Added `runMacs2` as a wrapper for `runMacs`

* Fixed error when using H2 in SimParam_setVarE
