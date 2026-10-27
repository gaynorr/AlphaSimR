# AlphaSimR 2.1.0.9023

## Changes to results from a seed

Seeded simulations run with earlier versions will not be reproduced exactly in these cases:

* `runMacs` and `runMacs2` draw random numbers in a different order, and now place each chromosome on absolute positions with its requested genetic length, so every MaCS founder population differs from earlier versions.
* Genetic maps whose first locus is not at position 0 now give different crossovers, because maps are no longer shifted to start at 0.
* `pedigreeCross` makes each generation's crosses in one call, which changes the order of random numbers.
* `addTrait` functions with `gamma=TRUE` now sample effects at the requested correlation.
* Traits with epistasis that are added, or rescaled by `SimParam$rescaleTraits`, with `useVarA=TRUE`, the default, are scaled to an additive variance that now accounts for linkage disequilibrium within epistatic pairs (see "Genome and traits"). Every founder population has some, so their effects, and error variances set from a heritability, differ from earlier versions.
* Genetic values, `genParam` output and new traits' effects can differ in the last few digits (around 1e-13), which can change results that depend on exact ties.
* The fast mode of `RRBLUP`, which replaces `fastRRBLUP`, estimates the variance components by REML unless they are supplied. With them supplied, it agrees with `fastRRBLUP` to about one part in a million.
* The reimplemented mixed model solvers (see below) give estimates that differ from earlier versions, so `RRBLUP`, `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA` can too.

## Breaking changes

* Removed `fastRRBLUP`. Use `RRBLUP(..., fastMode=TRUE)` in its place, which is faster and takes `Vu` and `Ve` as before. Unless they are supplied, fast mode estimates them by REML on `nSubsample` records, 3000 by default, where `fastRRBLUP` took them from the simulation's settings. `fastMode=TRUE` uses fast mode only when it is expected to be faster.
* Deprecated `RRBLUP2`, `RRBLUP_D2`, `RRBLUP_GCA2` and `RRBLUP_SCA2`. Their unnumbered counterparts now accept `Vu` and `Ve`, or an earlier fit as `varComp`, so `RRBLUP_D2(pop, Va=a, Vd=d, Ve=e, useEM=FALSE)` becomes `RRBLUP_D(pop, Vu=c(a, d), Ve=e)`.
* `RRBLUP_SCA` and `RRBLUP_SCA2` give an error for a population that is not diploid.
* The default `maxIter` of `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA` is now 100 REML iterations, up from 40.
* `SP` is now looked for only in the global environment, and a missing or invalid `simParam` gives a clear error.
* `pedigreeCross` no longer treats every parent missing from `id` as unknown. `NA`, and any values given in the new `unknownParent` argument, such as `"0"`, mean an unknown parent. Other missing parents are added as founders; with `matchID=FALSE`, only those named as a parent more than once.
* `pedigreeCross` with `matchID=TRUE` returns only the matched individuals and their descendants.
* `SimParam$importTrait` gives an error for duplicated `markerNames`, which were silently ignored.
* `addTrait` correlation arguments and `corE` must be valid correlation matrices; a covariance matrix is now refused.
* `addSegSite` gives an error for a site past the end of the chromosome unless `extendGenLen=TRUE`.
* Maps with negative, unsorted or non-finite positions are now an error instead of being shifted.
* A `SimParam` with recombination tracking that was saved by an earlier version can no longer track recombination or pull IBD haplotypes.
* The methods `addToRec`, `ibdHaplo`, `updateLastId` and `addToPed` are no longer documented, as they are for internal use.
* AlphaSimR now imports the Matrix package, a recommended package that ships with R.
* AlphaSimR now requires R 4.3.2 or later, the first version to ship Matrix 1.6, which pedigree BLUP needs.

## Pedigree analysis

* Added `pedigreeA` and `pedigreeInbreeding`, which calculate the numerator relationship matrix and inbreeding coefficients from the pedigree recorded with `setTrackPed`, for any even ploidy level.
* Added `pedBLUP`, which fits a single trait pedigree BLUP (animal) model with REML, or with supplied `Vu` and `Ve`. Its result, a new `PedSol` class, is used with `setEBV`, including for individuals created after the fit.
* Added `usePed` to `selectCross` and the new `selectOCS`, which measures expected fixation with the pedigree relationship matrix in place of genotypes.
* Added `nGen` to `getPed`, which returns the pedigree for more than one generation.

## Genomic selection

* Added `fastMode` to `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA`, which solves the model with much less memory on large populations.
* Added multiple trait models to the fast mode of `RRBLUP` and to `RRBLUP_D`, `RRBLUP_GCA`, `RRBLUP_SCA` and `solveRRBLUPMK`.
* Added `varComp` to the `RRBLUP` models, which reuses the variance components of an earlier fit. Fits now record their model in a new `model` slot of `RRsol`.
* `RRBLUPMemUse` takes `fastMode` and `nTraits`, and estimates the memory of the reimplemented solvers.
* Reimplemented the mixed model solvers `solveUVM`, `solveRRBLUP`, `solveMVM`, `solveRRBLUPMV`, `solveMKM` and `solveRRBLUPMK` from published methods, replacing code derived from packages with licenses incompatible with AlphaSimR's. Arguments and return values are unchanged, except that `LL` is now the restricted log-likelihood documented on each help page. The solvers drop linearly dependent columns of `X`, reject non-finite inputs, and are faster, in some cases by an order of magnitude or more.
* Added `onlyValues` to `setEBV`, which returns the estimates without changing the population.

## Selection and crossing

* Added `selectOCS`, and `restrInbr` to `selectCross`, for selection with a restriction on the expected increase in fixation.
* Added `onShortfall` to the selection functions, to warn, stop or stay silent when fewer candidates are available than requested.
* The selection functions give a clear error for an empty population.
* Added vector support for `nProgeny` in `self`, `randCross`, `randCross2`, `makeCross` and `makeCross2`.
* `pedigreeCross` adds a founder row for any individual named as a parent without a row of its own, and its default `maxCycle` is now the number of individuals in the pedigree.
* `pedigreeCross` and `hybridCross` accept `MapPop` or `NamedMapPop` objects to make founders before a simulation is set up.
* `hybridCross` accepts a `crossPlan` of IDs when `returnHybridPop=TRUE`.

## Genome and traits

* Chromosomes now have a genetic length that can extend past the last locus, stored in a new `genLen` slot of `MapPop` and read with `SimParam$genLen`, `$femaleGenLen` and `$maleGenLen`.
* Genetic maps now use absolute positions measured from the start of the chromosome.
* Added additive-by-additive epistasis to `SimParam$importTrait` through `epiPairs` and `epiEff`.
* Added `mutRate` to `runMacs`, and `predictSegSitesMoments` for theoretical predictions of the number of segregating sites.
* Added `asPoisson` and `asLogNormal`.
* Added `SimParam$finalizePheno`, a user-defined function applied to new phenotypes.
* Trait names can be used to choose traits in `setPheno`, `editGenomeTopQtl`, `SimParam$removeTrait` and `SimParam$switchTrait`.
* Added covariances between traits to the genic, Hardy-Weinberg and linkage disequilibrium parts of `genParam`, and the total non-additive deviations `nd`, `varN` and `gv_n`.
* `genParam` now accounts for linkage disequilibrium between the two loci of an epistatic pair, fitting the pair's epistatic value to their observed joint genotype frequencies. Under such disequilibrium `varA` and `varD` can be larger and `varAA` smaller than before, and a trait with epistasis but no dominance effects can have nonzero `dd` and `varD`, all of it in `covD_L`. Results under linkage equilibrium are unchanged.
* Added `mu_LE` to `genParam`, the mean genetic value with each locus's observed genotype frequencies but no linkage disequilibrium. With `mu_HW`, it splits the mean into the change from departures from Hardy-Weinberg genotype frequencies, `mu_LE - mu_HW`, and the change from linkage disequilibrium between epistatic loci, `mu - mu_LE`.
* `addTrait` functions with `gamma=TRUE` now give the requested correlation, or warn and use the attainable limit.
* Changed `mutate` into an S4 generic, so attaching AlphaSimR no longer masks `dplyr::mutate` for non-population arguments.

## MultiPop

* Added `[[`, `[`, `$` and `names` accessors and replacement methods, and a `show()` method.
* Added `mergeMultiPops`, `flattenMultiPop`, `unnameMultiPop`, `newEmptyMultiPop`, `splitPop`, `selectPop` and `nPop`.
* `mergePops`, `meanG`, `meanP`, `meanEBV`, `varG`, `varP`, `varEBV`, `gv`, `pheno` and `ebv` accept MultiPop objects. Added `meanGPop`, `meanPPop`, `varGPop` and `varPPop` for values between populations.
* `newPop` can create an empty population with a given ploidy.

## Threads and reproducibility

* Simulations, including `runMacs` and `runMacs2`, can now be reproduced with `set.seed()` regardless of the number of threads.
* Added `nThreads` arguments to the functions that use OpenMP. `SimParam$nThreads` now validates assignments, and `NULL` resets it to `getNumThreads()`.
* AlphaSimR reports OpenMP support and the default thread count when attached in interactive sessions; `options(AlphaSimR.quiet = TRUE)` silences this.

## Performance

* Faster crossing, genetic values, `genParam` and the functions built on it, often several times faster with many QTL.
* Faster `c()` and `mergePops`, and lower memory use when reading genotypes.
* `pullIbdHaplo` is much faster on large pedigrees, and recombination tracking uses about a sixteenth of the memory.
* MaCS samples the requested segregating sites while simulating, which makes `runMacs` faster and use less memory.

## Bug fixes

* Fixed the intercept of genetic values predicted by single trait `fastRRBLUP`, which shifted every predicted genetic value by a constant.
* Fixed `RRBLUP` failing when `traits` is a function returning a plain vector.
* Fixed multivariate `RRBLUP` and `solveMVM` failing with a rank deficient relationship matrix.
* Fixed `calcGCA` for crosses to a single tester with more than one record per cross.
* Fixed `selectFam` and `selectWithinFam` pooling every individual with an unknown parent into one family, and `selectWithinFam` ignoring `famType` for a MultiPop.
* Fixed `setMarkerHaplo` changing the genotypes of the population it was given.
* Fixed `mergePops` failing on a named list of populations whose `misc` slot holds a `Pop`.
* Fixed recombination tracking for quadrivalents in polyploids, which broke `pullIbdHaplo`.
* Fixed `SimParam$setRecombRatio` leaving the male centromeres empty and moving centromeres relative to their maps on repeated calls.
* Fixed `reduceGenome` always using the female centromeres.
* Fixed `SimParam$genMap` dropping chromosome names when sex-specific maps are in use.
* Fixed epistatic genic variances, which were wrong for unequal effects or allele frequencies, and the HWE substitution effect of a fixed locus in an epistatic pair.
* Fixed `getGvE` with multiple traits.
* Fixed several cases where invalid input could read or write memory outside a population's genotypes, including bad `haplo` values in the `pull*Haplo` functions and out-of-range indices in `self`, `hybridCross` and `runMacs`. These now give an error.
* Epistatic trait validity checks now reject malformed locus pairs.
* MaCS messages are now passed to R safely, once, and a MaCS failure reports what MaCS said.
* Input checks across the package now say which argument is wrong and what is expected.
* `parentAverage` and `mendelianSampling` with `use="bv"` are disabled because of a bug.

## Documentation

* Split the `SimParam` help page into one page per method, such as `?"SimParam$addTraitA"`.
* Reorganized the help pages of the `RRBLUP` models around the model each fits.
* Added vignettes on genomic selection and parallelization. The genome article is now a vignette.
* Added articles to the package website, <https://gaynorr.github.io/AlphaSimR/articles/>, on selection, pedigree BLUP, optimal contribution selection, custom traits and trait distributions, quantitative genetics, combining ability, demographic models, gene drop simulations, importing data, simulating a plant breeding program, and storing data in `misc` and `miscPop`, and revised the genetic recombination article.
* Corrected `citation("AlphaSimR")`, which split the authors' names into separate authors and gave the wrong page number for the AlphaSimR paper.
* Added references to the literature behind the claims in the vignettes and articles.

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
