# AlphaSimR 2.1.0.9022

* Shortened the `SimParam` help page. It now opens with a list of its methods grouped by task, linked to each method's section, and a single example that sets up a simulation step by step. That example replaces the separate examples of 17 methods, each of which repeated the same setup. The methods with examples that show something particular to them, `resetPed`, `altAddTraitAD`, `importTrait`, `setCorE` and `rescaleTraits`, keep their own.

* Moved the "Quantitative Genetics in AlphaSimR", "Genetic Recombination in AlphaSimR" and "General and Specific Combining Ability in AlphaSimR" vignettes out of the package and onto the package website, <https://gaynorr.github.io/AlphaSimR/articles/>, so that the package keeps the vignettes most useful when starting out. They can no longer be opened with `vignette()`, and the help pages and remaining vignettes that referred to them now link to the website. The help pages of `runMacs`, `runMacs2` and `gv` linked to articles that were already only on the website, and the package help page named a vignette that does not exist; these now point to the right place.

* Deprecated `RRBLUP2`, `RRBLUP_D2`, `RRBLUP_GCA2` and `RRBLUP_SCA2`, which now give a warning saying what replaces them. They still work as before.
  * `RRBLUP`, `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA` take `Vu` and `Ve`, or an earlier fit as `varComp`, with `fastMode=FALSE`, which used to be an error. For one trait they then solve the mixed model equations once with the variance components held fixed, with the solver of their own REML fit on whichever of the markers and the records is cheaper, and give the same result, to rounding, as the numbered model with `useEM=FALSE`. A variance of zero, which REML often estimates for a dominance or SCA term, leaves that term out. `Vu` holds one variance per random term, so `RRBLUP_D2(pop, Va=a, Vd=d, Ve=e, useEM=FALSE)` becomes `RRBLUP_D(pop, Vu=c(a, d), Ve=e)`, and likewise `VuF`, `VuM` and `VuD` of `RRBLUP_GCA2` and `RRBLUP_SCA2` become `Vu` in that order. Both `Vu` and `Ve` must be given. For several traits, supplied variance components still need `fastMode=TRUE`.
  * `RRBLUPMemUse` with `useEM=FALSE` estimates the memory of an unnumbered model given its variance components with `fastMode=FALSE`. Its estimate for `RRBLUP_D2`, `RRBLUP_GCA2` and `RRBLUP_SCA2` with `useEM=FALSE` now follows them to the model written on the records when there are fewer records than marker effects, which they already used.

* `RRBLUP_SCA` and `RRBLUP_SCA2` give an error for a population that is not diploid. The model is designed for diploid single crosses, and only fast mode refused other ploidy levels before.

* Reorganized the help pages of `RRBLUP`, `RRBLUP2`, `RRBLUP_D`, `RRBLUP_D2`, `RRBLUP_GCA`, `RRBLUP_GCA2`, `RRBLUP_SCA` and `RRBLUP_SCA2`. Each now opens with a short description and has a details section, which several arguments already referred to but `RRBLUP` lacked. The details state the model and say how it is fitted, with references, in separate parts for one trait and for several traits, followed by fast mode where the function has it. The pages of the four models ending in 2 say that they fit one trait only, by EM-REML or with the variance components taken as given. The `maxIter` argument of `RRBLUP_D2` is now described as the limit on EM iterations, which it is. Added Whittaker, Thompson and Denham (2000), Meuwissen, Hayes and Goddard (2001), and Dempster, Laird and Rubin (1977) to the package references.

* Added `fastMode` to `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA`, with `Vu`, `Ve` and `nSubsample` as in `RRBLUP`. Fast mode solves the mixed model equations of every random term and trait together by preconditioned conjugate gradient, iterating over the genotypes packed to one bit per allele, so the memory used stays close to the size of the genotypes. The variance components are estimated by the default mode's REML on `nSubsample` records, unless `Vu` and `Ve` are supplied, and the marker effects are solved for on every record. Given the variance components, the result is the BLUP of the same model that the default mode fits. For one trait, `Vu` is a vector with one value per random term, and for several traits a list of covariance matrices or the block diagonal matrix held in the `Vu` slot of a fit, so a fit's own variance components can be passed back. `RRBLUP_SCA` supports fast mode for diploids only. When every marker variance is zero, fast mode fits the fixed effects alone and returns marker effects of zero. The default `maxIter` of the three models is now `NULL`, which gives 100 iterations of the REML fit, as before, or 1000 iterations of the conjugate gradient solver in fast mode. Results in the default mode are exactly the same as before, including results from a seed.
  * The new arguments of `RRBLUP`, `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA` come after `...`, so they are given by name, and calls that give `simParam` and `nThreads` by position work as they did before `fastMode` was added to `RRBLUP`.
  * `fastMode=TRUE` in `RRBLUP`, `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA` asks for fast mode rather than requiring it. The cost of both modes is estimated from the numbers of records, markers, traits and `nSubsample`, and fast mode is used only when its cost is less than half that of the default mode, and never when `nSubsample` covers every record, which would repeat the default mode's REML fit. Otherwise the model is fitted exactly as with `fastMode=FALSE`. Supplying `Vu` and `Ve` always uses fast mode. Most of fast mode's time is its REML fit on the subset: with 20,000 records and 2,000 markers, the solve on every record took 1.5 to 3.5 seconds, and the REML fit on 5,000 records about as long as the default mode's fit on every record, because with fewer markers than records its cost is set by the number of marker effects. The default `nSubsample` is 3000, down from 5000, in all four models and in `RRBLUPMemUse`. In tests with 20,000 records, a subset of 3000 gave accuracies in a later generation within 0.001 on average, and 0.003 at worst, of a fit using every record or a subset of 8000, for all four models, and for two traits with `RRBLUP`, at heritabilities of 0.2 to 0.5. It took between three fifths and a twelfth of the time of a subset of 5000. Subsets of 2000 or fewer lost up to 0.011 for `RRBLUP_SCA`.
  * The fast mode of `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA` estimates its variance components from each term's cross product of the codes on the subset when there are at least as many markers as records in it, instead of first decomposing the cross product into a factor that the REML fit then multiplied back out. The estimates are the same to rounding. With 20,000 records, 5,000 markers and a subset of 5000, `RRBLUP_D` took 30 seconds in place of 80.
  * Added `varComp` to `RRBLUP`, `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA`, which takes an earlier fit and uses its `Vu` and `Ve` in fast mode, in place of supplying `Vu` and `Ve`. The fit must be of the same model, or of its counterpart ending in 2, for the same number of traits and on the same markers, or an error is given. Every fit now records the model it came from in a new `model` slot of `RRsol`, which is empty for an object created directly. Fits saved before this change have no model and are refused by `varComp`; their `Vu` and `Ve` slots can still be supplied directly.
  * `RRBLUPMemUse` takes `fastMode` for `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA`. For 100,000 records and 10,000 markers it estimates 3 GB for the fast mode of `RRBLUP_D` or `RRBLUP_GCA` against 49 GB for the default mode, and on 8 threads fast mode took 136 seconds for `RRBLUP_D` and 102 seconds for `RRBLUP_GCA`.

* Fixed `RRBLUP` failing with "argument is of length zero" when `traits` is a function returning a plain vector, which it now treats as one column as the other models do.

* Fixed `mergePops` failing with 'invalid class "Pop" object' when given a named list of populations whose `misc` slot holds a `Pop`.

* Removed `fastRRBLUP`, which is now the fast mode of `RRBLUP`. Call `RRBLUP` with `fastMode=TRUE` in its place; its arguments `Vu`, `Ve` and `nSubsample` are now arguments of `RRBLUP`, and `maxIter` limits the conjugate gradient iterations in fast mode. When fast mode is used, results, including results from a seed, are exactly the same as from `fastRRBLUP`, but `fastMode=TRUE` now uses fast mode only where it is expected to be faster, as described below. With the default `fastMode=FALSE`, `RRBLUP` behaves as before, and gives an error if `Vu` or `Ve` is supplied, because it always estimates them. `RRBLUPMemUse` takes `fastMode` for the fast mode of `RRBLUP`, and still accepts `model="fastRRBLUP"`, which it reads as `model="RRBLUP"` with `fastMode=TRUE`. Error messages from the fast mode now begin with "RRBLUP fastMode:" rather than "fastRRBLUP:".

* Fixed `selectFam` and `selectWithinFam` pooling every individual with an unknown parent into one family. An individual whose family is not known, because a parent defining it is unknown, is now treated as a family of its own, so that `selectWithinFam(founders, 1)` returns every founder rather than one. For `famType="B"` this applies when either parent is unknown. Parents of `"0"`, as `newPop` and the import functions record them, and of `NA`, as `pedigreeCross` records them, are both treated as unknown; with `NA`, `selectFam` previously dropped those individuals from selection altogether. Populations whose parents are all known give the same results as before.

* Added covariances between traits to the genic, Hardy-Weinberg and linkage disequilibrium parts of `genParam`. For each of `G`, `A`, `D`, `AA` and `N`, `genicCovX` is an nTraits by nTraits matrix of genic variances and covariances under Hardy-Weinberg and linkage equilibrium, `covX_HWMat` adjusts it for departures from Hardy-Weinberg proportions, and `covX_LMat` is the contribution of linkage disequilibrium, so that the three add up to `varX`. Their diagonals are the existing vectors `genicVarX`, `covX_HW` and `covX_L`, which are unchanged. Genic covariances come only from loci that both traits have as QTL, and for `AA`, only from pairs of loci that both traits treat as interacting, so traits on disjoint QTL have no genic covariance. All other outputs of `genParam`, and results from a seed, are exactly the same as before.

* Gave the checks on user input in `SimParam`, the founder and import functions, `editGenome`, `editGenomeTopQtl`, `setPhenoGCA`, `setPhenoProgTest`, the single trait `RRBLUP` models, `selectPop` and `writeRecords` messages that say which argument is wrong and what is expected, in place of R's generic "... is not TRUE". `SimParam$addSnpChipByName` and `SimParam$importTrait` now list the marker names that are not on the genetic map. The conditions checked are unchanged.

* Added `onShortfall` to `selectInd`, `selectFam`, `selectWithinFam`, `selectOP`, `selectOCS` and `selectCross`, for when fewer suitable candidates are available than requested. The default, `"warning"`, behaves as before: it warns and returns all of them. `"error"` stops instead, with a message giving how many were available, and `"none"` returns them without a warning, which suits `selectWithinFam` used as a cap on family size. The `selectWithinFam` warning now says how many families fall short.

* Changed the selection functions to give the error "pop contains no individuals to select from" for a population with no individuals, instead of failing later with an unrelated message. An empty population inside a `MultiPop` is returned unchanged. Fixed `selectFam` failing with "no rows to aggregate" when none of the individuals were eligible, which now returns no individuals, and `selectOCS` now says when there are no female or male candidates to select from. Fixed `selectWithinFam` ignoring `famType` for a `MultiPop`.

* Added an article on importing data, covering genetic maps, phased haplotypes, genotypes of inbred lines, PLINK text files, and starting a simulation from imported founders with `SimParam$importTrait` and `SimParam$addSnpChipByName`.

* Added `nGen` to `getPed`, which returns the pedigree for more than one generation. With `nGen` greater than 1, or `Inf` for all ancestors, the pedigree recorded in `SimParam` is followed, so pedigree tracking must be on. The result holds the individuals and their ancestors, each once, with parents before their progeny. Every meiosis counts as a generation, including each generation of selfing and each doubled haploid. The default `nGen=1` returns the parents in the population's slots, as before.

* Trait names can now be used to choose traits in `setPheno`, `editGenomeTopQtl`, `SimParam$removeTrait` and `SimParam$switchTrait`, as they already could in the selection and genomic selection functions. An unknown name gives an error. `getQtlMap`, `pullQtlGeno` and `pullQtlHaplo` already accepted trait names, and their documentation now says so.

* Added `onlyValues` to `setEBV`. When `TRUE`, the estimated values are returned as a matrix and the population is left unchanged, as `onlyPheno` does for `setPheno`, `setPhenoGCA` and `setPhenoProgTest`.

* Added `splitPop()` to split a Pop or MultiPop into nested subpopulations using grouping vectors, functions, or data frames.

* Expanded `calcPopValue()` to apply a function to each Pop, optionally flatten results to a selected nesting level, and retain source paths for simplified values.

* Added name-preservation options to `flattenMultiPop()` for retaining or constructing hierarchical names when flattening.

* Fixed the intercept of the genetic values predicted by single trait `fastRRBLUP`, which was found against the centered marker genotypes but used with the uncentered ones, shifting every genetic value predicted with `setEBV` by the same constant, the mean genotype times the marker effects. Predictions of genetic values now agree with `RRBLUP`. Breeding values and the ranking of individuals were not affected.

* Added multiple trait models to `fastRRBLUP`, `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA`. Give `traits` more than one trait and the traits are fitted together, with unstructured covariance matrices among the traits for the marker effects and for the errors, estimated by REML. Every trait must be recorded on every individual. Single trait fits are exactly the same as before.
  * `fastRRBLUP` estimates the covariance matrices by the multiple trait REML of `RRBLUP` on its subset of records, or takes them as matrices through `Vu` and `Ve`, and then separates the traits by a canonical transformation and solves them together by conjugate gradient. Given the covariance matrices, the result is the BLUP of the multiple trait model. A fit to several traits takes about as long as fitting them one at a time, and its memory stays close to the size of the genotypes.
  * `solveRRBLUPMK` fits a multiple trait model when `y` has more than one column, by average information REML in the Cholesky factors of the covariance matrices, on whichever of the variance matrix of the records and Henderson's equations is cheaper. `Vu` is then a list of covariance matrices, one per matrix in `Mlist`, and `Ve`, `beta` and each element of `u` have one column per trait. With one column of `y` its results are unchanged. Its arguments are no longer modified in place, and neither are those of `solveMKM`.
  * `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA` use it for several traits. The `Vu` slot of their result then holds the covariance matrices of the marker effects as a block diagonal matrix with rows and columns named "term:trait", such as "a:Trait1" and "d:Trait1" for `RRBLUP_D`, and `Ve` has the traits' names. The default `maxIter` of the three models is now 100 rather than 40, for the larger multiple trait fits; single trait fits that converged within 40 iterations give the same result.
  * `RRBLUPMemUse` estimates the memory of the multiple trait forms of these models through `nTraits`.

* Fixed `setMarkerHaplo` changing the genotypes of the population it was given, and of every R object sharing them, as well as those of the population it returns. The haplotypes are now written into a copy of each chromosome that has loci to set.

* Reimplemented the mixed model solvers `solveUVM`, `solveRRBLUP`, `solveMVM`, `solveRRBLUPMV`, `solveMKM` and `solveRRBLUPMK` from published methods. The earlier versions were translated from code in other R packages released under licenses that are not compatible with AlphaSimR's MIT license. The arguments and returned values are unchanged, and `LL`, where it is returned, is now the restricted log-likelihood documented on each help page, which differs from the earlier value by a constant. All six drop columns of `X` that are linearly dependent on the others and give them a fixed effect of zero, instead of failing or inverting a singular matrix. An intercept, a column whose values are all the same, is always kept, wherever it is in `X`. The fixed effects are fitted through an orthonormal basis for the columns of `X`, after a least squares fit has been taken off the response, so a response with a large mean or a covariate with a large offset costs no precision. `solveUVM`, `solveMVM` and `solveMKM` give an error for a `K` that is not positive semidefinite, whichever way the model is solved, and treat eigenvalues at rounding level the same way for the BLUPs as for the likelihood and on every solving path, keeping small eigenvalues of `K` that `Z` scales back up. One rule decides which eigenvalues are zero up to rounding, whichever of `M'M` and `MM'` is decomposed, so padding the markers with columns of zeros does not change the fit, and the parts of the response and of `X` outside the kernel's span are formed directly, so a covariate that lies within the span cannot give the fixed effects a negative precision when the residual variance is close to zero. All six give an error for missing or infinite values in their inputs. Multiplying `K`, or the markers, by a constant only rescales `Vu`, however large or small the constant. The search for the variance ratio in `solveUVM`, `solveRRBLUP` and the multiple trait solvers spans the full range of the kernel's eigenvalues, so a large eigenvalue whose direction is absorbed by the fixed effects no longer cuts it short. The random effects stay accurate when the residual variance is close to zero, as when a model is fitted to genetic values with the QTL as markers, where `RRBLUP` now reproduces the genetic values exactly. `RRBLUP`, `RRBLUP_D`, `RRBLUP_GCA`, `RRBLUP_SCA` and multiple trait `RRBLUP` use these solvers, so their estimates can differ from earlier versions.
  * `solveRRBLUP` and `solveUVM` use the spectral approach of FaST-LMM (Lippert et al., 2011) in place of EMMA. Estimates agree with earlier versions to about seven significant figures. `solveUVM` is two to thirteen times faster when there are at least as many records as levels, and about as fast otherwise, and `solveRRBLUP` takes about as long as before. When each record belongs to a single level, `solveUVM` and `solveMVM` read the covariance of the records straight out of `K` instead of multiplying by `Z`, which made fits with 1,000 records on 5,000 levels five to ten times faster. The FaST-LMM solver that was briefly available as `solveRRBLUP2` is now `solveRRBLUP`.
  * `solveMVM` and `solveRRBLUPMV` estimate the covariance matrices by REML, following Zhou and Stephens (2014), with EM steps refined by average information steps. Each trait is fitted on a standardized scale, so recording a trait in different units only rescales its results. They converge in tens of iterations rather than hundreds, and in tests the earlier versions sometimes stopped at `maxIter` before converging, so estimates can differ from earlier versions by several percent.
  * `solveMKM` and `solveRRBLUPMK` use average information REML (Gilmour et al., 1995) in place of Fisher scoring, and keep every variance component non-negative without reporting convergence while any component, held at zero or not, could still raise the likelihood. A random term that lies entirely within the fixed effects has its variance held at zero instead of stalling the fit. They now give an error when the fixed effects leave no residual degrees of freedom, as in a two-individual `RRBLUP_SCA` fit, where earlier versions returned estimates. Estimates agree with earlier versions to about three significant figures. When there are fewer random effects than records they now work with Henderson's mixed model equations instead of the variance matrix of the records, and choose between the two by cost, so `RRBLUP_D`, `RRBLUP_GCA` and `RRBLUP_SCA` with more individuals than markers are much faster: a two-kernel fit with 5,000 records and 500 markers per kernel took 0.6 seconds instead of 47. Fits with more effects than records are about one and a half to three times faster than in earlier versions. Trial points in the line search are evaluated from a Cholesky factor alone, and Henderson's equations are solved without forming their inverse, which made fits with fewer effects than records about 20 to 30 percent faster again.
  * `RRBLUPMemUse` now estimates the memory used by the new solvers, including the choice the dominance, GCA and SCA models make between the two ways of solving and the multiple trait optimizer, whose memory grows with the fourth power of the number of traits, and the multiple trait models' fixed effects, and the genomic selection vignette describes how each model is now fitted.

* Reduced the memory used by `fastRRBLUP`. The marker genotypes are now read directly from the population's packed genotypes, instead of from a copy of every chromosome made first, so the genotypes of loci that are not markers are no longer duplicated. Results are exactly the same as before.

* Reduced the memory used, and the time taken, by functions that create or combine populations or summarize their genotypes. Every function that reads a population's genotypes now reads the packed genotypes in place, as `fastRRBLUP` and the genetic value functions do. Results are exactly the same as before, including results from a seed, and do not depend on the number of threads.
  * Crossing functions, including `makeCross`, `randCross`, `self`, `makeDH` and `reduceGenome`, write the progeny's genotypes straight into the arrays they return instead of copying them on the way out. For 2,000 diploid progeny with 100,000 segregating sites, `randCross` took 0.04 seconds instead of 0.06. The time taken by many calls that each make a few individuals is little changed, because it is spent outside the crossing code.
  * `c()` and `mergePops` copy each population's genotypes once, straight into the merged population, instead of copying the merged genotypes a second time on the way out. For populations of 2,000 diploids with 100,000 segregating sites, `c()` of two populations is 1.9 times faster and `mergePops` of four is 2.2 times faster.
  * `RRBLUP_SCA` and `RRBLUP_SCA2` read each marker's maternal and paternal alleles once instead of reading the genotypes a third time, and no longer form an unused matrix the size of the marker matrix. `RRBLUP_D` and `RRBLUP_D2` free the genotypes before the model is fitted instead of holding them through it.
  * `setEBV` with a `targetPop` counts the target's allele frequencies without forming a matrix of its genotypes, which for 2,000 individuals and 500 markers took 27 milliseconds instead of 43.

* Sped up the iterations of `fastRRBLUP`. Once the variance components and the preconditioner have been worked out, the markers are packed to one bit per allele, a quarter of the memory for a diploid, and the products with the marker matrix are taken with lookup tables over four markers at a time. For a diploid with 3,000 to 10,000 individuals and 1,000 to 10,000 markers, each iteration is 4.4 to 5.2 times faster on one thread and 2.6 to 5.8 times faster on eight. The most memory used at once is unchanged. Results do not depend on the number of threads, and differ from earlier versions in the last few digits, around 1e-13, because sums are taken in a different order.

* Sped up `genParam`, and with it `varA`, `varD`, `varAA`, `bv`, `dd`, `aa`, `genicVarA` and the other functions built on it. Genotypes are no longer decoded into a matrix of dosages: the genotypes at each locus, or for epistatic traits the joint genotypes of each pair of loci, are counted directly from the packed genotypes, and breeding values, dominance deviations and the parts of genetic value are then worked out the same way as genetic values. Traits without epistasis that share QTL, such as correlated traits added together, have their genotypes counted and read once for all of them. For 10,000 individuals, `genParam` is about five times faster with 50,000 QTL, about seven times faster for five to twenty correlated traits, and 1.6 times faster with 10,000 QTL. For epistatic traits it is about three times faster with 50,000 QTL and 1.2 to 1.7 times faster with 10,000 QTL. Genic variances, means and average effects are exactly the same as before. Breeding values, dominance deviations, epistatic deviations and genetic values can differ from earlier versions in the last few digits, around 1e-13. Because every `SimParam$addTrait` function scales a new trait using these values from the founder population, new traits' effects and intercepts can also differ in the last few digits, around 1e-14, so seeded simulations can give results that differ from earlier versions in the last few digits and, where they depend on exact ties, differ more.

* Sped up calculating genetic values, which every function that creates individuals does for each trait. Genetic values are now worked out directly from the packed genotypes instead of first decoding them into a matrix of dosages, which for 10,000 individuals and 50,000 QTL was a 477 MB matrix. Additive traits, and dominance traits in diploids, use lookup tables covering four loci at a time, and traits added together with shared QTL, such as correlated traits, are worked out in a single pass. For 10,000 progeny `randCross` is about four times faster with 50,000 QTL, two to three times faster for five to twenty correlated traits on 10,000 shared QTL, and about 1.7 times faster for epistatic traits with 50,000 QTL. Traits with few QTL are about as fast as before. Results do not depend on the number of threads. Epistatic traits, polyploid traits with dominance, and groups of traits on sparse QTL give exactly the same genetic values as before. Other traits can differ from earlier versions in the last few digits, around 1e-13, because their effects are added in a different order, so seeded results that depend on exact ties can change. `resetPop` now calculates genetic values the same way as `newPop`, so the two always agree exactly, which they did not always do before.

* Fixed multivariate `RRBLUP` and `solveMVM` stopping with "inv_sympd(): matrix is singular or not positive definite" when the relationship matrix is rank deficient, as it is when there are fewer markers than individuals or, for `solveMVM`, more than one record per level. An eigenvalue of exactly zero made the genetic variance estimate `NaN`. The genetic variance update now leaves out eigenvalues at rounding level, so the estimates from a rank deficient fit converge in fewer iterations and differ slightly from earlier versions.

* Added additive-by-additive epistasis to `SimParam$importTrait`. The new `epiPairs` and `epiEff` arguments give the interacting pairs of loci by marker name and their epistatic effects, and the trait is created as `TraitAE`, or `TraitADE` when `domEff` is also given. Every marker must be in exactly one pair, and the two markers in a pair may be on different chromosomes. Markers do not need to be in map order. Calls without the new arguments behave as before.

* Changed `SimParam$importTrait` to give an error for duplicated `markerNames`. Previously only the first of the duplicates was used, and the effects given for the others were silently ignored.

* Fixed `calcGCA` for crosses to a single tester with more than one record per cross. With one male, `GCAf` returned the individual records, with parent IDs recycled onto them so some values were labeled with the wrong female, or stopped with an error about differing numbers of rows when the crosses had unequal numbers of records. It now returns one mean per female. `GCAm` had the same problem with one female and is fixed in the same way.

* Fixed building the vignettes and articles from a directory whose path contains Markdown punctuation, such as `--` or an apostrophe. Pandoc read the bibliography path as Markdown, so `--` became a dash and the bibliography could not be found.

* Sped up recombination tracking and reduced its memory use. `pullIbdHaplo` now builds IBD haplotypes in a single compiled pass instead of one R call per ancestor, so its cost grows roughly linearly with the pedigree instead of quadratically: on 40,000 diploid individuals it took 0.5 seconds instead of 44. Recombination histories and cached IBD haplotypes are now stored as one integer vector per individual, which cut the memory used by the history about sixteenfold, and crossing with `SimParam$setTrackRec(TRUE)` is faster. `SimParam$recHist` still returns the same nested lists, and `pullIbdHaplo` returns the same results. `SimParam$haplotypes` is now ordered by internal ID rather than by when each haplotype was computed. `SimParam$setFounderHap` now also accepts an unnamed list in internal ID order, and gives an error for names that are not internal IDs. A `SimParam` with recombination tracking that was saved by an earlier version can no longer be used to track recombination or to pull IBD haplotypes.

* Added chromosome genetic lengths, so a chromosome can extend past its last locus. A chromosome now runs from position 0 to its genetic length, crossovers are simulated over all of it, and those after the last locus affect no locus. The length is stored in a new `genLen` slot of `MapPop` and read from `SimParam` with `$genLen`, `$femaleGenLen` and `$maleGenLen`. `runMacs` and `runMacs2` keep the genetic length requested, `quickHaplo` keeps `genLen`, and `newMapPop`, `importHaplo`, `importInbredGeno`, `SimParam$switchGenMap`, `SimParam$switchFemaleMap` and `SimParam$switchMaleMap` take an optional `genLen` argument that defaults to the position of the last locus. The default centromere is placed at half of the genetic length, and `setRecombRatio` scales the lengths with the maps. `addSegSite` now gives an error for a site past the end of the chromosome unless `extendGenLen=TRUE`. A `MapPop` or `SimParam` saved by an earlier version ends each chromosome at its last locus, as before. Results from a given seed are unchanged when each chromosome ends at its last locus, and differ otherwise, including for every `runMacs` founder population.

* Fixed recombination tracking for quadrivalents in polyploids, which made `pullIbdHaplo` stop with an index out of bounds error, or return the wrong founder haplotypes, after crossing or `reduceGenome` when `SimParam$quadProb` was greater than zero. Simulated genotypes were not affected.

* Changed genetic maps to use absolute positions, measured from the start of the chromosome. The first locus is no longer moved to position 0 by `runMacs`, `runMacs2`, `importGenMap`, `importHaplo`, `importInbredGeno`, `addSegSite`, `SimParam$new`, `SimParam$switchGenMap`, `SimParam$switchFemaleMap` or `SimParam$switchMaleMap`. Meiosis now samples crossovers from position 0, so crossovers can fall before the first locus, where they change which parental chromosome a gamete starts on. Recombination between loci is unaffected for bivalents, but the default centromere and the quadrivalent exchange point are now placed relative to position 0. `runMacs` and `runMacs2` now place each centromere at half of the requested genetic length, for example 0.5 Morgans for a 1 Morgan chromosome, instead of half of the position of the last site returned. Maps with negative, unsorted or non-finite positions, and negative centromeres, are now an error instead of being shifted. Results from a given seed differ from earlier versions for any map whose first locus is not at 0, including every `runMacs` founder population.

* Added `restrInbr` to `selectCross`, which selects parents with a restriction on the expected increase in fixation, an approximation to optimal contribution selection with equal contributions. Fixation is measured at a SNP chip or at a trait's QTL, chosen with `snpChip` and `useQtl` as in `RRBLUP`. The target is set with `inbrTarget` and `inbrType`, either relative to the current population or as an absolute maximum. The default is `FALSE`, which leaves `selectCross` unchanged. The new function `selectOCS` makes the same selection without crossing, so the selected individuals can be crossed with a custom crossing plan, for example with `makeCross`.

* Fixed several places where invalid input could read or write memory outside a population's genotypes, or outside internal MaCS structures, instead of giving an error. These now stop with an error:
  * `haplo` in `pullSnpHaplo`, `pullQtlHaplo`, `pullSegSiteHaplo` and `pullMarkerHaplo` must be `"all"` or a single whole number from 1 to the ploidy level.
  * A `LociMap`, such as a trait passed to `SimParam$manAddTrait`, whose loci fall outside their chromosome.
  * `parents` in `self` and `crossPlan` in `hybridCross` with `returnHybridPop=TRUE` must index individuals in the population. `hybridCross` now also accepts a `crossPlan` of IDs for `returnHybridPop=TRUE`, as `makeCross2` does.
  * `segSites` in `runMacs` must have length 1 or `nChr`.
  * Population IDs in a `runMacs` `manualCommand` (the `-n`, `-g` and `-m` options and the `-en`, `-eg`, `-es`, `-ej` and `-em` events) must name a population that exists, and the number of populations given to `-I` must be a positive whole number.

* Changed MaCS to pass its messages to R after the chromosomes have been simulated, instead of writing them from worker threads, which is not safe. A warning about the command is now shown once rather than once per chromosome, and when MaCS fails the error includes what MaCS reported.

* Added a vignette on simulation of founder genomes from demographic history with MaCS presets, custom MaCS calls, and msprime through reticulate. Clarified parameter values in demographic models and documented segregating-site counts from simulations and theory, and reported laptop timings. We also provide a function `predictSegSitesMoments()` for theoretical predictions.

* Added `mutRate` to `runMacs()` to override a preset mutation rate. Expected segregating-site counts scale with the rate; reducing it reduces sequence diversity and memory usage, with runtime savings limited by ancestry simulation.

* Clarified `runMacs()` and `runMacs2()` documentation for site retention, mutation rates, historical population sizes, and the interaction between ploidy and inbred founders. Only `runMacs2()` rescales effective population sizes by ploidy.

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

* Update `meanG()`, `meanP()` and `meanEBV()` to support MultiPop objects. The output can be simplified via `simplify` and `level` arguments, which are passed to `calcPopValue()`. These functions calculate mean values within Pop objects.

* Add `meanGPop()` and `meanPPop()` to calculate mean genetic and phenotype values between Pops in a MultiPop.

* Update `varG()`, `varP()` and `varEBV()` to support MultiPop objects.

* Add `varGPop()` and `varPPop()` to calculate genetic and phenotype variances between populations in a MultiPop.

* Update `gv()`, `pheno()` and `ebv()` to support MultiPop objects.

* Add `nPop()` to count the number of Pop-class objects in a MultiPop.

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
