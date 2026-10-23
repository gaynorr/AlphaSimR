# Help pages for the methods of SimParam, one page per method.
#
# roxygen2 would otherwise put every method of an R6 class on the class's
# own page, which made that page very long and could only link to its
# methods by anchors within the page, which not every help viewer
# follows. DESCRIPTION therefore sets "Roxygen: list(r6 = FALSE)", and
# each method is documented here as a topic of its own. The methods
# themselves are in R/Class-SimParam.R, and their usage is generated from
# their formals by .simParamUsage, so it cannot drift from the code. The
# test in tests/testthat/test-SimParamDocs.R checks that every method a
# user can call has a page whose arguments match its formals.

# Returns an Rd usage section for a SimParam method. A call such as
# SP$addTraitA() is not one R CMD check can match to a function in the
# package, so the usage is wrapped in \special, which the check leaves
# alone. That keeps the usage in its usual place, ahead of the arguments.
.simParamUsage = function(method, call=paste0("SP$", method)){
  args = formals(SimParam$public_methods[[method]])
  txt = vapply(names(args), function(n){
    d = deparse(args[[n]])
    if(!nzchar(d[1L])) n else paste0(n, " = ", paste(d, collapse=" "))
  }, character(1))
  # Escape the characters that are special in Rd
  txt = gsub("([%{}\\\\])", "\\\\\\1", txt)
  if(length(txt)==0L){
    usage = paste0(call, "()")
  }else{
    usage = paste0(call, "(\n  ", paste(txt, collapse=",\n  "), "\n)")
  }
  paste0("\\usage{\\special{", usage, "}}")
}

#' Set how sexes are assigned
#'
#' @description Changes how sexes are determined in the simulation.
#' The default value is "no", indicating all individuals are hermaphrodites.
#' To add sexes to the simulation, run this function with "yes_sys" or
#' "yes_rand". The value "yes_sys" will systematically assign
#' sexes to newly created individuals as first male and then female.
#' Populations with an odd number of individuals will have one more male than
#' female. The value "yes_rand" will randomly assign a sex to each
#' individual.
#'
#' @evalRd .simParamUsage("setSexes")
#'
#' @param sexes acceptable values are "no", "yes_sys", or
#' "yes_rand"
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$setSexes("yes_sys")
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-setSexes
#' @aliases SimParam$setSexes
NULL

#' Turn pedigree tracking on or off
#'
#' @description Sets pedigree tracking for the simulation.
#' By default pedigree tracking is turned off. When turned on,
#' the pedigree of all individuals created will be tracked,
#' except those created by \code{\link{hybridCross}}. Turning
#' off pedigree tracking will turn off recombination tracking
#' if it is turned on.
#'
#' @evalRd .simParamUsage("setTrackPed")
#'
#' @param isTrackPed should pedigree tracking be on.
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$setTrackPed(TRUE)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-setTrackPed
#' @aliases SimParam$setTrackPed
NULL

#' Turn recombination tracking on or off
#'
#' @description Sets recombination tracking for the simulation.
#' By default recombination tracking is turned off. When turned
#' on recombination tracking will also turn on pedigree tracking.
#' Recombination tracking keeps records of all individuals created,
#' except those created by \code{\link{hybridCross}}, because their
#' pedigree is not tracked.
#'
#' @evalRd .simParamUsage("setTrackRec")
#'
#' @param isTrackRec should recombination tracking be on.
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$setTrackRec(TRUE)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-setTrackRec
#' @aliases SimParam$setTrackRec
NULL

#' Reset the pedigree
#'
#' @description Resets the internal lastId, the pedigree
#' and recombination tracking (if in use) to the
#' supplied lastId. Be careful using this function because
#' it may introduce a bug if you use individuals from
#' the deleted portion of the pedigree.
#'
#' @evalRd .simParamUsage("resetPed")
#'
#' @param lastId last ID to include in pedigree
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@id # 1:10
#'
#' #Create another population after reseting pedigree
#' SP$resetPed()
#' pop2 = newPop(founderPop, simParam=SP)
#' pop2@id # 1:10
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-resetPed
#' @aliases SimParam$resetPed
NULL

#' Add additive traits
#'
#' @description
#' Randomly assigns eligible QTLs for one or more additive traits.
#' If simulating more than one trait, all traits will be pleiotropic
#' with correlated additive genotypic effects.
#'
#' @evalRd .simParamUsage("addTraitA")
#'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of requested founder means (mean genetic values) for one or more traits
#' @param var a vector of desired genetic variances in the founder
#'   population for one or more traits
#' @param corA a matrix of correlations between additive genotypic effects
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#'   (the rate/scale parameter of the gamma distribution is accounted
#'   for via the desired level of genetic variance, the var argument)
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait(s)
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addTraitA
#' @aliases SimParam$addTraitA
NULL

#' Add additive and dominance traits
#'
#' @description
#' Randomly assigns eligible QTLs for one or more traits with dominance.
#' If simulating more than one trait, all traits will be pleiotropic
#' with correlated effects.
#'
#' @evalRd .simParamUsage("addTraitAD")
#'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of requested founder means (mean genetic values) for one or more traits
#' @param var a vector of desired genetic variances in the founder
#'   population for one or more traits
#' @param meanDD mean dominance degree
#' @param varDD variance of dominance degree
#' @param corA a matrix of correlations between additive genotypic effects
#' @param corDD a matrix of correlations between dominance degrees
#' @param useVarA if \code{TRUE}, scale genotypic effects to achieve
#'   the requested additive genetic variance (variance of breeding values)
#'   in the founder population.
#'   If \code{FALSE}, use total genetic variance (variance of genetic values).
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#'   (the rate/scale parameter of the gamma distribution is accounted
#'   for via the desired level of genetic variance, the var argument)
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait(s)
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAD(10, meanDD=0.5)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addTraitAD
#' @aliases SimParam$addTraitAD
NULL

#' Add additive and epistatic traits
#'
#' @description
#' Randomly assigns eligible QTLs for one or more additive and epistasis
#' traits. If simulating more than one trait, all traits will be pleiotropic
#' with correlated additive genotypic effects.
#'
#' @evalRd .simParamUsage("addTraitAE")
#'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of requested founder means (mean genetic values) for one or more traits
#' @param var a vector of desired genetic variances in the founder
#'   population for one or more traits
#' @param relAA ratio of additive-by-additive epistatic genetic variance
#'   to additive genetic variance used to sample epistatic genotypic effects,
#'   assuming a diploid population in Hardy-Weinberg and linkage equilibrium
#'   with allele frequencies of 0.5.
#'   The realized variance ratio in the founder population can differ.
#' @param corA a matrix of correlations between additive genotypic effects
#' @param corAA a matrix of correlations between additive-by-additive epistatic genotypic effects
#' @param useVarA if \code{TRUE}, scale genotypic effects to achieve
#'   the requested additive genetic variance (variance of breeding values)
#'   in the founder population.
#'   If \code{FALSE}, use total genetic variance (variance of genetic values).
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#'   (the rate/scale parameter of the gamma distribution is accounted
#'   for via the desired level of genetic variance, the var argument)
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait(s)
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAE(10, relAA=0.1)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addTraitAE
#' @aliases SimParam$addTraitAE
NULL

#' Add additive, dominance and epistatic traits
#'
#' @description
#' Randomly assigns eligible QTLs for one or more traits with dominance and
#' epistasis. If simulating more than one trait, all traits will be pleiotropic
#' with correlated effects.
#'
#' @evalRd .simParamUsage("addTraitADE")
#'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of requested founder means (mean genetic values) for one or more traits
#' @param var a vector of desired genetic variances in the founder
#'   population for one or more traits
#' @param meanDD mean dominance degree
#' @param varDD variance of dominance degree
#' @param relAA ratio of additive-by-additive epistatic genetic variance
#'   to additive genetic variance used to sample epistatic genotypic effects,
#'   assuming a diploid population in Hardy-Weinberg and linkage equilibrium
#'   with allele frequencies of 0.5.
#'   The realized variance ratio in the founder population can differ.
#' @param corA a matrix of correlations between additive genotypic effects
#' @param corDD a matrix of correlations between dominance degrees
#' @param corAA a matrix of correlations between additive-by-additive epistatic genotypic effects
#' @param useVarA if \code{TRUE}, scale genotypic effects to achieve
#'   the requested additive genetic variance (variance of breeding values)
#'   in the founder population.
#'   If \code{FALSE}, use total genetic variance (variance of genetic values).
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#'   (the rate/scale parameter of the gamma distribution is accounted
#'   for via the desired level of genetic variance, the var argument)
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait(s)
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitADE(10)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addTraitADE
#' @aliases SimParam$addTraitADE
NULL

#' Add additive traits with GxE
#'
#' @description
#' Randomly assigns eligible QTLs for one or more additive GxE traits.
#' If simulating more than one trait, all traits will be pleiotropic
#' with correlated effects.
#'
#' @evalRd .simParamUsage("addTraitAG")
#'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of requested founder means (mean genetic values) for one or more traits
#' @param var a vector of desired genetic variances in the founder
#'   population for one or more traits
#' @param varGxE a vector of total genotype-by-environment variances for the traits
#' @param varEnv a vector of environmental variances for one or more traits
#' @param corA a matrix of correlations between additive genotypic effects
#' @param corGxE a matrix of correlations between genotype-by-environment genotypic effects
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#'   (the rate/scale parameter of the gamma distribution is accounted
#'   for via the desired level of genetic variance, the var argument)
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait(s)
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAG(10, varGxE=2)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addTraitAG
#' @aliases SimParam$addTraitAG
NULL

#' Add additive and dominance traits with GxE
#'
#' @description
#' Randomly assigns eligible QTLs for a trait with dominance and GxE.
#'
#' @evalRd .simParamUsage("addTraitADG")
#'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single
#' value or nChr values.
#' @param mean a vector of requested founder means (mean genetic values) for one or more traits
#' @param var a vector of desired genetic variances in the founder
#'   population for one or more traits
#' @param varGxE a vector of total genotype-by-environment variances for the traits
#' @param varEnv a vector of environmental variances for one or more traits
#' @param meanDD mean dominance degree
#' @param varDD variance of dominance degree
#' @param corA a matrix of correlations between additive genotypic effects
#' @param corDD a matrix of correlations between dominance degrees
#' @param corGxE a matrix of correlations between genotype-by-environment genotypic effects
#' @param useVarA if \code{TRUE}, scale genotypic effects to achieve
#'   the requested additive genetic variance (variance of breeding values)
#'   in the founder population.
#'   If \code{FALSE}, use total genetic variance (variance of genetic values).
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#'   (the rate/scale parameter of the gamma distribution is accounted
#'   for via the desired level of genetic variance, the var argument)
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait(s)
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitADG(10, meanDD=0.5, varGxE=2)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addTraitADG
#' @aliases SimParam$addTraitADG
NULL

#' Add additive and epistatic traits with GxE
#'
#' @description
#' Randomly assigns eligible QTLs for one or more additive and epistasis
#' GxE traits. If simulating more than one trait, all traits will be pleiotropic
#' with correlated effects.
#'
#' @evalRd .simParamUsage("addTraitAEG")
#'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean a vector of requested founder means (mean genetic values) for one or more traits
#' @param var a vector of desired genetic variances in the founder
#'   population for one or more traits
#' @param relAA ratio of additive-by-additive epistatic genetic variance
#'   to additive genetic variance used to sample epistatic genotypic effects,
#'   assuming a diploid population in Hardy-Weinberg and linkage equilibrium
#'   with allele frequencies of 0.5.
#'   The realized variance ratio in the founder population can differ.
#' @param varGxE a vector of total genotype-by-environment variances for the traits
#' @param varEnv a vector of environmental variances for one or more traits
#' @param corA a matrix of correlations between additive genotypic effects
#' @param corAA a matrix of correlations between additive-by-additive epistatic genotypic effects
#' @param corGxE a matrix of correlations between genotype-by-environment genotypic effects
#' @param useVarA if \code{TRUE}, scale genotypic effects to achieve
#'   the requested additive genetic variance (variance of breeding values)
#'   in the founder population.
#'   If \code{FALSE}, use total genetic variance (variance of genetic values).
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#'   (the rate/scale parameter of the gamma distribution is accounted
#'   for via the desired level of genetic variance, the var argument)
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait(s)
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAEG(10, varGxE=2)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addTraitAEG
#' @aliases SimParam$addTraitAEG
NULL

#' Add additive, dominance and epistatic traits with GxE
#'
#' @description
#' Randomly assigns eligible QTLs for a trait with dominance,
#' epistasis and GxE.
#'
#' @evalRd .simParamUsage("addTraitADEG")
#'
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single
#' value or nChr values.
#' @param mean a vector of requested founder means (mean genetic values) for one or more traits
#' @param var a vector of desired genetic variances in the founder
#'   population for one or more traits
#' @param varGxE a vector of total genotype-by-environment variances for the traits
#' @param varEnv a vector of environmental variances for one or more traits
#' @param meanDD mean dominance degree
#' @param varDD variance of dominance degree
#' @param relAA ratio of additive-by-additive epistatic genetic variance
#'   to additive genetic variance used to sample epistatic genotypic effects,
#'   assuming a diploid population in Hardy-Weinberg and linkage equilibrium
#'   with allele frequencies of 0.5.
#'   The realized variance ratio in the founder population can differ.
#' @param corA a matrix of correlations between additive genotypic effects
#' @param corDD a matrix of correlations between dominance degrees
#' @param corAA a matrix of correlations between additive-by-additive epistatic genotypic effects
#' @param corGxE a matrix of correlations between genotype-by-environment genotypic effects
#' @param useVarA if \code{TRUE}, scale genotypic effects to achieve
#'   the requested additive genetic variance (variance of breeding values)
#'   in the founder population.
#'   If \code{FALSE}, use total genetic variance (variance of genetic values).
#' @param gamma should a gamma distribution be used instead of normal
#' @param shape the shape parameter for the gamma distribution
#'   (the rate/scale parameter of the gamma distribution is accounted
#'   for via the desired level of genetic variance, the var argument)
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait(s)
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitADEG(10, meanDD=0.5, varGxE=2)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addTraitADEG
#' @aliases SimParam$addTraitADEG
NULL

#' Add an additive and dominance trait by an alternative method
#'
#' @description
#' An alternative method for adding a trait with additive and dominance genotypic effects
#' to an AlphaSimR simulation. The function attempts to create a trait matching
#' user defined values for number of QTL, inbreeding depression, additive genetic
#' variance and dominance genetic variance.
#'
#' @evalRd .simParamUsage("altAddTraitAD")
#'
#' @param nQtlPerChr number of QTLs per chromosome.
#' Can be a single value or nChr values.
#' @param mean requested founder mean (mean genetic value) of the trait
#' @param varA desired additive variance
#' @param varD desired dominance variance
#' @param inbrDepr desired inbreeding depression, see details
#' @param limMeanDD limits for meanDD, see details
#' @param limVarDD limits for varDD, see details
#' @param silent should summary details be printed to the console
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#' @param name optional name for trait
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @details
#' This function will always add a trait to 'SimParam', unless an error occurs
#' with picking QTLs. The resulting trait will always have the requested founder mean and
#' additive genetic variance. However, it may not have the desired values for
#' inbreeding depression and dominance variance. Thus, it is strongly recommended
#' to check the output printed to the console to determine how close the trait's
#' parameters came to these desired values.
#'
#' The mean and additive genetic variance will always be achieved exactly. The
#' function attempts to achieve the desired dominance variance and inbreeding
#' depression while staying within the user supplied constraints for the
#' acceptable range of dominance degree mean and variance. If the desired values
#' are not being achieved, the acceptable range needs to be increased and/or the
#' number of QTL may need to be increased. There are no limits to setting the
#' range for dominance degree mean and variance, but care should be taken
#' with regard to the biological feasibility of the limits that are supplied.
#' The default limits were somewhat arbitrarily set, so I make no claim as to
#' how reasonable these limits are for routine use.
#'
#' Inbreeding depression in this function is defined as the difference in mean
#' genetic value between a population with the same allele frequency as the
#' reference population (population used to initialize SimParam) in
#' Hardy-Weinberg equilibrium compared to a population with the same allele
#' frequency that is fully inbred. This is equivalent to the amount the mean of
#' a population increases when going from an inbreeding coefficient of 1 (fully
#' inbred) to a population with an inbreeding coefficient of 0 (Hardy-Weinberg
#' equilibrium). Note that the sign of the value should (usually) be positive.
#' This corresponds to a detrimental effect of inbreeding when higher values of
#' the trait are considered biologically beneficial.
#'
#' Summary information on this trait is printed to the console when silent=FALSE.
#' The summary information reports the inbreeding depression and dominance
#' variance for the population as well as the dominance degree mean and variance
#' applied to the trait.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=100, nChr=1, segSites=100)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$altAddTraitAD(nQtlPerChr=10, mean=0, varA=1, varD=0.05, inbrDepr=0.2)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-altAddTraitAD
#' @aliases SimParam$altAddTraitAD
NULL

#' Add a trait from a LociMap
#'
#' @description
#' Manually add a new trait to the simulation. Trait must
#' be formatted as a \code{\link{LociMap-class}}. If the
#' trait is not already formatted, consider using importTrait.
#'
#' @evalRd .simParamUsage("manAddTrait")
#'
#' @param lociMap a new object descended from
#' \code{\link{LociMap-class}}
#' @param varE default error variance for phenotype, optional
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-manAddTrait
#' @aliases SimParam$manAddTrait
NULL

#' Import traits from QTL effects
#'
#' @description
#' Manually add a new trait(s) to the simulation. Unlike the
#' manAddTrait function, this function does not require
#' formatting the trait as a \code{\link{LociMap-class}}.
#' The formatting is performed automatically for the user,
#' with more user friendly data.frames or matrices taken as
#' inputs. This function works for A, AD, AE and ADE trait types.
#'
#' @evalRd .simParamUsage("importTrait")
#'
#' @param markerNames a vector of unique names for the QTL
#' @param addEff a matrix of additive genotypic effects (nLoci x nTraits).
#' Alternatively, a vector of length nLoci can be supplied for
#' a single trait.
#' @param domEff optional dominance genotypic effects for each locus
#' (nLoci x nTraits, or a vector of length nLoci for a single trait)
#' @param intercept optional trait intercepts in the genotypic parameterization
#' @param name optional name(s) for the trait(s)
#' @param varE default error variance for phenotype, optional
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#' @param epiPairs optional two-column matrix or data.frame of marker
#' names, with one row for each pair of loci that interact through
#' additive-by-additive epistasis. Must be supplied with \code{epiEff}.
#' @param epiEff optional additive-by-additive epistatic genotypic
#' effects, a matrix with one row for each row of \code{epiPairs} and
#' one column for each trait. Alternatively, a vector of length
#' \code{nrow(epiPairs)} can be supplied for a single trait.
#'
#' @details See \code{vignette("traits", package = "AlphaSimR")} for the
#' implemented trait model and genotypic effects.
#'
#' Supplying \code{epiPairs} and \code{epiEff} creates an AE trait,
#' or an ADE trait when \code{domEff} is also supplied. The epistatic
#' model only allows interactions between pairs of loci, and every
#' marker in \code{markerNames} must appear in exactly one pair, so
#' \code{nrow(epiPairs)} is half of \code{length(markerNames)}. The two
#' markers in a pair may be on different chromosomes. All traits
#' imported in one call share the same pairs, with effects taken from
#' the corresponding column of \code{epiEff}. Markers do not need to be
#' given in map order; the pairs are matched to the loci by name after
#' the loci are put in map order.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=2, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Import an AE trait, with one pair spanning both chromosomes
#' SP$importTrait(markerNames=c("1_1","1_5","2_3","2_8"),
#'                addEff=c(1,-1,0.5,0.5),
#'                epiPairs=rbind(c("1_1","1_5"),
#'                               c("2_3","2_8")),
#'                epiEff=c(0.2,-0.4),
#'                name="AE")
#'
#' #Import an ADE trait
#' SP$importTrait(markerNames=c("1_2","1_4","2_1","2_6"),
#'                addEff=c(1,-1,0.5,0.5),
#'                domEff=c(0.3,0.1,0,0.2),
#'                epiPairs=rbind(c("1_2","2_6"),
#'                               c("1_4","2_1")),
#'                epiEff=c(0.5,0.1),
#'                name="ADE")
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-importTrait
#' @aliases SimParam$importTrait
NULL

#' Set default error variances
#'
#' @description Defines default values for error
#' variances used in \code{\link{setPheno}}. These defaults
#' will be used to automatically generate phenotypes when new
#' populations are created. See the details section of \code{\link{setPheno}}
#' for more information about each argument and how they
#' should be used.
#'
#' @evalRd .simParamUsage("setVarE")
#'
#' @param h2 a vector of desired narrow-sense heritabilities
#' @param H2 a vector of desired broad-sense heritabilities
#' @param varE a vector or matrix of error variances
#' @param corE an optional matrix of error correlations
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-setVarE
#' @aliases SimParam$setVarE
NULL

#' Set default error correlations
#'
#' @description Defines a correlation structure for default
#' error variances. You must call \code{setVarE} first to define
#' the default error variances.
#'
#' @evalRd .simParamUsage("setCorE")
#'
#' @param corE a correlation matrix for the error variances
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10, mean=c(0,0), var=c(1,1), corA=diag(2))
#' SP$setVarE(varE=c(1,1))
#' E = 0.5*diag(2)+0.5 #Positively correlated error
#' SP$setCorE(E)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-setCorE
#' @aliases SimParam$setCorE
NULL

#' Restrict which segregating sites can be SNPs or QTL
#'
#' @description Sets restrictions on which segregating sites can
#' serve as a SNP and/or QTL. The default behavior of AlphaSimR is
#' to randomly sample QTL or SNP from all eligible sites and then
#' mark the sampled sites ineligible to be sampled as the other
#' type (e.g. if a site is sampled as a QTL it will be marked as
#' ineligible to be sampled as a SNP). This behavior is designed
#' to produce the most challenging scenario for genomic selection
#' when the markers used for prediction are not causal.
#'
#' Setting overlap=TRUE will prevent the addition of loci to the
#' ineligible list, but it won't remove sites already added to these
#' lists. Thus, timing of when restrSegSites is called matters. It
#' should be called before any addTrait or addSnpChip functions with
#' the overlap=TRUE argument to freely allow loci to overlap.
#'
#' The minQtlPerChr and minSnpPerChr arguments can be used with
#' overlap=FALSE to preallocate sites as QTL and SNP respectively. This
#' option is useful when simulating multiple traits and/or SNP chips,
#' because it can be used to guarantee that enough eligible sites are
#' available when running addTrait and or addSnpChip functions.
#'
#' @evalRd .simParamUsage("restrSegSites")
#'
#' @param minQtlPerChr the minimum number of segregating sites for
#' QTLs. Can be a single value or a vector of values for each chromosome.
#' @param minSnpPerChr the minimum number of segregating sites for SNPs.
#' Can be a single value or a vector of values for each chromosome.
#' @param excludeQtl an optional vector of segregating site names to
#' exclude from consideration as a viable QTL.
#' @param excludeSnp an optional vector of segregating site names to
#' exclude from consideration as a viable SNP.
#' @param overlap should SNP and QTL sites be allowed to overlap.
#' @param minSnpFreq minimum allowable frequency for SNP loci.
#' No minimum SNP frequency is used if value is NULL.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$restrSegSites(minQtlPerChr=5, minSnpPerChr=5)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-restrSegSites
#' @aliases SimParam$restrSegSites
NULL

#' Add a SNP chip
#'
#' @description
#' Randomly assigns eligible SNPs to a SNP chip
#'
#' @evalRd .simParamUsage("addSnpChip")
#'
#' @param nSnpPerChr number of SNPs per chromosome.
#' Can be a single value or nChr values.
#' @param minSnpFreq minimum allowable frequency for SNP loci.
#' If NULL, no minimum frequency is used.
#' @param refPop reference population for calculating SNP
#' frequency. If NULL, the founder population is used.
#' @param name optional name for chip
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addSnpChip(10)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addSnpChip
#' @aliases SimParam$addSnpChip
NULL

#' Add a SNP chip by marker name
#'
#' @description
#' Assigns SNPs to a SNP chip by supplying marker names. This function does
#' check against excluded SNPs and will not add the SNPs to the list of
#' excluded QTL for the purpose of avoiding overlap between SNPs and QTL.
#' Excluding these SNPs from being used as QTL can be accomplished using
#' the excludeQtl argument in SimParam's restrSegSites function.
#'
#' @evalRd .simParamUsage("addSnpChipByName")
#'
#' @param markers a vector of names for the markers
#' @param name optional name for chip
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addSnpChipByName(c("1_1","1_3"))
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addSnpChipByName
#' @aliases SimParam$addSnpChipByName
NULL

#' Add structured SNP chips
#'
#' @description
#' Randomly selects the number of snps in structure and then
#' assigns them to chips based on structure
#'
#' @evalRd .simParamUsage("addStructuredSnpChip")
#'
#' @param nSnpPerChr number of SNPs per chromosome.
#' Can be a single value or nChr values.
#' @param structure a matrix.  Rows are snp chips, columns are chips.
#' If value is true then that snp is on that chip.
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing.
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-addStructuredSnpChip
#' @aliases SimParam$addStructuredSnpChip
NULL

#' Set founder haplotypes
#'
#' @description
#' Allows for the manual setting of founder haplotypes. This functionality
#' is not fully documented, because it is still experimental.
#'
#' @evalRd .simParamUsage("setFounderHap")
#'
#' @param hapMap a list of founder haplotypes, one per individual in
#' the order of their internal IDs, or named by internal ID
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-setFounderHap
#' @aliases SimParam$setFounderHap
NULL

#' Set the female to male recombination ratio
#'
#' @description Set the relative recombination rates between males
#' and females. This allows for sex-specific recombination rates,
#' under the assumption of equivalent recombination landscapes.
#'
#' @evalRd .simParamUsage("setRecombRatio")
#'
#' @param femaleRatio relative ratio of recombination in females compared to
#' males. A value of 2 indicates twice as much recombination in females. The
#' value must be greater than 0. (default is 1)
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$setRecombRatio(2) #Twice as much recombination in females
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-setRecombRatio
#' @aliases SimParam$setRecombRatio
NULL

#' Replace the genetic map
#'
#' @description
#' Replaces existing genetic map.
#'
#' @evalRd .simParamUsage("switchGenMap")
#'
#' @param genMap a list of length nChr containing
#' numeric vectors for the absolute position of each segregating
#' site on a chromosome.
#' @param centromere a numeric vector of centromere
#' positions. If NULL, the centromeres are assumed to
#' be metacentric, sitting at half of each chromosome's
#' genetic length.
#' @param genLen an optional vector of chromosome genetic
#' lengths in Morgans, either one value or one per chromosome.
#' Each must be at least the position of the chromosome's
#' last locus. If NULL, each chromosome ends at its last locus.
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-switchGenMap
#' @aliases SimParam$switchGenMap
NULL

#' Replace the female genetic map
#'
#' @description
#' Replaces existing female genetic map.
#'
#' @evalRd .simParamUsage("switchFemaleMap")
#'
#' @param genMap a list of length nChr containing
#' numeric vectors for the absolute position of each segregating
#' site on a chromosome.
#' @param centromere a numeric vector of centromere
#' positions. If NULL, the centromeres are assumed to
#' be metacentric, sitting at half of each chromosome's
#' genetic length.
#' @param genLen an optional vector of chromosome genetic
#' lengths in Morgans, either one value or one per chromosome.
#' Each must be at least the position of the chromosome's
#' last locus. If NULL, each chromosome ends at its last locus.
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-switchFemaleMap
#' @aliases SimParam$switchFemaleMap
NULL

#' Replace the male genetic map
#'
#' @description
#' Replaces existing male genetic map.
#'
#' @evalRd .simParamUsage("switchMaleMap")
#'
#' @param genMap a list of length nChr containing
#' numeric vectors for the absolute position of each segregating
#' site on a chromosome.
#' @param centromere a numeric vector of centromere
#' positions. If NULL, the centromeres are assumed to
#' be metacentric, sitting at half of each chromosome's
#' genetic length.
#' @param genLen an optional vector of chromosome genetic
#' lengths in Morgans, either one value or one per chromosome.
#' Each must be at least the position of the chromosome's
#' last locus. If NULL, each chromosome ends at its last locus.
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-switchMaleMap
#' @aliases SimParam$switchMaleMap
NULL

#' Replace a trait
#'
#' @description
#' Switch a trait in the simulation.
#'
#' @evalRd .simParamUsage("switchTrait")
#'
#' @param traitPos an integer or a trait name indicating which trait
#' to switch
#' @param lociMap a new object descended from
#' \code{\link{LociMap-class}}
#' @param varE default error variance for phenotype, optional
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-switchTrait
#' @aliases SimParam$switchTrait
NULL

#' Remove a trait
#'
#' @description
#' Remove a trait from the simulation
#'
#' @evalRd .simParamUsage("removeTrait")
#'
#' @param traits an integer vector or a character vector of trait
#' names indicating which traits to remove
#' @param force should the check for a running simulation be
#' ignored. Only set to TRUE if you know what you are doing
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-removeTrait
#' @aliases SimParam$removeTrait
NULL

#' Rescale traits
#'
#' @description
#' Linearly scales all traits to achieve desired
#' values of means and variances in the founder population.
#'
#' @evalRd .simParamUsage("rescaleTraits")
#'
#' @param mean a vector of requested mean genetic values in the founder population
#' @param var a vector of requested genetic variances in the founder population
#' @param varEnv a vector of new environmental variances
#' @param varGxE a vector of new GxE variances
#' @param useVarA if \code{TRUE}, scale genotypic effects to achieve
#'   the requested additive genetic variance (variance of breeding values)
#'   in the founder population.
#'   If \code{FALSE}, use total genetic variance (variance of genetic values).
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{self$nThreads}.
#'
#' @note
#' By default the founder population is the population used to
#' initialize the SimParam object. This population can be changed by
#' replacing the population in the founderPop slot. You must run
#' \code{\link{resetPop}} on any existing populations to obtain the
#' new trait values.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' meanG(pop)
#'
#' #Change mean to 1
#' SP$rescaleTraits(mean=1)
#' #Run resetPop for change to take effect
#' pop = resetPop(pop, simParam=SP)
#' meanG(pop)
#'
#' @seealso \code{\link{SimParam}}, which lists every method
#'
#' @name SimParam-rescaleTraits
#' @aliases SimParam$rescaleTraits
NULL
