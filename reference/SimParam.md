# Simulation parameters

Container for global simulation parameters. Saving this object as SP
will allow it to be accessed by function defaults.

## Note

By default the founder population is the population used to initalize
the SimParam object. This population can be changed by replacing the
population in the founderPop slot. You must run
[`resetPop`](https://gaynorr.github.io/AlphaSimR/reference/resetPop.md)
on any existing populations to obtain the new trait values.

## Public fields

- `nThreads`:

  number of threads used on platforms with OpenMP support

- `snpChips`:

  list of SNP chips

- `invalidQtl`:

  list of segregating sites that aren't valid QTL

- `invalidSnp`:

  list of segregating sites that aren't valid SNP

- `founderPop`:

  founder population used for variance scaling

- `finalizePop`:

  function applied to newly created populations. Currently does nothing
  and should only be changed by expert users.

- `allowEmptyPop`:

  if true, population arguments with nInd=0 will return an empty
  population with a warning instead of an error.

- `v`:

  the crossover interference parameter for a gamma model of
  recombination. A value of 1 indicates no crossover interference (e.g.
  Haldane mapping function). A value of 2.6 approximates the degree of
  crossover interference implied by the Kosambi mapping function.
  (default is 2.6)

- `p`:

  the proportion of crossovers coming from a non-interfering pathway.
  (default is 0)

- `quadProb`:

  the probability of quadrivalent pairing in an autopolyploid. (default
  is 0)

## Active bindings

- `traitNames`:

  vector of trait names

- `snpChipNames`:

  vector of chip names

- `traits`:

  list of traits

- `nChr`:

  number of chromosomes

- `nTraits`:

  number of traits

- `nSnpChips`:

  number of SNP chips

- `segSites`:

  segregating sites per chromosome

- `sexes`:

  sexes used for mating

- `sepMap`:

  are there seperate genetic maps for males and females

- `genMap`:

  list of chromosome genetic maps

- `femaleMap`:

  list of chromosome genetic maps for females

- `maleMap`:

  list of chromosome genetic maps for males

- `centromere`:

  position of centromeres genetic map

- `femaleCentromere`:

  position of centromeres on female genetic map

- `maleCentromere`:

  position of centromeres on male genetic map

- `lastId`:

  last ID number assigned

- `isTrackPed`:

  is pedigree being tracked

- `pedigree`:

  pedigree matrix for all individuals

- `isTrackRec`:

  is recombination being tracked

- `recHist`:

  list of historic recombination events

- `haplotypes`:

  list of computed IBD haplotypes

- `varA`:

  additive genetic variance in founderPop

- `varG`:

  total genetic variance in founderPop

- `varE`:

  default error variance

- `version`:

  the version of AlphaSimR used to generate this object

- `activeQtl`:

  a LociMap representing all active QTL in simulation

- `qtlIndex`:

  a list of vectors giving trait specific QTL indices relative to all
  active QTL

## Methods

### Public methods

- [`SimParam$new()`](#method-SimParam-new)

- [`SimParam$setTrackPed()`](#method-SimParam-setTrackPed)

- [`SimParam$setTrackRec()`](#method-SimParam-setTrackRec)

- [`SimParam$resetPed()`](#method-SimParam-resetPed)

- [`SimParam$restrSegSites()`](#method-SimParam-restrSegSites)

- [`SimParam$setSexes()`](#method-SimParam-setSexes)

- [`SimParam$setFounderHap()`](#method-SimParam-setFounderHap)

- [`SimParam$addSnpChip()`](#method-SimParam-addSnpChip)

- [`SimParam$addSnpChipByName()`](#method-SimParam-addSnpChipByName)

- [`SimParam$addStructuredSnpChip()`](#method-SimParam-addStructuredSnpChip)

- [`SimParam$addTraitA()`](#method-SimParam-addTraitA)

- [`SimParam$addTraitAD()`](#method-SimParam-addTraitAD)

- [`SimParam$altAddTraitAD()`](#method-SimParam-altAddTraitAD)

- [`SimParam$addTraitAG()`](#method-SimParam-addTraitAG)

- [`SimParam$addTraitADG()`](#method-SimParam-addTraitADG)

- [`SimParam$addTraitAE()`](#method-SimParam-addTraitAE)

- [`SimParam$addTraitADE()`](#method-SimParam-addTraitADE)

- [`SimParam$addTraitAEG()`](#method-SimParam-addTraitAEG)

- [`SimParam$addTraitADEG()`](#method-SimParam-addTraitADEG)

- [`SimParam$manAddTrait()`](#method-SimParam-manAddTrait)

- [`SimParam$importTrait()`](#method-SimParam-importTrait)

- [`SimParam$switchTrait()`](#method-SimParam-switchTrait)

- [`SimParam$removeTrait()`](#method-SimParam-removeTrait)

- [`SimParam$setVarE()`](#method-SimParam-setVarE)

- [`SimParam$setCorE()`](#method-SimParam-setCorE)

- [`SimParam$rescaleTraits()`](#method-SimParam-rescaleTraits)

- [`SimParam$setRecombRatio()`](#method-SimParam-setRecombRatio)

- [`SimParam$switchGenMap()`](#method-SimParam-switchGenMap)

- [`SimParam$switchFemaleMap()`](#method-SimParam-switchFemaleMap)

- [`SimParam$switchMaleMap()`](#method-SimParam-switchMaleMap)

- [`SimParam$addToRec()`](#method-SimParam-addToRec)

- [`SimParam$ibdHaplo()`](#method-SimParam-ibdHaplo)

- [`SimParam$updateLastId()`](#method-SimParam-updateLastId)

- [`SimParam$addToPed()`](#method-SimParam-addToPed)

- [`SimParam$clone()`](#method-SimParam-clone)

------------------------------------------------------------------------

### Method `new()`

Starts the process of building a new simulation by creating a new
SimParam object and assigning a founder population to the class. It is
recommended that you save the object with the name "SP", because
subsequent functions will check your global environment for an object of
this name if their simParam arguments are NULL. This allows you to call
these functions without explicitly supplying a simParam argument with
every call.

#### Usage

    SimParam$new(founderPop)

#### Arguments

- `founderPop`:

  an object of
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)

------------------------------------------------------------------------

### Method `setTrackPed()`

Sets pedigree tracking for the simulation. By default pedigree tracking
is turned off. When turned on, the pedigree of all individuals created
will be tracked, except those created by
[`hybridCross`](https://gaynorr.github.io/AlphaSimR/reference/hybridCross.md).
Turning off pedigree tracking will turn off recombination tracking if it
is turned on.

#### Usage

    SimParam$setTrackPed(isTrackPed, force = FALSE)

#### Arguments

- `isTrackPed`:

  should pedigree tracking be on.

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$setTrackPed(TRUE)

------------------------------------------------------------------------

### Method `setTrackRec()`

Sets recombination tracking for the simulation. By default recombination
tracking is turned off. When turned on recombination tracking will also
turn on pedigree tracking. Recombination tracking keeps records of all
individuals created, except those created by
[`hybridCross`](https://gaynorr.github.io/AlphaSimR/reference/hybridCross.md),
because their pedigree is not tracked.

#### Usage

    SimParam$setTrackRec(isTrackRec, force = FALSE)

#### Arguments

- `isTrackRec`:

  should recombination tracking be on.

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$setTrackRec(TRUE)

------------------------------------------------------------------------

### Method `resetPed()`

Resets the internal lastId, the pedigree and recombination tracking (if
in use) to the supplied lastId. Be careful using this function because
it may introduce a bug if you use individuals from the deleted portion
of the pedigree.

#### Usage

    SimParam$resetPed(lastId = 0L)

#### Arguments

- `lastId`:

  last ID to include in pedigree

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}

    #Create population
    pop = newPop(founderPop, simParam=SP)
    pop@id # 1:10

    #Create another population after reseting pedigree
    SP$resetPed()
    pop2 = newPop(founderPop, simParam=SP)
    pop2@id # 1:10

------------------------------------------------------------------------

### Method `restrSegSites()`

Sets restrictions on which segregating sites can serve as a SNP and/or
QTL. The default behavior of AlphaSimR is to randomly sample QTL or SNP
from all eligible sites and then mark the sampled sites ineligible to be
sampled as the other type (e.g. if a site is sampled as a QTL it will be
marked as ineligible to be sampled as a SNP). This behavior is designed
to produce the most challenging scenario for genomic selection when the
markers used for prediction are not causal.

Setting overlap=TRUE will prevent the addition of loci to the ineligible
list, but it won't remove sites already added to these lists. Thus,
timing of when restrSegSites is called matters. It should be called
before any addTrait or addSnpChip functions with the overlap=TRUE
argument to freely allow loci to overlap.

The minQtlPerChr and minSnpPerChr arguments can be used with
overlap=FALSE to preallocate sites as QTL and SNP respectively. This
option is useful when simulating multiple traits and/or SNP chips,
because it can be used to guarantee that enough eligible sites are
available when running addTrait and or addSnpChip functions.

#### Usage

    SimParam$restrSegSites(
      minQtlPerChr = NULL,
      minSnpPerChr = NULL,
      excludeQtl = NULL,
      excludeSnp = NULL,
      overlap = FALSE,
      minSnpFreq = NULL
    )

#### Arguments

- `minQtlPerChr`:

  the minimum number of segregating sites for QTLs. Can be a single
  value or a vector values for each chromosome.

- `minSnpPerChr`:

  the minimum number of segregating sites for SNPs. Can be a single
  value or a vector values for each chromosome.

- `excludeQtl`:

  an optional vector of segregating site names to exclude from
  consideration as a viable QTL.

- `excludeSnp`:

  an optional vector of segregating site names to exclude from
  consideration as a viable SNP.

- `overlap`:

  should SNP and QTL sites be allowed to overlap.

- `minSnpFreq`:

  minimum allowable frequency for SNP loci. No minimum SNP frequency is
  used if value is NULL.

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$restrSegSites(minQtlPerChr=5, minSnpPerChr=5)

------------------------------------------------------------------------

### Method `setSexes()`

Changes how sexes are determined in the simulation. The default sexes is
"no", indicating all individuals are hermaphrodites. To add sexes to the
simulation, run this function with "yes_sys" or "yes_rand". The value
"yes_sys" will systematically assign sexes to newly created individuals
as first male and then female. Populations with an odd number of
individuals will have one more male than female. The value "yes_rand"
will randomly assign a sex to each individual.

#### Usage

    SimParam$setSexes(sexes, force = FALSE)

#### Arguments

- `sexes`:

  acceptable value are "no", "yes_sys", or "yes_rand"

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$setSexes("yes_sys")

------------------------------------------------------------------------

### Method `setFounderHap()`

Allows for the manual setting of founder haplotypes. This functionality
is not fully documented, because it is still experimental.

#### Usage

    SimParam$setFounderHap(hapMap)

#### Arguments

- `hapMap`:

  a list of founder haplotypes

------------------------------------------------------------------------

### Method `addSnpChip()`

Randomly assigns eligible SNPs to a SNP chip

#### Usage

    SimParam$addSnpChip(nSnpPerChr, minSnpFreq = NULL, refPop = NULL, name = NULL)

#### Arguments

- `nSnpPerChr`:

  number of SNPs per chromosome. Can be a single value or nChr values.

- `minSnpFreq`:

  minimum allowable frequency for SNP loci. If NULL, no minimum
  frequency is used.

- `refPop`:

  reference population for calculating SNP frequency. If NULL, the
  founder population is used.

- `name`:

  optional name for chip

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addSnpChip(10)

------------------------------------------------------------------------

### Method `addSnpChipByName()`

Assigns SNPs to a SNP chip by supplying marker names. This function does
check against excluded SNPs and will not add the SNPs to the list of
excluded QTL for the purpose of avoiding overlap between SNPs and QTL.
Excluding these SNPs from being used as QTL can be accomplished using
the excludeQtl argument in SimParam's restrSegSites function.

#### Usage

    SimParam$addSnpChipByName(markers, name = NULL)

#### Arguments

- `markers`:

  a vector of names for the markers

- `name`:

  optional name for chip

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    SP$addSnpChipByName(c("1_1","1_3"))

------------------------------------------------------------------------

### Method `addStructuredSnpChip()`

Randomly selects the number of snps in structure and then assigns them
to chips based on structure

#### Usage

    SimParam$addStructuredSnpChip(nSnpPerChr, structure, force = FALSE)

#### Arguments

- `nSnpPerChr`:

  number of SNPs per chromosome. Can be a single value or nChr values.

- `structure`:

  a matrix. Rows are snp chips, columns are chips. If value is true then
  that snp is on that chip.

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

------------------------------------------------------------------------

### Method `addTraitA()`

Randomly assigns eligible QTLs for one or more additive traits. If
simulating more than one trait, all traits will be pleiotropic with
correlated additive effects.

#### Usage

    SimParam$addTraitA(
      nQtlPerChr,
      mean = 0,
      var = 1,
      corA = NULL,
      gamma = FALSE,
      shape = 1,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  a vector of desired mean genetic values for one or more traits

- `var`:

  a vector of desired genetic variances for one or more traits

- `corA`:

  a matrix of correlations between additive effects

- `gamma`:

  should a gamma distribution be used instead of normal

- `shape`:

  the shape parameter for the gamma distribution (the rate/scale
  parameter of the gamma distribution is accounted for via the desired
  level of genetic variance, the var argument)

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait(s)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitA(10)

------------------------------------------------------------------------

### Method `addTraitAD()`

Randomly assigns eligible QTLs for one or more traits with dominance. If
simulating more than one trait, all traits will be pleiotropic with
correlated effects.

#### Usage

    SimParam$addTraitAD(
      nQtlPerChr,
      mean = 0,
      var = 1,
      meanDD = 0,
      varDD = 0,
      corA = NULL,
      corDD = NULL,
      useVarA = TRUE,
      gamma = FALSE,
      shape = 1,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  a vector of desired mean genetic values for one or more traits

- `var`:

  a vector of desired genetic variances for one or more traits

- `meanDD`:

  mean dominance degree

- `varDD`:

  variance of dominance degree

- `corA`:

  a matrix of correlations between additive effects

- `corDD`:

  a matrix of correlations between dominance degrees

- `useVarA`:

  tune according to additive genetic variance if true. If FALSE, tuning
  is performed according to total genetic variance.

- `gamma`:

  should a gamma distribution be used instead of normal

- `shape`:

  the shape parameter for the gamma distribution (the rate/scale
  parameter of the gamma distribution is accounted for via the desired
  level of genetic variance, the var argument)

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait(s)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitAD(10, meanDD=0.5)

------------------------------------------------------------------------

### Method `altAddTraitAD()`

An alternative method for adding a trait with additive and dominance
effects to an AlphaSimR simulation. The function attempts to create a
trait matching user defined values for number of QTL, inbreeding
depression, additive genetic variance and dominance genetic variance.

#### Usage

    SimParam$altAddTraitAD(
      nQtlPerChr,
      mean = 0,
      varA = 1,
      varD = 0,
      inbrDepr = 0,
      limMeanDD = c(0, 1.5),
      limVarDD = c(0, 0.5),
      silent = FALSE,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  desired mean of the trait

- `varA`:

  desired additive variance

- `varD`:

  desired dominance variance

- `inbrDepr`:

  desired inbreeding depression, see details

- `limMeanDD`:

  limits for meanDD, see details

- `limVarDD`:

  limits for varDD, see details

- `silent`:

  should summary details be printed to the console

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait

#### Details

This function will always add a trait to 'SimParam', unless an error
occurs with picking QTLs. The resulting trait will always have the
desired mean and additive genetic variance. However, it may not have the
desired values for inbreeding depression and dominance variance. Thus,
it is strongly recommended to check the output printed to the console to
determine how close the trait's parameters came to these desired values.

The mean and additive genetic variance will always be achieved exactly.
The function attempts to achieve the desired dominance variance and
inbreeding depression while staying within the user supplied constraints
for the acceptable range of dominance degree mean and variance. If the
desired values are not being achieved, the acceptable range need to be
increased and/or the number of QTL may need to be increased. There are
not limits to setting the range for dominance degree mean and variance,
but care should be taken to with regards to the biological feasibility
of the limits that are supplied. The default limits were somewhat
arbitrarily set, so I make not claim to how reasonable these limits are
for routine use.

Inbreeding depression in this function is defined as the difference in
mean genetic value between a population with the same allele frequency
as the reference population (population used to initialize SimParam) in
Hardy-Weinberg equilibrium compared to a population with the same allele
frequency that is fully inbred. This is equivalent to the amount the
mean of a population increases when going from an inbreeding coefficient
of 1 (fully inbred) to a population with an inbreeding coefficient of 0
(Hardy-Weinberg equilibrium). Note that the sign of the value should
(usually) be positive. This corresponds to a detrimental effect of
inbreeding when higher values of the trait are considered biologically
beneficial.

Summary information on this trait is printed to the console when
silent=FALSE. The summary information reports the inbreeding depression
and dominance variance for the population as well as the dominance
degree mean and variance applied to the trait.

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$altAddTraitAD(nQtlPerChr=10, mean=0, varA=1, varD=0.05, inbrDepr=0.2)

------------------------------------------------------------------------

### Method `addTraitAG()`

Randomly assigns eligible QTLs for one or more additive GxE traits. If
simulating more than one trait, all traits will be pleiotropic with
correlated effects.

#### Usage

    SimParam$addTraitAG(
      nQtlPerChr,
      mean = 0,
      var = 1,
      varGxE = 1e-06,
      varEnv = 0,
      corA = NULL,
      corGxE = NULL,
      gamma = FALSE,
      shape = 1,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  a vector of desired mean genetic values for one or more traits

- `var`:

  a vector of desired genetic variances for one or more traits

- `varGxE`:

  a vector of total genotype-by-environment variances for the traits

- `varEnv`:

  a vector of environmental variances for one or more traits

- `corA`:

  a matrix of correlations between additive effects

- `corGxE`:

  a matrix of correlations between GxE effects

- `gamma`:

  should a gamma distribution be used instead of normal

- `shape`:

  the shape parameter for the gamma distribution (the rate/scale
  parameter of the gamma distribution is accounted for via the desired
  level of genetic variance, the var argument)

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait(s)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitAG(10, varGxE=2)

------------------------------------------------------------------------

### Method `addTraitADG()`

Randomly assigns eligible QTLs for a trait with dominance and GxE.

#### Usage

    SimParam$addTraitADG(
      nQtlPerChr,
      mean = 0,
      var = 1,
      varEnv = 0,
      varGxE = 1e-06,
      meanDD = 0,
      varDD = 0,
      corA = NULL,
      corDD = NULL,
      corGxE = NULL,
      useVarA = TRUE,
      gamma = FALSE,
      shape = 1,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  a vector of desired mean genetic values for one or more traits

- `var`:

  a vector of desired genetic variances for one or more traits

- `varEnv`:

  a vector of environmental variances for one or more traits

- `varGxE`:

  a vector of total genotype-by-environment variances for the traits

- `meanDD`:

  mean dominance degree

- `varDD`:

  variance of dominance degree

- `corA`:

  a matrix of correlations between additive effects

- `corDD`:

  a matrix of correlations between dominance degrees

- `corGxE`:

  a matrix of correlations between GxE effects

- `useVarA`:

  tune according to additive genetic variance if true

- `gamma`:

  should a gamma distribution be used instead of normal

- `shape`:

  the shape parameter for the gamma distribution (the rate/scale
  parameter of the gamma distribution is accounted for via the desired
  level of genetic variance, the var argument)

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait(s)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitADG(10, meanDD=0.5, varGxE=2)

------------------------------------------------------------------------

### Method `addTraitAE()`

Randomly assigns eligible QTLs for one or more additive and epistasis
traits. If simulating more than one trait, all traits will be
pleiotropic with correlated additive effects.

#### Usage

    SimParam$addTraitAE(
      nQtlPerChr,
      mean = 0,
      var = 1,
      relAA = 0,
      corA = NULL,
      corAA = NULL,
      useVarA = TRUE,
      gamma = FALSE,
      shape = 1,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  a vector of desired mean genetic values for one or more traits

- `var`:

  a vector of desired genetic variances for one or more traits

- `relAA`:

  the relative value of additive-by-additive variance compared to
  additive variance in a diploid organism with allele frequency 0.5

- `corA`:

  a matrix of correlations between additive effects

- `corAA`:

  a matrix of correlations between additive-by-additive effects

- `useVarA`:

  tune according to additive genetic variance if true. If FALSE, tuning
  is performed according to total genetic variance.

- `gamma`:

  should a gamma distribution be used instead of normal

- `shape`:

  the shape parameter for the gamma distribution (the rate/scale
  parameter of the gamma distribution is accounted for via the desired
  level of genetic variance, the var argument)

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait(s)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitAE(10, relAA=0.1)

------------------------------------------------------------------------

### Method `addTraitADE()`

Randomly assigns eligible QTLs for one or more traits with dominance and
epistasis. If simulating more than one trait, all traits will be
pleiotropic with correlated effects.

#### Usage

    SimParam$addTraitADE(
      nQtlPerChr,
      mean = 0,
      var = 1,
      meanDD = 0,
      varDD = 0,
      relAA = 0,
      corA = NULL,
      corDD = NULL,
      corAA = NULL,
      useVarA = TRUE,
      gamma = FALSE,
      shape = 1,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  a vector of desired mean genetic values for one or more traits

- `var`:

  a vector of desired genetic variances for one or more traits

- `meanDD`:

  mean dominance degree

- `varDD`:

  variance of dominance degree

- `relAA`:

  the relative value of additive-by-additive variance compared to
  additive variance in a diploid organism with allele frequency 0.5

- `corA`:

  a matrix of correlations between additive effects

- `corDD`:

  a matrix of correlations between dominance degrees

- `corAA`:

  a matrix of correlations between additive-by-additive effects

- `useVarA`:

  tune according to additive genetic variance if true. If FALSE, tuning
  is performed according to total genetic variance.

- `gamma`:

  should a gamma distribution be used instead of normal

- `shape`:

  the shape parameter for the gamma distribution (the rate/scale
  parameter of the gamma distribution is accounted for via the desired
  level of genetic variance, the var argument)

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait(s)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitADE(10)

------------------------------------------------------------------------

### Method `addTraitAEG()`

Randomly assigns eligible QTLs for one or more additive and epistasis
GxE traits. If simulating more than one trait, all traits will be
pleiotropic with correlated effects.

#### Usage

    SimParam$addTraitAEG(
      nQtlPerChr,
      mean = 0,
      var = 1,
      relAA = 0,
      varGxE = 1e-06,
      varEnv = 0,
      corA = NULL,
      corAA = NULL,
      corGxE = NULL,
      useVarA = TRUE,
      gamma = FALSE,
      shape = 1,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  a vector of desired mean genetic values for one or more traits

- `var`:

  a vector of desired genetic variances for one or more traits

- `relAA`:

  the relative value of additive-by-additive variance compared to
  additive variance in a diploid organism with allele frequency 0.5

- `varGxE`:

  a vector of total genotype-by-environment variances for the traits

- `varEnv`:

  a vector of environmental variances for one or more traits

- `corA`:

  a matrix of correlations between additive effects

- `corAA`:

  a matrix of correlations between additive-by-additive effects

- `corGxE`:

  a matrix of correlations between GxE effects

- `useVarA`:

  tune according to additive genetic variance if true. If FALSE, tuning
  is performed according to total genetic variance.

- `gamma`:

  should a gamma distribution be used instead of normal

- `shape`:

  the shape parameter for the gamma distribution (the rate/scale
  parameter of the gamma distribution is accounted for via the desired
  level of genetic variance, the var argument)

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait(s)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitAEG(10, varGxE=2)

------------------------------------------------------------------------

### Method `addTraitADEG()`

Randomly assigns eligible QTLs for a trait with dominance, epistasis and
GxE.

#### Usage

    SimParam$addTraitADEG(
      nQtlPerChr,
      mean = 0,
      var = 1,
      varEnv = 0,
      varGxE = 1e-06,
      meanDD = 0,
      varDD = 0,
      relAA = 0,
      corA = NULL,
      corDD = NULL,
      corAA = NULL,
      corGxE = NULL,
      useVarA = TRUE,
      gamma = FALSE,
      shape = 1,
      force = FALSE,
      name = NULL
    )

#### Arguments

- `nQtlPerChr`:

  number of QTLs per chromosome. Can be a single value or nChr values.

- `mean`:

  a vector of desired mean genetic values for one or more traits

- `var`:

  a vector of desired genetic variances for one or more traits

- `varEnv`:

  a vector of environmental variances for one or more traits

- `varGxE`:

  a vector of total genotype-by-environment variances for the traits

- `meanDD`:

  mean dominance degree

- `varDD`:

  variance of dominance degree

- `relAA`:

  the relative value of additive-by-additive variance compared to
  additive variance in a diploid organism with allele frequency 0.5

- `corA`:

  a matrix of correlations between additive effects

- `corDD`:

  a matrix of correlations between dominance degrees

- `corAA`:

  a matrix of correlations between additive-by-additive effects

- `corGxE`:

  a matrix of correlations between GxE effects

- `useVarA`:

  tune according to additive genetic variance if true

- `gamma`:

  should a gamma distribution be used instead of normal

- `shape`:

  the shape parameter for the gamma distribution (the rate/scale
  parameter of the gamma distribution is accounted for via the desired
  level of genetic variance, the var argument)

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing.

- `name`:

  optional name for trait(s)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitADEG(10, meanDD=0.5, varGxE=2)

------------------------------------------------------------------------

### Method `manAddTrait()`

Manually add a new trait to the simulation. Trait must be formatted as a
[`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md).
If the trait is not already formatted, consider using importTrait.

#### Usage

    SimParam$manAddTrait(lociMap, varE = NA_real_, force = FALSE)

#### Arguments

- `lociMap`:

  a new object descended from
  [`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)

- `varE`:

  default error variance for phenotype, optional

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing

------------------------------------------------------------------------

### Method `importTrait()`

Manually add a new trait(s) to the simulation. Unlike the manAddTrait
function, this function does not require formatting the trait as a
[`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md).
The formatting is performed automatically for the user, with more user
friendly data.frames or matrices taken as inputs. This function only
works for A and AD trait types.

#### Usage

    SimParam$importTrait(
      markerNames,
      addEff,
      domEff = NULL,
      intercept = NULL,
      name = NULL,
      varE = NULL,
      force = FALSE
    )

#### Arguments

- `markerNames`:

  a vector of names for the QTL

- `addEff`:

  a matrix of additive effects (nLoci x nTraits). Alternatively, a
  vector of length nLoci can be supplied for a single trait.

- `domEff`:

  optional dominance effects for each locus

- `intercept`:

  optional intercepts for each trait

- `name`:

  optional name(s) for the trait(s)

- `varE`:

  default error variance for phenotype, optional

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing

------------------------------------------------------------------------

### Method `switchTrait()`

Switch a trait in the simulation.

#### Usage

    SimParam$switchTrait(traitPos, lociMap, varE = NA_real_, force = FALSE)

#### Arguments

- `traitPos`:

  an integer indicate which trait to switch

- `lociMap`:

  a new object descended from
  [`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)

- `varE`:

  default error variance for phenotype, optional

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing

------------------------------------------------------------------------

### Method `removeTrait()`

Remove a trait from the simulation

#### Usage

    SimParam$removeTrait(traits, force = FALSE)

#### Arguments

- `traits`:

  an integer vector indicating which traits to remove

- `force`:

  should the check for a running simulation be ignored. Only set to TRUE
  if you know what you are doing

------------------------------------------------------------------------

### Method `setVarE()`

Defines a default values for error variances used in
[`setPheno`](https://gaynorr.github.io/AlphaSimR/reference/setPheno.md).
These defaults will be used to automatically generate phenotypes when
new populations are created. See the details section of
[`setPheno`](https://gaynorr.github.io/AlphaSimR/reference/setPheno.md)
for more information about each arguments and how they should be used.

#### Usage

    SimParam$setVarE(h2 = NULL, H2 = NULL, varE = NULL, corE = NULL)

#### Arguments

- `h2`:

  a vector of desired narrow-sense heritabilities

- `H2`:

  a vector of desired broad-sense heritabilities

- `varE`:

  a vector or matrix of error variances

- `corE`:

  an optional matrix of error correlations

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitA(10)
    SP$setVarE(h2=0.5)

------------------------------------------------------------------------

### Method `setCorE()`

Defines a correlation structure for default error variances. You must
call `setVarE` first to define the default error variances.

#### Usage

    SimParam$setCorE(corE)

#### Arguments

- `corE`:

  a correlation matrix for the error variances

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$addTraitA(10, mean=c(0,0), var=c(1,1), corA=diag(2))
    SP$setVarE(varE=c(1,1))
    E = 0.5*diag(2)+0.5 #Positively correlated error
    SP$setCorE(E)

------------------------------------------------------------------------

### Method `rescaleTraits()`

Linearly scales all traits to achieve desired values of means and
variances in the founder population.

#### Usage

    SimParam$rescaleTraits(
      mean = 0,
      var = 1,
      varEnv = 0,
      varGxE = 1e-06,
      useVarA = TRUE
    )

#### Arguments

- `mean`:

  a vector of new trait means

- `var`:

  a vector of new trait variances

- `varEnv`:

  a vector of new environmental variances

- `varGxE`:

  a vector of new GxE variances

- `useVarA`:

  tune according to additive genetic variance if true

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    SP$addTraitA(10)

    #Create population
    pop = newPop(founderPop, simParam=SP)
    meanG(pop)

    #Change mean to 1
    SP$rescaleTraits(mean=1)
    \dontshow{SP$nThreads = 1L}
    #Run resetPop for change to take effect
    pop = resetPop(pop, simParam=SP)
    meanG(pop)

------------------------------------------------------------------------

### Method `setRecombRatio()`

Set the relative recombination rates between males and females. This
allows for sex-specific recombination rates, under the assumption of
equivalent recombination landscapes.

#### Usage

    SimParam$setRecombRatio(femaleRatio)

#### Arguments

- `femaleRatio`:

  relative ratio of recombination in females compared to males. A value
  of 2 indicate twice as much recombination in females. The value must
  be greater than 0. (default is 1)

#### Examples

    #Create founder haplotypes
    founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

    #Set simulation parameters
    SP = SimParam$new(founderPop)
    \dontshow{SP$nThreads = 1L}
    SP$setRecombRatio(2) #Twice as much recombination in females

------------------------------------------------------------------------

### Method `switchGenMap()`

Replaces existing genetic map.

#### Usage

    SimParam$switchGenMap(genMap, centromere = NULL)

#### Arguments

- `genMap`:

  a list of length nChr containing numeric vectors for the position of
  each segregating site on a chromosome.

- `centromere`:

  a numeric vector of centromere positions. If NULL, the centromere are
  assumed to be metacentric.

------------------------------------------------------------------------

### Method `switchFemaleMap()`

Replaces existing female genetic map.

#### Usage

    SimParam$switchFemaleMap(genMap, centromere = NULL)

#### Arguments

- `genMap`:

  a list of length nChr containing numeric vectors for the position of
  each segregating site on a chromosome.

- `centromere`:

  a numeric vector of centromere positions. If NULL, the centromere are
  assumed to be metacentric.

------------------------------------------------------------------------

### Method `switchMaleMap()`

Replaces existing male genetic map.

#### Usage

    SimParam$switchMaleMap(genMap, centromere = NULL)

#### Arguments

- `genMap`:

  a list of length nChr containing numeric vectors for the position of
  each segregating site on a chromosome.

- `centromere`:

  a numeric vector of centromere positions. If NULL, the centromere are
  assumed to be metacentric.

------------------------------------------------------------------------

### Method `addToRec()`

For internal use only.

#### Usage

    SimParam$addToRec(lastId, id, mother, father, isDH, hist, ploidy)

#### Arguments

- `lastId`:

  ID of last individual

- `id`:

  the name of each individual

- `mother`:

  vector of mother iids

- `father`:

  vector of father iids

- `isDH`:

  indicator for DH lines

- `hist`:

  new recombination history

- `ploidy`:

  ploidy level

------------------------------------------------------------------------

### Method `ibdHaplo()`

For internal use only.

#### Usage

    SimParam$ibdHaplo(iid)

#### Arguments

- `iid`:

  internal ID

------------------------------------------------------------------------

### Method `updateLastId()`

For internal use only.

#### Usage

    SimParam$updateLastId(lastId)

#### Arguments

- `lastId`:

  last ID assigned

------------------------------------------------------------------------

### Method `addToPed()`

For internal use only.

#### Usage

    SimParam$addToPed(lastId, id, mother, father, isDH)

#### Arguments

- `lastId`:

  ID of last individual

- `id`:

  the name of each individual

- `mother`:

  vector of mother iids

- `father`:

  vector of father iids

- `isDH`:

  indicator for DH lines

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimParam$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
## ------------------------------------------------
## Method `SimParam$new`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

## ------------------------------------------------
## Method `SimParam$setTrackPed`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$setTrackPed(TRUE)

## ------------------------------------------------
## Method `SimParam$setTrackRec`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$setTrackRec(TRUE)

## ------------------------------------------------
## Method `SimParam$resetPed`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create population
pop = newPop(founderPop, simParam=SP)
pop@id # 1:10
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"

#Create another population after reseting pedigree
SP$resetPed()
pop2 = newPop(founderPop, simParam=SP)
pop2@id # 1:10
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"

## ------------------------------------------------
## Method `SimParam$restrSegSites`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$restrSegSites(minQtlPerChr=5, minSnpPerChr=5)

## ------------------------------------------------
## Method `SimParam$setSexes`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$setSexes("yes_sys")

## ------------------------------------------------
## Method `SimParam$addSnpChip`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addSnpChip(10)

## ------------------------------------------------
## Method `SimParam$addSnpChipByName`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addSnpChipByName(c("1_1","1_3"))

## ------------------------------------------------
## Method `SimParam$addTraitA`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

## ------------------------------------------------
## Method `SimParam$addTraitAD`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)

## ------------------------------------------------
## Method `SimParam$altAddTraitAD`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$altAddTraitAD(nQtlPerChr=10, mean=0, varA=1, varD=0.05, inbrDepr=0.2)
#> A new trait called Trait1 was added. 
#>    varD = 0.05000002 
#>    inbrDepr = 0.2000004 
#>    meanDD = 0.04665813 
#>    varDD = 0.1382312 

## ------------------------------------------------
## Method `SimParam$addTraitAG`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAG(10, varGxE=2)

## ------------------------------------------------
## Method `SimParam$addTraitADG`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitADG(10, meanDD=0.5, varGxE=2)

## ------------------------------------------------
## Method `SimParam$addTraitAE`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAE(10, relAA=0.1)

## ------------------------------------------------
## Method `SimParam$addTraitADE`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitADE(10)

## ------------------------------------------------
## Method `SimParam$addTraitAEG`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAEG(10, varGxE=2)

## ------------------------------------------------
## Method `SimParam$addTraitADEG`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitADEG(10, meanDD=0.5, varGxE=2)

## ------------------------------------------------
## Method `SimParam$setVarE`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)

## ------------------------------------------------
## Method `SimParam$setCorE`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10, mean=c(0,0), var=c(1,1), corA=diag(2))
SP$setVarE(varE=c(1,1))
E = 0.5*diag(2)+0.5 #Positively correlated error
SP$setCorE(E)
#> Warning: This function has been deprecated. Use simParam$setVarE instead.

## ------------------------------------------------
## Method `SimParam$rescaleTraits`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)
meanG(pop)
#>        Trait1 
#> -6.661338e-17 

#Change mean to 1
SP$rescaleTraits(mean=1)
#Run resetPop for change to take effect
pop = resetPop(pop, simParam=SP)
meanG(pop)
#> Trait1 
#>      1 

## ------------------------------------------------
## Method `SimParam$setRecombRatio`
## ------------------------------------------------

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$setRecombRatio(2) #Twice as much recombination in females
```
