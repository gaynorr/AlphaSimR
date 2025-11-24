# Create founder haplotypes using MaCS

Uses the MaCS software to produce founder haplotypes (Chen et al. 2009)
.

## Usage

``` r
runMacs(
  nInd,
  nChr = 1,
  segSites = NULL,
  inbred = FALSE,
  species = "GENERIC",
  split = NULL,
  ploidy = 2L,
  manualCommand = NULL,
  manualGenLen = NULL,
  nThreads = NULL
)
```

## Arguments

- nInd:

  number of individuals to simulate

- nChr:

  number of chromosomes to simulate

- segSites:

  number of segregating sites to keep per chromosome. A value of NULL
  results in all sites being retained.

- inbred:

  should founder individuals be inbred

- species:

  species history to simulate. See details.

- split:

  an optional historic population split in terms of generations ago.

- ploidy:

  ploidy level of organism

- manualCommand:

  user provided MaCS options. For advanced users only.

- manualGenLen:

  user provided genetic length. This must be supplied if using
  manualCommand. If not using manualCommand, this value will replace the
  predefined genetic length for the species. However, this the genetic
  length is only used by AlphaSimR and is not passed to MaCS, so MaCS
  still uses the predefined genetic length. For advanced users only.

- nThreads:

  if OpenMP is available, this will allow for simulating chromosomes in
  parallel. If the value is NULL, the number of threads is automatically
  detected.

## Value

an object of
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

## Details

There are currently three species histories available: GENERIC, CATTLE,
WHEAT, and MAIZE.

The GENERIC history is meant to be a reasonable all-purpose choice. It
runs quickly and models a population with an effective populations size
that has gone through several historic bottlenecks. This species history
is used as the default arguments in the
[`runMacs2`](https://gaynorr.github.io/AlphaSimR/reference/runMacs2.md)
function, so the user should examine this function for the details of
how the species is modeled.

The CATTLE history is based off of real genome sequence data (MacLeod et
al. 2013) .

The WHEAT (Gaynor et al. 2017) and MAIZE (Hickey et al. 2014) histories
have been included due to their use in previous simulations. However, it
should be noted that neither faithfully simulates its respective
species. This is apparent by the low number of segregating sites
simulated by each history relative to their real-world analogs.
Adjusting these histories to better represent their real-world analogs
would result in a drastic increase to runtime.

## References

Chen GK, Marjoram P, Wall JD (2009). “Fast and Flexible Simulation of
DNA Sequence Data.” *Genome Research*, **19**, 136-142.
<https://genome.cshlp.org/content/19/1/136>.  
  
Gaynor RC, Gorjanc G, Bentley AR, Ober ES, Howell P, Jackson R, Mackay
IJ, Hickey JM (2017). “A Two-Part Strategy for Using Genomic Selection
to Develop Inbred Lines.” *Crop Science*, **57**(5), 2372–2386. ISSN
0011-183X,
[doi:10.2135/cropsci2016.09.0742](https://doi.org/10.2135/cropsci2016.09.0742)
,
<https://acsess.onlinelibrary.wiley.com/doi/full/10.2135/cropsci2016.09.0742>.  
  
Hickey JMDS, Crossa J, Hearne S, Babu R, Prasanna BM, Grondona M,
Zambelli A, Windhausen VS, Mathews K, Gorjanc G (2014). “Evaluation of
Genomic Selection Training Population Designs and Genotyping Strategies
in Plant Breeding Programs Using Simulation.” *Crop Science*, **54**(4),
1476-1488.
[doi:10.2135/cropsci2013.03.0195](https://doi.org/10.2135/cropsci2013.03.0195)
.  
  
MacLeod IM, Larkin DM, Lewin HAHBJ, Goddard ME (2013). “Inferring
Demography from Runs of Homozygosity in Whole-Genome Sequence, with
Correction for Sequence Errors.” *Molecular Biology and Evolution*,
**30**(9), 2209–2223.
[doi:10.1093/molbev/mst125](https://doi.org/10.1093/molbev/mst125) .

## Examples

``` r
# Creates a populations of 10 outbred individuals
# Their genome consists of 1 chromosome and 100 segregating sites
if (FALSE) { # \dontrun{
founderPop = runMacs(nInd=10,nChr=1,segSites=100)
} # }
```
