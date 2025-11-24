# Alternative wrapper for MaCS

A wrapper function for
[`runMacs`](https://gaynorr.github.io/AlphaSimR/reference/runMacs.md).
This wrapper is designed to provide a more intuitive interface for
writing custom commands in MaCS (Chen et al. 2009) . It effectively
automates the creation of an appropriate line for the manualCommand
argument in
[`runMacs`](https://gaynorr.github.io/AlphaSimR/reference/runMacs.md)
using user supplied variables, but only allows for a subset of the
functionality offered by this argument. The default arguments of this
function were chosen to match species="GENERIC" in
[`runMacs`](https://gaynorr.github.io/AlphaSimR/reference/runMacs.md).

## Usage

``` r
runMacs2(
  nInd,
  nChr = 1,
  segSites = NULL,
  Ne = 100,
  bp = 1e+08,
  genLen = 1,
  mutRate = 2.5e-08,
  histNe = c(500, 1500, 6000, 12000, 1e+05),
  histGen = c(100, 1000, 10000, 1e+05, 1e+06),
  inbred = FALSE,
  split = NULL,
  ploidy = 2L,
  returnCommand = FALSE,
  nThreads = NULL
)
```

## Arguments

- nInd:

  number of individuals to simulate

- nChr:

  number of chromosomes to simulate

- segSites:

  number of segregating sites to keep per chromosome

- Ne:

  effective population size

- bp:

  base pair length of chromosome

- genLen:

  genetic length of chromosome in Morgans

- mutRate:

  per base pair mutation rate

- histNe:

  effective population size in previous generations

- histGen:

  number of generations ago for effective population sizes given in
  histNe

- inbred:

  should founder individuals be inbred

- split:

  an optional historic population split in terms of generations ago

- ploidy:

  ploidy level of organism

- returnCommand:

  should the command passed to manualCommand in
  [`runMacs`](https://gaynorr.github.io/AlphaSimR/reference/runMacs.md)
  be returned. If TRUE, MaCS will not be called and the command is
  returned instead.

- nThreads:

  if OpenMP is available, this will allow for simulating chromosomes in
  parallel. If the value is NULL, the number of threads is automatically
  detected.

## Value

an object of
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
or if returnCommand is true a string giving the MaCS command passed to
the manualCommand argument of
[`runMacs`](https://gaynorr.github.io/AlphaSimR/reference/runMacs.md).

## References

Chen GK, Marjoram P, Wall JD (2009). “Fast and Flexible Simulation of
DNA Sequence Data.” *Genome Research*, **19**, 136-142.
<https://genome.cshlp.org/content/19/1/136>.

## Examples

``` r
# Creates a populations of 10 outbred individuals
# Their genome consists of 1 chromosome and 100 segregating sites
# The command is equivalent to using species="GENERIC" in runMacs
if (FALSE) { # \dontrun{
founderPop = runMacs2(nInd=10,nChr=1,segSites=100)

# runMacs() Implementation of the cattle demography following
#  Macleod et al. (2013) https://doi.org/10.1093/molbev/mst125
cattleChrSum = 2.8e9 # https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/
(cattleChrBp = cattleChrSum / 30)
recRate = 9.26e-09
(cattleGenLen = recRate * cattleChrBp)
mutRate = 1.20e-08
runMacs2(nInd = 10, nChr = 1, Ne = 90, bp = cattleChrBp,
         genLen = cattleGenLen, mutRate = 1.20e-08,
         histNe  = c(120, 250, 350, 1000, 1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000),
         histGen = c(  3,   6,  12,   18,   24,  154,  454,  654, 1754,  2354,  3354, 33154),
         returnCommand = TRUE)
} # }
```
