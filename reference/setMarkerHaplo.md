# Set marker haplotypes

Manually sets the haplotypes in a population for all individuals at one
or more loci.

## Usage

``` r
setMarkerHaplo(pop, haplo, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`RawPop-class`](https://gaynorr.github.io/AlphaSimR/reference/RawPop-class.md)
  or
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

- haplo:

  a matrix of haplotypes, see details

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md),
  not used if pop is
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

## Value

an object of the same class as the "pop" input

## Details

The format of the haplotype matrix should match the format of the output
from
[`pullMarkerHaplo`](https://gaynorr.github.io/AlphaSimR/reference/pullMarkerHaplo.md)
with the option haplo="all". Thus, it is recommended that this function
is first used to extract the haplotypes and that any desired changes be
made to the output of pullMarkerHaplo before passing the matrix to
setMarkerHaplo. Any changes made to QTL may potentially result in
changes to an individuals genetic value. These changes will be reflected
in the gv and/or gxe slot. All other slots will remain unchanged, so the
ebv and pheno slots will not reflect the new genotypes.

## Examples

``` r
# Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)

# Extract haplotypes for marker "1_1"
H = pullMarkerHaplo(founderPop, markers="1_1")

# Set the first haplotype to 1
H[1,1] = 1L

# Set marker haplotypes
founderPop = setMarkerHaplo(founderPop, haplo=H)
```
