# Finds positions of loci by marker name

Used to generate lociPerChr and lociLoc objects for a set of markers.
These objects can be passed other functions for pulling genotypes or
haplotypes.

## Usage

``` r
mapLoci(markers, genMap)
```

## Arguments

- markers:

  a vector of marker names

- genMap:

  a genetic map in AlphaSimR's internal genetic map format

## Value

A list containing lociPerChr and lociLoc that can be
