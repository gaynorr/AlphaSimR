# Raw Population

The raw population class contains only genotype data.

## Usage

``` r
# S4 method for class 'RawPop'
x[i]

# S4 method for class 'RawPop'
c(x, ...)

# S4 method for class 'RawPop'
show(object)

isRawPop(x)
```

## Arguments

- x:

  a 'RawPop' object

- i:

  index of individuals

- ...:

  additional 'RawPop' objects

- object:

  a 'RawPop' object

## Methods (by generic)

- `[`: Extract RawPop by index

- `c(RawPop)`: Combine multiple RawPops

- `show(RawPop)`: Show population summary

## Functions

- `isRawPop()`: Test if object is of a RawPop class

## Slots

- `nInd`:

  number of individuals

- `nChr`:

  number of chromosomes

- `ploidy`:

  level of ploidy

- `nLoci`:

  number of loci per chromosome

- `geno`:

  list of nChr length containing chromosome genotypes. Each element is a
  three dimensional array of raw values. The array dimensions are nLoci
  by ploidy by nInd.
