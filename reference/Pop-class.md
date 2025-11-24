# Population

Extends
[`RawPop-class`](https://gaynorr.github.io/AlphaSimR/reference/RawPop-class.md)
to add sex, genetic values, phenotypes, and pedigrees.

## Usage

``` r
# S4 method for class 'Pop'
x[i]

# S4 method for class 'Pop'
c(x, ...)

# S4 method for class 'Pop'
show(object)

# S4 method for class 'Pop'
length(x)
```

## Arguments

- x:

  a 'Pop' object

- i:

  index of individuals

- ...:

  additional 'Pop' objects

- object:

  a 'Pop' object

## Methods (by generic)

- `[`: Extract Pop by index or id

- `c(Pop)`: Combine multiple Pops

- `show(Pop)`: Show population summary

- `length(Pop)`: Number of individuals in Pop (the same as nInd())

## Slots

- `id`:

  an individual's identifier

- `iid`:

  an individual's internal identifier

- `mother`:

  the identifier of the individual's mother

- `father`:

  the identifier of the individual's father

- `sex`:

  sex of individuals: "M" for males, "F" for females, and "H" for
  hermaphrodites

- `nTraits`:

  number of traits

- `gv`:

  matrix of genetic values. When using GxE traits, gv reflects gv when
  p=0.5. Dimensions are nInd by nTraits.

- `pheno`:

  matrix of phenotypic values. Dimensions are nInd by nTraits.

- `ebv`:

  matrix of estimated breeding values. Dimensions are nInd rows and a
  variable number of columns.

- `gxe`:

  list containing GxE slopes for GxE traits

- `fixEff`:

  a fixed effect relating to the phenotype. Used by genomic selection
  models but otherwise ignored.

- `misc`:

  a list whose elements correspond to additional miscellaneous nodes
  with the items for individuals in the population (see example in
  [`newPop`](https://gaynorr.github.io/AlphaSimR/reference/newPop.md)) -
  we support vectors and matrices or objects that have a generic length
  and subset method. This list is normally empty and exists solely as an
  open slot available for uses to store extra information about
  individuals.

- `miscPop`:

  a list of any length containing optional meta data for the population
  (see example in
  [`newPop`](https://gaynorr.github.io/AlphaSimR/reference/newPop.md)).
  This list is empty unless information is supplied by the user. Note
  that the list is emptied every time the population is subsetted or
  combined because the meta data for old population might not be valid
  anymore.

## See also

[`newPop`](https://gaynorr.github.io/AlphaSimR/reference/newPop.md),
[`newEmptyPop`](https://gaynorr.github.io/AlphaSimR/reference/newEmptyPop.md),
[`resetPop`](https://gaynorr.github.io/AlphaSimR/reference/resetPop.md)
