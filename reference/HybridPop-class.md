# Hybrid population

A lightweight version of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
for hybrid lines. Memory is saved by not storing genotypic data.

## Usage

``` r
# S4 method for class 'HybridPop'
x[i]

# S4 method for class 'HybridPop'
c(x, ...)

isHybridPop(x)
```

## Arguments

- x:

  a 'HybridPop'

- i:

  index of individuals

- ...:

  additional 'HybridPop' objects

## Methods (by generic)

- `[`: Extract HybridPop using index or id

- `c(HybridPop)`: Combine multiple HybridPops

## Functions

- `isHybridPop()`: Test if object is of a HybridPop class

## Slots

- `nInd`:

  number of individuals

- `id`:

  an individual's identifier

- `mother`:

  the identifier of the individual's mother

- `father`:

  the identifier of the individual's father

- `nTraits`:

  number of traits

- `gv`:

  matrix of genetic values. When using GxE traits, gv reflects gv when
  p=0.5. Dimensions are nInd by nTraits.

- `pheno`:

  matrix of phenotypic values. Dimensions are nInd by nTraits.

- `gxe`:

  list containing GxE slopes for GxE traits
