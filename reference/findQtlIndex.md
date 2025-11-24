# Find trait QTL index

Compares to a
[`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)
objects to determine if the first one is a superset of the second. If it
is, the function returns NULL. If it is not, the function return a
[`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)
object that is a superset of both a
[`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)
objects.

## Usage

``` r
findQtlIndex(activeQtl, traitQtl)
```

## Arguments

- activeQtl:

  a
  [`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)
  representing all active QTL

- traitQtl:

  a
  [`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)
  representing QTL for a trait of interest

## Value

an integer vector for relative positions of traitQtl on activeQtl
