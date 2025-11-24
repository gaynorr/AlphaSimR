# Find LociMap superset

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
findLociMapSuperset(lociMap1, lociMap2)
```

## Arguments

- lociMap1:

  a
  [`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)
  that is tested to determine if it is a superset

- lociMap2:

  a second
  [`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)
  that is tested

## Value

NULL if locMap1 is a superset, or a
[`LociMap-class`](https://gaynorr.github.io/AlphaSimR/reference/LociMap-class.md)
if it is not
