# Test if individuals of a population are female or male

Test if individuals of a population are female or male

## Usage

``` r
isFemale(x)

isMale(x)
```

## Arguments

- x:

  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Value

logical

## Functions

- `isMale()`: Test if individuals of a population are female or male

## Examples

``` r
founderGenomes <- quickHaplo(nInd = 3, nChr = 1, segSites = 100)
SP <- SimParam$new(founderGenomes)
SP$setSexes(sexes = "yes_sys")
pop <- newPop(founderGenomes)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found

isFemale(pop)
#> Error: object 'pop' not found
isMale(pop)
#> Error: object 'pop' not found

pop[isFemale(pop)]
#> Error: object 'pop' not found
pop[isFemale(pop)]@sex
#> Error: object 'pop' not found
```
