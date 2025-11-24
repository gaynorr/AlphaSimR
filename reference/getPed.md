# Get pedigree

Returns the population's pedigree as stored in the id, mother and father
slots. NULL is returned if the input population lacks the required.

## Usage

``` r
getPed(pop)
```

## Arguments

- pop:

  a population

## Examples

``` r
# Create a founder population
founderPop = quickHaplo(2,1,2)

# Set simulation parameters
SP = SimParam$new(founderPop)

# Create a population
pop = newPop(founderPop, simParam=SP)

# Get the pedigree
getPed(pop)
#>   id mother father
#> 1  1      0      0
#> 2  2      0      0

# Returns NULL when a population lacks a pedigree
getPed(founderPop)
#> NULL
```
