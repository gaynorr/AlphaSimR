# Selection index

Calculates values of a selection index given trait values and weights.
This function is intended to be used in combination with selection
functions working on populations such as
[`selectInd`](https://gaynorr.github.io/AlphaSimR/reference/selectInd.md).

## Usage

``` r
selIndex(Y, b, scale = FALSE)
```

## Arguments

- Y:

  a matrix of trait values

- b:

  a vector of weights

- scale:

  should Y be scaled and centered

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
#Model two genetically correlated traits
G = 1.5*diag(2)-0.5 #Genetic correlation matrix
SP$addTraitA(10, mean=c(0,0), var=c(1,1), corA=G)
SP$setVarE(h2=c(0.5,0.5))

#Create population
pop = newPop(founderPop, simParam=SP)

#Calculate Smith-Hazel weights
econWt = c(1, 1)
b = smithHazel(econWt, varG(pop), varP(pop))

#Selection 2 best individuals using Smith-Hazel index
#selIndex is used as a trait
pop2 = selectInd(pop, nInd=2, trait=selIndex,
                 simParam=SP, b=b)
```
