# Create new population

Creates a new
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
from an object of of the Pop superclass.

## Usage

``` r
.newPop(
  rawPop,
  id = NULL,
  mother = NULL,
  father = NULL,
  iMother = NULL,
  iFather = NULL,
  isDH = NULL,
  femaleParentPop = NULL,
  maleParentPop = NULL,
  hist = NULL,
  simParam = NULL,
  ...
)
```

## Arguments

- rawPop:

  an object of the pop superclass

- id:

  optional id for new individuals

- mother:

  optional id for mothers

- father:

  optional id for fathers

- iMother:

  optional internal id for mothers

- iFather:

  optional internal id for fathers

- isDH:

  optional indicator for DH/inbred individuals

- femaleParentPop:

  optional population of female parents

- maleParentPop:

  optional population of male parents

- hist:

  optional recombination history

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

- ...:

  additional arguments passed to the finalizePop function in simParam

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
