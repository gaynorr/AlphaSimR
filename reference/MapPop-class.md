# Raw population with genetic map

Extends
[`RawPop-class`](https://gaynorr.github.io/AlphaSimR/reference/RawPop-class.md)
to add a genetic map. This is the first object created in a simulation.
It is used for creating initial populations and setting traits in the
[`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md).

## Usage

``` r
# S4 method for class 'MapPop'
x[i]

# S4 method for class 'MapPop'
c(x, ...)

isMapPop(x)
```

## Arguments

- x:

  a 'MapPop' object

- i:

  index of individuals

- ...:

  additional 'MapPop' objects

## Methods (by generic)

- `[`: Extract MapPop by index

- `c(MapPop)`: Combine multiple MapPops

## Functions

- `isMapPop()`: Test if object is of a MapPop class

## Slots

- `genMap`:

  list of chromosome genetic maps

- `centromere`:

  vector of centromere positions

- `inbred`:

  indicates whether the individuals are fully inbred
