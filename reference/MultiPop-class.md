# Multi-Population

The mega-population represents a population of populations. It is
designed to behave like a list of populations.

## Usage

``` r
# S4 method for class 'MultiPop'
x[i]

# S4 method for class 'MultiPop'
x[[i]]

# S4 method for class 'MultiPop'
c(x, ...)

# S4 method for class 'MultiPop'
length(x)

isMultiPop(x)
```

## Arguments

- x:

  a 'MultiPop' object

- i:

  index of populations or mega-populations

- ...:

  additional 'MultiPop' or 'Pop' objects

## Methods (by generic)

- `[`: Extract MultiPop by index

- `[[`: Extract Pop by index

- `c(MultiPop)`: Combine multiple MultiPops

- `length(MultiPop)`: Number of pops in MultiPop

## Functions

- `isMultiPop()`: Test if object is of a MultiPop class

## Slots

- `pops`:

  list of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  and/or `MultiPop-class`
