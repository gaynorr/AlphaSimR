# Raw population with genetic map and id

Extends
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
to add id, mother and father.

## Usage

``` r
# S4 method for class 'NamedMapPop'
x[i]

# S4 method for class 'NamedMapPop'
c(x, ...)

isNamedMapPop(x)
```

## Arguments

- x:

  a 'NamedMapPop' object

- i:

  index of individuals

- ...:

  additional 'NamedMapPop' objects

## Methods (by generic)

- `[`: Extract NamedMapPop by index

- `c(NamedMapPop)`: Combine multiple NamedMapPops

## Functions

- `isNamedMapPop()`: Test if object is a NamedMapPop class

## Slots

- `id`:

  an individual's identifier

- `mother`:

  the identifier of the individual's mother

- `father`:

  the identifier of the individual's father
