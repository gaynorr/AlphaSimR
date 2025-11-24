# Import haplotypes

Formats haplotype in a matrix format to an AlphaSimR population that can
be used to initialize a simulation. This function serves as wrapper for
[`newMapPop`](https://gaynorr.github.io/AlphaSimR/reference/newMapPop.md)
that utilizes a more user friendly input format.

## Usage

``` r
importHaplo(haplo, genMap, ploidy = 2L, ped = NULL)
```

## Arguments

- haplo:

  a matrix of haplotypes

- genMap:

  genetic map as a data.frame. The first three columns must be: marker
  name, chromosome, and map position (Morgans). Marker name and
  chromosome are coerced using as.character. See
  [`importGenMap`](https://gaynorr.github.io/AlphaSimR/reference/importGenMap.md)

- ploidy:

  ploidy level of the organism

- ped:

  an optional pedigree for the supplied genotypes. See details.

## Value

a
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
if ped is NULL, otherwise a
[`NamedMapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/NamedMapPop-class.md)

## Details

The optional pedigree can be a data.frame, matrix or a vector. If the
object is a data.frame or matrix, the first three columns must include
information in the following order: id, mother, and father. All values
are coerced using as.character. If the object is a vector, it is assumed
to only include the id. In this case, the mother and father will be set
to "0" for all individuals.

## Examples

``` r
haplo = rbind(c(1,1,0,1,0),
              c(1,1,0,1,0),
              c(0,1,1,0,0),
              c(0,1,1,0,0))
colnames(haplo) = letters[1:5]

genMap = data.frame(markerName=letters[1:5],
                    chromosome=c(1,1,1,2,2),
                    position=c(0,0.5,1,0.15,0.4))

ped = data.frame(id=c("a","b"),
                 mother=c(0,0),
                 father=c(0,0))

founderPop = importHaplo(haplo=haplo, 
                         genMap=genMap,
                         ploidy=2L,
                         ped=ped)
```
