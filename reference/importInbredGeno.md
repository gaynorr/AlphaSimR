# Import inbred, diploid genotypes

Formats the genotypes from inbred, diploid lines to an AlphaSimR
population that can be used to initialize a simulation. An attempt is
made to automatically detect 0,1,2 or -1,0,1 genotype coding.
Heterozygotes or probabilistic genotypes are allowed, but will be
coerced to the nearest homozygote. Pedigree information is optional and
when provided will be passed to the population for easier identification
in the simulation.

## Usage

``` r
importInbredGeno(geno, genMap, ped = NULL)
```

## Arguments

- geno:

  a matrix of genotypes

- genMap:

  genetic map as a data.frame. The first three columns must be: marker
  name, chromosome, and map position (Morgans). Marker name and
  chromosome are coerced using as.character. See
  [importGenMap](https://gaynorr.github.io/AlphaSimR/reference/importGenMap.md)

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
geno = rbind(c(2,2,0,2,0),
             c(0,2,2,0,0))
colnames(geno) = letters[1:5]

genMap = data.frame(markerName=letters[1:5],
                    chromosome=c(1,1,1,2,2),
                    position=c(0,0.5,1,0.15,0.4))

ped = data.frame(id=c("a","b"),
                 mother=c(0,0),
                 father=c(0,0))

founderPop = importInbredGeno(geno=geno,
                              genMap=genMap,
                              ped=ped)
```
