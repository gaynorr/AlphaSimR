# Import genetic map

Formats a genetic map stored in a data.frame to AlphaSimR's internal
format. Map positions must be in Morgans.

## Usage

``` r
importGenMap(genMap)
```

## Arguments

- genMap:

  genetic map as a data.frame. The first three columns must be: marker
  name, chromosome, and map position (Morgans). Marker name and
  chromosome are coerced using as.character.

## Value

a list of named vectors

## Examples

``` r
genMap = data.frame(markerName=letters[1:5],
                    chromosome=c(1,1,1,2,2),
                    position=c(0,0.5,1,0.15,0.4))

asrMap = importGenMap(genMap=genMap)

str(asrMap)
#> List of 2
#>  $ 1: Named num [1:3] 0 0.5 1
#>   ..- attr(*, "names")= chr [1:3] "a" "b" "c"
#>  $ 2: Named num [1:2] 0 0.25
#>   ..- attr(*, "names")= chr [1:2] "d" "e"
```
