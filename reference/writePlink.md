# Writes a Pop-class as PLINK files

Writes a Pop-class to PLINK PED and MAP files. The arguments for this
function were chosen for consistency with
[`RRBLUP2`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP2.md).
The base pair coordinate will the locus position as stored in AlphaSimR
and not an actual base pair position. This is because AlphaSimR doesn't
track base pair positions, only relative positions for the loci used in
the simulation.

## Usage

``` r
writePlink(
  pop,
  baseName,
  traits = 1,
  use = "pheno",
  snpChip = 1,
  useQtl = FALSE,
  simParam = NULL,
  ...
)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- baseName:

  basename for PED and MAP files.

- traits:

  an integer indicating the trait to write, a trait name, or a function
  of the traits returning a single value.

- use:

  what to use for PLINK's phenotype field. Either phenotypes "pheno",
  genetic values "gv", estimated breeding values "ebv", breeding values
  "bv", or random values "rand".

- snpChip:

  an integer indicating which SNP chip genotype to use

- useQtl:

  should QTL genotypes be used instead of a SNP chip. If TRUE, snpChip
  specifies which trait's QTL to use, and thus these QTL may not match
  the QTL underlying the phenotype supplied in traits.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

- ...:

  additional arguments if using a function for traits

## Examples
