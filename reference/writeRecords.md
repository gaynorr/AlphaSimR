# Write data records

Saves a population's phenotypic and marker data to a directory.

## Usage

``` r
writeRecords(
  pop,
  dir,
  snpChip = 1,
  useQtl = FALSE,
  includeHaplo = FALSE,
  append = TRUE,
  simParam = NULL
)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- dir:

  path to a directory for saving output

- snpChip:

  which SNP chip genotype to save. If useQtl=TRUE, this value will
  indicate which trait's QTL genotype to save. A value of 0 will skip
  writing a snpChip.

- useQtl:

  should QTL genotype be written instead of SNP chip genotypes.

- includeHaplo:

  should markers be separated by female and male haplotypes.

- append:

  if true, new records are added to any existing records. If false, any
  existing records are deleted before writing new records. Note that
  this will delete all files in the 'dir' directory.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)
