# Contributing to AlphaSimR

Thank you for helping to improve AlphaSimR. This page explains how changes
are made to the package. Please also read the
[code of conduct](CODE_OF_CONDUCT.md).

## Reporting bugs and suggesting features

Open an [issue](https://github.com/gaynorr/AlphaSimR/issues) using one of
the templates. For a bug, include a short script that reproduces it,
starting from `library(AlphaSimR)` and setting a seed with `set.seed()`.
Questions about using the package are welcome in
[Discussions](https://github.com/gaynorr/AlphaSimR/discussions).

For anything larger than a small fix, please open an issue to discuss it
before starting work, so that effort is not spent on a change that cannot
be accepted.

## Branches

- `master` is the released version, the one on CRAN.
- `devel` is the development version. All changes go to `devel` first.

Start each change from an up-to-date `devel` on a branch of its own, named
for the change, and open a pull request into `devel`. Pull requests into
`master` are made by the maintainer for releases. Both branches require a
pull request and passing checks.

## Making a change

AlphaSimR is an R package with a large amount of compiled C++ code, so you
need a working build toolchain: Rtools on Windows, the Xcode command line
tools on macOS, or a compiler on Linux.

[`AGENTS.md`](../AGENTS.md) in the root of the repository describes the
package's structure, conventions, and the changes that break things without
an error, such as changes to seeded results. It is written for AI coding
agents but is the most complete guide for anyone changing the code. In
particular:

- **Do not add dependencies** without discussing it first.
- **Follow the style of the function you are changing.** Most of the code
  uses `=` for assignment and no space before an opening parenthesis or
  brace, as in `if(x){`. There is no automatic formatter; do not reformat
  lines your change does not need.
- **Use American English** in code, comments, documentation and tests.
- **Write comments that explain why**, in full sentences.

### Documentation

Help pages are written as roxygen comments in `R/`. You do not need to
regenerate `man/` or `NAMESPACE` yourself: a workflow does it on every push
that changes `R/`, and commits the result to your branch. If you do
document locally, use the roxygen2 version named in `DESCRIPTION`.

A new exported function needs a place in the `reference:` index of
`_pkgdown.yml`, or the website build fails.

### Tests

Add tests for new behavior in `tests/testthat/`, following the file you are
adding to. Run the tests with:

```r
devtools::test()
```

Slow tests are skipped unless `NOT_CRAN` is set, as it is by
`devtools::test()` and in continuous integration. Before opening a pull
request, check the package with:

```r
devtools::check()
```

### Seeded results

A simulation run from a fixed seed must give the same results. Changing
the order in which random numbers are drawn, or how many individuals one
call produces, changes every seeded result even though no test of behavior
fails. `tests/testthat/test-seededResults.R` stores expected results from
seeded simulations and fails when they change. If your change moves them
on purpose, update the stored values and describe the change under
"Changes to results from a seed" in `NEWS.md`.

### NEWS

Add an entry to the top section of `NEWS.md` for anything a user would
notice, in a sentence or two, under the heading it belongs to. Do not
change `Version` or `Date` in `DESCRIPTION`; the maintainer sets those.

## Continuous integration

Every pull request into `devel` or `master` runs `R CMD check` on Windows,
macOS and Linux, on the current, previous, development and oldest
supported versions of R, builds the package website, and measures test
coverage. The website is published only from `master`. Sanitizer and
valgrind checks of the compiled code can be run by hand with the R-hub
workflow.
