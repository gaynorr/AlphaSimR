# AGENTS.md

Guidance for AI agents working on AlphaSimR. It covers what is not visible
in the code, and the mistakes this package punishes quietly. Read it before
changing anything.

## Ground rules

**Do not change any files without express consent.** Propose the change,
say what you would edit, and wait. This applies to every file, including
tests and documentation. The package has compiled code, so without an R
toolchain nothing you write can be verified: once a change has been agreed
and made, say plainly what has not been built or run and what the
maintainer needs to do.

**Do not add dependencies.** `Depends` is R (>= 4.0.0); `Imports` is
`Rcpp`, `Rdpack`, `methods` and `R6`; `LinkingTo` is `Rcpp`,
`RcppArmadillo`, `BH` and `dqrng`; `Suggests` is `knitr`, `rmarkdown` and
`testthat`. That list is a decision rather than an accident, and adding to
it needs discussing first.

## Map of the package

The S4 population classes branch rather than form a chain. `RawPop` is the
base. One branch adds the genome: `MapPop` extends `RawPop` with a genetic
map, and `NamedMapPop` extends `MapPop` with IDs. The other branch is
`Pop`, which extends `RawPop` directly with traits, phenotypes and a
simulation's bookkeeping. A `Pop` therefore has no genetic map of its own;
the map lives in the `SimParam`. `MultiPop` is a container of populations,
and `HybridPop`, in `R/Class-HybridPop.R`, is separate from all of them.

- `R/Class-SimParam.R` is the largest file in the package at around three
  thousand lines.
- `R/Class-Pop.R` holds the classes above along with their subsetting and
  combining methods, which is where several of the traps below live.
- `R/mergePops.R` holds `mergePops`, which is how populations are combined
  and is what `c()` on a `Pop` routes through.
- `R/crossing.R` holds the crossing functions and is the hot path. Per-call
  overhead there is large relative to the per-individual work, so a loop
  calling `makeCross` once per individual is far slower than one call
  making the same individuals together. `pedigreeCross` batches by
  generation for this reason. Batching is not a free speed-up: it changes
  results, for the reason given under reproducibility below.
- `src/rng.cpp` holds the random number design. A C++ generator is seeded
  from R's state, and randomness inside a call comes from substreams of it
  rather than from R.
- `.github/workflows/` holds `R-CMD-check.yaml`, `pkgdown.yaml`,
  `test-coverage.yaml` and `document.yaml`. Only the last of these writes
  to the repository; see Generated files.
- `vignettes/articles/` is built by pkgdown but not woven by `R CMD build`,
  because vignette subdirectories are not searched. Its code still has to
  run when the site is built.

## Code style

The prevailing style is `=` for assignment and no space before an opening
paren or brace, as in `function(x){`, `if(x){` and `}else{`. `R/crossing.R`
is entirely in it: ninety-seven `if(`, no `if (`, no `<-`. Newer additions
use `if (x) {` and occasionally `<-`, and they sit inside otherwise
old-style files; in `R/phenotypes.R` and `R/popSummary.R` the two styles
are at about parity. Follow the function you are working in rather than the
file as a whole.

Six files carry `# fmt: skip file` on their first line: `R/AlphaSimR.R`,
`R/Class-Pop.R`, `R/mergePops.R`, `R/popSummary.R`,
`tests/testthat/test-mergePops.R` and `tests/testthat/test-popSummary.R`.
Do not add or remove the marker.

Comments explain *why*, in full prose sentences, in `src/` as well as `R/`.

Use American English spelling in code, comments, documentation, tests
and `NEWS.md`: "penalize", "behavior", "toward", not "penalise",
"behaviour", "towards".

## Changes that break things silently

None of these produces an error, and most will not fail an existing test.

**Reproducibility depends on the order of operations, and no test defends
it.** A simulation run from a fixed seed must give the same answer. Random
numbers are consumed as individuals are made, so changing the order in
which crosses happen, or how many individuals one call produces, changes
every seeded result. A crossing call builds `nChr * countBlocks(n)` RNG
substreams and every individual within a block shares one, so the batch
size feeds directly into which numbers an individual gets. The `n` there is
the number of parents in most paths but the number of progeny in the
`nProgeny` and reduce-genome paths, so read the call site rather than
assuming.

Separately, `.newPop` draws from R's own stream: `sample()` for sex when
`sexes` is `"yes_rand"`, and whatever `setPheno` uses. So R's stream
advances once per population built, on top of the C++ substreams.

There is no stored expected output anywhere in the package: no
`tests/testthat/_snaps`, no `expect_snapshot`, no fixtures.
`tests/testthat/test-reproducibility.R` only compares one run against
another run in the same build — the same seed twice, and one thread against
two — and the thread comparisons skip themselves when only one thread is
available. A change that moves the RNG stream therefore passes the entire
suite. Review and a `NEWS.md` entry are the only defense.

**Results must not depend on `nThreads`.** The invariant is the constant in
`src/misc.h`:

    const arma::uword nWorkBlocks = 64;

`countBlocks` in `src/misc.cpp` returns one block for an empty input, one
block per item below sixty-four, and sixty-four above that. The number of
blocks is therefore fixed by the size of the problem and never by the
thread count, which fixes the split of the work, the number of accumulators
a function allocates, and the order those accumulators are summed in, on
every machine. `blockStart` hands out the ranges, and `makeWorkRngs` in
`src/meiosis.cpp` gives each of the `nChr * nBlocks` work items its own RNG
substream. The dependency runs one way only: the thread count is capped by
the work, never the reverse. Deriving the block count from `nThreads` is
the obvious optimization, it compiles, and it passes every test.

Most other block-partitioned loops draw no random numbers, so for them the
invariant is about summation order alone; that is true of `MME.cpp`,
`getGv.cpp`, `getHybridGv.cpp` and `calcGenParam.cpp`.
`src/altAddTraitAD.cpp` is the exception: it draws `2 * nLoci` normal
deviates outside its block loop, so editing it can move seeded results even
though its loop is a plain reduction.

**`SimParam` is an R6 class, so it has reference semantics.** Assigning to
a field of a `simParam` inside a function changes the caller's simulation
permanently, not just for the call. Never overwrite a user's settings this
way: `pedigreeCross` refuses the recombination settings `v`, `p` and
`quadProb` when it is given a `Pop` for exactly this reason. Some fields
are plain entries in the R6 `public` list, including those three, and
others are active bindings over `private` ones with validation attached.
Check which kind a neighbouring field is before adding one.

**Build new individuals with `.newPop`, not with `new("Pop", ...)`.**
`.newPop` derives each individual's `iid` from `simParam$lastId` and then
calls exactly one of `addToRec`, `addToPed` or `updateLastId`, depending on
what the `SimParam` is tracking, to record the new individuals and move the
counter on. Constructing a `Pop` directly for individuals that do not yet
exist skips all of that and leaves the ID counter and the recorded pedigree
quietly out of step. Constructing one directly is fine where no new
individuals are involved, which is why `newEmptyPop` and `mergePops` both
do it.

Calling the recorders by hand is riskier than it looks. `addToRec` and
`addToPed` infer the count from `lastId - simParam$lastId` and reject a
non-positive result, so a repeat raises an error. `updateLastId`, which is
the path taken when the pedigree is not being tracked, only checks that the
new value is greater than *or equal to* the old one: calling it twice does
nothing and calling it with too large a value silently skips IDs. Nothing
will tell you.

`.newPop` also guards the genetic value and phenotype work behind
`if(simParam$nTraits>=1)`, so a `SimParam` with no traits is legal.
`pedigreeCross` relies on this when it builds a temporary `SimParam` for a
map population. It ends by calling `simParam$finalizePop`, a user-supplied
hook that runs on every population the package creates, so a change to what
`.newPop` returns has to leave that hook working.

**Subsetting and merging lose the miscellaneous slots.** There is no
replacement method for `RawPop`, `MapPop`, `NamedMapPop` or `Pop`, so
changing part of a population means subsetting it, rebuilding with
`mergePops()`, and index-subsetting the result with a permutation to
restore the original order. Subsetting empties `miscPop` outright, and
`mergePops` drops `misc` with a warning when the populations do not carry
the same named elements. Check what is in them before taking that route.
`MultiPop` is the exception to all of this: it has `[<-`, `[[<-`, `$<-` and
`names<-`.

**`c()` on a `Pop` can return something other than a `Pop`.** If any
argument is a `MultiPop`, the method delegates and the result is a
`MultiPop`. Otherwise it routes through `mergePops()`.

**Three of the six `c()` methods test the class strictly, on purpose.**
`RawPop` and `MapPop` compare with `identical()` against
`as.character(class(y))` rather than using `is()`, because both classes are
extended and `is()` would let a subclass be combined into its parent and
silently drop the extra slots. The code comment above the `RawPop` method
explains why the comparison is written the way it is; read it before
simplifying. `HybridPop` tests strictly too, in a different form.
`NamedMapPop`, `Pop` and `MultiPop` use `is()` or no class test at all, and
are right to, because nothing extends those classes. Do not make the six
consistent with each other.

## Generated files

`man/*.Rd` and `NAMESPACE` are generated from roxygen comments, and
`.github/workflows/document.yaml` regenerates them on every push that
touches `R/**`, committing `man/`, `NAMESPACE` and `DESCRIPTION` back to
the branch. So `man/` in a working tree may lag behind the roxygen above
it, and that is normal rather than something to fix by hand. The workflow
installs whatever roxygen is current, so it can rewrite the
`Config/roxygen2/version` pin in `DESCRIPTION` even though nothing else
should edit that file. If you document locally, use the pinned version; a
different one rewrites every `.Rd` and buries the real change in the diff.

`src/RcppExports.cpp` and `R/RcppExports.R` come from
`Rcpp::compileAttributes()`, which no workflow runs. Changing a
`// [[Rcpp::export]]` signature means regenerating them and rebuilding.

## Third-party code in src/

Part of `src/` is vendored from the MaCS simulator. `LICENSE.note` lists
which files, and each of them says so in its own header; check one or the
other rather than guessing from names, because several files that look
third-party are not. Those files have been modified by the AlphaSimR
authors and may be modified further, but `LICENSE.note` carries a summary
of the modifications and any change to one of them needs a matching change
there.

## Tests

- testthat, edition 2. This is the default by omission: there is no
  `Config/testthat/edition` field in `DESCRIPTION`, and adding one would
  retire `context()` and change `expect_equal()` semantics package-wide.
  Do not add it.
- Most files open with `context()`; a few do not. Follow the surrounding
  file.
- `skip_on_cran()` guards two kinds of test: slow ones, and ones that
  compare a Monte Carlo estimate against a tolerance and so fail
  occasionally. `skip_if(getNumThreads() < 2)` or `skip_if_not` guards
  anything that compares thread counts.
- Internal functions are tested through `AlphaSimR:::`.
- Larger test files divide blocks into sections with numbered banner
  comments, and most block names carry an uppercase tag naming the
  behavior under test, as in `test_that("EXTEND a parent used once ...")`.
  `tests/testthat/test-pedigreeCross.R` is the current example.
- A genetic test should assert something that must hold rather than a value
  that happened to come out. The Mendelian bound holds for a direct cross;
  it does not survive selfing or doubled haploid production, where only the
  weaker "a locus fixed in both ancestors stays fixed" bound does.
- Roxygen examples that build a `SimParam` nearly always include
  `\dontshow{SP$nThreads = 1L}` so that checks run single threaded. New
  examples need it too.

### Local test coverage

For behavior-changing work, use local test coverage to check whether tests
exercise the changed lines and to investigate gaps reported by CI.
From the package root (the directory containing `DESCRIPTION`), run:

```sh
Rscript -e 'cov <- covr::package_coverage(clean = TRUE); print(cov); covr::report(cov)'
```

This requires `covr` and a working package build toolchain. Inspect the
changed functions in the report and add tests for relevant uncovered paths.
Coverage shows which lines ran; assertions still need to verify the expected
behavior. This supplements the existing tests and R CMD check.

## Finishing a change

Add an entry to `NEWS.md` under its top heading for anything a user would
notice. Do not edit `Version` or `Date` in `DESCRIPTION`; the maintainer
sets those. `Version` in `DESCRIPTION` and the top heading of `NEWS.md` are
kept in step for releases and may drift apart on the development branch, so
a difference between those two is not a bug to fix.

A new file at the top level of the package needs an entry in
`.Rbuildignore` unless it is meant to ship. This file has such an entry and
does not ship.
