## What this changes

<!-- Describe the change and why it is needed. Link any issue it closes, as in "Closes #123". -->

## Checklist

- [ ] The branch starts from an up-to-date `devel`, and this pull request targets `devel`
- [ ] Tests cover the new or changed behavior, and `devtools::check()` passes
- [ ] `NEWS.md` has an entry for anything a user would notice
- [ ] Seeded results are unchanged, or the stored values in `tests/testthat/test-seededResults.R` are updated and the change is listed under "Changes to results from a seed" in `NEWS.md`
- [ ] No dependency is added to `DESCRIPTION`, or it has been agreed in an issue
- [ ] A new exported function is in the `reference:` index of `_pkgdown.yml`
- [ ] A change to `Depends: R` in `DESCRIPTION` is matched by the oldest-R job in `.github/workflows/R-CMD-check.yaml`, and by the required checks
- [ ] Text written by an AI agent and posted to GitHub is marked as such (see `AGENTS.md`)
