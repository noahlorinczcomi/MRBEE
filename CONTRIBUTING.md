# Contributing to MRBEE

Thank you for your interest in contributing. This document covers how to set up a development environment, run tests, and submit changes.

## Development setup

**Requirements:** R >= 4.1.0

```r
# Install development dependencies
install.packages(c("devtools", "testthat", "roxygen2", "renv"))

# Clone the repo, then from the project root:
renv::restore()          # restore exact dependency versions
devtools::document()     # regenerate docs from roxygen2 comments
devtools::load_all()     # load the package for interactive use
```

## Running tests

```r
devtools::test()         # run the full test suite
devtools::check()        # run R CMD check (what CI runs)
```

Tests live in `tests/testthat/`. Each test file covers one module:

| File | Covers |
|------|--------|
| `test-allele-harmonise.R` | `allele_harmonise()` |
| `test-MRBEE-IMRP.R` | `MRBEE.IMRP()`, `MRBEE.IMRP.Egger()` |
| `test-MRBEE-IMRP-UV.R` | `MRBEE.IMRP.UV()` |
| `test-Joint-test.R` | `Joint.test()` |
| `test-errorCov.R` | `errorCov()` |

## Making changes

1. Fork the repository and create a branch from `main`.
2. Make your changes. Add or update tests for any changed behavior.
3. Run `devtools::document()` if you edit roxygen2 comments.
4. Run `devtools::check()` and resolve any ERRORs or WARNINGs before opening a PR.
5. Open a pull request against `main` and fill in the PR template.

## Code style

- Use `<-` for assignment (not `=` outside function arguments).
- Use `TRUE`/`FALSE` (not `T`/`F`).
- 2-space indentation.
- Keep lines under 100 characters where practical.

## Updating the lockfile

After adding or upgrading a dependency, regenerate `renv.lock`:

```r
renv::snapshot(type = "explicit")
```

Commit `renv.lock` alongside the `DESCRIPTION` change.

## Reporting bugs

Please open a GitHub issue using the **Bug Report** template. Include the output of `sessionInfo()`.

## Questions

Reach out to the maintainers listed in `DESCRIPTION` or open an issue.
