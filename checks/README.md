# checks/

Scripts for verifying numerical consistency across branches.

## compare_branches.R

Sources the `R/` files from two branches into isolated R environments and checks
that every exported function returns numerically identical results.

### Prerequisites

- **MASS** — already in `renv.lock`
- **optparse** — install separately; not a package dependency:
  ```r
  install.packages("optparse")
  ```

### Usage

```
Rscript checks/compare_branches.R [options]

Options:
  -a BRANCH, --branch-a=BRANCH  First git branch or ref  [default: HEAD]
  -b BRANCH, --branch-b=BRANCH  Second git branch or ref [default: main]
  --label-a=LABEL               Display name for branch A [default: -a value]
  --label-b=LABEL               Display name for branch B [default: -b value]
  -t TOL, --tol=TOL             Numerical tolerance for all.equal() [default: 1e-8]
  -h, --help                    Show this help message
```

### Examples

Compare the current branch against `main`:

```bash
Rscript checks/compare_branches.R
```

Specify branches and human-readable labels explicitly:

```bash
Rscript checks/compare_branches.R \
  --branch-a admin \
  --branch-b main \
  --label-a admin \
  --label-b main
```

Use a stricter tolerance:

```bash
Rscript checks/compare_branches.R --tol 1e-12
```
