# Branch Comparison: `admin` vs `main`

**Date:** 2026-04-25
**Script:** `checks/compare_branches.R`

## Method

Both branches are sourced into isolated R environments using `sys.source()` (no package installation required beyond what is in `renv.lock`). The same synthetic dataset (`set.seed(1234)`, n = 200, p = 3) is passed to each function and outputs are compared with `all.equal(..., tolerance = 1e-8)`.

The script creates and cleans up git worktrees automatically. Pass branch names directly:

```bash
Rscript checks/compare_branches.R --branch-a admin --branch-b main
```

## Results

| Function | Output | Result |
|---|---|---|
| `MRBEE.IMRP` | `theta` | PASS |
| `MRBEE.IMRP` | `covtheta` | PASS |
| `MRBEE.IMRP` | `delta` | PASS |
| `MRBEE.IMRP.Egger` | `theta` | PASS |
| `MRBEE.IMRP.Egger` | `covtheta` | PASS |
| `MRBEE.IMRP.Egger` | `delta` | PASS |
| `MRBEE.IMRP.UV` | `theta` | PASS |
| `MRBEE.IMRP.UV` | `vartheta` | PASS |
| `MRBEE.IMRP.UV` | `delta` | PASS |

**9 / 9 checks passed.** Both branches produce numerically identical results.
