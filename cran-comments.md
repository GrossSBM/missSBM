
# missSBM 1.1.1

Patch release fixing check failures found by CRAN's own machines after 1.1.0 was accepted:

- r-oldrel-macos-arm64 segfault during checks: caused by `blockmodels`'s
  internal use of fork-based `parallel::mclapply` for parallel
  initialization exploration, which is unsafe on macOS in combination with
  the Accelerate/vecLib BLAS-LAPACK framework. Fixed by forcing
  `nbCores = 1` in the two test files that build a `blockmodels`-based
  reference model for comparison; no change to missSBM's own code.
  (https://www.r-project.org/nosvn/R.check/r-oldrel-macos-arm64/missSBM-00check.html)

No other changes since 1.1.0. See NEWS.md for details.

## Tested environments

* tested locally on Ubuntu Linux 24.04.4 LTS, R-release, GCC

* win-builder: Windows Server 2022, R-oldrelease/R-release/R-devel, 64 bit
* GitHub Actions (R-CMD-check.yaml): Ubuntu 24.04 (R-devel/R-release/R-oldrel-1),
    Windows Server 2022 (R-release), macOS (R-release)

## R CMD check results

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mno-omit-leaf-frame-pointer’

0 errors ✔ | 0 warnings ✔ | 1 note ✖
