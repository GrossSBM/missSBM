
# missSBM 1.1.1

Patch release fixing two install/check failures found by CRAN's own
pretest/check machines after 1.1.0 was accepted:

- Debian pretest install failure (`undefined symbol: dgelsd_`): caused by
  `.Rbuildignore` erroneously excluding `src/Makevars` from the built
  tarball, so `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` was never
  shipped and LAPACK's `dgelsd` (used internally by RcppArmadillo's
  `arma::solve()`) was left unresolved on that machine.
  (https://win-builder.r-project.org/incoming_pretest/missSBM_1.1.0_20260721_165706/Debian/00install.out)

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
