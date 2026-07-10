
# missSBM 1.1.0

Fixes the WARNING currently shown on the CRAN check page for missSBM
(https://cran.r-project.org/web/checks/check_results_missSBM.html) on
r-devel-linux-x86_64-debian-gcc, r-devel-linux-x86_64-fedora-gcc and
r-patched-linux-x86_64: a `-Wmismatched-new-delete` warning raised in
RcppArmadillo's `memory.hpp`, triggered while compiling `src/packing.cpp`.
This release replaces the NLopt/CCSAQ optimizer previously used to fit the
covariate connectivity parameters with a builtin Newton-Raphson solver;
`src/packing.cpp` and `src/nlopt_wrapper.cpp` are removed and `nloptr` is no
longer a dependency, so the warning has nothing left to be raised from.

Also in this release:

* `missSBM_fit` now exposes `split()`, `merge()`, `candidates_split()` and
  `candidates_merge()` as instance methods (previously inlined in
  `missSBM_collection`)
* several profiling-driven performance fixes to `estimateMissSBM()`'s hot
  paths (roughly 3x faster on our benchmark, bit-identical results)
* bug fixes: a biased fill value in `partlyObservedNetwork$imputation()`, a
  crash and a closed-form algebra bug in the "degree" sampling design, and a
  step-back consistency bug in `missSBM_fit$doVEM()`

See NEWS.md for the full changelog.

## Tested environments

* tested locally on Ubuntu Linux 24.04.2 LTS, R-release, GCC

* tested remotely with win-builder
  - Windows Server 2022, R-oldrelease, 64 bit
  - Windows Server 2022, R-release, 64 bit
  - Windows Server 2022, R-devel, 64 bit

* tested remotely via GitHub Actions (R-CMD-check.yaml)
  - Ubuntu 24.04, R-devel
  - Ubuntu 24.04, R-release
  - Ubuntu 24.04, R-oldrel-1
  - Windows Server 2022, R-release, 64 bit
  - macOS, R-release

## R CMD check results

0 errors | 0 warnings | 2 NOTEs

* compilation used a non-portable compiler flag (`-mno-omit-leaf-frame-pointer`);
  this comes from the local R/toolchain configuration, not from the package's
  own build flags
* two examples (`estimateMissSBM`, `missSBM_fit`) take marginally more than
  5s to run
