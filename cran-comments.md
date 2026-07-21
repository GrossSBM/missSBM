
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

## Breaking change

`estimateMissSBM()`'s `control` argument must now be built with the new
`missSBM_param()` helper (one named, documented, defaulted argument per
option) instead of a raw `list(...)`; passing a plain list now errors with a
message pointing at the replacement. This also means unknown/typo'd option
names now fail immediately instead of being silently ignored.

We checked the single reverse dependency, `gsbm` (which has missSBM in
`Suggests`, used only in a vignette): it calls
`missSBM::estimateMissSBM(A, vBlocks, "node")` with no `control` argument, so
it is unaffected by this change.

Also in this release: several new opt-in/on-by-default features to make
model selection across numbers of blocks more robust when the requested
number of blocks exceeds what the network actually supports (VEM class
collapse recovery, a node-swap polishing step, an alternative warm-start
initialization scheme), a refactor exposing `missSBM_fit`'s split/merge
exploration as instance methods, and several profiling-driven performance
fixes to `estimateMissSBM()`'s hot paths (roughly 3x faster on our
benchmark, bit-identical results). See NEWS.md for the full changelog.

## Tested environments

* tested locally on Ubuntu Linux 24.04.4 LTS, R-release, GCC

* win-builder: Windows Server 2022, R-oldrelease/R-release/R-devel, 64 bit
* GitHub Actions (R-CMD-check.yaml): Ubuntu 24.04 (R-devel/R-release/R-oldrel-1),
    Windows Server 2022 (R-release), macOS (R-release)

## R CMD check results

Duration: 2m 27.6s

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mno-omit-leaf-frame-pointer’

0 errors ✔ | 0 warnings ✔ | 1 note ✖
