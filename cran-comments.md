# missSBM 1.2.0

Update of missSBM 1.1.0 (1.1.1 was prepared but not submitted; its fix is included).

- The estimation of the missing-not-at-random sampling designs and the model selection by the
  ICL change (see NEWS.md): the results of some designs, and the selected models, can differ from
  those of 1.1.0.
- Memory and speed: no n x n matrix for the designs without covariates, several computations in
  C++.
- Bug fixes, including the check failure found by CRAN on r-oldrel-macos-arm64 after 1.1.0 was
  accepted (a segfault in two test files building a reference model with `blockmodels`, through
  `sbm`, in parallel with fork-based `mclapply`): these tests now run sequentially
  (`nbCores = 1`).
  (https://www.r-project.org/nosvn/R.check/r-oldrel-macos-arm64/missSBM-00check.html)
- A new vignette.

## Tested environments

* local: Ubuntu Linux 24.04, R-release, GCC
* win-builder: R-oldrelease, R-release, R-devel
* GitHub Actions (R-CMD-check.yaml): Ubuntu (R-devel, R-release, R-oldrel-1), Windows
  (R-release), macOS (R-release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s): '-mno-omit-leaf-frame-pointer'

  This flag comes from the local R configuration, not from the package.

## Reverse dependencies

gsbm (Suggests) uses missSBM in its vignette (`estimateMissSBM(A, 1:10, "node")`): this call is
unchanged in missSBM 1.2.0.
