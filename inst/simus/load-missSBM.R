## Loads missSBM for a campaign (sourced by the scripts, which run from the root of the package).
##
##  - SIMU_LIB set: the version installed in this library, compiled with optimization as for a
##    user. To compare two versions, install each one in its own library, e.g.
##      git worktree add ../missSBM-v1.1.0 v1.1.0
##      R CMD INSTALL --preclean -l <lib-v1.1.0> ../missSBM-v1.1.0
##      R CMD INSTALL --preclean -l <lib-v1.2.0> .
##  - otherwise: the sources of the working directory, with devtools::load_all(). The C++ code is
##    then compiled without optimization: fine for statistical results, not for timings.
##
## The internal functions are attached in both cases. `simu_version` ("v1.2.0") labels the
## results by default; `simu_source` says where the package comes from.
lib <- Sys.getenv("SIMU_LIB", "")
if (nzchar(lib)) {
  suppressMessages(library(missSBM, lib.loc = lib))
  attach(asNamespace("missSBM"), name = "missSBM_internals", warn.conflicts = FALSE)
  simu_version <- paste0("v", utils::packageVersion("missSBM", lib.loc = lib))
  simu_source  <- "installed"
} else {
  suppressMessages(devtools::load_all(quiet = TRUE))
  simu_version <- paste0("v", read.dcf("DESCRIPTION")[, "Version"])
  simu_source  <- paste("sources at", system("git rev-parse --short HEAD", intern = TRUE))
}
future::plan("sequential") # the campaigns run their jobs in parallel instead
