# CRAN submission: RTsampler 1.0.0

## Test environments
* Local: Rocky Linux 8.10, R 4.4.1 (GCC 8.5.0), x86_64
* win-builder: R 4.5.1 (ucrt) — release & devel — OK
* R-hub: windows-x86_64, macOS (arm64), linux-x86_64 — OK

## R CMD check results
0 errors | 0 warnings | 2 notes

### Note 1: Installed package size
* Installed size is ~7.0 MB (libs ~6.8 MB).
* Reason: compiled C++ with Armadillo templates (Rcpp/RcppArmadillo backends) increases the size of the shared library.
* Mitigation: no bundled large assets; no debug symbols; examples and tests use small data and run quickly.

### Note 2: Non-standard top-level file
* A `cran-comments.md` file is included at the top level for the reviewers’ convenience.

## Policies & portability
* No network access; no writing outside `tempdir()`; parallelism disabled by default on CRAN.
* Uses BLAS/LAPACK via `$(BLAS_LIBS) $(LAPACK_LIBS) $(FLIBS)`; no hard-coded toolchain flags.
* C++ standard: `CXX_STD = CXX17` (declared in `src/Makevars{.win}`); `SystemRequirements: C++17` is noted in `DESCRIPTION`.


