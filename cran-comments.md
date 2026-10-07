## Submission

spCF 0.2.2.1 fixes the test failure reported by CRAN for 0.2.2 under the
alternative BLAS checks (BLIS, OpenBLAS), to be corrected before 2026-10-27.

* The failing test (test-predict.R) compared predictions made in different
  batches for exact identity. With BLIS/OpenBLAS the matrix products of
  different sizes are evaluated in a different order, and the results differed
  by about 1e-16. The predict() tests now compare with a tolerance of 1e-10.
  The package code is unchanged.

## Test environments

* local macOS (aarch64-apple-darwin), R 4.6.0

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.
