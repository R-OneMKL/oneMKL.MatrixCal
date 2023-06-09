% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneMKL.MatrixCal-package.R
\docType{package}
\name{onemkl.MatrixCal-package}
\alias{onemkl.MatrixCal-package}
\title{`oneMKL.MatrixCal` Package}
\description{
The `oneMKL.MatrixCal` package provides an efficient and powerful tool for
conducting linear algebra computations in the R environment. In the package,
matrix operations, such as factorization, inversion, and multiplication,
are performed using the optimized algorithms embedded in the
**Intel oneAPI Math Kernel Library** (`oneMKL`). Because of this feature,
the computational performance can be significantly improved. To make the package
function work correctly with Intel `oneMKL`, users must first install the package,
`oneMKL`, to import Intel `oneMKL` into R. The package also supports multithreaded
computations, enabling users to utilize the benefit of multiple CPU cores,
i.e., parallel computing, to handle large-scale matrix computations.
The package is not only compatible with Linux, but aslo Windows operating systems,
which can benefit a wide range of users.
}
