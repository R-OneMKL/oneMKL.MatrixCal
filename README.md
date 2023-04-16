## oneMKL.MatrixCal

The `oneMKL.MatrixCal` package provide matrix operation functions based on Intel® oneAPI Math Kernel Library (oneMKL).

### Installation

1. To build the package from source, Windows users will need to have [Rtools](http://cran.csie.ntu.edu.tw/bin/windows/Rtools/) installed.

2.  "oneMKL.MatrixCal" package is only available for Windows and Linux operating systems due to the availability of the Anaconda mkl package on these platforms. It is important to note that this package is not supported on Mac due to the lack of support from Intel MKL for Mac M1/M2 CPUs. 

3. You can install this package through our `drat` repository:

```r
# for Windows users
install.packages(c("oneMKL"), repos="https://R-OneMKL.github.io/drat", type="source")
install.packages(c("oneMKL.MatrixCal"), repos="https://R-OneMKL.github.io/drat")

# for Linux users
install.packages(c("oneMKL", "oneMKL.MatrixCal"), repos="https://R-OneMKL.github.io/drat")
```

Or, to get this package from github:

```r
# install.packages('remotes')
remotes::install_github("R-OneMKL/oneMKL") # users must first install the "oneMKL" package to construct the connection between R and oneMKL
remotes::install_github("R-OneMKL/oneMKL.MatrixCal")
```

### Notice

1. The function 'solve` in vanilla R will be broken for solving inverse matrices (large size) after loading MKL in the UNIX system. MKL uses INT64 ipiv in `dgesv`, but R uses INT32. Hence, it causes the issue. Therefore, we will replace `base::solve` with `fMatSolve` by employing `rlang` to avoid incorrect results.


### License

The `oneMKL` package is made available under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license.

The Intel oneMKL Library is licensed under the (Intel Simplified Software License)[https://www.intel.com/en-us/license/intel-simplified-software-license], as described at (oneMKL License FAQ)[https://www.intel.com/content/www/us/en/developer/articles/license/onemkl-license-faq.html].
