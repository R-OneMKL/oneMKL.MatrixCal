## oneMKLUtil

The `oneMKLUtil` package provide some utiliy functions based on `oneMKL` package.

### Installation

To build the package from source, Windows users will need to have [Rtools](http://cran.csie.ntu.edu.tw/bin/windows/Rtools/) installed.

Because the availability of Anaconda `mkl` package, we only support Windows and Linux.
Also, note that this package does not support Mac because Intel MKL does not support Mac M1/M2 CPUs.

You can install this package through our `drat` repository:

```r
# for windows 
install.packages(c("oneMKL"), repos="https://R-OneMKL.github.io/drat", type="source")
install.packages(c("oneMKLUtil"), repos="https://R-OneMKL.github.io/drat")

# for Linux
install.packages(c("oneMKL", "oneMKLUtil"), repos="https://R-OneMKL.github.io/drat")
```

Or, to get this package from github:

```r
# install.packages('remotes')
remotes::install_github("R-OneMKL/oneMKL") # must install oneMKL first
remotes::install_github("R-OneMKL/oneMKLUtil")
```

### Notice

1. `solve` in vanilla R will be broken for solving inverse matrices (large size) after loading MKL in UNIX system. MKL uses INT64 ipiv in `dgesv`, but R uses INT32. Hence, it causes the issue. Therefore, we will replace `base::solve` with `fMatSolve` by employing `rlang` to avoid the incorrect results.


### License

The `oneMKL` package is made available under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license.

The Intel MKL Library is licensed under the (Intel Simplified Software License)[https://www.intel.com/en-us/license/intel-simplified-software-license], as described at (Intel MKL License FAQ)[https://www.intel.com/content/www/us/en/developer/articles/license/onemkl-license-faq.html].
