## oneMKL.MatrixCal

The `oneMKL.MatrixCal` package provides an efficient and powerful tool for 
conducting linear algebra computations in the R environment. It achieves this by 
utilizing the optimized algorithms embedded in the **Intel oneAPI Math Kernel Library** 
(`oneMKL`) to accelerate matrix operations, such as factorization, inversion, and 
multiplication. To ensure proper installation, users must first install the `oneMKL` 
package to establish the connection between R and Intel `oneMKL`. The package also supports 
multithreaded computations, enabling users to utilize multiple CPU cores for a faster 
performance on large-scale matrix computations. Its compatibility with Windows and 
Linux operating systems makes it a versatile option for many users.

### Installation

1. To build the package from source, Windows users will need to have [Rtools](http://cran.csie.ntu.edu.tw/bin/windows/Rtools/) installed.

2. You can install this package through our `drat` repository:

```r
# for Windows users (Because we don't provide the binary package of oneMKL)
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


### License

The `oneMKL` package is made available under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license.

The Intel oneMKL Library is licensed under the (Intel Simplified Software License)[https://www.intel.com/en-us/license/intel-simplified-software-license], as described at (oneMKL License FAQ)[https://www.intel.com/content/www/us/en/developer/articles/license/onemkl-license-faq.html].
