## oneMKL.MatrixCal

The `oneMKL.MatrixCal` package provides an efficient and powerful tool for
conducting linear algebra computations in the R environment. In the package,
matrix operations, such as decomposition, inversion, and multiplication,
are performed using the optimized algorithms embedded in the
**Intel oneAPI Math Kernel Library** (`oneMKL`). Because of this feature,
the computational performance can be significantly improved. To make the package
function work correctly with Intel `oneMKL`, users must first install the package,
`oneMKL`, to import Intel `oneMKL` into R. The package also supports multithreaded
computations, enabling users to utilize the benefit of multiple CPU cores,
i.e., parallel computing, to handle large-scale matrix computations.
The package is not only compatible with Linux, but aslo Windows operating systems,
which can benefit a wide range of users.


### Installation

1. To build the package from source, Windows users will need to have [Rtools](http://cran.csie.ntu.edu.tw/bin/windows/Rtools/) installed.

2. You can install this package through our `drat` repository:

```r
install.packages(c("oneMKL", "oneMKL.MatrixCal"), repos=c("https://cloud.r-project.org/", "https://R-OneMKL.github.io/drat"), type = "source")
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
