## oneMKLUtil

The `oneMKLUtil` package provide some utiliy functions based on `oneMKL` package.

### Installation

To build the package from source, Windows users will need to have [Rtools](http://cran.csie.ntu.edu.tw/bin/windows/Rtools/) installed.

Because the availability of Anaconda `mkl` package, we only support Windows and Linux.
Also, note that this package does not support Mac because Intel MKL does not support Mac M1/M2 CPUs.

To get this package from github:

```r
install.packages('remotes')
remotes::install_github("R-OneMKL/oneMKLUtil")
```

### License

The `oneMKL` package is made available under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license.

The Intel MKL Library is licensed under the (Intel Simplified Software License)[https://www.intel.com/en-us/license/intel-simplified-software-license], as described at (Intel MKL License FAQ)[https://www.intel.com/content/www/us/en/developer/articles/license/onemkl-license-faq.html].
