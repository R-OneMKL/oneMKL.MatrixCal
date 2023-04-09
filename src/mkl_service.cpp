// Copyright (C) 2022             Ching-Chuan Chen
//
// This file is part of oneMKL.
//
// oneMKL is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// oneMKL is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with oneMKL. If not, see <http://www.gnu.org/licenses/>.

// [[Rcpp::depends(oneMKL)]]

#include <oneMKL.h>
#include <string>

//' Function to get the version of Intel MKL
//'
//' @return The version of Intel MKL
//' @examples
//' getMKLVersion()
//' @export
// [[Rcpp::export]]
std::string getMKLVersion() {
  int len=198;
  char buf[198];
  mkl_get_version_string(buf, len);
  std::string mklVersionString(buf);
  return(mklVersionString);
}

//' Function to get/set the number of threads used in Intel MKL
//'
//' @param nThreads The number of threads you want to use in Intel MKL.
//' @return The number of threads.
//'
//' @examples
//' getMKLThreads() # Default is the number of CPUs cores on your PC
//'
//' \dontrun{
//' setMKLThreads(1)
//' getMKLThreads() # 1
//' }
//' @name mkl_threads
//' @export
// [[Rcpp::export]]
int setMKLThreads(int nThreads) {
  mkl_set_num_threads(nThreads);
  return nThreads;
}

//' @name mkl_threads
//' @export
// [[Rcpp::export]]
int getMKLThreads() {
  return mkl_get_max_threads();
}

