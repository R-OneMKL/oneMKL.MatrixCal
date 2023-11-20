library(oneMKL.MatrixCal)
library(dqrng)

hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }

testMatOps <- function() {
  dqset.seed(100)
  x <- matrix(dqrnorm(3e4), 300)
  z <- matrix(dqrnorm(3e4), 300)

  # fMatProd,fMatTransProd
  checkEquals(t(x) %*% x, fMatProd(t(x), x))
  checkEquals(t(x) %*% x, fMatTransProd(x, x))

  checkException(fMatProd(x, x))
  checkException(fMatTransProd(t(x), x))

  # test is_X_symmetric
  s <- hilbert(300)
  checkEquals(s %*% x, fMatProd(s, x, TRUE))

  b <- matrix(dqrnorm(300), 300)
  checkEquals(s %*% b, fMatProd(s, b, TRUE))

  # fMatDet
  XtX <- fMatTransProd(x, x)
  checkEquals(fMatDet(XtX), det(XtX))

  checkException(fMatDet(x))
}

testMatInverse <- function() {
  # check condition is allowed
  checkException(fMatInv(matrix(1:2, 1)))
  checkException(fMatSolve(matrix(1:2, 1), 1:3))
  checkException(fMatLeastSquare(matrix(1:2, 1), 1:3))

  # fMatInv 2*2
  z <- matrix(c(4, 3, 3, 5), 2)
  checkEquals(fMatInv(z), solve(z))
  checkEquals(fMatInv(z) %*% z, diag(1, 2))
  checkEquals(fMatInv(z, TRUE) %*% z, diag(1, 2))

  # fMatInv 9*9
  set.seed(100)
  z2 <- matrix(dqrnorm(81), 9)
  checkEquals(fMatInv(z2), solve(z2))
  checkEquals(fMatInv(z2) %*% z2, diag(1, 9))

  # test least square
  set.seed(100)
  x <- matrix(dqrnorm(3e4), 300, 100)
  beta <- matrix(dqrnorm(1e2), 100)
  y <- fMatProd(x, beta) + dqrnorm(300)
  checkEquals(solve(t(x) %*% x, t(x) %*% y), fMatSolve(fMatTransProd(x, x), fMatTransProd(x, y)))
  checkEquals(solve(t(x) %*% x, t(x) %*% y), fMatLeastSquare(x, y))
  checkEquals(solve(t(x) %*% x, t(x) %*% y), fMatLeastSquare(x, y, FALSE))

  # test if X is not full-rank
  if (require("pracma", quietly = TRUE)) {
    set.seed(100)
    x <- matrix(dqrnorm(3e3), 30, 100)
    beta <- matrix(dqrnorm(1e2), 100)
    y <- fMatProd(x, beta) + dqrnorm(30)
    checkEquals(pracma::pinv(t(x) %*% x) %*% t(x) %*% y, fMatLeastSquare(x, y, is_X_full_rank = FALSE))
  } else {
    cat("R package 'pracma' cannot be loaded -- the least square for X which is not full-rank will be skipped.\n")
  }
}
