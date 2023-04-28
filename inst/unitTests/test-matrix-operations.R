library(oneMKL.MatrixCal)
testMatOps <- function() {
  x <- matrix(rnorm(3e4), 300)
  z <- matrix(rnorm(3e4), 300)

  # fMatProd,fMatTransProd
  checkEquals(t(x) %*% x, fMatProd(t(x), x))
  checkEquals(t(x) %*% x, fMatTransProd(x, x))

  # fMatAdd, fMatSubtract
  checkEquals(fMatAdd(x, z), x + z)
  checkEquals(fMatSubtract(x, z), x - z)

  # fMatDet
  XtX <- fMatTransProd(x, x)
  checkEquals(fMatDet(XtX), det(XtX))
}

testMatInverse <- function() {
  x <- matrix(rnorm(3e4), 300, 100)
  y <- fMatTransProd(x, matrix(rnorm(3e2), 300)) + rnorm(100)
  XtX <- fMatTransProd(x, x)

  # fMatSolve (linear model)
  checkEquals(solve(t(x) %*% x, y), fMatSolve(XtX, y))

  # fMatInv 2*2
  z <- matrix(c(4, 3, 3, 5), 2)
  checkEquals(fMatInv(z), solve(z))
  checkEquals(fMatInv(z) %*% z, diag(1, 2))
  checkEquals(fMatInv(z, TRUE) %*% z, diag(1, 2))
}
