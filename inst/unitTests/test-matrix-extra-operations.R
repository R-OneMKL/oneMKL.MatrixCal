library(oneMKL.MatrixCal)
library(dqrng)

testMatRowWiseOps <- function() {
  dqset.seed(100)
  x <- matrix(dqrnorm(3e4), 300)
  z <- matrix(dqrnorm(3e4), 300)

  # fMatAdd, fMatSubtract
  checkEquals(fMatAdd(x, z), x + z)
  checkEquals(fMatSubtract(x, z), x - z)

  checkException(fMatAdd(x, b))
  checkException(fMatSubtract(x, b))
}

testMatRowWiseOps <- function() {
  dqset.seed(100)
  x1 <- matrix(dqrnorm(3e4), 300)
  x2 <- matrix(as.integer(round(dqrnorm(3e4))), 300)

  # fMatRowSum, fMatRowMin, fMatRowMax
  checkTrue(all.equal(fMatRowSum(x1), rowSums(x1)))
  checkTrue(all.equal(fMatRowMin(x1), apply(x1, 1, min)))
  checkTrue(all.equal(fMatRowMax(x1), apply(x1, 1, max)))

  checkTrue(all.equal(fMatRowSum(x2), rowSums(x2)))
  checkTrue(all.equal(fMatRowMin(x2), apply(x2, 1, min)))
  checkTrue(all.equal(fMatRowMax(x2), apply(x2, 1, max)))
  checkEquals(typeof(fMatRowSum(x2)), "integer")
  checkEquals(typeof(fMatRowMin(x2)), "integer")
  checkEquals(typeof(fMatRowMax(x2)), "integer")
}

testMatColWiseOps <- function() {
  dqset.seed(100)
  x1 <- matrix(dqrnorm(3e4), 300)
  x2 <- matrix(as.integer(round(dqrnorm(3e4))), 300)

  # fMatColSum, fMatColMin, fMatColMax
  checkTrue(all.equal(as.vector(fMatColSum(x1)), colSums(x1)))
  checkTrue(all.equal(as.vector(fMatColMin(x1)), apply(x1, 2, min)))
  checkTrue(all.equal(as.vector(fMatColMax(x1)), apply(x1, 2, max)))

  checkTrue(all.equal(as.vector(fMatColSum(x2)), colSums(x2)))
  checkTrue(all.equal(as.vector(fMatColMin(x2)), apply(x2, 2, min)))
  checkTrue(all.equal(as.vector(fMatColMax(x2)), apply(x2, 2, max)))
  checkEquals(typeof(fMatColSum(x2)), "integer")
  checkEquals(typeof(fMatColMin(x2)), "integer")
  checkEquals(typeof(fMatColMax(x2)), "integer")
}

testMatElementWiseOps <- function() {
  dqset.seed(100)
  x1 <- matrix(dqrnorm(3e4), 300)
  x2 <- matrix(as.integer(round(dqrnorm(3e4))), 300)
  x3 <- matrix(dqrnorm(3e4), 300)

  # fMatElementWiseProduct
  checkTrue(all.equal(fMatElementWiseProduct(x1, x1), x1*x1))
  checkTrue(all.equal(fMatElementWiseProduct(x2, x2), x2*x2))
  checkEquals(typeof(fMatElementWiseProduct(x2, x2)), "integer")

  # fMatElementWiseDivide
  checkTrue(all.equal(fMatElementWiseDivide(x1, x3), x1/x3))

  # check condition is validated
  checkException(fMatElementWiseProduct(x1, x1[1:100, ]))
  checkException(fMatElementWiseDivide(x1, x1[1:100, ]))
}

testMatSpecialOps <- function() {
  dqset.seed(100)
  x1 <- matrix(dqrnorm(3e4), 300)
  x2 <- matrix(dqrnorm(3e4), 300)

  checkTrue(all.equal(fMatSumDiffSquared(x1, x2), sum((x1-x2)^2)))
  checkException(fMatSumDiffSquared(x1, x2[1:100, ]))
}
