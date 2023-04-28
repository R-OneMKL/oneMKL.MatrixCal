library(oneMKL.MatrixCal)
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }

testMatChol <- function() {
  m <- matrix(c(5,1,1,3),2,2)
  checkEquals(fMatChol(m), chol(m))
}

testMatSvd <- function() {
  X <- hilbert(9)[, 1:6]
  s <- fMatSvd(X)
  D <- diag(as.vector(s$d))
  checkEquals(s$u[ , 1:6] %*% D %*% t(s$v), X)
  checkEquals(t(s$u[ , 1:6]) %*% X %*% s$v, D)
}

testMatEigen <- function() {
  X <- cbind(c(1,-1), c(-1,1))
  eigenResVanillaR <- eigen(X)
  eigenRes <- fMatEigen(X)
  eigenVals <- Re(eigenRes$values)
  Z1 <- abs(Re(eigenRes$vectors) %*% eigenResVanillaR$vectors)
  checkEquals(
    as.vector(eigenVals[order(eigenVals)]),
    eigenResVanillaR$values[order(eigenResVanillaR$values)]
  )
  checkTrue(all.equal(matrix(1, 2, 2) - diag(1, 2, 2), Z1) || all.equal(Z1, diag(1, 2, 2)))

  X <- hilbert(16)
  eigenRes <- fMatEigen(X)
  checkTrue(all(abs(Re(eigenRes$vectors %*% diag(eigenRes$values) %*% solve(eigenRes$vectors)) - X) < 1e-6))
}

testOtherDecom <- function() {
  # fMatLu
  X <- matrix(rnorm(9), 3, 3)
  luRes <- fMatLu(X)
  checkEquals(solve(luRes$P) %*% luRes$L %*% luRes$U, X)

  # fMatQr
  qrRes <- fMatQr(X)
  checkEquals(qrRes$Q %*% qrRes$R, X)
}
