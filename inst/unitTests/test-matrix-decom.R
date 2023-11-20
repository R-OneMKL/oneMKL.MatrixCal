library(oneMKL.MatrixCal)
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }

testMatChol <- function() {
  X <- hilbert(9)
  checkEquals(fMatChol(X), chol(X))
  checkException(fMatChol(x[ , 1:2]))
}

testMatQr <- function() {
  X <- hilbert(16)
  qrRes <- fMatQR(X)
  checkEquals(qrRes$Q %*% qrRes$R, X)

  qrRes2 <- fMatQR(X, TRUE)
  checkEquals(qrRes2$Q %*% qrRes2$R %*% fMatInv(qrRes2$P), X)
}

testMatSvd <- function() {
  X <- hilbert(9)[, 1:6]
  s <- fMatSVD(X)
  D <- diag(as.vector(s$d))
  checkEquals(s$U[ , 1:6] %*% D %*% t(s$V), X)
  checkEquals(t(s$U[ , 1:6]) %*% X %*% s$V, D)
}

testMatEigen <- function() {
  X <- cbind(c(1, -1), c(-1, 1))
  eigenResVanillaR <- eigen(X)
  eigenRes <- fMatEigen(X, TRUE)
  Z1 <- abs(eigenRes$vectors %*% eigenResVanillaR$vectors)
  checkEquals(
    as.vector(eigenRes$values[order(eigenRes$values)]),
    eigenResVanillaR$values[order(eigenResVanillaR$values)]
  )
  checkTrue(all.equal(matrix(1, 2, 2) - diag(1, 2, 2), Z1) || all.equal(Z1, diag(1, 2, 2)))

  X <- hilbert(16)
  eigenRes1 <- fMatEigen(X)
  checkTrue(all(abs(Re(eigenRes1$vectors %*% diag(eigenRes1$values) %*% solve(eigenRes1$vectors)) - X) < 1e-6))
  checkException(fMatChol(x[ , 1:2]))

  eigenRes2 <- fMatEigen(X, TRUE)
  checkTrue(all(abs(eigenRes2$vectors %*% diag(eigenRes2$values) %*% fMatInv(eigenRes2$vectors) - X) < 1e-6))
}

testMatLu <- function() {
  X <- hilbert(16)
  luRes <- fMatLU(X)
  checkEquals(fMatInv(luRes$P) %*% luRes$L %*% luRes$U, X)
}
