library(oneMKL.MatrixCal)
testMatChol <- function() {
  m <- matrix(c(5,1,1,3),2,2)
  checkEquals(fMatChol(m), chol(m))
  checkEquals(fMatChol(m, FALSE), t(chol(m)))
}

testMatSvd <- function() {
  hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
  X <- hilbert(9)[, 1:6]
  s <- fMatSvd(X)
  D <- diag(as.vector(s$d))
  checkEquals(s$u[ , 1:6] %*% D %*% t(s$v), X)
  checkEquals(t(s$u[ , 1:6]) %*% X %*% s$v, D)
}

checkMatEigen <- function() {
  X <- cbind(c(1,-1), c(-1,1))
  eigenResVanillaR <- eigen(X)
  eigenRes1 <- fMatEigen(X, TRUE)
  Z1 <- abs(eigenRes1$vectors %*% eigenResVanillaR$vectors)
  checkEquals(
    as.vector(eigenRes1$values[order(eigenRes1$values)]),
    eigenResVanillaR$values[order(eigenResVanillaR$values)]
  )
  checkTrue(all.equal(matrix(1, 2, 2) - diag(1, 2, 2), Z1) || all.equal(Z1, diag(1, 2, 2)))

  eigenRes2 <- fMatEigen(X, FALSE)
  Z2 <- abs(Re(eigenRes2$vectors) %*% eigenResVanillaR$vectors)

  stopifnot(all(Im(eigenRes2$values) == 0))
  stopifnot(all(Im(eigenRes2$vectors) == 0))
  checkEquals(
    Re(eigenRes2$values)[order(Re(eigenRes2$values))],
    eigenResVanillaR$values[order(eigenResVanillaR$values)]
  )
  checkTrue(all.equal(matrix(1, 2, 2) - diag(1, 2, 2), Z2) || all.equal(Z2, diag(1, 2, 2)))
}

checkOtherDecom <- function() {
  # fMatLu
  X <- matrix(rnorm(9), 3, 3)
  luRes <- fMatLu(X, TRUE)
  checkEquals(t(luRes$P) %*% luRes$L %*% luRes$U, X)

  # fMatSchur
  schurRes <- fMatSchur(X)
  checkEquals(schurRes$U %*% schurRes$S %*% t(schurRes$U), X)

  # fMatQr
  qrRes <- fMatQr(X)
  checkEquals(qrRes$Q %*% qrRes$R, X)
}
