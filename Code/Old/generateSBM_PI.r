generateSBM_PI <- function( B_kxk, PI_1xk, n, undirected = TRUE, hollow = TRUE) {
  Z_nxk = generateZ_nxk(PI_1xk, n);
  P <- Z_nxk %*% B_kxk %*% t(Z_nxk);
  d <- dim(P)
  P_vect <- as.vector(P)
  bino <- rbinom(length(P_vect), 1, P_vect);
  bino <- matrix(bino,d);
  bino <- upper.tri(bino)*bino;
  A <- bino+t(bino);
  myGraph <- graph.adjacency( A, mode = "undirected", diag = hollow);
  mySBM <- list( P =P, B=B_kxk, A=A, Z=Z_nxk, PI = PI_1xk, myGraph=myGraph, k = ncol(B_kxk), n = n);
  return(mySBM)
}