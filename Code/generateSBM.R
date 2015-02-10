#1/30/2015: MDK
# functions to create model for an exact Stochastic Block Model

#Parameters: B  -> Block Matrix
#            PI -> Class Likelihood
#            n  -> Number of vertices
# Return: P-> Edge Probability Matrix
#         A -> Instance of an Adjacency Matrix
#         Z -> Vertex Membership
#         myGraph-> iGraph
#         k -> number of classess
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

generateSBM_Z <- function( B_kxk, Z_nxk, undirected = TRUE, hollow = TRUE) {
  P <- Z_nxk %*% B_kxk %*% t(Z_nxk);
  d <- dim(P)
  P_vect <- as.vector(P)
  bino <- rbinom(length(P_vect), 1, P_vect);
  bino <- matrix(bino,d);
  bino <- upper.tri(bino)*bino;
  A <- bino+t(bino);
  myGraph <- graph.adjacency( A, mode = "undirected", diag = hollow);
  mySBM <- list( P =P, B=B_kxk, A=A, Z=Z_nxk, myGraph=myGraph, k = ncol(B_kxk), n = nrow(Z_nxk));
  return(mySBM)
}