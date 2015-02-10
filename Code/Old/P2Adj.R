# 2/6/2015 MDK
# Construct an adjacency matrix, given an edge Probabilty matrix
P2Adj <- function( P, undirected = TRUE, hollow = TRUE) {
  d = dim(P);
  P_vect <- as.vector(P)
  bino <- rbinom(length(P_vect), 1, P_vect);
  bino <- matrix(bino,d);
  if (!undirected){
    if(hollow){
      diag(bino) <- 0;
    }
    return(bino)
  }
  bino <- upper.tri(bino, TRUE)*bino;
  A <- bino+t(bino);
  if (hollow){
    diag(A) <- 0;
  } else {
    diag(A) <- diag(A)/2;
  }
  return(A)
}