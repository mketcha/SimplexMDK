# 2/6/2015 MDK
# Generic Functions for Dealing with Graphs

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

# P2Adj <- function( P, undirected = TRUE, hollow = TRUE) {
#   d = dim(P);
#   if (!undirected){
#     P_vect <- as.vector(P)
#     bino <- rbinom(length(P_vect), 1, P_vect);
#     bino <- matrix(bino,d);
#     if(hollow){
#       diag(bino) <- 0;
#     }
#     return(bino)
#   }
#   #bino <- upper.tri(bino, TRUE)*bino;
#   P_vect <- P[upper.tri(P,TRUE)];
#   bino <- rbinom(length(P_vect), 1, P_vect);
#   A <- matrix(0,d,d);
#   A[upper.tri(A,TRUE)] <- bino;
#   A <- A+t(A);
#   #A <- symMatrix(bino,d[1], byrow=TRUE)
#   if (hollow){
#     diag(A) <- 0;
#   } else {
#     diag(A) <- diag(A)/2;
#   }
#   return(A)
# }

# ASE of a symmetric matrix, given number of eigenvalues to calculate
spectEmbed <- function(M, no, DAdjust = FALSE, opt = igraph.arpack.default) {
  if (DAdjust){
    diag(M) <- rowSums(M)/(dim(M)[1]-1);
  }
  if (length(M[,1])==2){
    sp <- eigen(M);
    if (no == 1){
      X = sp$vectors[,1]*sqrt(abs(sp$values[1]))
      lambda = sp$values[1];
    } else if (no==2){
      X = sp$vectors %*% diag(sqrt(abs(sp$values)));
      lambda = sp$values
    } else {
      stop
    }
    ret = list(X=X, lamda = lambda)
    return(ret)
    
  } else{
    sp = eigs(M, no, options = opt, which = "LM")
    if (no == 1){
      X = sp$vectors*sqrt(abs(sp$values))
    } else {
      isn <- -1*(sp$values<0);
      X = sp$vectors %*% diag(isn*sqrt(abs(sp$values)));
    }
    ret = list(X=X, lamda = sp$values)
    return(ret)
  }
}

# generate random class counts given probabilities and number of vertices
VertexClassCount <- function(PI_1xk, n){
  r <- runif(n);
  Z_k <- rep(0,length(PI_1xk))
  PI_c = c(0,cumsum(PI_1xk));
  for (i in 1:(length(PI_c)-1)){
    Z_k[i] <- sum(r>=PI_c[i] & r<PI_c[i+1]);
    if (Z_k[i] ==0)
      return(VertexClassCount(PI_1xk,n))
  }
  return(Z_k)
}

# generate random class vertex assignments given probabilities and number of vertices
generateZ_nxk <- function(PI_1xk, n){
#   r <- runif(n);
#   Z_n <- rep(0,n)
#   PI_c = c(0,cumsum(PI_1xk));
#   for (i in 1:(length(PI_c)-1)){
#     Z_n[r>=PI_c[i] & r<PI_c[i+1]] = i;
#   }
#   Z_nxk <- matrix(0, n, length(PI_1xk));
#   for (i in 1:length(Z_n)) {
#     Z_nxk[i,Z_n[i]] = 1;
#   }
#   return(Z_nxk)
  VCC <- VertexClassCount(PI_1xk, n);
  Z <- matrix(0, n, length(VCC));
  count = 1;
  for (i in 1:length(VCC)){
    Z[count:(count+VCC[i]-1),i] <-1;
    count <- count+VCC[i]
  }
  return(Z);
}