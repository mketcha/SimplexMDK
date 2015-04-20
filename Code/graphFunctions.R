# 2/6/2015 MDK
# Generic Functions for Dealing with Graphs
require("igraph")
require("rARPACK")
require("doParallel")
require("foreach")
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

breakScalar <- function(N,d){
  n = round(N/d);
  ret = c();
  for (i in 1:(d-1)){
    ret[i] = n;
  }
  ret[i+1] <- n+ N%%d;
  return(ret)
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
      isn <- 2*((-1*(sp$values<0))+.5);
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
      isn <- 2*((-1*(sp$values<0))+.5);
      X = sp$vectors %*% diag(sqrt(abs(sp$values)));
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

generateZ_nxk_exact <- function(PI_1xk, n){
  VCC <- PI_1xk*n;
  Z <- matrix(0, n, length(VCC));
  count = 1;
  for (i in 1:length(VCC)){
    Z[count:(count+VCC[i]-1),i] <-1;
    count <- count+VCC[i]
  }
  return(Z);
}

calcCov <- function(X, x){
  d <- dim(X)[2]
  n <- dim(X)[1]
  Delta = (t(X) %*% X)/n;
  DInv = ginv(Delta);
  inside = matrix(0,d,d);
  for (i in 1:n){
    inside = inside + (t(X[i,,drop = FALSE]) %*% X[i,,drop = FALSE])*as.vector((x %*% t(X[i,,drop = FALSE]))*(1-(x %*% t(X[i,,drop = FALSE]))));
  }
  inside  = inside/n;
  Cov = DInv %*% inside %*% DInv;
  return(Cov)
}

makeP_Bar <- function(P,m,n){
P_Bar <- matrix(0,n,n);
if (m*n < 2000){
  for (h in 1:m) {
    g <- P2Adj(P);
    P_Bar <- P_Bar +g;
  }
} else {
  P_Bar <- foreach(mb = breakScalar(m, 3), .export = c("P2Adj")) %dopar%{
    Ps_i = matrix(0,n,n);
    for (h in 1:mb) {
      g <- P2Adj(P);
      Ps_i <- Ps_i +g;
    }
    return(Ps_i)
  }
  P_Bar <- Reduce('+',P_Bar);
}
P_Bar <- P_Bar/m;
return(P_Bar)
}

Schein <- function(M, rank, iter, myoptions){
  for (i in 1:iter){
    ASE <- spectEmbed(M, rank, DAdjust=FALSE, opt = myoptions)
    diag(M) <- diag(ASE$X %*% t(ASE$X));
  }
  return(ASE)
}