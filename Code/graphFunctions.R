# 5/25/2015 MDK
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

# Helper function for make P_bar
breakScalar <- function(N,d){
  n = floor(N/d);
  ret = c();
  for (i in 1:(d-1)){
    ret[i] = n;
  }
  ret[i+1] <- n+ N%%d;
  return(ret)
}

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

# Assigns vertices with exact membership probabilities
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

# Embedded vector covariance matrix calculation
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

# Creates matrix of m averaged nxn graphs generated from the bernoulli of P
# If Binomial TRUE draws from Binomial(P_ij,m).  If set to false, then m nxn Bernoulli matrices will be summed -> slow
# Parallelization implemented if Binomial = False and Serial = False
makeP_Bar <- function(P,m,n, hollow = TRUE, Binomial = TRUE, serial = FALSE){
  P_Bar <- matrix(0,n,n);
  if (Binomial = True) {
    P_vect <- as.vector(P)
    P_Bar <- rbinom(length(P_vect), m, P_vect);
    P_Bar <- matrix(P_Bar,n);
  }else if ((m*n < 2000)||serial == TRUE){
    for (h in 1:m) {
      g <- P2Adj(P, hollow = hollow);
      P_Bar <- P_Bar +g;
    }
  } else {
    P_Bar <- foreach(mb = breakScalar(m, 3), .export = c("P2Adj")) %dopar%{
      Ps_i = matrix(0,n,n);
      for (h in 1:mb) {
        g <- P2Adj(P, hollow = hollow);
        Ps_i <- Ps_i +g;
      }
      return(Ps_i)
    }
    P_Bar <- Reduce('+',P_Bar);
  }
  P_Bar <- P_Bar/m;
  if (hollow){
    diag(P_Bar) <- 0;
  }
  return(P_Bar)
}

# Implementation of Scheinerman iterarations for given iteration number
# Convergence cut-off not yet implemented
Schein <- function(M, rank, iter, myoptions){
  for (i in 1:iter){
    ASE <- spectEmbed(M, rank, DAdjust=FALSE, opt = myoptions)
    diag(M) <- diag(ASE$X %*% t(ASE$X));
  }
  return(ASE)
}

#Estimation of variance for zero mean random matrix defined by (P_hat-A_Bar)
# Output is approximation for sigma^2 = (lambda_{k+1}/n)^2, where lambda_{k+1} is the (k+1)th eigenvalue of A_Bar
calcSig= function(B, PI,m,n){
  k = dim(B)[1];
  s =0;
  for (i in 1:k){
    for (j in 1:k){
      sig = B[i,j]*(1-B[i,j]);
      s = s+(sig*PI[i]*PI[j]);
    }
  }
  return(s/m/n*4)
}

saveFig <- function(folder, filename){
  dev.copy(png, paste(folder, filename, "_", format(Sys.Date(), format ="%m%d%y"), ".png"))
  dev.off();
  dev.copy(pdf, paste(folder, filename, "_", format(Sys.Date(), format ="%m%d%y"), ".pdf"))
  dev.off();
}