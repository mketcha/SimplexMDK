# 2/9/2015: MDK
# functions to generate generalized SBM using Dirilecht Random Variable

#function: Generate dSBM
#Parameters: B  -> Block Matrix
#            PI -> Class Likelihood
#            n  -> Number of vertices
#            r  -> Dirilecht Parameter
#            opt -> ASE options
# Return: P-> Edge Probability Matrix
generateDirilechtSBM <- function(B, PI, n, r, opt = igraph.arpack.default) {
  k <- length(PI);
  nVec <- VertexClassCount(PI, n)
  BASE <- spectEmbed(B,k, opt)
  optimOpt  = list( outer.iter= 10000, inner.iter = 10000)
  #rot <- solnp(as.vector(diag(k)), fun = diriOpti, eqfun = diriEq, control = optimOpt, ineqfun= diriIneq, LB = as.vector(-matrix(1,k,k)), UB =as.vector(matrix(1,k,k)), ineqLB = as.vector(matrix(-Inf,k,k)), ineqUB =as.vector(matrix(0,k,k)), nVertex = k, dimLatentPosition = k, xHatTmp = as.vector(BASE$X))
  rot <- nloptr(as.vector(diag(k)), eval_f = diriOpti, lb = as.vector(-matrix(1,k,k)), ub = as.vector(matrix(1,k,k)), eval_g_ineq = diriIneq, eval_g_eq = diriEq, opts = list(algorithm = "NLOPT_GN_ISRES", maxeval = 10000), nVertex = k, dimLatentPosition = k, xHatTmp = as.vector(BASE$X));
  BASE$X <- BASE$X %*% matrix(rot$solution,k,k);
  xDir <- list();
  for (i in 1:k){
    tmp <- r*c(BASE$X[i,], 1-sum(BASE$X[i,])) + rep(1,k+1);
    xTmp <- drchrnd(tmp, nVec[i]);
    xDir[[i]] <- xTmp[,1:k];
  }
  xDir <- do.call(rbind,xDir);
  P <- xDir %*% t(xDir);
  diag(P) <-0;
  return(P)
}

#function: Generate ASE of B matrix
#Parameters: B  -> Block Matrix
#            opt -> ASE options
# Return: BASE-> ASE of B matrix
generateDirilechtBASE <- function(B, opt = igraph.arpack.default){
  k <- dim(B)[1];
  repeat{
    BASE <- spectEmbed(B,k, opt)
    optimOpt  = list( outer.iter= 10000, inner.iter = 10000)
    rot <- nloptr(as.vector(diag(k)), eval_f = diriOpti, lb = as.vector(-matrix(1,k,k)), ub = as.vector(matrix(1,k,k)), eval_g_ineq = diriIneq, eval_g_eq = diriEq, opts = list(algorithm = "NLOPT_GN_ISRES", maxeval = 10000), nVertex = k, dimLatentPosition = k, xHatTmp = as.vector(BASE$X));
    BASE$X <- BASE$X %*% matrix(rot$solution,k,k);
    
    d_graph <- generateDirilechtSBM_BASE( BASE, rep(1/k,k), 100, 500);
    if (sum(is.na(d_graph$xDir)) == 0) {
      break
    }
    print("err_BASE")
  }
  return(BASE)
}

#function: Generate dSBM given ASE of B matrix
#Parameters: BASE  -> ASE of Block Matrix
#            PI -> Class Likelihood
#            n  -> Number of vertices
#            r  -> Dirilecht Parameter
# Return: P   -> Edge Probability Matrix
#        nVec -> number of vertices in each class
generateDirilechtSBM_BASE <- function(BASE, PI, n, r) {
  k <- length(PI);
  nVec <- VertexClassCount(PI, n)
  xDir <- list();
  for (i in 1:k){
    tmp <- r*c(BASE$X[i,], 1-sum(BASE$X[i,])) + rep(1,k+1);
    xTmp <- drchrnd(tmp, nVec[i]);
    xDir[[i]] <- xTmp[,1:k];
  }
  xDir <- do.call(rbind,xDir);
  P <- xDir %*% t(xDir);
  diag(P) <-0;
  return(list(P = P, nVec = nVec, xDir = xDir))
}

# function: generate Dirilecht RV
# Parameters: a-> shape parameter
#             n -> number to be generated
drchrnd <- function(a,n){
  p <- length(a);
  d <- matrix(0,n,p);
  for (i in 1:length(a)){
    d[,i] <- rgamma(n,a[i]);
  }
  rs <-rowSums(d);
  n <- replicate(p,rs);
  return(d/n)
}

# Cost Function for restriction of latent space to first orthant
diriOpti <- function(x, nVertex, dimLatentPosition, xHatTmp){
  rotation <- matrix(x, dimLatentPosition, dimLatentPosition);
  xHat <- matrix(xHatTmp,dimLatentPosition,dimLatentPosition) %*% rotation;
  val = -min(abs(xHat));
  return(val)
}

# Inequality Function for restriction of latent space to first orthant
diriIneq <- function(x, nVertex, dimLatentPosition, xHatTmp){
  rotation <- matrix(x, dimLatentPosition, dimLatentPosition);
  xHat <- matrix(xHatTmp,dimLatentPosition,dimLatentPosition) %*% rotation;
  Ineq <- c();
  for (i in 1:nVertex){
    for (j in 1:dimLatentPosition){
      Ineq <- c(Ineq,-xHat[j,i]);
    }
  }
  return(Ineq)
}

# Equality Function for restriction of latent space to first orthant
diriEq <- function(x, nVertex, dimLatentPosition, xHatTmp){
  rotation <- matrix(x, dimLatentPosition, dimLatentPosition);
  xHat <- matrix(xHatTmp,dimLatentPosition,dimLatentPosition) %*% rotation;
  Eq <- c();
  for (i in 1:dimLatentPosition){
    for (j in 1:dimLatentPosition){
      if (i==j){
        Eq <- c(Eq, (rotation[j,, drop=FALSE] %*% t(rotation[i,,drop=FALSE])) -1);
      } else {
        Eq <- c(Eq, rotation[j,,drop=FALSE] %*% t(rotation[i,,drop=FALSE]));
      }
    }
  }
  return(Eq)
}