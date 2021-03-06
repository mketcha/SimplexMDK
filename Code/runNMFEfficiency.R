setwd("~/Simplex/SimplexMDK/Code")
library("igraph")
library("NMF")
library("vioplot")
library("rARPACK")
library("mclust")
library("nloptr")
library("fields")
source("graphFunctions.R")
source("generateDirilechtSBM.R");

MN_M = matrix(c(100,5,100,10, 100, 25, 100,50, 100,75,100,100,100,250, 100,500),8,2,byrow=TRUE)
MN_N = matrix(c(1,100,10,100,50,100,100, 100,500,100,1000,100,2000,100, 5000,100),8,2,byrow=TRUE)

B <- matrix(c(.42, .2, .3, .2,.5, .15, .3, .15, .7),3);
#B <- matrix(c(.42, .2, .2, .7),2);
PI <- c(.33, .33,.34)
#PI <- c(.5,.5);
trials <- 1000;
max_k <- 4;
nclust <- 3;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
eff_M = c();
eff_N = c();


for (run in 1:8) {
  n <- MN_M[run,2];
  m <- MN_M[run,1];
  error_Phat_noClust <- matrix(0, trials, max_k);
  error_Phat <- matrix(0, trials, max_k);
  error_mean <- matrix(0, trials, 1);
  for (i in 1:trials) {
    if (i%%100 ==0)
      print(i)
    Z <-generateZ_nxk(PI,n)
    P = Z %*% B %*% t(Z);
    P_Bar <- matrix(0,n,n);
    for (h in 1:m) {
      g <- P2Adj(P);
      P_Bar <- P_Bar +g;
    }
    P_Bar <- P_Bar/m;
    err_v <- P-P_Bar;
    err_v <- err_v[upper.tri(err_v)];
    error_mean[i,1] <- mean(err_v^2);
    j = 3;
    diag(P_Bar) <- rowSums(P_Bar)/(dim(P_Bar)[1]-1);
    ASE <- nmf(P_Bar,j, "brunet")
    P_hat_noClust <- basis(ASE) %*% coef(ASE);
    
    err_v_noClust <- P-P_hat_noClust;
    err_v_noClust <- err_v_noClust[upper.tri(err_v_noClust)];
    error_Phat_noClust[i,j] <- mean(err_v_noClust^2);
    
    
    #}
  }
  print(mean(error_mean));
  print(mean(error_Phat_noClust[,j]))
  print(mean(error_Phat[,j]))
  eff_M = cbind(eff_M,mean(error_Phat_noClust[,j])/mean(error_mean))
}

for (run in 1:8) {
  n <- MN_N[run,2];
  m <- MN_N[run,1];
  error_Phat_noClust <- matrix(0, trials, max_k);
  error_Phat <- matrix(0, trials, max_k);
  error_mean <- matrix(0, trials, 1);
  for (i in 1:trials) {
    if (i%%100 ==0)
      print(i)
    Z <-generateZ_nxk(PI,n)
    P = Z %*% B %*% t(Z);
    P_Bar <- matrix(0,n,n);
    for (h in 1:m) {
      g <- P2Adj(P);
      P_Bar <- P_Bar +g;
    }
    P_Bar <- P_Bar/m;
    err_v <- P-P_Bar;
    err_v <- err_v[upper.tri(err_v)];
    error_mean[i,1] <- mean(err_v^2);
    j = 3;
    diag(P_Bar) <- rowSums(P_Bar)/(dim(P_Bar)[1]-1);
    ASE <- nmf(P_Bar,j, "brunet")
    P_hat_noClust <- basis(ASE) %*% coef(ASE);
    
    err_v_noClust <- P-P_hat_noClust;
    err_v_noClust <- err_v_noClust[upper.tri(err_v_noClust)];
    error_Phat_noClust[i,j] <- mean(err_v_noClust^2);
    
  }
  print(mean(error_mean));
  print(mean(error_Phat_noClust[,j]))
  print(mean(error_Phat[,j]))
  eff_N = cbind(eff_N,mean(error_Phat_noClust[,j])/mean(error_mean))
}
