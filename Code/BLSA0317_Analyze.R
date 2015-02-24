setwd("~/Simplex/SimplexMDK/Code")
library("igraph")
library("vioplot")
library("rARPACK")
library("mclust")
library("nloptr")
library("fields")
library("cvTools")
library("R.matlab")
source("graphFunctions.R")
source("generateDirilechtSBM.R");

Mat <- readMat("../Data/BLSA0317.mat");
A_all <- Mat$Abinary;
n <- 70;
M <- 49;
nTrials= 500;
groups = c(30, 15, 5, 4, 3, 2, 1)
max_k = 4;
nclust <- 2;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
attach(mtcars)
par(mfrow=c(4,2));
error_Phat_noClust <- matrix(0, M, max_k);
error_Phat <- matrix(0, M, max_k);
error_mean <- matrix(0, M, 1);

for (i in 1:M){
  A_bar <- matrix(0,n,n);
  
  for (j in 1:M) {
    if (i!=j)
      A_bar <- A_bar + A_all[,,j] + t(A_all[,,j]);
  }
  A_i <- A_all[,,i] + t(A_all[,,i]);
  A_bar <- A_bar/(M-1)
  
#   err_v <- A_i-A_bar;
#   err_v <- err_v[upper.tri(err_v)];
#   error_mean[i,1] <- mean(err_v^2);
  
  error_mean[i,1] <- norm(A_bar-A_i, type = "F");
  for (j in 1:max_k) {
    repeat{
      ASE <- try(spectEmbed(A_bar,j, myoptions))
      if (class(ASE) != "try-error") {
        break
      }
      print("err")
    }
    cl <- try(Mclust(ASE$X,nclust))
    if (class(cl) == "try-error") {
      i <- i-1;
      break
    }
    
    clust = matrix(0,n,j)
    if (is.vector(cl$parameters$mean))
      cl$parameters$mean <- matrix(cl$parameters$mean, ncol = nclust)
    for (cent in 1:nclust){
      if (sum(cl$classification==cent)!=0)
        clust[cl$classification==cent,] <- matrix(rep(t(cl$parameters$mean[,cent, drop = FALSE]),sum(cl$classification==cent)), nrow = sum(cl$classification==cent), byrow=TRUE);
    }
    P_hat_noClust <- ASE$X %*% t(ASE$X);
    diag(P_hat_noClust) <-0;
#     err_v_noClust <- A_i-P_hat_noClust;
#     err_v_noClust <- err_v_noClust[upper.tri(err_v_noClust)];
#     error_Phat_noClust[i,j] <- mean(err_v_noClust^2);
    
    error_Phat_noClust[i,j] <- norm(P_hat_noClust-A_i, type = "F");
    
    P_hat <- clust %*% t(clust)
    diag(P_hat) <-0;
#     err_v <- A_i-P_hat;
#     err_v <- err_v[upper.tri(err_v)];
#     error_Phat[i,j] <- mean(err_v^2);
    
    error_Phat[i,j] <- norm(P_hat-A_i, type = "F");
    
  }
}

  vioplot( error_mean, error_Phat_noClust[,1], error_Phat[,1], error_Phat_noClust[,2], error_Phat[,2], error_Phat_noClust[,3], error_Phat[,3], error_Phat_noClust[,4], error_Phat[,4], names = c("Naive", "k=1", "k=1", "k=2", "k=2", "k=3", "k=3","k=4", "k=4"), col =c("blue"));
  vioplot(error_mean, at = 1, col = c("magenta"), add= TRUE)
  vioplot(error_Phat_noClust[,1], at = 2, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,2], at = 4, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,3], at = 6, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,4], at = 8, col = c("green"), add= TRUE)
  #vioplot(error_Phat_noClust[,5], at = 10, col = c("green"), add= TRUE)
  title(main=paste("Frobenius Norm of Error:  M = ", M-1, ", N = ", n,", clusters = ", nclust))
  #legend(locator(1), legend= c("Naive", "X*X'", "Up to 10 clusters"), fill = c("magenta","green", "blue"))

for (run in 1:length(groups)){
  Ncv <- groups[run];
  error_Phat_noClust <- matrix(0, nTrials*(M-Ncv), max_k);
  error_Phat <- matrix(0, nTrials*(M-Ncv), max_k);
  error_mean <- matrix(0, nTrials*(M-Ncv), 1);
  for (i in 1:nTrials) {
    randV <- runif(49);
    minval <- sort(randV,partial=M-Ncv+1)[M-Ncv+1];
    ingrp <- (randV >= minval);
    A_bar <- matrix(0,n,n);
    
    for (j in 1:M) {
      if (ingrp[j])
        A_bar <- A_bar + A_all[,,j] + t(A_all[,,j]);
    }
    
    A_i <- A_all[,,!ingrp]
    A_bar <- A_bar/(sum(ingrp))
    
    #   err_v <- A_i-A_bar;
    #   err_v <- err_v[upper.tri(err_v)];
    #   error_mean[i,1] <- mean(err_v^2);
    
    for (h in 1:(M-Ncv)){
      error_mean[(i-1)*(M-Ncv)+h,1] <- norm(A_bar-(A_i[,,h]+t(A_i[,,h])), type = "F");
    }
    for (j in 1:max_k) {
      repeat{
        ASE <- try(spectEmbed(A_bar,j, myoptions))
        if (class(ASE) != "try-error") {
          break
        }
        print("err")
      }
      cl <- try(Mclust(ASE$X,nclust))
      if (class(cl) == "try-error") {
        i <- i-1;
        break
      }
      
      clust = matrix(0,n,j)
      if (is.vector(cl$parameters$mean))
        cl$parameters$mean <- matrix(cl$parameters$mean, ncol = nclust)
      for (cent in 1:nclust){
        if (sum(cl$classification==cent)!=0)
          clust[cl$classification==cent,] <- matrix(rep(t(cl$parameters$mean[,cent, drop = FALSE]),sum(cl$classification==cent)), nrow = sum(cl$classification==cent), byrow=TRUE);
      }
      P_hat_noClust <- ASE$X %*% t(ASE$X);
      diag(P_hat_noClust) <-0;
      #     err_v_noClust <- A_i-P_hat_noClust;
      #     err_v_noClust <- err_v_noClust[upper.tri(err_v_noClust)];
      #     error_Phat_noClust[i,j] <- mean(err_v_noClust^2);
      
      for (h in 1:(M-Ncv)){
        error_Phat_noClust[(i-1)*(M-Ncv)+h,j] <- norm(P_hat_noClust-(A_i[,,h]+t(A_i[,,h])), type = "F");
      }
      
      P_hat <- clust %*% t(clust)
      diag(P_hat) <-0;
      #     err_v <- A_i-P_hat;
      #     err_v <- err_v[upper.tri(err_v)];
      #     error_Phat[i,j] <- mean(err_v^2);
      for (h in 1:(M-Ncv)){
        error_Phat[(i-1)*(M-Ncv)+h,j] <- norm(P_hat-(A_i[,,h]+t(A_i[,,h])), type = "F");
      }
      
    }
  }
  
  vioplot( error_mean, error_Phat_noClust[,1], error_Phat[,1], error_Phat_noClust[,2], error_Phat[,2], error_Phat_noClust[,3], error_Phat[,3], error_Phat_noClust[,4], error_Phat[,4], names = c("Naive", "k=1", "k=1", "k=2", "k=2", "k=3", "k=3","k=4", "k=4"), col =c("blue"));
  vioplot(error_mean, at = 1, col = c("magenta"), add= TRUE)
  vioplot(error_Phat_noClust[,1], at = 2, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,2], at = 4, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,3], at = 6, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,4], at = 8, col = c("green"), add= TRUE)
  #vioplot(error_Phat_noClust[,5], at = 10, col = c("green"), add= TRUE)
  title(main=paste("Frobenius Norm of Error:  M = ", Ncv, ", N = ", n,", clusters = ", nclust))
  #legend(locator(1), legend= c("Naive", "X*X'", "Up to 10 clusters"), fill = c("magenta","green", "blue"))
}