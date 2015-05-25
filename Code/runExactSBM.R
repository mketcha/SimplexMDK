setwd("~/Simplex/SimplexMDK/Code")
library("igraph")
library("vioplot")
library("rARPACK")
library("mclust")
library("nloptr")
library("fields")
source("graphFunctions.R")
source("generateDirilechtSBM.R");


MN = matrix(c(1,10,10,100,100,10, 100, 100),4,2,byrow=TRUE)
#MN = matrix(c(1,100,100,100),4,2,byrow=TRUE)
#B <- matrix(c(.22, .22, .22,.6),2);
#B <- matrix(c(.42, .42, .42,.5),2);
B <- matrix(c(.42, .2, .2, .7),2);
#B <- matrix(c(.42, .2, .3, .2,.5, .15, .3, .15, .7),3);
#B <- matrix(c(.42, .23, .35, .23,.5, .15, .35, .15, .7),3);
#B <- matrix(c(.42, .2, .3 ,.05, .2, .49, .35, .15, .3, .35, .6, .25, .05 , .15,.25, .7),4);
#PI <- c( .25, .25, .25, .25);
#PI <- c(.33, .33,.34)
PI <- c(.5,.5);
trials <- 1000;
max_k <- 4;
nclust <- 2;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000

attach(mtcars)
par(mfrow=c(3,2));

for (run in 1:4) {
  n <- MN[run,2];
  m <- MN[run,1];
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
    #     err_v <- err_v[upper.tri(err_v)];
    #     error_mean[i,1] <- mean(err_v^2);
    error_mean[i,1] <- norm(err_v, type="2");
    #diag(P_Bar) <- diag(P);
    for (j in 1:max_k) {
      repeat{
        ASE <- try(spectEmbed(P_Bar,j, DAdjust=TRUE, opt = myoptions))
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
      err_v_noClust <- P-P_hat_noClust;
      #       err_v_noClust <- err_v_noClust[upper.tri(err_v_noClust)];
      #       error_Phat_noClust[i,j] <- mean(err_v_noClust^2);
      error_Phat_noClust[i,j] <- norm(err_v_noClust, type = "2");
      
      P_hat <- clust %*% t(clust)
      err_v <- P-P_hat;
      #       err_v <- err_v[upper.tri(err_v)];
      #       error_Phat[i,j] <- mean(err_v^2);
      error_Phat[i,j] <- norm(err_v, type = "2");
      
    }
  }
  print(mean(error_mean));
  print(mean(error_Phat_noClust[,2]))
  print(mean(error_Phat[,2]))
  vioplot( error_mean, error_Phat_noClust[,1], error_Phat[,1], error_Phat_noClust[,2], error_Phat[,2], error_Phat_noClust[,3], error_Phat[,3], error_Phat_noClust[,4], error_Phat[,4], names = c("Naive", "d=1", "d=1", "d=2", "d=2", "d=3", "d=3","d=4", "d=4"), col =c("blue"));
  vioplot(error_mean, at = 1, col = c("magenta"), add= TRUE)
  vioplot(error_Phat_noClust[,1], at = 2, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,2], at = 4, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,3], at = 6, col = c("green"), add= TRUE)
  vioplot(error_Phat_noClust[,4], at = 8, col = c("green"), add= TRUE)
  #vioplot(error_Phat_noClust[,5], at = 10, col = c("green"), add= TRUE)
  title(main=paste("M = ", m, ", N = ", n,", Exact SBM, clusters = ", nclust))
  #legend(locator(1), legend= c("Naive", "X*X'", "Up to 10 clusters"), fill = c("magenta","green", "blue"))
  
  # vioplot( error_mean, error_Phat_noClust[,1], error_Phat[,1], error_Phat_noClust[,2], error_Phat[,2], error_Phat_noClust[,3], error_Phat[,3], error_Phat_noClust[,4], error_Phat[,4], error_Phat_noClust[,5], error_Phat[,5], names = c("Naive", "k=1", "k=1", "k=2", "k=2", "k=3", "k=3","k=4", "k=4", "k=5", "k=5"), col =c("blue"));
  # vioplot(error_mean, at = 1, col = c("magenta"), add= TRUE)
  # vioplot(error_Phat_noClust[,1], at = 2, col = c("green"), add= TRUE)
  # vioplot(error_Phat_noClust[,2], at = 4, col = c("green"), add= TRUE)
  # vioplot(error_Phat_noClust[,3], at = 6, col = c("green"), add= TRUE)
  # vioplot(error_Phat_noClust[,4], at = 8, col = c("green"), add= TRUE)
  # vioplot(error_Phat_noClust[,5], at = 10, col = c("green"), add= TRUE)
  # #vioplot(error_Phat_noClust[,5], at = 10, col = c("green"), add= TRUE)
  # title(main=paste("M = ", m, ", N = ", n,", r = ", r,", clusters = ", nclust))
  # #legend(locator(1), legend= c("Naive", "X*X'", "Up to 10 clusters"), fill = c("magenta","green", "blue"))
  
  #   vioplot( error_mean, error_Phat_noClust[,2], error_Phat[,2], error_Phat_noClust[,3], error_Phat[,3], error_Phat_noClust[,4], error_Phat[,4], names = c("Naive", "k=2", "k=2", "k=3", "k=3","k=4", "k=4"), col =c("blue"));
  #   vioplot(error_mean, at = 1, col = c("magenta"), add= TRUE)
  #   vioplot(error_Phat_noClust[,2], at = 2, col = c("green"), add= TRUE)
  #   vioplot(error_Phat_noClust[,3], at = 4, col = c("green"), add= TRUE)
  #   vioplot(error_Phat_noClust[,4], at = 6, col = c("green"), add= TRUE)
  #   #vioplot(error_Phat_noClust[,5], at = 10, col = c("green"), add= TRUE)
  #   title(main=paste("M = ", m, ", N = ", n,", r = ", r,", clusters = ", nclust))
}
legend(locator(1), legend= c("Naive", "X_hat*X_hat'", "2 clusters"), fill = c("magenta","green", "blue"))
image.plot(P[,ncol(P):1], col = gray((0:100)/100))
title("P")