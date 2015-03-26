setwd("~/Simplex/SimplexMDK/Code")
library("igraph")
library("vioplot")
library("rARPACK")
library("mclust")
library("nloptr")
library("fields")
source("graphFunctions.R")
source("generateDirilechtSBM.R");

B <- matrix(c(.42, .42, .42,.5),2);
PI = c(.6,.4);
k <-2;
n = 1000;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000;
M <- c(1, 5, 10, 50, 100, 1000)
attach(mtcars)
par(mfrow=c(3,2));
#BASE <- generateDirilechtBASE(B, myoptions)

Z <-generateZ_nxk(PI,n)
P = Z %*% B %*% t(Z);

for (i in 1:length(M)){
  A_Bar <- matrix(0,n,n);
  for (h in 1:M[i]) {
    g <- P2Adj(P);
    A_Bar <- A_Bar +g;
  }
  A_Bar <- A_Bar/M[i];
  ASE <- spectEmbed(A_Bar, k, DAdjust = TRUE, opt = myoptions)
  
  plot(ASE$X[,1],ASE$X[,2], ylim = c(-1,1), xlim = c(-1,1))
  points(ASE$X[(Z[,1]==1),1],ASE$X[(Z[,1]==1),2], col ="blue")
  points(ASE$X[(Z[,2]==1),1],ASE$X[(Z[,2]==1),2], col ="red")
  title(main=paste("Estimated Latent Space:  M = ", M[i], ", N = ", n))
}
