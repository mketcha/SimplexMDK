# Cross Validation Look at BLSA data
# MDK
rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/Paper")
library("igraph")
library("rARPACK")
library("fields")
library("R.matlab")
source("../graphFunctions.R")

Mat <- readMat("../../Data/BLSA0317.mat");
A_all <- Mat$Abinary;

n <- 70;
M <- 49;
nTrials= 500;
test_size = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20)
k = 68;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
attach(mtcars)
par(mfrow=c(4,2));
A_Sum = rowSums(A_all[], dims = 2)
A_Sum = A_Sum +t(A_Sum)
err_P_Hat = matrix(0, length(test_size), nTrials)
err_A_Bar = matrix(0, length(test_size), nTrials)

err_P_Hat_l = c()
err_A_Bar_l = c()

for (i in 1:length(test_size)){
  print(i)
  s = test_size[i];
  for (j in 1:nTrials){
    samp = sample.int(M,s);
    A_Bar = rowSums(A_all[,,samp], dims = 2);
    A_Bar = A_Bar+t(A_Bar)
    P_Bar = (A_Sum-A_Bar)/(M-s);
    A_Bar = (A_Bar)/s;
    ASE = Schein(A_Bar, k, 15, myoptions)
    P_Hat = ASE$X %*% t(ASE$X);
    P_Hat[P_Hat>1] =1;
    P_Hat[P_Hat<0] = 0;
    err_v <- P_Bar-P_Hat;
    err_v <- err_v[upper.tri(err_v)];
    err_P_Hat[i,j] <- mean(err_v^2);
    
    err_v <- P_Bar-A_Bar;
    err_v <- err_v[upper.tri(err_v)];
    err_A_Bar[i,j] <- mean(err_v^2);
  }
}

for (i in 1:M){
  samp = i;
  s =1;
  A_Bar = A_all[,,samp];
  A_Bar = A_Bar+t(A_Bar)
  P_Bar = (A_Sum-A_Bar)/(M-s);
  A_Bar = (A_Bar)/s;
  ASE = Schein(A_Bar, k, 15, myoptions)
  P_Hat = ASE$X %*% t(ASE$X);
  P_Hat[P_Hat>1] =1;
  P_Hat[P_Hat<0] = 0;
  err_v <- P_Bar-P_Hat;
  err_v <- err_v[upper.tri(err_v)];
  err_P_Hat_l[i] <- mean(err_v^2);
  
  err_v <- P_Bar-A_Bar;
  err_v <- err_v[upper.tri(err_v)];
  err_A_Bar_l[i] <- mean(err_v^2);
}

rowSums(err_P_Hat)
rowSums(err_A_Bar)
