# Cross Validation Look at BLSA data
# MDK
rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/HCP")
library("igraph")
library("rARPACK")
library("fields")
library("R.matlab")
source("../CoRR/ZhuGhodsi.R")
source("../graphFunctions.R")


ns = c(25, 50, 100, 200, 300);
nTrials= 500;
test_size = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50, 100)
k = 10;
scheins = c(15, 15, 10, 5, 5);
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
attach(mtcars)
par(mfrow=c(1,1));

MSE_ABar = matrix(0, length(ns), length(test_size));
MSE_PHat = matrix(0, length(ns), length(test_size));

for (ni in 1:length(ns)){
  n = ns[ni];
  M = 451*1;
  A_all <- array(rep(0, n*n*M), dim=c(n, n, M))
  for (sub in 1:451) {
    for (scan in 1:1) {
      Mat <- readMat(paste("../../../Data/HCP/subject", sub, "/scan", scan ,"/sg_", n,"_fmri_thresholded_0.2.mat", sep=""));
      A_all[,, (sub-1)*1 +scan] <- Mat$G;
    }
  }
  
  schein = scheins[ni];
  
  A_Sum = rowSums(A_all[], dims = 2)
  
  #A_Sum = A_Sum +t(A_Sum)
  err_P_Hat = matrix(0, length(test_size), nTrials)
  err_A_Bar = matrix(0, length(test_size), nTrials)
  
  err_P_Hat_l = c()
  err_A_Bar_l = c()
  
  for (i in 1:length(test_size)){
    print(c(ni, i))
    s = test_size[i];
    for (j in 1:nTrials){
      samp = sample.int(M,s);
      A_Bar = rowSums(A_all[,,samp], dims = 2);
      #   A_Bar = A_Bar+t(A_Bar)
      P_Bar = (A_Sum-A_Bar)/(M-s);
      A_Bar = (A_Bar)/s;
      k = findCutoff(A_Bar, 3);
      ASE = Schein(A_Bar, k, schein, myoptions)
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
  
  MSE_PHat[ni,] = rowSums(err_P_Hat)/nTrials;
  MSE_ABar[ni,] = rowSums(err_A_Bar)/nTrials;
  
  # for (i in 1:M){
  #   samp = i;
  #   s =1;
  #   A_Bar = A_all[,,samp];
  #   A_Bar = A_Bar+t(A_Bar)
  #   P_Bar = (A_Sum-A_Bar)/(M-s);
  #   A_Bar = (A_Bar)/s;
  #   ASE = Schein(A_Bar, k, 15, myoptions)
  #   P_Hat = ASE$X %*% t(ASE$X);
  #   P_Hat[P_Hat>1] =1;
  #   P_Hat[P_Hat<0] = 0;
  #   err_v <- P_Bar-P_Hat;
  #   err_v <- err_v[upper.tri(err_v)];
  #   err_P_Hat_l[i] <- mean(err_v^2);
  #   
  #   err_v <- P_Bar-A_Bar;
  #   err_v <- err_v[upper.tri(err_v)];
  #   err_A_Bar_l[i] <- mean(err_v^2);
  # }

}