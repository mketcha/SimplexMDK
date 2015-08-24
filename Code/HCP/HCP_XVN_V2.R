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
# n = 50;
# M = 451*1;
# A_all <- array(rep(0, n*n*M), dim=c(n, n, M))
# for (sub in 1:451) {
#   for (scan in 1:1) {
#     Mat <- readMat(paste("../../../Data/HCP/subject", sub, "/scan", scan ,"/sg_", n,"_fmri_thresholded_0.2.mat", sep=""));
#     A_all[,, (sub-1)*1 +scan] <- Mat$G;
#   }
# }


nTrials= 200;
test_size = 2:100;

schein = 5;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
attach(mtcars)
par(mfrow=c(1,1));
# A_Sum = rowSums(A_all[], dims = 2)
# A_Sum = A_Sum +t(A_Sum)
# err_P_Hat = matrix(0, length(test_size), nTrials)
# err_A_Bar = matrix(0, length(test_size), nTrials)
err_P_Hat = rep(0,nTrials);
err_A_Bar = rep(0,nTrials);
prev = 1;
CutOff = c();
for (N in 1:length(ns)){
  n = ns[N];
  M = 451*1;
  A_all <- array(rep(0, n*n*M), dim=c(n, n, M))
  for (sub in 1:451) {
    for (scan in 1:1) {
      Mat <- readMat(paste("../../../Data/HCP/subject", sub, "/scan", scan ,"/sg_", n,"_fmri_thresholded_0.2.mat", sep=""));
      A_all[,, sub] <- Mat$G;
    }
  }
  A_Sum = rowSums(A_all[], dims = 2);
  for (i in prev:length(test_size)){
    s = test_size[i];
    print(c(n,s))
    for (j in 1:nTrials){
      samp = sample.int(M,s);
      A_Bar = rowSums(A_all[,,samp], dims = 2);
      #   A_Bar = A_Bar+t(A_Bar)
      P_Bar = (A_Sum-A_Bar)/(M-s);
      A_Bar = (A_Bar)/s;
      
      k = findCutoff(A_Bar,3);
      ASE = Schein(A_Bar, k, schein, myoptions)
      P_Hat = ASE$X %*% t(ASE$X);
      P_Hat[P_Hat>1] =1;
      P_Hat[P_Hat<0] = 0;
      err_v <- P_Bar-P_Hat;
      err_v <- err_v[upper.tri(err_v)];
      err_P_Hat[j] <- mean(err_v^2);
      
      err_v <- P_Bar-A_Bar;
      err_v <- err_v[upper.tri(err_v)];
      err_A_Bar[j] <- mean(err_v^2);
    }
    if (mean(err_P_Hat)>mean(err_A_Bar)){
      prev = i;
      CutOff[N] = test_size[i];
      break;
    }
  }
}
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

CutOff