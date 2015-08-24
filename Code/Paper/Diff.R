# Relative Efficiency Analysis
rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/Paper")
library("igraph")
library("rARPACK")
library("mclust")
library("fields")
source("../graphFunctions.R")
library("doParallel")
library("foreach")

N = c(10, 20, 50, 100, 250, 500, 1000,1500)


B <- matrix(c(.42, .2, .2, .7),2);
PI <- c(.5,.5);
k = 2;
trials <- 1000;

diff_v = matrix(0,trials, length(N) );


for (run in 1:length(N)){
  print(run)
  n <- N[run];
  m <- 100
  Z <-generateZ_nxk_exact(PI,n);
  P = Z %*% B %*% t(Z);
  var_P_temp = matrix(0,3,trials);
  eff_temp = matrix(0,3,trials);
  diag(P)<-0;
  for (i in 1:trials){
    A_Bar = makeP_Bar(P,m,n)
    ASE = Schein(A_Bar, k, 5, myoptions);
    P_hat = ASE$X %*% t(ASE$X);
    err_v <- A_Bar-P_hat;
    err_v <- err_v[upper.tri(err_v)];
    diff_v[i,run] <- mean(err_v^2);
  }
}
dff = colMeans(diff_v)
plot(N, dff, type= 'b', pch = 19)
