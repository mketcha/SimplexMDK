rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/EF_NMF_PlotFunc/")
library("doParallel")
library("foreach")
source("../graphFunctions.R")
mcl = makeCluster(3);
registerDoParallel(mcl);
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
B <- matrix(c(.3, .05, .05, .9),2);
K <- dim(B)[1];
#B <- matrix(c(.1, .05, .05, .9),2);
B_Assign = matrix(c(1, 0, 2, 3),2);
PIs = matrix(c(.99,.01, .95,.05,.9,.1,.75,.25,.5,.5,.25,.75,.1,.9,.05,.95,.01,.99),9,2,byrow=TRUE)
ratio = matrix(0,9,3);
VarPij = matrix(0,9,3);
dan = matrix(rep(1/PIs[,1]+1/PIs[,2],3), ncol = 3)
m = 100;
max_j = 4;
err_block <- rep(list(matrix(0,3,9)),max_j)
eff_block <- rep(list(matrix(0,3,9)),max_j)
eff_block_sub <- rep(list(matrix(0,3,9)),max_j)
err_Naive_block <- matrix(0,3,9)
error_mean <- matrix(0, 9, 1);
P_err_mean <- matrix(max_j,9)
dan = matrix(0,9,K*(K+1)/2);
trials = 10;
lamdas = c();
lamda_est <- c();
for (run in 1:9){
  comb = matrix(0, K*(K+1)/2, 2)
  count = 1;
  for (i in 1:K){
    for (j in i:K){
      comb[count,] = c(i,j);
      count =count+1
    }
  }
  print(run)
  PI <- PIs[run,];
  n = 2000;
  for (i in 1:dim(comb)[1]) {
    #dan[run,i] = (1/PI[comb[i,1]]+1/PI[comb[i,2]])*(B[B_Assign==i])*(1-B[B_Assign==i])/n
    dan[run,i] = (1/PI[comb[i,1]]+1/PI[comb[i,2]])/n
  }
  Z <-generateZ_nxk_exact(PI,n);
  P = Z %*% B %*% t(Z);
  lamda_est[run] = limLamda(P,m,n)
  P_Assign = Z %*% B_Assign %*% t(Z);
  for (i in 1:trials){
    P_Bar = makeP_Bar(P,m,n, hollow = FALSE, serial = FALSE)
    Err = P_Bar-P;
    err_v <- P-P_Bar;
    err_v <- err_v[upper.tri(err_v)];
    error_mean[run,1] <- mean(err_v^2);
    ASE_Err = spectEmbed(Err, 1, DAdjust=FALSE, opt = myoptions);
    lamdas[run] = ASE_Err$lamda[1]
    for (j in 1:max_j){
      ASE = spectEmbed(P_Bar, j, DAdjust=FALSE, opt = myoptions)
      P_hat = ASE$X %*% t(ASE$X);
      P_hat2 = (P-P_hat)^2;
      err_P_hat = P_hat2[upper.tri(P_hat2)]
      
      err_Assign = P_Assign[upper.tri(P_Assign)];
      for (k in 1:3) {
        err_block[[j]][k, run] = err_block[[j]][k, run]+mean((err_P_hat[(err_Assign==k)]))/trials;
        err_Naive_block[k, run] = err_Naive_block[k, run] +mean((err_v[(err_Assign==k)])^2)/trials;
        eff_block[[j]][k,run] = eff_block[[j]][k,run]+mean(( err_P_hat[(err_Assign==k)]))/mean((err_v[(err_Assign==k)])^2)/trials;
        eff_block_sub[[j]][k,run] = eff_block[[j]][k,run] - dan[run,k]/trials;
      }
    }
  }
}