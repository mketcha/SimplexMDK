rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/EF_NMF_PlotFunc/")
library("doParallel")
library("foreach")
source("../graphFunctions.R")
mcl = makeCluster(3);
registerDoParallel(mcl);
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
B <- matrix(c(.2, .1, .1, .4),2);
#B <- matrix(c(.3, .4, .4, .3),2);
B_Assign = matrix(c(1, 0, 2, 3),2);
PI <- c(.5,.5)
k <- dim(B)[1];
n = 1500;
m = 10000;
Z <-generateZ_nxk_exact(PI,n);
P = Z %*% B %*% t(Z);
P = .25*matrix(1,n,n);
P_Assign = Z %*% B_Assign %*% t(Z);

P_Bar = makeP_Bar(P,m,n, hollow = FALSE, serial = FALSE)
Err = P_Bar-P;

ASE= spectEmbed(P_Bar, n, DAdjust=FALSE, opt = myoptions);
ASE_err= spectEmbed(P-P_Bar, n, DAdjust=FALSE, opt = myoptions);
#ASE3= spectEmbed(P_Bar, 3, DAdjust=FALSE, opt = myoptions);
P_hat = ASE$X[,3] %*% t(ASE$X[,3]);
P_hat2 = (P_hat)^2;
var_k3 = c()
for (i in 1:3){
  var_k3[i] = mean(P_hat2[P==B[B_Assign==i]])
}

P_hat_err = ASE_err$X[,1] %*% t(ASE_err$X[,1]);
P_hat2_err = (P_hat_err)^2;
var_k3_err = c()
for (i in 1:3){
  var_k3_err[i] = mean(P_hat2_err[P==B[B_Assign==i]])
}
# var_i = c();
# for(i in 1:n){
#   P_hat = ASE$X[,i] %*% t(ASE$X[,i]);
#   var_i[i] = mean(as.vector(P_hat)^2);
# }

# for (i in 1:1000){
#   for (j in 1:i){
#     if (sum(breakScalar(i,j)) != i)
#       print(c(i,j))
#   }
# }