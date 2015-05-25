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
B_Assign = matrix(c(1, 0, 2, 3),2);
PI <- c(.5,.5)
k <- dim(B)[1];
n = 500;
m = 100;
Z <-generateZ_nxk_exact(PI,n);
P = Z %*% B %*% t(Z);
P = .4*matrix(1,n,n);
P_Assign = Z %*% B_Assign %*% t(Z);

P_Bar = makeP_Bar(P,m,n, hollow = FALSE, serial = FALSE)
Err = P_Bar-P;

ASE= spectEmbed(P_Bar, n, DAdjust=FALSE, opt = myoptions);
P_hat = ASE$X[,1] %*% t(ASE$X[,1]);

P_hat_err2 = ASE$X[,2] %*% t(ASE$X[,2]);
P_hat2_err2 = (P_hat_err2)^2;

ASE_err= spectEmbed(P_Bar-P_hat, n, DAdjust=FALSE, opt = myoptions);
#ASE3= spectEmbed(P_Bar, 3, DAdjust=FALSE, opt = myoptions);

P_hat2 = (P_hat)^2;

P_hat_err = ASE_err$X[,1] %*% t(ASE_err$X[,1]);
P_hat2_err = (P_hat_err)^2;
pred = calcSig(matrix(.4,1,1), c(1), m,n)
print(norm(P-P_Bar,'2')/sqrt(4*.1*.9/m*n))
