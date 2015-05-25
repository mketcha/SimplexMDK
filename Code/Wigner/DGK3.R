rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/EF_NMF_PlotFunc/")
library("doParallel")
library("foreach")
source("../graphFunctions.R")
mcl = makeCluster(3);
registerDoParallel(mcl);
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
B <- matrix(c(.4, .5, .5, .4),2);
B1 <- matrix(c(.4, .4, .4, .4),2);
B2 <- matrix(c(.4, .2, .2, .4),2);
B3 <- matrix(c(.4, .1, .1, .4),2);
B4 <- matrix(c(.4, .08, .08, .4),2);
B5 <- matrix(c(.4, .07, .07, .4),2);
B6 <- matrix(c(.4, .06, .06, .4),2);
B7 <- matrix(c(.4, .05, .05, .4),2);
B8 <- matrix(c(.4, .04, .04, .4),2);
B9 <- matrix(c(.4, .03, .03, .4),2);
B10 <- matrix(c(.4, .02, .02, .4),2);
B11 <- matrix(c(.4, .01, .01, .4),2);
Bs <- list(B,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11);
X = c(.5,.4,.2,.1,.08,.07,.06,.05,.04,.03,.02,.01)
out = matrix(0,12,3)
Pout = rep(0,12);
for (run in 1:12){
  print(run)
  B <- Bs[[run]]
  B_Assign = matrix(c(1, 0, 2, 3),2);
  PI <- c(.5,.5)
  k <- dim(B)[1];
  n = 1000;
  m = 1000;
  Z <-generateZ_nxk_exact(PI,n);
  P = Z %*% B %*% t(Z);
  #P = .25*matrix(1,n,n);
  P_Assign = Z %*% B_Assign %*% t(Z);
  var_k3_err = c(0,0,0)
  for (i in 1:10){
    P_Bar = makeP_Bar(P,m,n, hollow = FALSE, serial = FALSE)
  #   ASE= spectEmbed(P_Bar, n, DAdjust=FALSE, opt = myoptions);
    ASE_err= spectEmbed(P-P_Bar, n, DAdjust=FALSE, opt = myoptions);
    #ASE3= spectEmbed(P_Bar, 3, DAdjust=FALSE, opt = myoptions);
  #   P_hat = ASE$X[,3] %*% t(ASE$X[,3]);
  #   P_hat2 = (P_hat)^2;
  #   var_k3 = c()
  #   for (i in 1:3){
  #     var_k3[i] = mean(P_hat2[P==B[B_Assign==i]])
  #   }
    
    P_hat_err = ASE_err$X[,1] %*% t(ASE_err$X[,1]);
    P_hat2_err = (P_hat_err)^2;
    
    for (i in 1:3){
      var_k3_err[i] = var_k3_err[i] + mean(P_hat2_err[P==B[B_Assign==i]])/10
      Pout[run] = Pout[run]+mean(P_hat2_err)/10
    }
  }
  out[run,] = var_k3_err
}