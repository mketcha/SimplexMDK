rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/EF_NMF_PlotFunc/")
library("doParallel")
library("foreach")
source("../graphFunctions.R")
mcl = makeCluster(3);
registerDoParallel(mcl);
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
B_1 <- matrix(c(.6, .2, .2, .6),2);
B_2 <- matrix(c(.7, .2, .2, .7),2);
B_3 <- matrix(c(.8, .2, .2, .8),2);
B_4 <- matrix(c(.9, .2, .2, .9),2);
B_5 <- matrix(c(.93, .2, .2, .93),2);
B_6 <- matrix(c(.94, .2, .2, .94),2);
B_7 <- matrix(c(.95, .2, .2, .95),2);
B_8 <- matrix(c(.96, .2, .2, .96),2);
B_9 <- matrix(c(.97, .2, .2, .97),2);
B_10 <- matrix(c(.98, .2, .2, .98),2);
B_11 <- matrix(c(.99, .2, .2, .99),2);
B <- matrix(c(.5, .2, .2, .5),2);
B1 <- matrix(c(.4, .2, .2, .4),2);
B2 <- matrix(c(.3, .2, .2, .3),2);
B3 <- matrix(c(.2, .2, .2, .2),2);
B4 <- matrix(c(.1, .2, .2, .1),2);
B5 <- matrix(c(.07, .2, .2, .07),2);
B6 <- matrix(c(.06, .2, .2, .06),2);
B7 <- matrix(c(.05, .2, .2, .05),2);
B8 <- matrix(c(.04, .2, .2, .04),2);
B9 <- matrix(c(.03, .2, .2, .03),2);
B10 <- matrix(c(.02, .2, .2, .02),2);
B11 <- matrix(c(.01, .2, .2, .01),2);
Bs <- list(B_11,B_10,B_9,B_8,B_7,B_6,B_5,B_4,B_3,B_2,B_1,B,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11);
X = c(.99,.98,.97,.96,.95,.94,.93,.9,.8,.7,.6,.5,.4,.3,.2,.1,.07,.06,.05,.04,.03,.02,.01)
out = matrix(0,length(X),3)
Pout = c();
est = c();
for (run in 1:length(X)){
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
  est[run] = calcSig(B, PI,m,n)
  Pout[run] = mean(P_hat2_err)
  out[run,] = var_k3_err
}