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

m = 100;
n = 500;
k = 2;
trials <- 1000;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000

B <- matrix(c(.42, .2, .2, .7),2);
PIs <- matrix(c(.1,.9,.2,.8,.3,.7,.4,.6,.5,.5,.6,.4,.7,.3,.8,.2,.9,.1),9,2,byrow= TRUE);
B_Assign = matrix(c(1, 0, 2, 3),2);
modes = dim(PIs)[1];

eff_pred = c()
var_P_hat_pred = matrix(0,3,modes);

for (run in 1:modes){
  eff_pred[run] = (1/PIs[run,1]+1/PIs[run,2])/n
  for (i in 1:3){
    var_P_hat_pred[i,run] =  (1/PIs[run,1]+1/PIs[run,2])*B[B_Assign==i]*(1-B[B_Assign==i])/n/m
  }
}

eff = matrix(0,3, modes);
var_eff = matrix(0,3,modes);

var_P_hat = matrix(0,3,modes);
var_var_P_hat  = matrix(0,3,modes);

for (run in 1:modes){
  print(run)
  PI = PIs[run,]
  Z <-generateZ_nxk_exact(PI,n);
  P = Z %*% B %*% t(Z);
  var_P_temp = matrix(0,3,trials);
  eff_temp = matrix(0,3,trials);
  diag(P) = 0;
  for (i in 1:trials){
    A_Bar = makeP_Bar(P,m,n)
    ASE = Schein(A_Bar, k, 1, myoptions);
    P_hat = ASE$X %*% t(ASE$X);
    P_hat_var = (P-P_hat)^2;
    A_Bar_var = (P-A_Bar)^2;
    for (j in 1:3){
      var_P_temp[j,i] = mean(P_hat_var[P==B[B_Assign==j]]);
      eff_temp[j,i] =  mean(P_hat_var[P==B[B_Assign==j]])/mean(A_Bar_var[P==B[B_Assign==j]])
    }
  }
  for (j in 1:3){
    eff[j,run] = mean(eff_temp[j,])
    var_eff[j,run] = var(eff_temp[j,])
    var_P_hat[j,run] = mean(var_P_temp[j,])
    var_var_P_hat[j,run]  = var(var_P_temp[j,]);
  }
}

eff_pred = matrix(0,3,modes);
var_P_hat_pred = matrix(0,3,modes);
list = matrix(c(1,1,1,2,2,2), 3,2, byrow=TRUE)
for (run in 1:modes){
  for (i in 1:3){
    eff_pred[i,run] = (1/PIs[run,list[i,1]]+1/PIs[run,list[i,2]])/n
    var_P_hat_pred[i,run] =  (1/PIs[run,list[i,1]]+1/PIs[run,list[i,2]])*B[B_Assign==i]*(1-B[B_Assign==i])/n/m
  }
}

par(mfrow=c(2,2));
matplot( PIs[,1], t(var_P_hat_pred),type = c("b","b","b"), xlab=  expression(paste(rho[i], " , ", 1-rho[j])), ylim = c(0,.0001),ylab = "Predicted Variance",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))

matplot( PIs[,1], t(var_P_hat),type = c("b","b","b"), xlab=  expression(paste(rho[i], " , ", 1-rho[j])),ylim = c(0,.0001), ylab = "Simulation Variance",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))
matplot( PIs[,1], t(eff_pred),type = c("b","b","b"), xlab= expression(paste(rho[i], " , ", 1-rho[j])), ylim =c(0,.04), ylab = "Predicted RE",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))

matplot( PIs[,1], t(eff),type = c("b","b","b"), xlab=  expression(paste(rho[i], " , ", 1-rho[j])), ylim=c(0,.04), ylab = "Simulation RE",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))
# matplot( N_mat, cbind(t(eff_M[,1:7]))*N_mat,type = c("b","b","b"), xlab= "N", ylab = "N*RE",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))
# title("M = 100")
# matplot( c(-200,2200), c(eff_pred_N, eff_pred_N)*1000,type = c("l","l","l"), xlab= "N", ylab = "RE*N",lty = c(2,2,2), lwd = 3, col = c("black"), add= TRUE)
# matplot( M_mat, cbind(t(eff_N[,1:7])*1000),type = c("b","b","b"), xlab= "M", xlim = c(-100,2100), ylim = c(3.8,4.2), ylab = "N*RE",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))
# title("N = 1000")
# matplot( c(-200,2200), c(eff_pred_N, eff_pred_N)*1000,type = c("l","l","l"), xlab= "N", ylab = "RE*N",lty = c(2,2,2), lwd = 4, col = c("black"), add= TRUE)

