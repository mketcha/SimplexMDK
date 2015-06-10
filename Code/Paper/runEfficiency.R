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

MN_M = matrix(c(100,10, 100, 20, 100,50,100,100,100,250,100, 500,100,1000),7,2,byrow=TRUE)
MN_N = matrix(c(1,1000,10,1000,50,1000,100, 1000,500,1000,1000,1000,2000,1000),7,2,byrow=TRUE)

modes = dim(MN_M)[1];

B <- matrix(c(.42, .2, .2, .7),2);
PI <- c(.5,.5);
k = 2;
trials <- 1000;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000

B_Assign = matrix(c(1, 0, 2, 3),2);

eff_M = matrix(0,3, modes);
var_eff_M = matrix(0,3,modes);

var_P_hat_M = matrix(0,3,modes);
var_var_P_hat_M  = matrix(0,3,modes);

eff_N = matrix(0,3, modes);
var_eff_N = matrix(0,3,modes);

var_P_hat_N = matrix(0,3,modes);
var_var_P_hat_N  = matrix(0,3,modes);

for (run in 1:modes){
  print(run)
  n <- MN_M[run,2];
  m <- MN_M[run,1];
  Z <-generateZ_nxk_exact(PI,n);
  P = Z %*% B %*% t(Z);
  var_P_temp = matrix(0,3,trials);
  eff_temp = matrix(0,3,trials);
  diag(P)<-0;
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
    eff_M[j,run] = mean(eff_temp[j,])
    var_eff_M[j,run] = var(eff_temp[j,])
    var_P_hat_M[j,run] = mean(var_P_temp[j,])
    var_var_P_hat_M[j,run]  = var(var_P_temp[j,]);
  }
}

for (run in 1:modes){
  print(run)
  n <- MN_N[run,2];
  m <- MN_N[run,1];
  Z <-generateZ_nxk_exact(PI,n);
  P = Z %*% B %*% t(Z);
  diag(P)<-0;
  var_P_temp = matrix(0,3,trials);
  eff_temp = matrix(0,3,trials);
  for (i in 1:trials){
    A_Bar = makeP_Bar(P,m,n)
    ASE = try(Schein(A_Bar, k, 1, myoptions));
    if (class(ASE) == "try-error") {
      i = i-1;
      next
    }
    P_hat = ASE$X %*% t(ASE$X);
    P_hat_var = (P-P_hat)^2;
    A_Bar_var = (P-A_Bar)^2;
    for (j in 1:3){
      var_P_temp[j,i] = mean(P_hat_var[P==B[B_Assign==j]]);
      eff_temp[j,i] =  mean(P_hat_var[P==B[B_Assign==j]])/mean(A_Bar_var[P==B[B_Assign==j]])
    }
  }
  for (j in 1:3){
    eff_N[j,run] = mean(eff_temp[j,])
    var_eff_N[j,run] = var(eff_temp[j,])
    var_P_hat_N[j,run] = mean(var_P_temp[j,])
    var_var_P_hat_N[j,run]  = var(var_P_temp[j,]);
  }
}

var_P_hat_M_pred = matrix(0,3,modes);
var_P_hat_N_pred = matrix(0,3,modes);
for (run in 1:modes){
  for (i in 1:3){
    var_P_hat_M_pred[i,run]= 4*B[B_Assign==i]*(1-B[B_Assign==i])/MN_M[run,1]/MN_M[run,2];
    var_P_hat_N_pred[i,run]= 4*B[B_Assign==i]*(1-B[B_Assign==i])/MN_N[run,1]/MN_N[run,2];
  }
}
eff_pred_N = 2*2/1000;
eff_pred_M = 2*2/MN_M[,2];


par(mfrow=c(2,2));
M_mat = matrix(rep(MN_N[,1],3), ncol =3)
N_mat = matrix(rep(MN_M[,2],3), ncol =3)
Nvar_P_hat_M_pred = var_P_hat_M_pred*t(N_mat);
Mvar_P_hat_N_pred = var_P_hat_N_pred*t(M_mat);
matplot( N_mat, cbind(t(var_P_hat_M[,1:7]))*N_mat,type = c("b","b","b"), xlab= "N", ylab = "N*Variance",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))
title("M = 100")
matplot( c(-50,1050), t(Nvar_P_hat_M_pred[,1:2]),type = c("l","l","l"), xlab= "N", ylab = "RE*N",lty = c(2,2,2), lwd = 1.5, col = c("red", "blue", "green"), add= TRUE)
matplot( M_mat, cbind(t(var_P_hat_N[,1:7])*M_mat),type = c("b","b","b"), xlab= "M", ylab = "N*Variance",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))
title("N = 1000")
matplot( c(-50,2050), t(Mvar_P_hat_N_pred[,1:2]),type = c("l","l","l"), xlab= "N", ylab = "RE*N",lty = c(2,2,2), lwd = 1.5, col = c("red", "blue", "green"), add= TRUE)

matplot( N_mat, cbind(t(eff_M[,1:7]))*N_mat,type = c("b","b","b"), xlab= "N", ylab = "N*RE",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))
title("M = 100")
matplot( c(-200,2200), c(eff_pred_N, eff_pred_N)*1000,type = c("l","l","l"), xlab= "N", ylab = "RE*N",lty = c(2,2,2), lwd = 3, col = c("black"), add= TRUE)
matplot( M_mat, cbind(t(eff_N[,1:7])*1000),type = c("b","b","b"), xlab= "M", xlim = c(-100,2100), ylim = c(3.8,4.2), ylab = "N*RE",lty = c(1,1,1,2,2,2), pch = c(15, 16, 17), lwd = 2, col = c("red", "blue", "green"))
title("N = 1000")
matplot( c(-200,2200), c(eff_pred_N, eff_pred_N)*1000,type = c("l","l","l"), xlab= "N", ylab = "RE*N",lty = c(2,2,2), lwd = 4, col = c("black"), add= TRUE)

# plotCI(x = rep(MN_N[,1],6), y = as.vector(cbind(t(eff_M[,1:7]), t(eff_M_NMF[,1:7]))/(matrix(rep(invN,6), ncol =6))),uiw = as.vector(sqrt(cbind(t(eff_M_var[,1:7]), t(eff_M_var_NMF[,1:7])))/matrix(rep(invN,6), ncol =6)),col = rep(c("red", "blue", "green","orange","cyan","darkgreen"),each = 7),add=T)
legend(locator(1), c("ASE: k, d = 2", "ASE: k, d = 3", "ASE: k, d = 4","NMF: k, d = 2", "NMF: k, d = 3", "NMF: k, d = 4"), lty = c(1,1,1,2,2,2), col = c("red", "blue", "green","orange","cyan","darkgreen"), lwd = 2)
title("Relative Efficiency*N, N = 500")