# Examination of Block Size for BLSA data
setwd("~/Simplex/SimplexMDK/Code")
library("igraph")
library("vioplot")
library("rARPACK")
library("mclust")
library("nloptr")
library("fields")
library("cvTools")
library("R.matlab")
source("graphFunctions.R")
source("generateDirilechtSBM.R");

Mat <- readMat("../Data/BLSA0317.mat");
A_all <- Mat$Abinary;
A_w <- Mat$Aweighted;
n <- 70;
M <- 49;
max_k = 4;
nclust <- 2;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 2000000
# attach(mtcars)
# par(mfrow=c(3,2));
error_Phat_noClust <- matrix(0, M, max_k);
error_Phat <- matrix(0, M, max_k);
error_mean <- matrix(0, M, 1);
A_bar <- matrix(0,n,n);
A_wbar<- matrix(0,n,n);
A_HL = matrix(0,n,n)
for (i in 1:n){
  for (j in i:n){
    meds = rep(0, M*(M+1)/2)
    a = 1
    for (k in 1:M){
      for (l in k:M){
         meds[a] <- (A_w[i,j,k]+A_w[i,j,l])/2
         a = a+1;
      }
    }
    A_HL[i,j] = median(meds);
    A_HL[j,i] = median(meds);
  }
}
for (j in 1:M) {
    A_bar <- A_bar + A_all[,,j] + t(A_all[,,j]);
    A_wbar <- A_wbar + A_w[,,j] + t(A_w[,,j]);
}
A_bar <- A_bar/(M)
A_wbar <- A_wbar/(M)
# image.plot(A_bar[,ncol(A_bar):1], col = gray((0:100)/100))
# title("A_bar")
# d <- c(2,3,4,5,6,7,8)
# for (i in 1:length(d)){
#   ASE <- try(spectEmbed(A_bar,d[i], DAdjust=TRUE, opt = myoptions))
#   cl <- Mclust(ASE$X,d[i])
#   ASE_Xn = matrix(0, dim(ASE$X)[1],dim(ASE$X)[2])
#   count = 1;
#   for (j in 1:d[i]){
#     for (k in 1:n){
#       if (cl$classification[k]==j){
#         ASE_Xn[count,] = ASE$X[k,];
#         count = count +1;
#       }
#     }
#   }
#   P_hat = ASE_Xn %*% t(ASE_Xn)
#   image.plot(P_hat[,ncol(P_hat):1], col = gray((0:100)/100))
#   title(main=paste("P_hat, d,k = ",d[i]))
# }
d = 5
ASE_HL <- try(spectEmbed(A_HL,d, DAdjust=TRUE, opt = myoptions))
ASE_wbar <- try(spectEmbed(A_wbar,d, DAdjust=TRUE, opt = myoptions))
cl_HL <- Mclust(ASE_HL$X,d)
cl_wbar <- Mclust(ASE_wbar$X,d)
ASE_wn = eigen(A_wbar);
# X = ASE_n$vectors %*% diag(sqrt(abs(ASE_n$values)));
ASE_Xn_HL = matrix(0, dim(ASE_HL$X)[1],dim(ASE_HL$X)[2])
ASE_Xn_wbar = matrix(0, dim(ASE_wbar$X)[1],dim(ASE_wbar$X)[2])
A_barX_wn = matrix(0, dim(ASE_wn$vectors)[1],dim(ASE_wn$vectors)[2])
count1 = 1;
count2 = 1;
for (j in 1:d){
  for (k in 1:n){
    if (cl_HL$classification[k]==j){
      ASE_Xn_HL[count1,] = ASE_HL$X[k,];
      A_barX_wn[count1,] = ASE_wn$vectors[k,]
      ASE_Xn_wbar[count1,] = ASE_wbar$X[k,];
      count1 = count1 +1;
    }
    if (cl_wbar$classification[k]==j){
      
      count2 = count2+1;
    }
  }
}



ASE <- try(spectEmbed(A_bar,d, DAdjust=TRUE, opt = myoptions))
cl <- Mclust(ASE$X,d)
ASE_n = eigen(A_bar);
# X = ASE_n$vectors %*% diag(sqrt(abs(ASE_n$values)));
ASE_Xn = matrix(0, dim(ASE$X)[1],dim(ASE$X)[2])
A_barX_n = matrix(0, dim(ASE_n$vectors)[1],dim(ASE_n$vectors)[2])

count = 1;
for (j in 1:d){
  for (k in 1:n){
    if (cl$classification[k]==j){
      ASE_Xn[count,] = ASE$X[k,];
      A_barX_n[count,] = ASE_n$vectors[k,]
      count = count +1;
    }
  }
}

s= 10
separ = matrix(1, s,n);
A_bar_wn = A_barX_wn %*% diag(ASE_wn$values) %*%t(A_barX_wn)
A_bar_wn[A_bar_wn>1] = 1;
A_bar_wn[A_bar_wn<0] = 0;

P_hat_HL = ASE_Xn_HL %*% t(ASE_Xn_HL)
P_hat_HL[P_hat_HL>1] = 1;
P_hat_HL[P_hat_HL<0] = 0;

P_hat_wbar = ASE_Xn_wbar %*% t(ASE_Xn_wbar)
P_hat_wbar[P_hat_wbar>1] = 1;
P_hat_wbar[P_hat_wbar<0] = 0;
Weighted = rbind(A_bar_wn, separ,P_hat_HL, separ,P_hat_wbar)


A_bar_n = A_barX_n %*% diag(ASE_n$values) %*%t(A_barX_n)
A_bar_n[A_bar_n>1] = 1;
A_bar_n[A_bar_n<0] = 0;

P_hat = ASE_Xn %*% t(ASE_Xn)
P_hat[P_hat>1] = 1;
P_hat[P_hat<0] = 0;
UnWeighted = rbind(A_bar_n, separ,P_hat)
image.plot(UnWeighted[,ncol(UnWeighted):1], col = gray((0:100)/100), asp =70/(2*70+1*s))
title("Weighted Graphs")
# image.plot(A_bar_n[,ncol(A_bar_n):1], col = gray((0:100)/100), asp =1)
# title("A_bar_unweighted")


# image.plot(P_hat_HL[,ncol(P_hat_HL):1], col = gray((0:100)/100), asp =1)
# title(main=paste("P_hat_HL, d,k = ",d))


# image.plot(P_hat[,ncol(P_hat):1], col = gray((0:100)/100), asp =1)
# title(main=paste("P_hat, d,k = ",d))

# image.plot(P_hat_wbar[,ncol(P_hat_wbar):1], col = gray((0:100)/100), asp =1)
# title(main=paste("P_hat_wbar, d,k = ",d))
