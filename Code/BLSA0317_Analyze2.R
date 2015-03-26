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
n <- 70;
M <- 49;
max_k = 4;
nclust <- 2;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 2000000
attach(mtcars)
par(mfrow=c(2,1));
error_Phat_noClust <- matrix(0, M, max_k);
error_Phat <- matrix(0, M, max_k);
error_mean <- matrix(0, M, 1);
A_bar <- matrix(0,n,n);

for (j in 1:M) {
    A_bar <- A_bar + A_all[,,j] + t(A_all[,,j]);
}
A_bar <- A_bar/(M)
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
A_bar_n = A_barX_n %*% diag(ASE_n$values) %*%t(A_barX_n)
image.plot(A_bar_n[,ncol(A_bar_n):1], col = gray((0:100)/100))
title("A_bar")
P_hat = ASE_Xn %*% t(ASE_Xn)
image.plot(P_hat[,ncol(P_hat):1], col = gray((0:100)/100))
title(main=paste("P_hat, d,k = ",d))
