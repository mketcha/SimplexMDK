# Cross Validation Look at BLSA data
# MDK
rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/CoRR")
library("igraph")
library("rARPACK")
library("fields")
library("R.matlab")
source("ZhuGhodsi.R")
source("../graphFunctions.R")


subjects = readLines("../../../Data/CoRR/subnames.txt")
n = 788;
M = 232*2;
A_all <- array(rep(0, n*n*M), dim=c(n, n, M))
for (sub in 1:232) {
  for (session in 1:2) {
    Mat <- readMat(paste("../../../Data/CoRR/", subjects[sub], "/session_", session,"/", "smg_thresholded_0.125_788roi_vR_subject_", subjects[sub], "_session_",session,".mat", sep=""));
    Mat$threshSMG[Mat$threshSMG>0]<-1;
    A_all[,, (sub-1)*2 +session] <- Mat$threshSMG;
  }
}


nTrials= 500;
test_size = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 50, 100, 150, 200)

schein = 5;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
attach(mtcars)
par(mfrow=c(4,2));
A_Sum = rowSums(A_all[], dims = 2)
k = findCutoff(A_Sum/M,3);
#A_Sum = A_Sum +t(A_Sum)

err_P_Hat = array(0, dim=c(length(test_size),nTrials,6));
err_A_Bar = matrix(0, length(test_size), nTrials)

err_P_Hat_l = c()
err_A_Bar_l = c()
CU = c();
for (Z in 1:6){
  k = findCutoff(A_Sum/M,Z);
  CU[Z] = k;
}
  
MSE_PHat = matrix(0,length(test_size),6);
MSE_ABar = c();
for (i in 1:length(test_size)){
  print(i)
  s = test_size[i];
  for (j in 1:nTrials){
    if ((j %% 100)==0){
      print(c(i, j))
    }
    samp = sample.int(M,s);
    A_Bar = rowSums(A_all[,,samp], dims = 2);
 #   A_Bar = A_Bar+t(A_Bar)
    P_Bar = (A_Sum-A_Bar)/(M-s);
    A_Bar = (A_Bar)/s;
    ASE = Schein(A_Bar, 788, schein, myoptions)
    for (Z in 1:6){
      k = CU[Z];
      indz = which(abs(ASE$lamda) >= sort( abs(ASE$lamda), decreasing=T)[k], arr.ind=TRUE)
      P_Hat = ASE$X[,indz] %*% t(ASE$X[,indz]);
      P_Hat[P_Hat>1] =1;
      P_Hat[P_Hat<0] = 0;
      err_v <- P_Bar-P_Hat;
      err_v <- err_v[upper.tri(err_v)];
      err_P_Hat[i,j,Z] <- mean(err_v^2);
    }
    err_v <- P_Bar-A_Bar;
    err_v <- err_v[upper.tri(err_v)];
    err_A_Bar[i,j] <- mean(err_v^2);
  }
  MSE_ABar[i] = sum(err_A_Bar[i,])/nTrials;
  for (z2 in 1:6){
    MSE_PHat[i,z2] = sum(err_P_Hat[i,,z2])/nTrials;
  }
}

# for (i in 1:M){
#   samp = i;
#   s =1;
#   A_Bar = A_all[,,samp];
#   A_Bar = A_Bar+t(A_Bar)
#   P_Bar = (A_Sum-A_Bar)/(M-s);
#   A_Bar = (A_Bar)/s;
#   ASE = Schein(A_Bar, k, 15, myoptions)
#   P_Hat = ASE$X %*% t(ASE$X);
#   P_Hat[P_Hat>1] =1;
#   P_Hat[P_Hat<0] = 0;
#   err_v <- P_Bar-P_Hat;
#   err_v <- err_v[upper.tri(err_v)];
#   err_P_Hat_l[i] <- mean(err_v^2);
#   
#   err_v <- P_Bar-A_Bar;
#   err_v <- err_v[upper.tri(err_v)];
#   err_A_Bar_l[i] <- mean(err_v^2);
# }

# rowSums(err_P_Hat)
# rowSums(err_A_Bar)
# 
# par(mfrow=c(1,1));
# plot(test_size, rowSums(err_P_Hat)/500, ylab = "MSE", xlab = "M", type ='b', pch =19, col = 'red')
# lines(test_size, rowSums(err_A_Bar)/500, type ='b', pch =17, col = 'blue')
# legend(locator(1),legend =c("A_Bar", "P_Hat"),pch = c(17,19), col = c("blue", "red"))
par(mfrow=c(1,1));
matplot(test_size, MSE_PHat[,1:6], ylab = "MSE", xlab = "M", type ='b', pch = c(15,18,19,25,8), col = c('red','blue','green','purple','cyan', 'coral'))
lines(test_size, MSE_ABar, type ='b', pch =17, col = 'black')

legend(locator(1),legend =c("A_Bar", "P_Hat: d=1", "P_Hat: d=7","P_Hat: d=26","P_Hat: d=285","P_Hat: d=587","P_Hat: d=725"),pch = c(17,15,18,19,25,8), col =  c('black','red','blue','green','purple','cyan', 'coral'))


###############################################
library("mclust")
A_Bar_n = (A_Sum)/M;
k = findCutoff(A_Bar_n,2);
ASE <- Schein(A_Bar_n, k, schein, myoptions)
ASE_n = eigen(A_Bar_n);
# X = ASE_n$vectors %*% diag(sqrt(abs(ASE_n$values)));
ASE_Xn = matrix(0, dim(ASE$X)[1],dim(ASE$X)[2])
A_barX_n = matrix(0, dim(ASE_n$vectors)[1],dim(ASE_n$vectors)[2])
cl <- Mclust(ASE$X,7);

count = 1;
for (j in 1:7){
  for (k in 1:n){
    if (cl$classification[k]==j){
      ASE_Xn[count,] = ASE$X[k,];
      A_barX_n[count,] = ASE_n$vectors[k,]
      count = count +1;
    }
  }
}


A_bar_n = A_barX_n %*% diag(ASE_n$values) %*%t(A_barX_n)
A_bar_n[A_bar_n>1] = 1;
A_bar_n[A_bar_n<0] = 0;

P_hat = ASE_Xn %*% t(ASE_Xn)
P_hat[P_hat>1] = 1;
P_hat[P_hat<0] = 0;
diag(A_bar_n)=0;
diag(P_hat)=0;
se= n/7;
separ = matrix(1, round(se),n);

UnWeighted = rbind(A_bar_n, separ,P_hat)
par(mfrow=c(1,1));
image.plot(UnWeighted[,ncol(UnWeighted):1], col = gray((0:100)/100), asp =n/(2*n+1*se), axes=FALSE)
title("A_Bar                                            P_Hat")
# title("A_Bar                                                                              P_Hat")
difff = A_bar_n-P_hat;
image.plot(difff[,ncol(difff):1], col = gray((0:100)/100), asp =1, axes=FALSE)
title("A_Bar - P_Hat")
