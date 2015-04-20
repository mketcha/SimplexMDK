# Michael Ketcha
# 

source("../graphFunctions.R")
library("combinat")
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
B <- matrix(c(.42, .2, .3, .2,.5, .15, .3, .15, .7),3);
k <- dim(B)[1];
PIs = matrix(c(1/3, 1/3, 1/3,
               1/4, 1/4, 2/4,
               1/6, 1/3, 1/2,
               1/10, 3/10, 3/5),
               4, 3,byrow=TRUE)
runs = dim(PIs)[1]
ratio = matrix(0,runs,k*(k+1)/2);
VarPij = matrix(0,runs,k*(k+1)/2);
dan = matrix(0,runs,k*(k+1)/2);
for (run in 1:runs){
  PI <- PIs[run,];
  n = 10000;
  Z <-generateZ_nxk_exact(PI,n);
  P = Z %*% B %*% t(Z);
  ind  = c(1:k);
  for (i in 1:k){
    ind[i] = match(1, Z[,i])
  } 
  ASE <- spectEmbed(P, 5, DAdjust=FALSE, opt = myoptions)
  
  Xi = ASE$X[ind,]
  
  Sig = rep(list(diag(k),k))
  for (i in 1:k){
    Sig[[i]] = calcCov(ASE$X, Xi[i, ,drop = FALSE])
  }
  
  comb = matrix(0, k*(k+1)/2, 2)
  count = 1;
  for (i in 1:k){
    for (j in i:k){
      comb[count,] = c(i,j);
      count =count+1
    }
  }
  
  NMVarP = c()
  MVarA = c()
  
  for (i in 1:dim(comb)[1]) {
    MVarA[i] = (t(Xi[comb[i,1],]) %*% Xi[comb[i,2],])*(1-t(Xi[comb[i,1],]) %*% Xi[comb[i,2],])
    NMVarP[i] = t(Xi[comb[i,1],]) %*% Sig[[comb[i,2]]] %*% Xi[comb[i,1],] + t(Xi[comb[i,2],]) %*% Sig[[comb[i,1]]] %*% Xi[comb[i,2],]
    dan[run,i] = 1/PI[comb[i,1]]+1/PI[comb[i,2]]
  }

  VarPij[run,] = NMVarP
  ratio[run,] = NMVarP/MVarA;
}