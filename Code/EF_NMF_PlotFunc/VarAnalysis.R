source("../graphFunctions.R")
library("combinat")
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
B <- matrix(c(.4, .2, .2, .7),2);
PIs = matrix(c(.99,.01, .95,.05,.9,.1,.75,.25,.5,.5,.25,.75,.1,.9,.05,.95,.01,.99),9,2,byrow=TRUE)
ratio = matrix(0,9,3);
VarPij = matrix(0,9,3);
dan = matrix(0,9,3);
for (run in 1:9){
  PI <- PIs[run,];
  k <- dim(B)[1];
  n = 10000;
  Z <-generateZ_nxk_exact(PI,n);
  P = Z %*% B %*% t(Z);
  ind  = c(1:k);
  for (i in 1:k){
    ind[i] = match(1, Z[,i])
  } 
  ASE <- spectEmbed(P, k, DAdjust=FALSE, opt = myoptions)
  
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
  }
  dan[run,] = c(1/PI[1] +1/PI[1], 1/PI[1] +1/PI[2], 1/PI[2] +1/PI[2])
  VarPij[run,] = NMVarP
  ratio[run,] = NMVarP/MVarA;
}