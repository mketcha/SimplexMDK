# Cross Validation Look at BLSA data
# MDK
rm(list = ls())
setwd("~/Simplex/SimplexMDK/Code/CoRR")
library("igraph")
library("rARPACK")
library("fields")
library("R.matlab")
library("doSNOW")
library("foreach")
source("ZhuGhodsi.R")
source("../graphFunctions.R")

cl<-makeCluster(5) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

subjects = readLines("../../../Data/CoRR/subnames.txt")
pheno = read.csv("../../../Data/CoRR/Phenotype.txt")
n = 788;
M = 232*2;
NPerm = 10;
popSizes = c(20, 30, 50, 100, 200)
sex = c();
Males = c();
Females = c();
A_all <- array(rep(0, n*n*M), dim=c(n, n, M))
for (sub in 1:232) {
  for (session in 1:2) {
    sex[(sub-1)*2 +session] = as.numeric(pheno$SEX[match(as.numeric(subjects[sub]), pheno$SUBID)]) -1;
    if (sex[(sub-1)*2 +session] == 1){
      Females = c(Females, (sub-1)*2 +session)
    } else {
      Males = c(Males, (sub-1)*2 +session)
    }
    Mat <- readMat(paste("../../../Data/CoRR/", subjects[sub], "/session_", session,"/", "smg_thresholded_0.125_788roi_vR_subject_", subjects[sub], "_session_",session,".mat", sep=""));
    Mat$threshSMG[Mat$threshSMG>0]<-1;
    A_all[,, (sub-1)*2 +session] <- Mat$threshSMG;
  }
}

schein = 5;
k = 7;
# permDistHat = matrix(0, length(popSizes), NPerm);
# permDistBar = matrix(0, length(popSizes), NPerm);


# trueDistHat = c();
# trueDistBar = c();
# pValsHat = c();
# pValsBar = c();

#for (i in 1:length(popSizes)){
permResults <- foreach(i=1:NPerm, .export = c("eigs")) %dopar% { 
  print(i);
  s = popSizes[i];
  permDistHat = rep(0, NPerm);
  permDistBar = rep(0, NPerm);
  indM = sample(length(Males),s/2);
  indF = sample(length(Females), s/2);
  trueM = rowSums(A_all[,,Males[indM]], dims = 2)/length(indM);
  trueF = rowSums(A_all[,,Females[indF]], dims = 2)/length(indF);
  pop = c(Males[indM], Females[indF]);
  for (j in 1:NPerm){
    if ((j %% 100)==0){
      print(c(i, j))
    }
    sM = sample(s, s/2)
    permM = rowSums(A_all[,,pop[sM]], dims = 2)/length(sM);
    permF = rowSums(A_all[,,pop[-sM]], dims = 2)/length(sM);
    ASEM = Schein(permM, k, schein, myoptions)
    P_HatM = ASEM$X %*% t(ASEM$X);
    ASEF = Schein(permF, k, schein, myoptions)
    P_HatF = ASEF$X %*% t(ASEF$X);
    
    err_v_hat <- P_HatF-P_HatM;
    err_v_hat <- err_v_hat[upper.tri(err_v_hat)];
    #permDistHat[i,j] <- mean(err_v_hat^2);
    permDistHat[j] <- mean(err_v_hat^2);
    
    err_v_bar <- permF-permM;
    err_v_bar <- err_v_bar[upper.tri(err_v_bar)];
    #permDistBar[i,j] <- mean(err_v_hat^2);
    permDistBar[j] <- mean(err_v_hat^2);
  }
  ASEM = Schein(trueM, k, schein, myoptions)
  P_HatM = ASEM$X %*% t(ASEM$X);
  ASEF = Schein(trueF, k, schein, myoptions)
  P_HatF = ASEF$X %*% t(ASEF$X);
  
  err_v_hat <- P_HatF-P_HatM;
  err_v_hat <- err_v_hat[upper.tri(err_v_hat)];
  #trueDistHat[i] <- mean(err_v_hat^2);
  trueDistHat <- mean(err_v_hat^2);
  
  err_v_bar <- trueF-trueM;
  err_v_bar <- err_v_bar[upper.tri(err_v_bar)];
  #trueDistBar[i] <- mean(err_v_hat^2);
  trueDistBar <- mean(err_v_hat^2);
  
  #   pValsHat = sum(trueDistHat[i]<permDistHat[i,])/NPerm;
  #   pValsBar = sum(trueDistBar[i]<permDistBar[i,])/NPerm;
  pValsHat = sum(trueDistHat<permDistHat)/NPerm;
  pValsBar = sum(trueDistBar<permDistBar)/NPerm;
  
  permResults <- list(pValsHat = pValsHat, pValsBar = pValsBar, trueDistHat = trueDistHat, trueDistBar = trueDistBar, permDistHat = permDistHat, permDistBar = permDistBar)
  return( permResults)
  
}

stopCluster(cl)