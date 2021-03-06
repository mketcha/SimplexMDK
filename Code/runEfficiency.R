# Relative Efficiency Analysis
setwd("~/Simplex/SimplexMDK/Code")
library("igraph")
library("vioplot")
library("rARPACK")
library("mclust")
library("nloptr")
library("fields")
source("graphFunctions.R")
source("generateDirilechtSBM.R");
library("doParallel")
library("foreach")


mcl = makeCluster(3);
registerDoParallel(mcl);

MN_M = matrix(c(100,5,100,10, 100, 25, 100,50, 100,75,100,100,100,250,100, 500),8,2,byrow=TRUE)
MN_N = matrix(c(1,100,10,100,50,100,100, 100,500,100,1000,100,2000,100, 5000,100),8,2,byrow=TRUE)
#MN = matrix(c(1,100,100,100),4,2,byrow=TRUE)
#B <- matrix(c(.22, .22, .22,.6),2);
#B <- matrix(c(.42, .42, .42,.5),2);
B <- matrix(c(.42, .2, .2, .7),2);
###############B <- matrix(c(.42, .2, .3, .2,.5, .15, .3, .15, .7),3);
#B <- matrix(c(.42, .23, .35, .23,.5, .15, .35, .15, .7),3);
#B <- matrix(c(.42, .2, .3 ,.05, .2, .49, .35, .15, .3, .35, .6, .25, .05 , .15,.25, .7),4);
#PI <- c( .25, .25, .25, .25);
###############PI <- c(.33, .33,.34)
PI <- c(.05,.95);
trials <- 1000;
max_k <- 4;
nclust <- 3;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
eff_M = c();
eff_N = c();
eff_M_S = c();
eff_N_S = c();

#Varying N
for (run in 1:8) {
  n <- MN_M[run,2];
  m <- MN_M[run,1];
  error_Phat_noClust <- matrix(0, trials, max_k);
  error_Phat <- matrix(0, trials, max_k);
  error_mean <- matrix(0, trials, 1);
  for (i in 1:trials) {
    if (i%%100 ==0)
      print(i)
    Z <-generateZ_nxk(PI,n)
    P = Z %*% B %*% t(Z);
    P_Bar <- matrix(0,n,n);
    if (m < 500 && n < 100){
      for (h in 1:m) {
        g <- P2Adj(P);
        P_Bar <- P_Bar +g;
      }
    } else {
      P_Bar <- foreach(mb = breakScalar(m, 3)) %dopar%{
            Ps_i = matrix(0,n,n);
            for (h in 1:mb) {
              g <- P2Adj(P);
              Ps_i <- Ps_i +g;
            }
            return(Ps_i)
      }
      P_Bar <- Reduce('+',P_Bar);
    }

    P_Bar <- P_Bar/m;
    err_v <- P-P_Bar;
    err_v <- err_v[upper.tri(err_v)];
    error_mean[i,1] <- mean(err_v^2);
    j = 2;
    #error_mean[i,1] <- norm(err_v, type="2");
    diag(P_Bar) <- diag(P);
    #for (j in 1:max_k) {
     repeat{
        ASE <- try(spectEmbed(P_Bar,j, DAdjust=FALSE, opt = myoptions))
        if (class(ASE) != "try-error") {
          break
        }
        print("err")
      }
      P_hat_noClust <- ASE$X %*% t(ASE$X);
      err_v_noClust <- P-P_hat_noClust;
      err_v_noClust <- err_v_noClust[upper.tri(err_v_noClust)];
      error_Phat_noClust[i,j] <- mean(err_v_noClust^2);

      
    #}
  }
  print(mean(error_mean));
  print(mean(error_Phat_noClust[,j]))
  print(mean(error_Phat[,j]))
  eff_M = cbind(eff_M,mean(error_Phat_noClust[,j])/mean(error_mean))
}

# Varying M
for (run in 1:8) {
  n <- MN_N[run,2];
  m <- MN_N[run,1];
  error_Phat_noClust <- matrix(0, trials, max_k);
  error_Phat <- matrix(0, trials, max_k);
  error_mean <- matrix(0, trials, 1);
  eff_N_bySample <- matrix(0, trials, 1);
  for (i in 1:trials) {
    if (i%%100 ==0)
      print(i)
    Z <-generateZ_nxk(PI,n)
    P = Z %*% B %*% t(Z);
    P_Bar <- makeP_Bar(P,m,n)
    err_v <- err_v[upper.tri(err_v)];
    error_mean[i,1] <- mean(err_v^2);
    #error_mean[i,1] <- norm(err_v, type="2");
    diag(P_Bar) <- diag(P);
    #for (j in 1:max_k) {
    j =3;
      repeat{
        ASE <- try(spectEmbed(P_Bar,j, DAdjust=TRUE, opt = myoptions))
        if (class(ASE) != "try-error") {
          break
        }
       # print("err")
      }
      P_hat_noClust <- ASE$X %*% t(ASE$X);
      err_v_noClust <- P-P_hat_noClust;
      err_v_noClust <- err_v_noClust[upper.tri(err_v_noClust)];
      error_Phat_noClust[i,j] <- mean(err_v_noClust^2);
      eff_N_bySample[i,1] <- mean(err_v_noClust^2)/error_mean[i,1];
      
    #}
  }
  print(mean(error_mean));
  print(mean(error_Phat_noClust[,j]))
  print(mean(error_Phat[,j]))
  eff_N_S = cbind
  eff_N = cbind(eff_N,mean(error_Phat_noClust[,j])/mean(error_mean))
}

d_M = data.frame(
  N = MN_M[,2])
 ,ASE_M = eff_M
 ,NMF_M = eff_NMF_M
)