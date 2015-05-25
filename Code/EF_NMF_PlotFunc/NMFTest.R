setwd("~/Simplex/SimplexMDK/Code/EF_NMF_PlotFunc/")
library("NMF")
source("../graphFunctions.R")
B = matrix(c(.42, .2, .2, .7),2);
PI = c(.5,.5)
n = 100;
m = 50;
myoptions <- igraph.arpack.default
myoptions$maxiter <- 20000
Z <-generateZ_nxk(PI,n)
P = Z %*% B %*% t(Z);
P_Bar = makeP_Bar(P,m,n)
ASE = nmf(P_Bar,2, "brunet")

Orig <- basis(ASE) %*% coef(ASE);

diag(P_Bar) <- 20000000000000;# rowSums(P_Bar)/(dim(P_Bar)[1]-1);

ASE_aug = nmf(P_Bar,2, '.M#brunet', seed  = "nndsvd")
Aug <- basis(ASE_aug) %*% coef(ASE_aug);

ASE_aug2 = nmf(P_Bar,2, "brunet")

Orig_Aug <- basis(ASE_aug2) %*% coef(ASE_aug2);
# err = c()
# for (i in 1:20){
#   ASE = Schein(P_Bar, 2, i, myoptions)
#   P_hat_noClust <- ASE$X %*% t(ASE$X);
#   err_v_noClust <- P-P_hat_noClust;
#   err_v_noClust <- err_v_noClust[upper.tri(err_v_noClust)];
#   err[i] = mean(err_v_noClust^2);
# }
