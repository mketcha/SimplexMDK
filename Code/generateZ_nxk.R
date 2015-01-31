generateZ_nxk <- function(PI_1xk, n){
  r <- runif(n);
  Z_n <- rep(0,n)
  PI_c = c(0,cumsum(PI_1xk));
  for (i in 1:(length(PI_c)-1)){
    Z_n[r>=PI_c[i] & r<PI_c[i+1]] = i;
  }
  Z_nxk <- matrix(0, n, length(PI_1xk));
  for (i in 1:length(Z_n)) {
    Z_nxk[i,Z_n[i]] = 1;
  }
  return(Z_nxk)
}