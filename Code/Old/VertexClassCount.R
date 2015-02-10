VertexClassCount <- function(PI_1xk, n){
  r <- runif(n);
  Z_k <- rep(0,length(PI_1xk))
  PI_c = c(0,cumsum(PI_1xk));
  for (i in 1:(length(PI_c)-1)){
    Z_k[i] <- sum(r>=PI_c[i] & r<PI_c[i+1]);
    if (Z_k[i] ==0)
      return(VertexClassCount(PI_1xk,n))
  }
  return(Z_k)
}