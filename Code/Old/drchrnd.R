drchrnd <- function(a,n){
  p <- length(a);
  d <- matrix(0,n,p);
  for (i in 1:length(a)){
    d[,i] <- rgamma(n,a[i]);
  }
  rs <-rowSums(d);
  n <- replicate(p,rs);
  return(d/n)
}