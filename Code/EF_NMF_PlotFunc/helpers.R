source('../graphFunctions.R')
makeP_Bar <- function(P, n, m){
  P_Bar <- matrix(0,n,n);
  if (m < 50 && n < 25){
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
}