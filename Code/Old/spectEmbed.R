spectEmbed <- function(M,no, opt = igraph.arpack.default) {
  if (length(M[,1])==2){
    sp <- eigen(M);
    if (no == 1){
      X = sp$vectors[,1]*sqrt(abs(sp$values[1]))
      lambda = sp$values[1];
    } else if (no==2){
      X = sp$vectors %*% diag(sqrt(abs(sp$values)));
      lambda = sp$values
    } else {
      stop
    }
    ret = list(X=X, lamda = lambda)
    return(ret)
    
  } else{
    sp = eigs(M, no, options = opt)
    if (no == 1){
      X = sp$vectors*sqrt(abs(sp$values))
    } else {
      X = sp$vectors %*% diag(sqrt(abs(sp$values)));
    }
    ret = list(X=X, lamda = sp$values)
    return(ret)
  }
}