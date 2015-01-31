spectEmbedUndirected = function(M,no, opt = igraph.arpack.default) {
  sp = eigs_sym(M, no, options = opt)
  if (no == 1){
    X = sp$vectors*sqrt(abs(sp$values))
  } else {
    X = sp$vectors %*% diag(sqrt(abs(sp$values)));
  }
  ret = list(X=X, lamda = sp$values)
  return(ret)
}