logdet <- function(A) {
  y <- try(2*sum(log(diag(chol(A)))));
  if (class(y) == "try-error") {
    y = 0;
  }
  return(y);
}

approxeq <- function(a,b, tol = 1e-2, rel = FALSE){
  a = as.vector(a);
  b = as.vector(b);
  if (length(a)!=length(b)){
    return(FALSE);
  }
  d = abs(a-b);
  if (rel){
    p = sum((d/(abs(a)+2.2204e-16))>tol);
  } else {
    p = sum(d>tol);
  }
  return(p==0);
}

gaussLogprob <- function(mu, Sigma, X){
  mu = as.vector(mu);
  if (length(mu)==1){
    X = as.vector(X);
  }
  
  if (!is.vector(X)) {
    N = dim(X)[1];
    d = dim(X)[2];
  } else {
    N = length(X)
    d = 1;
  }

  if (d == 1){
    X = as.vector(X) - mu;
  } else {
    X = sweep(X,2,mu)
  }
  if (!is.vector(Sigma)){
    NS = dim(Sigma)[1];
    dS = dim(Sigma)[2];
  } else {
    NS = length(Sigma)
    dS = 1;
  }

  if ((dS == 1) && (NS >1)){
    sig2 = matrix(rep(as.vector(Sigma),N), N, length(Sigma), byrow=TRUE);
    tmp = -(X^2)/(2*sig2)-.5*log(2*pi*sig2);
    logp = rowSums(tmp);
  } else {
    R = chol(Sigma);
    logp = -.5*rowSums((X %*% solve(R))*X);
    logZ = .5*d*log(2*pi)+sum(log(diag(R)));
    logp = logp-logZ;
  }
  return(logp);
}

screeplotChooseDim <- function( lambdas){
  Lmax = length(lambdas);
  ndx = 1:Lmax;
  ll = c();
  for (q in 1:(Lmax-1)){
    group1 = ndx<=q;
    group2 = ndx >q;
    mu1 = mean(lambdas[group1]);
    v1 = max(0,var(lambdas[group1]), na.rm=TRUE);
    len1 = sum(group1);
    v1 = v1*(len1-1)/len1;
    mu2 = mean(lambdas[group2]);
    v2 = max(0,var(lambdas[group2]), na.rm=TRUE);
    len2 = sum(group2);
    v2 = v2*(len2-1)/len2;
    v = (len1*v1+len2*v2)/Lmax;
    ll[q] = sum(gaussLogprob(mu1,v,lambdas[group1])) + sum(gaussLogprob(mu2,v,lambdas[group2])) 
  }
  return(which.max(ll));
}

findCutoff <- function(M, k){
  lambdas = eigen(M,only.values=TRUE);
  lambdas = sort(abs(lambdas$values),decreasing=TRUE);
  cutoff = 0;
  for (i in 1:k){
    cut = screeplotChooseDim(lambdas[(cutoff+1):length(lambdas)]);
    cutoff = cutoff+cut;
  }
  return(cutoff);
}