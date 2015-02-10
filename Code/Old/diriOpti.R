diriOpti <- function(x, nVertex, dimLatentPosition, xHatTmp){
  rotation <- matrix(x, dimLatentPosition, dimLatentPosition);
  xHat <- matrix(xHatTmp,dimLatentPosition,dimLatentPosition) %*% rotation;
  val = -min(abs(xHat));
  return(val)
}

diriIneq <- function(x, nVertex, dimLatentPosition, xHatTmp){
  rotation <- matrix(x, dimLatentPosition, dimLatentPosition);
  xHat <- matrix(xHatTmp,dimLatentPosition,dimLatentPosition) %*% rotation;
  Ineq <- c();
  for (i in 1:nVertex){
    for (j in 1:dimLatentPosition){
      Ineq <- c(Ineq,-xHat[j,i]);
    }
  }
  return(Ineq)
}

diriEq <- function(x, nVertex, dimLatentPosition, xHatTmp){
  rotation <- matrix(x, dimLatentPosition, dimLatentPosition);
  xHat <- matrix(xHatTmp,dimLatentPosition,dimLatentPosition) %*% rotation;
  Eq <- c();
  for (i in 1:dimLatentPosition){
    for (j in 1:dimLatentPosition){
      if (i==j){
        Eq <- c(Eq, (rotation[j,, drop=FALSE] %*% t(rotation[i,,drop=FALSE])) -1);
      } else {
        Eq <- c(Eq, rotation[j,,drop=FALSE] %*% t(rotation[i,,drop=FALSE]));
      }
    }
  }
  return(Eq)
}