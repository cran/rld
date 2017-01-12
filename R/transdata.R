
#transfer input data to required data structure

transdata <- function(data, ndlevel, nexposure){
  
  sdata <- data
  n <- nrow(sdata)
  varnames <- names(sdata)
  
  indxtime <- match("time", varnames, nomatch = 0)
  if(indxtime==0) stop("Variable 'time' is required in the dataset")
  
  indxdelta <- match("delta", varnames, nomatch = 0)
  if(indxdelta==0) stop("Variable 'delta' is required in the dataset")
  
  time <- sdata[,indxtime]
  delta <- sdata[,indxdelta]
  
  Xmat <- sdata[,-c(indxtime, indxdelta)]
  ncolXmat <- ifelse(is.null(ncol(Xmat)), 1, ncol(Xmat))
  
  id <- c(1:n)
  
  repdata <- matrix(NA, nrow = sum(time), ncol = 4+ncolXmat)
  colnames(repdata) <- c("id", "time", "delta", "dose", varnames[-c(indxtime,indxdelta)])
  
  for (i in 1:n){
    Ci <- time[i]
    deltai <- delta[i]
    dose <- rep(1:ndlevel, nexposure)[1:Ci]
    
    if (ncolXmat==1){
      X <- Xmat[i]
      augX <- rep(X, times = Ci)
    } else if (ncolXmat>1){
      X <- Xmat[i,]
      augX <- X[rep(seq(nrow(X)), each = Ci),]
    }
    
    rowstart <- min(which(is.na(repdata[,1])==TRUE))
    rowend <- rowstart+Ci-1
    
    repdata[rowstart:rowend, 1] <- id[i]
    repdata[rowstart:rowend, 2] <- c(1:Ci)
    
    if(deltai>0){
      repdata[rowstart:(rowend-1), 3] <- 0
      repdata[rowend, 3] <- 1
    }else{
      repdata[rowstart:rowend, 3] <- 0
    }
    
    repdata[rowstart:rowend, 4] <- dose
    repdata[rowstart:rowend, 5:(4+ncolXmat)] <- as.matrix(augX)
  }
  
  #correction
  for(i in 1:n){
    deltai <- delta[i]
    lastipoi <- sum(sdata$time[1:i])
    if(repdata[lastipoi,3]!=deltai) {repdata[lastipoi,3] <- deltai}
  }
  
  repdata <- data.frame(repdata)
  return(repdata)
}
